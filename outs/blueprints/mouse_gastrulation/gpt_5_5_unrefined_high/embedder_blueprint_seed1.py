#!/usr/bin/env python3
"""Declarative count-matrix embedder blueprint.

This module exposes EmbedderBlueprint, a small declarative builder for
mean-expression cell-type embeddings, and a CLI that reads an AnnData
count matrix and writes a parquet of axis values per cell.

AGENTS: EDIT ONLY `build_blueprint` at the bottom of this file. The
EmbedderBlueprint API is documented in per-method docstrings on the
class. Verification:

    python3 embedder_blueprint.py --describe

This must exit cleanly before the task is done.
"""
from __future__ import annotations

import argparse
import difflib
import re

import numpy as np
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class ModuleSpec:
    """A named set of genes (e.g. marker genes, or genes from a shared pathway)."""

    name: str
    genes: tuple[str, ...]
    resource_keys: tuple[str, ...] = ()


@dataclass(frozen=True)
class Resource:
    """A literature resource (paper, atlas, or curated database).

    Resources are declared once per blueprint via add_resource() and
    cited by key from gene modules. DOI is required for paper/atlas
    (the canonical handle) and optional for database (the gene-set
    identifier in `citation` is the canonical handle there).
    """

    key: str
    kind: str           # "paper" | "atlas" | "database"
    citation: str       # "Author Year Journal" or "DB gene-set ID"
    doi: str | None     # "10.xxxx/yyyy"; may be None when kind="database"


@dataclass(frozen=True)
class AxisSpec:
    """One dimension of the output embedding.

    A plain axis has empty `negatives` and no `gate`.
    A contrast axis has exactly one negative and `negative_weight=1.0`.
    A gated axis has a `gate` module name and typically
    `negative_weight=0.5`.
    """

    name: str
    positive: str
    negatives: tuple[str, ...] = ()
    negative_weight: float = 1.0
    gate: str | None = None


# Crossref's recommended DOI validator: prefix "10." + 4-9 digits,
# slash, then suffix. Bare URLs and "doi:" prefixes are rejected.
_DOI_PATTERN = re.compile(r"^10\.\d{4,9}/[-._;()/:A-Za-z0-9]+$")
_VALID_RESOURCE_KINDS = frozenset({"paper", "atlas", "database"})

# Sentinels in the skeleton's commented `<<EXAMPLE_*>>` / `_DELETE_ME`
# scaffolding. Any of these in a resource / module / gene / axis name
# means the agent left placeholder text in place.
_PLACEHOLDER_MARKERS = ("<<", ">>", "EXAMPLE_", "_DELETE_ME", "PLACEHOLDER")


class EmbedderBlueprint:
    """Declarative builder for a count-matrix cell embedding."""

    def __init__(self, name: str) -> None:
        self._name = name
        self._modules: dict[str, ModuleSpec] = {}
        self._axes: list[AxisSpec] = []
        self._clip: float | None = 4.0
        self._resources: dict[str, Resource] = {}

    @property
    def name(self) -> str:
        return self._name

    @property
    def modules(self) -> dict[str, ModuleSpec]:
        return dict(self._modules)

    @property
    def axes(self) -> tuple[AxisSpec, ...]:
        return tuple(self._axes)

    @property
    def resources(self) -> dict[str, Resource]:
        return dict(self._resources)

    def add_resource(
        self,
        key: str,
        *,
        kind: str,
        citation: str,
        doi: str | None = None,
    ) -> None:
        """Register a literature resource that gene modules can cite.

        `kind` must be one of "atlas", "database", or "paper". `doi`
        is required for `paper` / `atlas` (must match `10.xxxx/yyyy`)
        and optional for `database` — for databases the gene-set
        identifier inside `citation` (e.g. `HALLMARK_GLYCOLYSIS`,
        `KEGG hsa04060`, `GO:0007049`) is the canonical handle.

        Example:
            blueprint.add_resource(
                "key1",
                kind="paper",
                citation="Author1 Year1 Journal1",
                doi="10.xxxx/yyyy",
            )
        """
        # Dedupe-first matches the ordering in add_gene_module: identifier
        # collision wins over field validation so callers get a focused
        # error when re-registering an existing key.
        if not isinstance(key, str) or not key:
            raise ValueError("key must be non-empty for resource")
        if key in self._resources:
            raise ValueError(f"Resource already registered: {key}")
        if kind not in _VALID_RESOURCE_KINDS:
            raise ValueError(
                f"invalid kind: {kind!r}; must be one of "
                f"{sorted(_VALID_RESOURCE_KINDS)}"
            )
        if not isinstance(citation, str) or not citation:
            raise ValueError(f"citation must be non-empty for resource: {key}")
        if kind in ("paper", "atlas"):
            if not isinstance(doi, str) or not doi:
                raise ValueError(
                    f"doi is required for kind={kind!r} resource: {key}"
                )
        if doi is not None:
            if not isinstance(doi, str) or not _DOI_PATTERN.match(doi):
                raise ValueError(f"invalid DOI for resource {key}: {doi!r}")
        self._resources[key] = Resource(
            key=key, kind=kind, citation=citation, doi=doi,
        )

    def add_gene_module(
        self,
        name: str,
        genes: Iterable[str],
        *,
        resources: Iterable[str],
    ) -> None:
        """Register a named set of genes (a "module").

        Modules are the atoms that axes combine. Use any name you like;
        the name is what `add_axis` and friends refer to.

        `resources` is a non-empty iterable of resource keys previously
        registered via `add_resource(...)`. The registry refuses to
        register a module without at least one resource, and refuses
        any unknown resource key.

        Example:
            blueprint.add_gene_module(
                "module1",
                ["gene1", "gene2", "gene3"],
                resources=["key1", "key2"],
            )
        """
        if name in self._modules:
            raise ValueError(f"Module already registered: {name}")
        keys = tuple(resources)
        if not keys:
            raise ValueError(f"resources must be non-empty for module: {name}")
        for key in keys:
            if key not in self._resources:
                raise ValueError(f"unknown resource: {key}")
        self._modules[name] = ModuleSpec(
            name=name, genes=tuple(genes), resource_keys=keys,
        )

    def add_axis(
        self,
        name: str,
        module: str,
    ) -> None:
        """Register a plain axis whose value is the module's own score.

        Use when a single program cleanly separates one population.

        Example:
            blueprint.add_axis("axis1", module="module1")
        """
        self._require_module(module)
        self._store_axis(AxisSpec(name=name, positive=module))

    def add_contrast_axis(
        self,
        name: str,
        positive_module: str,
        negative_module: str,
    ) -> None:
        """Register a contrast axis: value = score(pos) - score(neg).

        Use to contrast different groups of cells, for example when two
        programs occupy opposite ends of one biological dimension (e.g.
        excitatory vs inhibitory).

        Example:
            blueprint.add_contrast_axis(
                "module1_vs_module2",
                positive_module="module1",
                negative_module="module2",
            )
        """
        self._require_module(positive_module)
        self._require_module(negative_module)
        self._store_axis(
            AxisSpec(
                name=name,
                positive=positive_module,
                negatives=(negative_module,),
                negative_weight=1.0,
            )
        )

    def add_gated_axis(
        self,
        name: str,
        gate_module: str,
        positive_module: str,
        negative_modules: Iterable[str] = (),
        negative_weight: float = 0.5,
    ) -> None:
        """Register a gated axis: value = score(gate) * (score(pos) - w * mean(negs)).

        The gate suppresses the axis in cells that do not express the
        gating program. Two valid uses:

        1. Subtype refinement within a class — pass `negative_modules`
           listing the sibling subtypes you want to discriminate against.
           Example: distinguishing cells of type A among cells of type B.

        2. State scoping — leave `negative_modules` empty. Use this when
           the positive program is a state (e.g. activity-dependent IEGs)
           that is not exclusive to one cell class, and the gate restricts
           the signal to that class. With no negatives, the axis reduces
           to score(gate) * score(pos).

        Example:
            blueprint.add_gated_axis(
                "module1_subtype1",
                gate_module="module1",
                positive_module="module3",
                negative_modules=["module4", "module5"],
                negative_weight=0.3,
            )
        """
        self._require_module(gate_module)
        self._require_module(positive_module)
        negs = tuple(negative_modules)
        for n in negs:
            self._require_module(n)
        self._store_axis(
            AxisSpec(
                name=name,
                positive=positive_module,
                negatives=negs,
                negative_weight=negative_weight,
                gate=gate_module,
            )
        )

    def set_clip(self, value: float | None) -> None:
        """Set the post-zscore clip threshold.

        Z-scored axis values are clipped to [-value, +value]. Pass None
        to disable clipping. Default is 4.0.

        Example:
            blueprint.set_clip(None)  # disable clipping
        """
        self._clip = value

    def _require_module(self, module_name: str) -> None:
        if module_name not in self._modules:
            raise ValueError(f"unknown module: {module_name}")

    def _store_axis(self, axis: AxisSpec) -> None:
        if any(existing.name == axis.name for existing in self._axes):
            raise ValueError(f"Axis already registered: {axis.name}")
        self._axes.append(axis)

    def _axis_score(
        self,
        axis: AxisSpec,
        score_by_module: dict[str, float],
    ) -> float:
        pos = score_by_module[axis.positive]
        if axis.negatives:
            neg_mean = sum(
                score_by_module[n] for n in axis.negatives
            ) / len(axis.negatives)
        else:
            neg_mean = 0.0
        inner = pos - axis.negative_weight * neg_mean
        if axis.gate is None:
            return inner
        gate = score_by_module[axis.gate]
        return gate * inner

    def _normalize(self, matrix: np.ndarray) -> np.ndarray:
        """Z-score columns; constant columns become all-zero (no divide by zero)."""
        mean = matrix.mean(axis=0)
        std = matrix.std(axis=0)
        safe_std = np.where(std == 0.0, 1.0, std)
        z = (matrix - mean) / safe_std
        if self._clip is not None:
            return np.clip(z, -self._clip, self._clip)
        return z

    def build_embedding(
        self,
        h5ad_path: Path | str,
        output_path: Path | str,
        *,
        layer: str | None = None,
        n_min_genes: int = 0,
    ) -> None:
        """Read an AnnData h5ad, score every cell against every
        registered axis using mean expression of module genes,
        normalize, and write a parquet with columns
        cell_id + dim_0..dim_{N-1}.

        If `layer` is given, expression values come from
        `adata.layers[layer]`; otherwise from `adata.X`. Values are
        used as-is (no library-size normalization, no log1p).

        Axes whose raw score is constant across all cells are preserved
        as all-zero columns (z-scoring a constant column yields zero;
        the column is still written so dim indices stay aligned with
        the order returned by ``axes``).

        `n_min_genes` (default 0) excludes any axis whose *positive*
        module has fewer than this many genes present in the dataset's
        var_names (case-insensitive). Surviving axes are compactly
        renumbered dim_0..dim_{M-1} in declaration order; the default
        of 0 is a no-op (present >= 0 always holds). The kept/dropped
        decision per axis is recorded in the coverage sidecar.

        Called by the CLI; agents do not call it directly.
        """
        h5ad_path = Path(h5ad_path)
        output_path = Path(output_path)
        if not h5ad_path.exists():
            raise FileNotFoundError(h5ad_path)
        if not self._axes:
            raise ValueError("no axes registered; build_blueprint must add at least one axis")

        import anndata as ad
        import pandas as pd
        import scipy.sparse as sp
        import time as _time

        # Fence the input read to report load_data_s separately.
        _t_read0 = _time.perf_counter()
        adata = ad.read_h5ad(h5ad_path)
        self._last_load_data_s = _time.perf_counter() - _t_read0
        try:
            import resource as _res, sys as _s
            _rss = (_res.getrusage(_res.RUSAGE_SELF).ru_maxrss
                    + _res.getrusage(_res.RUSAGE_CHILDREN).ru_maxrss)
            self._last_rss_after_load_mb = _rss / 1024.0**2 if _s.platform == "darwin" else _rss / 1024.0
        except Exception:
            self._last_rss_after_load_mb = None
        gene_names = adata.var_names.to_numpy()
        # Case-insensitive gene matching: the engine never writes a gene name,
        # it only indexes into the matrix, so dataset display-case is preserved
        # while the agent's casing is tolerated. Two genes that differ only in
        # case collapse to the first one's column (setdefault), matching the
        # de-facto behavior of `adata.var_names_make_unique`-style pipelines.
        gene_index: dict[str, int] = {}
        for i, g in enumerate(gene_names):
            gene_index.setdefault(g.lower(), i)

        X = adata.layers[layer] if layer is not None else adata.X
        is_sparse = sp.issparse(X)
        if is_sparse:
            X = X.tocsr()

        cell_ids = adata.obs_names.tolist()
        n_cells = adata.n_obs

        # Per-module mean expression, vectorised over all cells: subset the
        # matrix to the union of module genes once, densify that small slice,
        # then take per-module column means. Axis scores are then plain numpy
        # over length-n_cells vectors (see _axis_score), so there is no
        # per-cell Python loop. Missing genes are dropped from the average.
        union = sorted({gene_index[g.lower()] for m in self._modules.values()
                        for g in m.genes if g.lower() in gene_index})
        local = {orig: j for j, orig in enumerate(union)}
        Xd = np.asarray(X[:, union].toarray() if is_sparse else X[:, union], dtype=float)
        score_by_module: dict[str, np.ndarray] = {}
        for name, module in self._modules.items():
            cols = [local[gene_index[g.lower()]] for g in module.genes if g.lower() in gene_index]
            score_by_module[name] = Xd[:, cols].mean(axis=1) if cols else np.zeros(n_cells)

        # Optional axis exclusion: drop any axis whose positive module has
        # fewer than n_min_genes present in this dataset (case-insensitive
        # ∩ var_names). Default 0 keeps every axis. Surviving axes are
        # compactly renumbered dim_0..dim_{M-1} below; self._axes is left
        # intact so describe() and re-runs are unaffected.
        present_by_module = {
            name: sum(1 for g in module.genes if g.lower() in gene_index)
            for name, module in self._modules.items()
        }
        active_axes = [a for a in self._axes
                       if present_by_module[a.positive] >= n_min_genes]
        dropped_axes = [a for a in self._axes if a not in active_axes]
        if not active_axes:
            raise ValueError(
                f"all {len(self._axes)} axes excluded by n_min_genes={n_min_genes}; "
                "no axis's positive module has enough present genes for this dataset"
            )
        if dropped_axes:
            print(
                f"[n_min_genes] excluded {len(dropped_axes)}/{len(self._axes)} axes "
                f"(positive module <{n_min_genes} present genes): "
                + ", ".join(a.name for a in dropped_axes)
            )
        matrix = np.column_stack([self._axis_score(axis, score_by_module) for axis in active_axes])

        # Coverage telemetry: print one line per module that has >=1
        # gene missing from the matrix's var_names (compared case-insensitively),
        # with did-you-mean suggestions for typos. Never raises. The same
        # gene-presence counts and per-module score-nonzero fractions are also
        # collected for the structured sidecar written below.
        # Maintain a lower→original mapping so close-match suggestions can
        # come back in the dataset's display case.
        vocab_lower_to_orig: dict[str, str] = {}
        for g in gene_names.tolist():
            vocab_lower_to_orig.setdefault(g.lower(), g)
        sorted_vocab_lower = sorted(vocab_lower_to_orig.keys())
        coverage_modules: list[dict] = []
        for module in self._modules.values():
            missing = sorted(g for g in module.genes if g.lower() not in vocab_lower_to_orig)
            n_total = len(module.genes)
            n_present = n_total - len(missing)
            score = score_by_module[module.name]
            frac_nonzero = float((score > 0).mean()) if n_cells > 0 else 0.0
            coverage_modules.append({
                "name": module.name,
                "n_genes": n_total,
                "n_genes_present": n_present,
                "n_genes_missing": len(missing),
                "missing_genes": missing,
                "fraction_nonzero_cells": frac_nonzero,
            })
            if not missing:
                continue
            parts: list[str] = []
            for gene in missing:
                # Case-blind close-match search: compare lowercased gene
                # against lowercased vocab; then resolve back to the dataset's
                # original display case for the suggestion string.
                suggestions_lower = difflib.get_close_matches(
                    gene.lower(), sorted_vocab_lower, n=2, cutoff=0.7
                )
                suggestions = [vocab_lower_to_orig[s] for s in suggestions_lower]
                if suggestions:
                    sugg_str = ", ".join(f"'{s}'" for s in suggestions)
                    parts.append(f"{gene} (did you mean: {sugg_str})")
                else:
                    parts.append(gene)
            print(
                f"[coverage] module '{module.name}': "
                f"{n_present}/{n_total} genes present; "
                f"missing: {', '.join(parts)}"
            )

        normalized = self._normalize(matrix)

        frame = pd.DataFrame({"cell_id": cell_ids})
        for idx in range(normalized.shape[1]):
            frame[f"dim_{idx}"] = normalized[:, idx]

        output_path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_parquet(output_path, index=False)

        # Structured dataset-coverage sidecar next to the parquet. The
        # embedder is the source of truth for which declared genes are
        # actually used (case-sensitive ∩ var_names) and for the per-module
        # mean's non-zero rate across cells; extract_run.py reads this
        # back via its --parquet flag and folds it into summary.json.
        import json as _json
        coverage_record = {
            "h5ad_path": str(h5ad_path),
            "layer": layer,
            "n_cells": int(n_cells),
            "n_vars": int(len(gene_names)),
            "n_min_genes": int(n_min_genes),
            "clip": self._clip,
            "modules": coverage_modules,
            "axes": [
                {
                    "name": a.name,
                    "positive_module": a.positive,
                    "n_present": present_by_module[a.positive],
                    "kept": a in active_axes,
                    "dim": active_axes.index(a) if a in active_axes else None,
                }
                for a in self._axes
            ],
        }
        sidecar = output_path.parent / f"{output_path.stem}_coverage.json"
        sidecar.write_text(_json.dumps(coverage_record, indent=2) + "\n")

    def describe(self) -> str:
        """Return a human-readable listing of modules, axes, and config.

        Used by the `--describe` CLI flag.
        """
        lines = [f"Blueprint: {self._name}"]
        lines.append(f"  clip={self._clip}")
        lines.append(f"  modules: {len(self._modules)} registered")
        for module in self._modules.values():
            lines.append(
                f"    - {module.name} "
                f"({len(module.genes)} genes, "
                f"{len(module.resource_keys)} resource_keys)"
            )
        lines.append(f"  axes: {len(self._axes)} registered")
        for idx, axis in enumerate(self._axes):
            kind = self._axis_kind(axis)
            lines.append(
                f"    - dim_{idx}: {axis.name} [{kind}] "
                f"{self._axis_formula(axis)}"
            )
        return "\n".join(lines)

    @staticmethod
    def _axis_kind(axis: AxisSpec) -> str:
        if axis.gate is not None:
            return "gated"
        if axis.negatives:
            return "contrast"
        return "plain"

    @staticmethod
    def _axis_formula(axis: AxisSpec) -> str:
        if axis.gate is None and not axis.negatives:
            return axis.positive
        inner = axis.positive
        if axis.negatives:
            neg_expr = " + ".join(axis.negatives)
            if len(axis.negatives) > 1:
                neg_expr = f"mean[{neg_expr}]"
            inner = f"{inner} - {axis.negative_weight}*{neg_expr}"
        if axis.gate is not None:
            return f"{axis.gate} * ({inner})"
        return inner


def _check_placeholders(blueprint: EmbedderBlueprint) -> None:
    """Raise if any <<EXAMPLE_*>> / _DELETE_ME placeholder survived
    into a resource / module / gene / axis name or citation."""
    bad: set[str] = set()

    def scan(label: str, value: str) -> None:
        if any(m in value for m in _PLACEHOLDER_MARKERS):
            bad.add(f"{label}: {value!r}")

    for r in blueprint._resources.values():
        scan("resource key", r.key)
        scan("citation", r.citation)
    for m in blueprint._modules.values():
        scan("module name", m.name)
        for g in m.genes:
            scan(f"gene in module {m.name!r}", g)
    for a in blueprint._axes:
        scan("axis name", a.name)

    if bad:
        raise RuntimeError(
            "Placeholder text still present (replace every "
            "<<EXAMPLE_*>> / _DELETE_ME):\n  - "
            + "\n  - ".join(sorted(bad))
        )


# ============================================================
# AGENT EDIT POINT
# Edit only inside build_blueprint() below. Do not add new
# helpers, classes, or methods. Use only the EmbedderBlueprint
# API (see per-method docstrings on the class).
# ============================================================

def build_blueprint(name: str) -> EmbedderBlueprint:
    """Agent: define resources, gene modules, and axes here.

    Use the EmbedderBlueprint API (see per-method docstrings on the
    class for full signatures and examples):
      - add_resource(key, *, kind, citation, doi)            # declare a resource once
      - add_gene_module(name, genes, *, resources)           # cite resources by key
      - add_axis(name, module)
      - add_contrast_axis(name, positive_module, negative_module)
      - add_gated_axis(name, gate_module, positive_module,
                        negative_modules=(), negative_weight=0.5)
      - set_clip(value)  # optional; default 4.0

    The number of axes you register determines the embedding dimension.
    """
    blueprint = EmbedderBlueprint(name)

    blueprint.add_resource(
        "pijuan2019_gastrulation_atlas",
        kind="atlas",
        citation="Pijuan-Sala 2019 Nature mouse gastrulation and early organogenesis atlas",
        doi="10.1038/s41586-019-0933-9",
    )
    blueprint.add_resource(
        "argelaguet2019_gastrulation_multiome",
        kind="atlas",
        citation="Argelaguet 2019 Nature mouse gastrulation single-cell multi-omics atlas",
        doi="10.1038/s41586-019-1825-8",
    )
    blueprint.add_resource(
        "nowotschin2019_gut_endoderm",
        kind="paper",
        citation="Nowotschin 2019 Nature mouse gut endoderm single-cell landscape",
        doi="10.1038/s41586-019-1127-1",
    )
    blueprint.add_resource(
        "scialdone2016_mesoderm",
        kind="paper",
        citation="Scialdone 2016 Nature early mouse mesoderm diversification single-cell profiling",
        doi="10.1038/nature18633",
    )
    blueprint.add_resource(
        "magnusdottir2013_pgc",
        kind="paper",
        citation="Magnusdottir 2013 Nature Cell Biology PGC specification transcription-factor network",
        doi="10.1038/ncb2798",
    )
    blueprint.add_resource(
        "reactome_cell_cycle",
        kind="database",
        citation="Reactome R-MMU-1640170 Cell Cycle",
        doi=None,
    )

    blueprint.add_gene_module(
        "epiblast_pluripotent",
        ["Pou5f1", "Nanog", "Sox2", "Otx2", "Fgf5", "Utf1", "Dnmt3b"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "primitive_streak",
        ["T", "Eomes", "Mixl1", "Wnt3", "Fgf8", "Evx1", "Mesp1", "Tbx6"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome", "scialdone2016_mesoderm"],
    )
    blueprint.add_gene_module(
        "embryonic_ectoderm",
        ["Sox2", "Otx2", "Pou3f1", "Zic2", "Zic3", "Cdh1", "Epcam"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "primordial_germ_cell",
        ["Prdm1", "Prdm14", "Tfap2c", "Dppa3", "Nanos3", "Kit", "Pou5f1"],
        resources=["pijuan2019_gastrulation_atlas", "magnusdottir2013_pgc"],
    )
    blueprint.add_gene_module(
        "neuromesodermal_progenitor",
        ["T", "Sox2", "Cdx2", "Nkx1-2", "Fgf8", "Wnt3a", "Cyp26a1"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "definitive_endoderm",
        ["Sox17", "Foxa2", "Cxcr4", "Gsc", "Cer1", "Lhx1", "Hhex"],
        resources=["pijuan2019_gastrulation_atlas", "nowotschin2019_gut_endoderm"],
    )
    blueprint.add_gene_module(
        "node_notochord",
        ["Noto", "Foxa2", "Shh", "T", "Chrd", "Nog", "Dnah11", "Lefty1"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "nascent_mesoderm",
        ["Tbx6", "Msgn1", "Mesp1", "Mixl1", "Fgf8", "Cdx1", "Cdx2", "Snai1"],
        resources=["pijuan2019_gastrulation_atlas", "scialdone2016_mesoderm"],
    )
    blueprint.add_gene_module(
        "extraembryonic_mesoderm",
        ["Tbx4", "Hand1", "Hoxa10", "Bmp4", "Msx1", "Cdx2", "Pdgfra"],
        resources=["pijuan2019_gastrulation_atlas", "scialdone2016_mesoderm"],
    )
    blueprint.add_gene_module(
        "paraxial_presomitic_mesoderm",
        ["Tbx6", "Msgn1", "Meox1", "Mesp2", "Dll1", "Ripply2", "Dmrt2"],
        resources=["pijuan2019_gastrulation_atlas", "scialdone2016_mesoderm"],
    )
    blueprint.add_gene_module(
        "somitic_mesoderm",
        ["Pax3", "Meox1", "Tcf15", "Uncx", "Dmrt2", "Myf5", "Pax7"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "intermediate_mesoderm",
        ["Osr1", "Pax2", "Lhx1", "Wt1", "Gata3"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "lateral_plate_mesoderm",
        ["Foxf1", "Hand1", "Hand2", "Prrx1", "Bmp4", "Wt1", "Kdr"],
        resources=["pijuan2019_gastrulation_atlas", "scialdone2016_mesoderm"],
    )
    blueprint.add_gene_module(
        "cardiac_mesoderm",
        ["Mesp1", "Nkx2-5", "Isl1", "Tbx5", "Gata4", "Mef2c", "Hand2"],
        resources=["pijuan2019_gastrulation_atlas", "scialdone2016_mesoderm"],
    )
    blueprint.add_gene_module(
        "pharyngeal_mesoderm",
        ["Tbx1", "Pitx2", "Alx3", "Alx4", "Twist1", "Tcf21"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "hematoendothelial_mesoderm",
        ["Kdr", "Tal1", "Lmo2", "Gata2", "Fli1", "Runx1", "Etv2"],
        resources=["pijuan2019_gastrulation_atlas", "scialdone2016_mesoderm"],
    )
    blueprint.add_gene_module(
        "endothelial",
        ["Pecam1", "Cdh5", "Tek", "Kdr", "Sox7", "Sox18", "Esam", "Eng"],
        resources=["pijuan2019_gastrulation_atlas", "scialdone2016_mesoderm"],
    )
    blueprint.add_gene_module(
        "primitive_erythroid",
        ["Gata1", "Klf1", "Hbb-bt", "Hbb-y", "Hba-a1", "Alas2", "Car2", "Epor"],
        resources=["pijuan2019_gastrulation_atlas", "scialdone2016_mesoderm"],
    )
    blueprint.add_gene_module(
        "primitive_myeloid",
        ["Spi1", "Csf1r", "Lyz2", "Itgam", "Adgre1", "Mpeg1", "C1qa"],
        resources=["pijuan2019_gastrulation_atlas", "scialdone2016_mesoderm"],
    )
    blueprint.add_gene_module(
        "neural_plate",
        ["Sox1", "Sox2", "Pax6", "Zic1", "Zic2", "Hesx1", "Otx2"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "anterior_neural",
        ["Six3", "Hesx1", "Foxg1", "Lhx5", "Otx2", "Rax"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "posterior_neural",
        ["Sox2", "Sox1", "Pax6", "Cdx2", "Hoxb1", "Hoxa2", "Hoxb2", "Nkx1-2"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "neural_crest",
        ["Sox10", "Foxd3", "Tfap2a", "Snai2", "Pax3", "Msx1", "Ngfr"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "surface_ectoderm",
        ["Krt8", "Krt18", "Krt19", "Tfap2a", "Grhl2", "Dlx5", "Perp"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "otic_placode",
        ["Pax2", "Pax8", "Eya1", "Six1", "Dlx5", "Sox9"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "visceral_endoderm",
        ["Apoa1", "Apoa2", "Afp", "Ttr", "Amn", "Cubn", "Hnf4a", "Krt8"],
        resources=["pijuan2019_gastrulation_atlas", "nowotschin2019_gut_endoderm"],
    )
    blueprint.add_gene_module(
        "anterior_visceral_endoderm",
        ["Cer1", "Lefty1", "Dkk1", "Hhex", "Hesx1", "Gsc", "Otx2"],
        resources=["pijuan2019_gastrulation_atlas", "nowotschin2019_gut_endoderm"],
    )
    blueprint.add_gene_module(
        "parietal_endoderm",
        ["Pdgfra", "Sparc", "Lama1", "Col4a1", "Col4a2", "Plat", "Dab2"],
        resources=["pijuan2019_gastrulation_atlas", "nowotschin2019_gut_endoderm"],
    )
    blueprint.add_gene_module(
        "yolk_sac_endoderm",
        ["Apoa1", "Apoa2", "Afp", "Ttr", "Apob", "Hnf4a", "Rbp4"],
        resources=["pijuan2019_gastrulation_atlas", "nowotschin2019_gut_endoderm"],
    )
    blueprint.add_gene_module(
        "gut_endoderm",
        ["Sox17", "Foxa2", "Hnf1b", "Onecut1", "Hhex", "Cldn6", "Cldn7"],
        resources=["pijuan2019_gastrulation_atlas", "nowotschin2019_gut_endoderm"],
    )
    blueprint.add_gene_module(
        "anterior_foregut",
        ["Hhex", "Foxa2", "Sox2", "Onecut1", "Irx3", "Shh"],
        resources=["pijuan2019_gastrulation_atlas", "nowotschin2019_gut_endoderm"],
    )
    blueprint.add_gene_module(
        "posterior_hindgut",
        ["Cdx2", "Cdx1", "Hoxb9", "Hoxa9", "Foxa2", "Gata4"],
        resources=["pijuan2019_gastrulation_atlas", "nowotschin2019_gut_endoderm"],
    )
    blueprint.add_gene_module(
        "extraembryonic_ectoderm",
        ["Cdx2", "Elf5", "Eomes", "Tfap2c", "Tead4", "Hand1", "Krt8"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "chorion_trophoblast",
        ["Gcm1", "Esx1", "Cdx2", "Elf5", "Eomes", "Tfap2c"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "trophoblast_giant_cell",
        ["Hand1", "Prl3d1", "Prl3b1", "Tpbpa", "Ctsq", "Ascl2"],
        resources=["pijuan2019_gastrulation_atlas", "argelaguet2019_gastrulation_multiome"],
    )
    blueprint.add_gene_module(
        "proliferation",
        ["Mki67", "Top2a", "Cdk1", "Ccnb1", "Ccna2", "Pcna", "Mcm5", "Mcm6", "Ube2c", "Cenpf"],
        resources=["pijuan2019_gastrulation_atlas", "reactome_cell_cycle"],
    )

    blueprint.add_axis("epi_pluripotentRoot", module="epiblast_pluripotent")
    blueprint.add_axis("ve_visceralEndodermRoot", module="visceral_endoderm")
    blueprint.add_axis("exe_extraembryonicEctodermRoot", module="extraembryonic_ectoderm")
    blueprint.add_axis("cycle_proliferationRoot", module="proliferation")

    blueprint.add_gated_axis(
        "streak_primitiveStreak",
        gate_module="epiblast_pluripotent",
        positive_module="primitive_streak",
        negative_modules=["embryonic_ectoderm", "primordial_germ_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "ect_embryonicEctoderm",
        gate_module="epiblast_pluripotent",
        positive_module="embryonic_ectoderm",
        negative_modules=["primitive_streak", "primordial_germ_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "pgc_primordialGermCell",
        gate_module="epiblast_pluripotent",
        positive_module="primordial_germ_cell",
        negative_modules=["primitive_streak", "embryonic_ectoderm"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "nmp_neuromesodermalProgenitor",
        gate_module="epiblast_pluripotent",
        positive_module="neuromesodermal_progenitor",
        negative_modules=["primitive_streak", "embryonic_ectoderm"],
        negative_weight=0.5,
    )

    blueprint.add_gated_axis(
        "mes_nascentMesoderm",
        gate_module="primitive_streak",
        positive_module="nascent_mesoderm",
        negative_modules=["definitive_endoderm", "node_notochord", "extraembryonic_mesoderm"],
        negative_weight=0.35,
    )
    blueprint.add_gated_axis(
        "end_definitiveEndoderm",
        gate_module="primitive_streak",
        positive_module="definitive_endoderm",
        negative_modules=["nascent_mesoderm", "node_notochord", "extraembryonic_mesoderm"],
        negative_weight=0.35,
    )
    blueprint.add_gated_axis(
        "axial_nodeNotochord",
        gate_module="primitive_streak",
        positive_module="node_notochord",
        negative_modules=["nascent_mesoderm", "definitive_endoderm", "extraembryonic_mesoderm"],
        negative_weight=0.35,
    )
    blueprint.add_gated_axis(
        "extraMes_extraembryonicMesoderm",
        gate_module="primitive_streak",
        positive_module="extraembryonic_mesoderm",
        negative_modules=["nascent_mesoderm", "definitive_endoderm", "node_notochord"],
        negative_weight=0.35,
    )

    blueprint.add_gated_axis(
        "mes_paraxialPresomitic",
        gate_module="nascent_mesoderm",
        positive_module="paraxial_presomitic_mesoderm",
        negative_modules=["lateral_plate_mesoderm", "intermediate_mesoderm", "cardiac_mesoderm", "pharyngeal_mesoderm", "hematoendothelial_mesoderm"],
        negative_weight=0.3,
    )
    blueprint.add_gated_axis(
        "mes_lateralPlate",
        gate_module="nascent_mesoderm",
        positive_module="lateral_plate_mesoderm",
        negative_modules=["paraxial_presomitic_mesoderm", "intermediate_mesoderm", "cardiac_mesoderm", "pharyngeal_mesoderm", "hematoendothelial_mesoderm"],
        negative_weight=0.3,
    )
    blueprint.add_gated_axis(
        "mes_intermediateMesoderm",
        gate_module="nascent_mesoderm",
        positive_module="intermediate_mesoderm",
        negative_modules=["paraxial_presomitic_mesoderm", "lateral_plate_mesoderm", "cardiac_mesoderm", "pharyngeal_mesoderm", "hematoendothelial_mesoderm"],
        negative_weight=0.3,
    )
    blueprint.add_gated_axis(
        "mes_cardiacMesoderm",
        gate_module="nascent_mesoderm",
        positive_module="cardiac_mesoderm",
        negative_modules=["paraxial_presomitic_mesoderm", "lateral_plate_mesoderm", "intermediate_mesoderm", "pharyngeal_mesoderm", "hematoendothelial_mesoderm"],
        negative_weight=0.3,
    )
    blueprint.add_gated_axis(
        "mes_pharyngealMesoderm",
        gate_module="nascent_mesoderm",
        positive_module="pharyngeal_mesoderm",
        negative_modules=["paraxial_presomitic_mesoderm", "lateral_plate_mesoderm", "intermediate_mesoderm", "cardiac_mesoderm", "hematoendothelial_mesoderm"],
        negative_weight=0.3,
    )
    blueprint.add_gated_axis(
        "mes_hematoendothelial",
        gate_module="nascent_mesoderm",
        positive_module="hematoendothelial_mesoderm",
        negative_modules=["paraxial_presomitic_mesoderm", "lateral_plate_mesoderm", "intermediate_mesoderm", "cardiac_mesoderm", "pharyngeal_mesoderm"],
        negative_weight=0.3,
    )
    blueprint.add_gated_axis(
        "mes_somiticMesoderm",
        gate_module="paraxial_presomitic_mesoderm",
        positive_module="somitic_mesoderm",
        negative_modules=[],
        negative_weight=0.5,
    )

    blueprint.add_gated_axis(
        "hemo_endothelial",
        gate_module="hematoendothelial_mesoderm",
        positive_module="endothelial",
        negative_modules=["primitive_erythroid", "primitive_myeloid"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "hemo_primitiveErythroid",
        gate_module="hematoendothelial_mesoderm",
        positive_module="primitive_erythroid",
        negative_modules=["endothelial", "primitive_myeloid"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "hemo_primitiveMyeloid",
        gate_module="hematoendothelial_mesoderm",
        positive_module="primitive_myeloid",
        negative_modules=["endothelial", "primitive_erythroid"],
        negative_weight=0.5,
    )

    blueprint.add_gated_axis(
        "ect_neuralPlate",
        gate_module="embryonic_ectoderm",
        positive_module="neural_plate",
        negative_modules=["surface_ectoderm", "otic_placode"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "ect_surfaceEctoderm",
        gate_module="embryonic_ectoderm",
        positive_module="surface_ectoderm",
        negative_modules=["neural_plate", "otic_placode"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "ect_oticPlacode",
        gate_module="embryonic_ectoderm",
        positive_module="otic_placode",
        negative_modules=["neural_plate", "surface_ectoderm"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "neural_anteriorNeural",
        gate_module="neural_plate",
        positive_module="anterior_neural",
        negative_modules=["posterior_neural", "neural_crest"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "neural_posteriorNeural",
        gate_module="neural_plate",
        positive_module="posterior_neural",
        negative_modules=["anterior_neural", "neural_crest"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "neural_neuralCrest",
        gate_module="neural_plate",
        positive_module="neural_crest",
        negative_modules=["anterior_neural", "posterior_neural"],
        negative_weight=0.5,
    )

    blueprint.add_gated_axis(
        "end_gutEndoderm",
        gate_module="definitive_endoderm",
        positive_module="gut_endoderm",
        negative_modules=["visceral_endoderm"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "gut_anteriorForegut",
        gate_module="gut_endoderm",
        positive_module="anterior_foregut",
        negative_modules=["posterior_hindgut"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "gut_posteriorHindgut",
        gate_module="gut_endoderm",
        positive_module="posterior_hindgut",
        negative_modules=["anterior_foregut"],
        negative_weight=0.5,
    )

    blueprint.add_gated_axis(
        "ve_anteriorVisceralEndoderm",
        gate_module="visceral_endoderm",
        positive_module="anterior_visceral_endoderm",
        negative_modules=["parietal_endoderm", "yolk_sac_endoderm"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "ve_parietalEndoderm",
        gate_module="visceral_endoderm",
        positive_module="parietal_endoderm",
        negative_modules=["anterior_visceral_endoderm", "yolk_sac_endoderm"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "ve_yolkSacEndoderm",
        gate_module="visceral_endoderm",
        positive_module="yolk_sac_endoderm",
        negative_modules=["anterior_visceral_endoderm", "parietal_endoderm"],
        negative_weight=0.5,
    )

    blueprint.add_gated_axis(
        "exe_chorionTrophoblast",
        gate_module="extraembryonic_ectoderm",
        positive_module="chorion_trophoblast",
        negative_modules=["trophoblast_giant_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "exe_trophoblastGiantCell",
        gate_module="extraembryonic_ectoderm",
        positive_module="trophoblast_giant_cell",
        negative_modules=["chorion_trophoblast"],
        negative_weight=0.5,
    )

    blueprint.set_clip(4.0)

    return blueprint


def _clip_value(s: str):
    """argparse type for --clip: a positive float, or 'none' to disable.

    Returns the string sentinel "none" for disable (distinct from the
    flag's default of None, which means "not passed; keep the blueprint's
    value").
    """
    if s.strip().lower() in ("none", "off", "disable"):
        return "none"
    try:
        v = float(s)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"--clip must be a number or 'none', got {s!r}"
        )
    if v <= 0:
        raise argparse.ArgumentTypeError(
            "--clip must be positive (or 'none' to disable)"
        )
    return v


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a count-matrix cell embedding from a declarative blueprint."
        )
    )
    parser.add_argument(
        "--input",
        type=Path,
        default=Path("data/processed/tasic2016_qc.h5ad"),
        help="Path to input AnnData .h5ad with a cells x genes matrix.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("embedding.parquet"),
        help="Path to output parquet with cell_id + dim_0..dim_{N-1}.",
    )
    parser.add_argument(
        "--layer",
        type=str,
        default=None,
        help=(
            "Optional AnnData layer to read expression from "
            "(e.g. 'counts'). If omitted, adata.X is used."
        ),
    )
    parser.add_argument(
        "--clip",
        type=_clip_value,
        default=None,
        help=(
            "Override the z-score clip threshold. If omitted, the blueprint's "
            "value is used (default 4.0). Pass a positive number, or 'none' to "
            "disable clipping."
        ),
    )
    parser.add_argument(
        "--n-min-genes",
        type=int,
        default=0,
        dest="n_min_genes",
        help=(
            "Exclude any axis whose positive module has fewer than this many "
            "genes present in the dataset's var_names (default 0 = keep all axes)."
        ),
    )
    parser.add_argument(
        "--describe",
        action="store_true",
        help="Print the blueprint's axis list and exit without embedding.",
    )
    args = parser.parse_args()
    if args.n_min_genes < 0:
        parser.error("--n-min-genes must be >= 0")
    return args


def main() -> None:
    # Imported here (not at module level) so importing this module to call
    # build_blueprint() stays side-effect-free; only the CLI entrypoint needs these.
    import platform
    import resource
    import sys
    import time

    def _rss_mb() -> float:
        r = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
             + resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
        return r / 1024.0**2 if sys.platform == "darwin" else r / 1024.0

    args = parse_args()
    _t0 = time.perf_counter()
    blueprint = build_blueprint("default")          # ~ms — instantiate the declarative model
    _t1 = time.perf_counter()
    _rss_after_obtain = _rss_mb()
    _check_placeholders(blueprint)                  # fail-fast on un-replaced <<EXAMPLE_*>> sentinels
    if args.clip is not None:                        # CLI override wins over set_clip
        blueprint.set_clip(None if args.clip == "none" else args.clip)
    if args.describe:
        print(blueprint.describe())
        return
    blueprint.build_embedding(
        args.input, args.output, layer=args.layer, n_min_genes=args.n_min_genes
    )
    _t2 = time.perf_counter()
    _peak = _rss_mb()
    _walltime = _t2 - _t0
    _obtain = _t1 - _t0                                              # build_blueprint()
    _load = float(getattr(blueprint, "_last_load_data_s", 0.0) or 0.0)  # read inside build_embedding
    _build = max(0.0, (_t2 - _t1) - _load)                          # scoring + parquet write
    _other = _walltime - _load - _obtain - _build
    _rss_load = float(getattr(blueprint, "_last_rss_after_load_mb", 0.0) or _peak)
    _emb_bytes = args.output.stat().st_size if args.output.exists() else 0
    _host = platform.node().replace(" ", "_") or "?"
    print(f"Wrote {args.output}")
    # Time taxonomy (sums to walltime_s): load_data (read the input) + obtain_model
    # (build_blueprint) + build_embedding (score every cell + write parquet) + other.
    # model_source=authored: the "model" is this ~tens-of-KB script, written by the LLM.
    print(f"[timing] walltime_s={_walltime:.4f} load_data_s={_load:.4f} "
          f"obtain_model_s={_obtain:.4f} build_embedding_s={_build:.4f} other_s={_other:.4f} "
          f"model_source=authored model_download_s=0.0 "
          f"peak_rss_mb={_peak:.1f} rss_after_load_mb={_rss_load:.1f} "
          f"rss_after_obtain_mb={_rss_after_obtain:.1f} rss_after_build_mb={_peak:.1f} "
          f"model_storage_bytes={Path(__file__).stat().st_size} "
          f"embedding_bytes={_emb_bytes} device=cpu host={_host}")


if __name__ == "__main__":
    main()
