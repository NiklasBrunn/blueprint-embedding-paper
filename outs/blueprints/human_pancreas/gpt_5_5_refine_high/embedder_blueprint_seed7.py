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
        "segerstolpe_2016_islet_sct",
        kind="paper",
        citation="Segerstolpe 2016 Cell Metabolism",
        doi="10.1016/j.cmet.2016.08.020",
    )
    blueprint.add_resource(
        "muraro_2016_pancreas_atlas",
        kind="atlas",
        citation="Muraro 2016 Cell Systems",
        doi="10.1016/j.cels.2016.09.002",
    )
    blueprint.add_resource(
        "baron_2016_pancreas_map",
        kind="atlas",
        citation="Baron 2016 Cell Systems",
        doi="10.1016/j.cels.2016.08.011",
    )
    blueprint.add_resource(
        "tabula_sapiens_2022",
        kind="atlas",
        citation="Tabula Sapiens Consortium 2022 Science",
        doi="10.1126/science.abl4896",
    )
    blueprint.add_resource(
        "craig_schapiro_2025_pancreas_endothelium",
        kind="paper",
        citation="Craig-Schapiro 2025 Nature Communications",
        doi="10.1038/s41467-024-55415-3",
    )
    blueprint.add_resource(
        "kegg_pancreatic_secretion",
        kind="database",
        citation="KEGG hsa04972 Pancreatic secretion",
        doi=None,
    )
    blueprint.add_resource(
        "kegg_insulin_secretion",
        kind="database",
        citation="KEGG hsa04911 Insulin secretion",
        doi=None,
    )
    blueprint.add_resource(
        "msigdb_emt",
        kind="database",
        citation="MSigDB HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
        doi=None,
    )
    blueprint.add_resource(
        "go_angiogenesis",
        kind="database",
        citation="GO:0001525 angiogenesis",
        doi=None,
    )
    blueprint.add_resource(
        "go_immune_response",
        kind="database",
        citation="GO:0006955 immune response",
        doi=None,
    )
    blueprint.add_resource(
        "go_t_cell_activation",
        kind="database",
        citation="GO:0042110 T cell activation",
        doi=None,
    )
    blueprint.add_resource(
        "go_b_cell_activation",
        kind="database",
        citation="GO:0042113 B cell activation",
        doi=None,
    )
    blueprint.add_resource(
        "go_nk_cell_cytotoxicity",
        kind="database",
        citation="GO:0042267 natural killer cell mediated cytotoxicity",
        doi=None,
    )
    blueprint.add_resource(
        "go_antigen_processing_presentation",
        kind="database",
        citation="GO:0019882 antigen processing and presentation",
        doi=None,
    )
    blueprint.add_resource(
        "go_humoral_immune_response",
        kind="database",
        citation="GO:0006959 humoral immune response",
        doi=None,
    )
    blueprint.add_resource(
        "go_extracellular_matrix_organization",
        kind="database",
        citation="GO:0030198 extracellular matrix organization",
        doi=None,
    )

    blueprint.add_gene_module(
        "endocrine_core",
        ["CHGA", "CHGB", "SCG2", "PCSK1", "PCSK2", "ISL1", "PAX6", "NEUROD1"],
        resources=[
            "segerstolpe_2016_islet_sct",
            "muraro_2016_pancreas_atlas",
            "tabula_sapiens_2022",
        ],
    )
    blueprint.add_gene_module(
        "beta_cell",
        ["INS", "IAPP", "MAFA", "PDX1", "NKX6-1", "SLC30A8", "G6PC2", "UCN3", "PCSK1"],
        resources=[
            "segerstolpe_2016_islet_sct",
            "muraro_2016_pancreas_atlas",
            "craig_schapiro_2025_pancreas_endothelium",
            "kegg_insulin_secretion",
        ],
    )
    blueprint.add_gene_module(
        "alpha_cell",
        ["GCG", "GC", "TTR", "ARX", "IRX2", "MAFB", "TM4SF4"],
        resources=[
            "segerstolpe_2016_islet_sct",
            "muraro_2016_pancreas_atlas",
            "craig_schapiro_2025_pancreas_endothelium",
        ],
    )
    blueprint.add_gene_module(
        "delta_cell",
        ["SST", "LEPR", "PRG4", "HHEX", "RBP4", "PCSK1N"],
        resources=[
            "segerstolpe_2016_islet_sct",
            "muraro_2016_pancreas_atlas",
            "craig_schapiro_2025_pancreas_endothelium",
        ],
    )
    blueprint.add_gene_module(
        "gamma_pp_cell",
        ["PPY", "PYY", "SLC6A4", "FEV"],
        resources=[
            "segerstolpe_2016_islet_sct",
            "muraro_2016_pancreas_atlas",
            "craig_schapiro_2025_pancreas_endothelium",
        ],
    )
    blueprint.add_gene_module(
        "epsilon_cell",
        ["GHRL"],
        resources=[
            "segerstolpe_2016_islet_sct",
            "baron_2016_pancreas_map",
            "craig_schapiro_2025_pancreas_endothelium",
        ],
    )
    blueprint.add_gene_module(
        "acinar_secretory",
        ["PRSS1", "PRSS2", "CPA1", "CPA2", "CTRB1", "CTRB2", "CELA3A", "CPB1", "PNLIP", "CEL", "AMY2A"],
        resources=[
            "muraro_2016_pancreas_atlas",
            "baron_2016_pancreas_map",
            "craig_schapiro_2025_pancreas_endothelium",
            "kegg_pancreatic_secretion",
        ],
    )
    blueprint.add_gene_module(
        "regenerating_acinar",
        ["REG1A", "REG1B", "REG3A", "REG3G", "SPINK1"],
        resources=[
            "muraro_2016_pancreas_atlas",
            "baron_2016_pancreas_map",
            "kegg_pancreatic_secretion",
        ],
    )
    blueprint.add_gene_module(
        "ductal_epithelial",
        ["KRT19", "KRT7", "SOX9", "CFTR", "SPP1", "MMP7", "SFRP5", "EPCAM", "CLDN10", "ANXA4", "SLC4A4"],
        resources=[
            "muraro_2016_pancreas_atlas",
            "baron_2016_pancreas_map",
            "craig_schapiro_2025_pancreas_endothelium",
            "kegg_pancreatic_secretion",
        ],
    )
    blueprint.add_gene_module(
        "endothelial_core",
        ["PECAM1", "CDH5", "KDR", "VWF", "ESAM", "ENG", "PLVAP", "RAMP2", "ERG", "FLI1"],
        resources=[
            "baron_2016_pancreas_map",
            "tabula_sapiens_2022",
            "craig_schapiro_2025_pancreas_endothelium",
            "go_angiogenesis",
        ],
    )
    blueprint.add_gene_module(
        "islet_endothelium",
        ["ACE", "PASK", "F2RL3", "ESM1", "CXCR4", "ACKR3", "UNC5B", "LAMA4", "CREM", "COL13A1", "ANGPTL2", "THBS1"],
        resources=[
            "craig_schapiro_2025_pancreas_endothelium",
            "go_angiogenesis",
        ],
    )
    blueprint.add_gene_module(
        "acinar_endothelium",
        ["AQP1", "CCL14", "GPIHBP1", "ATOH8", "CA4", "RGCC", "FABP4", "EMCN", "CLEC14A", "TM4SF1"],
        resources=[
            "craig_schapiro_2025_pancreas_endothelium",
            "go_angiogenesis",
        ],
    )
    blueprint.add_gene_module(
        "arterial_endothelium",
        ["SERPINE2", "SEMA3G", "FBLN2", "FBLN5", "HEY1", "GJA5"],
        resources=[
            "craig_schapiro_2025_pancreas_endothelium",
            "go_angiogenesis",
        ],
    )
    blueprint.add_gene_module(
        "venous_endothelium",
        ["ACKR1", "MADCAM1", "CPE", "FAM155A", "NR2F2"],
        resources=[
            "craig_schapiro_2025_pancreas_endothelium",
            "go_angiogenesis",
        ],
    )
    blueprint.add_gene_module(
        "lymphatic_endothelium",
        ["PROX1", "PDPN", "LYVE1", "FLT4", "CCL21", "MMRN1"],
        resources=[
            "craig_schapiro_2025_pancreas_endothelium",
            "go_angiogenesis",
        ],
    )
    blueprint.add_gene_module(
        "stromal_matrix",
        ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "MMP2", "BGN", "SPARC", "PDGFRA"],
        resources=[
            "muraro_2016_pancreas_atlas",
            "baron_2016_pancreas_map",
            "craig_schapiro_2025_pancreas_endothelium",
            "msigdb_emt",
        ],
    )
    blueprint.add_gene_module(
        "pericyte_smooth_muscle",
        ["PDGFRB", "RGS5", "MCAM", "CSPG4", "NOTCH3", "ABCC9", "KCNJ8", "MYH11"],
        resources=[
            "baron_2016_pancreas_map",
            "craig_schapiro_2025_pancreas_endothelium",
            "msigdb_emt",
            "go_angiogenesis",
        ],
    )
    blueprint.add_gene_module(
        "quiescent_stellate",
        ["RBP1", "CYGB", "DES", "GFAP", "PDGFRB", "RGS5"],
        resources=[
            "baron_2016_pancreas_map",
            "craig_schapiro_2025_pancreas_endothelium",
            "msigdb_emt",
        ],
    )
    blueprint.add_gene_module(
        "activated_stellate_myofibroblast",
        ["ACTA2", "TAGLN", "COL1A1", "COL3A1", "POSTN", "THY1", "MMP2", "COL11A1"],
        resources=[
            "baron_2016_pancreas_map",
            "craig_schapiro_2025_pancreas_endothelium",
            "msigdb_emt",
            "go_extracellular_matrix_organization",
        ],
    )
    blueprint.add_gene_module(
        "immune_leukocyte",
        ["PTPRC", "LST1", "TYROBP", "CD74", "HLA-DRA", "HLA-DPA1", "CXCL8"],
        resources=[
            "baron_2016_pancreas_map",
            "tabula_sapiens_2022",
            "craig_schapiro_2025_pancreas_endothelium",
            "go_immune_response",
        ],
    )
    blueprint.add_gene_module(
        "macrophage_myeloid",
        ["CD68", "LYZ", "LST1", "C1QA", "C1QB", "TYROBP", "FCGR3A", "AIF1"],
        resources=[
            "baron_2016_pancreas_map",
            "craig_schapiro_2025_pancreas_endothelium",
            "go_immune_response",
        ],
    )
    blueprint.add_gene_module(
        "t_cell",
        ["CD3D", "CD3E", "CD2", "TRAC", "CD8A", "IL7R"],
        resources=[
            "baron_2016_pancreas_map",
            "craig_schapiro_2025_pancreas_endothelium",
            "go_t_cell_activation",
        ],
    )
    blueprint.add_gene_module(
        "cytotoxic_nk_cell",
        ["NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "FCGR3A", "CTSW"],
        resources=[
            "baron_2016_pancreas_map",
            "tabula_sapiens_2022",
            "craig_schapiro_2025_pancreas_endothelium",
            "go_t_cell_activation",
            "go_nk_cell_cytotoxicity",
        ],
    )
    blueprint.add_gene_module(
        "dendritic_apc",
        ["FCER1A", "CD1C", "CLEC10A", "CLEC9A", "LILRA4", "IRF8", "HLA-DRA", "HLA-DPA1"],
        resources=[
            "baron_2016_pancreas_map",
            "tabula_sapiens_2022",
            "craig_schapiro_2025_pancreas_endothelium",
            "go_immune_response",
            "go_antigen_processing_presentation",
        ],
    )
    blueprint.add_gene_module(
        "b_cell",
        ["CD37", "MS4A1", "CD79A", "CD79B", "BANK1", "CD74"],
        resources=[
            "baron_2016_pancreas_map",
            "craig_schapiro_2025_pancreas_endothelium",
            "go_b_cell_activation",
        ],
    )
    blueprint.add_gene_module(
        "plasma_cell",
        ["MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1", "IGKC", "IGLC2"],
        resources=[
            "tabula_sapiens_2022",
            "craig_schapiro_2025_pancreas_endothelium",
            "go_humoral_immune_response",
        ],
    )
    blueprint.add_gene_module(
        "mast_cell",
        ["KIT", "TPSAB1", "TPSB2", "CPA3", "HDC", "MS4A2"],
        resources=[
            "baron_2016_pancreas_map",
            "craig_schapiro_2025_pancreas_endothelium",
            "go_immune_response",
        ],
    )
    blueprint.add_gene_module(
        "neural_ganglial_schwann",
        ["KIRREL3", "VAT1L", "NRG3", "SOX10", "PLP1", "MPZ", "S100B"],
        resources=[
            "baron_2016_pancreas_map",
            "craig_schapiro_2025_pancreas_endothelium",
        ],
    )

    blueprint.add_axis("endocrine_islet_compartment", module="endocrine_core")
    blueprint.add_axis("acinar_exocrine_compartment", module="acinar_secretory")
    blueprint.add_axis("ductal_epithelial_compartment", module="ductal_epithelial")
    blueprint.add_axis("stromal_matrix_compartment", module="stromal_matrix")
    blueprint.add_axis("vascular_endothelial_compartment", module="endothelial_core")
    blueprint.add_axis("resident_immune_compartment", module="immune_leukocyte")
    blueprint.add_axis("neural_ganglial_schwann_compartment", module="neural_ganglial_schwann")
    blueprint.add_gated_axis(
        "endocrine_beta_identity",
        gate_module="endocrine_core",
        positive_module="beta_cell",
        negative_modules=["alpha_cell", "delta_cell", "gamma_pp_cell", "epsilon_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "endocrine_alpha_identity",
        gate_module="endocrine_core",
        positive_module="alpha_cell",
        negative_modules=["beta_cell", "delta_cell", "gamma_pp_cell", "epsilon_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "endocrine_delta_identity",
        gate_module="endocrine_core",
        positive_module="delta_cell",
        negative_modules=["beta_cell", "alpha_cell", "gamma_pp_cell", "epsilon_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "endocrine_gamma_pp_identity",
        gate_module="endocrine_core",
        positive_module="gamma_pp_cell",
        negative_modules=["beta_cell", "alpha_cell", "delta_cell", "epsilon_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "endocrine_epsilon_identity",
        gate_module="endocrine_core",
        positive_module="epsilon_cell",
        negative_modules=["beta_cell", "alpha_cell", "delta_cell", "gamma_pp_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "acinar_regenerating_substate",
        gate_module="acinar_secretory",
        positive_module="regenerating_acinar",
        negative_modules=[],
        negative_weight=0.4,
    )
    blueprint.add_gated_axis(
        "islet_endothelial_specialization",
        gate_module="endothelial_core",
        positive_module="islet_endothelium",
        negative_modules=["acinar_endothelium", "arterial_endothelium", "venous_endothelium", "lymphatic_endothelium"],
        negative_weight=0.4,
    )
    blueprint.add_gated_axis(
        "acinar_endothelial_specialization",
        gate_module="endothelial_core",
        positive_module="acinar_endothelium",
        negative_modules=["islet_endothelium", "arterial_endothelium", "venous_endothelium", "lymphatic_endothelium"],
        negative_weight=0.4,
    )
    blueprint.add_gated_axis(
        "arterial_endothelial_identity",
        gate_module="endothelial_core",
        positive_module="arterial_endothelium",
        negative_modules=["venous_endothelium", "lymphatic_endothelium"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "venous_endothelial_identity",
        gate_module="endothelial_core",
        positive_module="venous_endothelium",
        negative_modules=["arterial_endothelium", "lymphatic_endothelium"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "lymphatic_endothelial_identity",
        gate_module="endothelial_core",
        positive_module="lymphatic_endothelium",
        negative_modules=["arterial_endothelium", "venous_endothelium"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "stromal_pericyte_smooth_muscle_identity",
        gate_module="stromal_matrix",
        positive_module="pericyte_smooth_muscle",
        negative_modules=["quiescent_stellate", "activated_stellate_myofibroblast"],
        negative_weight=0.4,
    )
    blueprint.add_gated_axis(
        "stromal_quiescent_stellate_identity",
        gate_module="stromal_matrix",
        positive_module="quiescent_stellate",
        negative_modules=["pericyte_smooth_muscle", "activated_stellate_myofibroblast"],
        negative_weight=0.4,
    )
    blueprint.add_gated_axis(
        "stromal_activated_stellate_myofibroblast_state",
        gate_module="stromal_matrix",
        positive_module="activated_stellate_myofibroblast",
        negative_modules=["pericyte_smooth_muscle", "quiescent_stellate"],
        negative_weight=0.4,
    )
    blueprint.add_gated_axis(
        "immune_macrophage_identity",
        gate_module="immune_leukocyte",
        positive_module="macrophage_myeloid",
        negative_modules=["t_cell", "cytotoxic_nk_cell", "dendritic_apc", "b_cell", "plasma_cell", "mast_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "immune_t_cell_identity",
        gate_module="immune_leukocyte",
        positive_module="t_cell",
        negative_modules=["macrophage_myeloid", "cytotoxic_nk_cell", "dendritic_apc", "b_cell", "plasma_cell", "mast_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "immune_cytotoxic_nk_identity",
        gate_module="immune_leukocyte",
        positive_module="cytotoxic_nk_cell",
        negative_modules=["macrophage_myeloid", "dendritic_apc", "t_cell", "b_cell", "plasma_cell", "mast_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "immune_dendritic_apc_identity",
        gate_module="immune_leukocyte",
        positive_module="dendritic_apc",
        negative_modules=["macrophage_myeloid", "cytotoxic_nk_cell", "t_cell", "b_cell", "plasma_cell", "mast_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "immune_b_cell_identity",
        gate_module="immune_leukocyte",
        positive_module="b_cell",
        negative_modules=["macrophage_myeloid", "dendritic_apc", "t_cell", "cytotoxic_nk_cell", "plasma_cell", "mast_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "immune_plasma_cell_identity",
        gate_module="immune_leukocyte",
        positive_module="plasma_cell",
        negative_modules=["macrophage_myeloid", "dendritic_apc", "t_cell", "cytotoxic_nk_cell", "b_cell", "mast_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "immune_mast_cell_identity",
        gate_module="immune_leukocyte",
        positive_module="mast_cell",
        negative_modules=["macrophage_myeloid", "dendritic_apc", "t_cell", "cytotoxic_nk_cell", "b_cell", "plasma_cell"],
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
