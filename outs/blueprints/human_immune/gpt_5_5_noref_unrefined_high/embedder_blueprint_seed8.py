#!/usr/bin/env python3
"""Declarative count-matrix embedder blueprint.

This module exposes EmbedderBlueprint, a small declarative builder for
mean-expression cell-type embeddings, and a CLI that reads an AnnData
count matrix and writes a parquet of axis values per cell. Gene modules
are defined by their gene lists alone; there is no literature-resource
or citation mechanism.

AGENTS: EDIT ONLY `build_blueprint` at the bottom of this file. The
EmbedderBlueprint API is documented in per-method docstrings on the
class. Verification:

    python3 embedder_blueprint_noref.py --describe

This must exit cleanly before the task is done.
"""
from __future__ import annotations

import argparse
import difflib

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


# Sentinels in the skeleton's commented `<<EXAMPLE_*>>` / `_DELETE_ME`
# scaffolding. Any of these in a module / gene / axis name means the
# agent left placeholder text in place.
_PLACEHOLDER_MARKERS = ("<<", ">>", "EXAMPLE_", "_DELETE_ME", "PLACEHOLDER")


class EmbedderBlueprint:
    """Declarative builder for a count-matrix cell embedding."""

    def __init__(self, name: str) -> None:
        self._name = name
        self._modules: dict[str, ModuleSpec] = {}
        self._axes: list[AxisSpec] = []
        self._clip: float | None = 4.0

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
    def resources(self) -> dict:
        # No literature-resource mechanism in this variant; always empty.
        # Retained so blueprint introspection (extract_run / regression_guard)
        # sees a uniform, resource-free shape across all runs.
        return {}

    def add_gene_module(
        self,
        name: str,
        genes: Iterable[str],
    ) -> None:
        """Register a named set of genes (a "module").

        Modules are the atoms that axes combine. Use any name you like;
        the name is what `add_axis` and friends refer to.

        Example:
            blueprint.add_gene_module(
                "module1",
                ["gene1", "gene2", "gene3"],
            )
        """
        if name in self._modules:
            raise ValueError(f"Module already registered: {name}")
        self._modules[name] = ModuleSpec(name=name, genes=tuple(genes))

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
                f"    - {module.name} ({len(module.genes)} genes)"
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
    into a module / gene / axis name."""
    bad: set[str] = set()

    def scan(label: str, value: str) -> None:
        if any(m in value for m in _PLACEHOLDER_MARKERS):
            bad.add(f"{label}: {value!r}")

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
    """Agent: define gene modules and axes here.

    Use the EmbedderBlueprint API (see per-method docstrings on the
    class for full signatures and examples):
      - add_gene_module(name, genes)
      - add_axis(name, module)
      - add_contrast_axis(name, positive_module, negative_module)
      - add_gated_axis(name, gate_module, positive_module,
                        negative_modules=(), negative_weight=0.5)
      - set_clip(value)  # optional; default 4.0

    The number of axes you register determines the embedding dimension.
    """
    blueprint = EmbedderBlueprint(name)

    blueprint.add_gene_module(
        "leukocyte_core",
        ["PTPRC", "LCP1", "CORO1A", "ARHGAP15", "IKZF1", "TYROBP", "SPI1"],
    )
    blueprint.add_gene_module(
        "hspc_primitive",
        ["CD34", "SPINK2", "AVP", "CRHBP", "HLF", "MEIS1", "MLLT3", "PROCR", "GATA2", "MPL"],
    )
    blueprint.add_gene_module(
        "lymphoid_progenitor",
        ["IL7R", "DNTT", "RAG1", "RAG2", "VPREB1", "IGLL1", "MME", "EBF1", "PAX5"],
    )
    blueprint.add_gene_module(
        "cell_cycle",
        ["MKI67", "TOP2A", "UBE2C", "HMGB2", "TYMS", "PCNA", "STMN1", "TUBB", "CENPF"],
    )
    blueprint.add_gene_module(
        "erythroid_early",
        ["GATA1", "KLF1", "TAL1", "ALAS2", "GYPA", "SLC4A1", "EPOR", "HBA1", "HBA2"],
    )
    blueprint.add_gene_module(
        "erythroid_late",
        ["HBB", "HBA1", "HBA2", "AHSP", "CA1", "BLVRB", "KLF1", "GYPA", "SLC4A1"],
    )
    blueprint.add_gene_module(
        "megakaryocyte_platelet",
        ["PPBP", "PF4", "ITGA2B", "GP1BA", "GP9", "TUBB1", "NFE2", "VWF", "MPL", "FERMT3"],
    )
    blueprint.add_gene_module(
        "myeloid_core",
        ["SPI1", "LYZ", "CTSS", "TYROBP", "FCER1G", "CST3", "AIF1", "MNDA", "LST1"],
    )
    blueprint.add_gene_module(
        "granulocyte_immature",
        ["MPO", "ELANE", "AZU1", "PRTN3", "CTSG", "CEBPE", "DEFA4", "BPI", "LTF"],
    )
    blueprint.add_gene_module(
        "neutrophil_mature",
        ["S100A8", "S100A9", "FCGR3B", "CSF3R", "CXCR2", "MNDA", "VCAN", "FPR1", "NAMPT"],
    )
    blueprint.add_gene_module(
        "eosinophil",
        ["CLC", "PRG2", "EPX", "RNASE2", "RNASE3", "CCR3", "SIGLEC8", "ALOX15", "GATA1"],
    )
    blueprint.add_gene_module(
        "basophil_mast",
        ["MS4A2", "TPSAB1", "TPSB2", "CPA3", "HDC", "KIT", "GATA2", "CLC", "FCER1A"],
    )
    blueprint.add_gene_module(
        "monocyte_macrophage_core",
        ["LYZ", "LST1", "CTSS", "S100A8", "S100A9", "FCN1", "AIF1", "TYROBP", "CD14"],
    )
    blueprint.add_gene_module(
        "classical_monocyte",
        ["CD14", "FCN1", "VCAN", "S100A8", "S100A9", "FPR1", "MNDA", "CTSS", "LYZ"],
    )
    blueprint.add_gene_module(
        "nonclassical_monocyte",
        ["FCGR3A", "MS4A7", "CX3CR1", "LILRB1", "IFITM3", "ABI3", "RHOB", "TCF7L2", "HLA-DPA1"],
    )
    blueprint.add_gene_module(
        "marrow_macrophage",
        ["C1QA", "C1QB", "C1QC", "CD68", "APOE", "MRC1", "MARCO", "MSR1", "LIPA"],
    )
    blueprint.add_gene_module(
        "dendritic_core",
        ["HLA-DRA", "HLA-DPA1", "HLA-DPB1", "HLA-DRB1", "ITGAX", "FLT3", "IRF8", "CD74"],
    )
    blueprint.add_gene_module(
        "cdc2",
        ["CD1C", "FCER1A", "CLEC10A", "CD1E", "ENHO", "HLA-DQA1", "HLA-DQB1", "IRF4"],
    )
    blueprint.add_gene_module(
        "cdc1",
        ["CLEC9A", "XCR1", "BATF3", "CADM1", "IRF8", "THBD", "IDO1", "CLNK"],
    )
    blueprint.add_gene_module(
        "pdc",
        ["CLEC4C", "IL3RA", "LILRA4", "TCF4", "GZMB", "IRF7", "PLD4", "SCT", "JCHAIN"],
    )
    blueprint.add_gene_module(
        "lymphocyte_core",
        ["CD3D", "CD3E", "TRAC", "NKG7", "GNLY", "CD79A", "MS4A1", "IL7R", "LTB"],
    )
    blueprint.add_gene_module(
        "b_lineage",
        ["CD79A", "CD79B", "MS4A1", "CD19", "BANK1", "PAX5", "CD74", "HLA-DRA", "BLK"],
    )
    blueprint.add_gene_module(
        "developing_b",
        ["DNTT", "RAG1", "RAG2", "VPREB1", "VPREB3", "IGLL1", "MME", "EBF1", "PAX5"],
    )
    blueprint.add_gene_module(
        "plasma_cell",
        ["MZB1", "XBP1", "PRDM1", "SDC1", "JCHAIN", "IGKC", "IGHG1", "DERL3", "FKBP11"],
    )
    blueprint.add_gene_module(
        "t_cell_core",
        ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "IL7R", "LTB", "TCF7"],
    )
    blueprint.add_gene_module(
        "cd4_t",
        ["CD4", "IL7R", "CCR7", "TCF7", "LEF1", "MAL", "ANXA1", "LTB"],
    )
    blueprint.add_gene_module(
        "cd8_cytotoxic_t",
        ["CD8A", "CD8B", "GZMK", "GZMH", "NKG7", "PRF1", "CCL5", "CST7", "KLRG1"],
    )
    blueprint.add_gene_module(
        "treg",
        ["FOXP3", "IL2RA", "CTLA4", "IKZF2", "TIGIT", "TNFRSF18", "LAYN", "RTN3"],
    )
    blueprint.add_gene_module(
        "innate_like_t",
        ["TRDC", "TRGC1", "TRGC2", "KLRB1", "TRAV1-2", "SLC4A10", "ZBTB16", "NCR3", "CXCR6"],
    )
    blueprint.add_gene_module(
        "nk_cell",
        ["NKG7", "GNLY", "PRF1", "KLRD1", "KLRF1", "FGFBP2", "FCGR3A", "TYROBP", "GZMB"],
    )
    blueprint.add_gene_module(
        "endothelial",
        ["PECAM1", "VWF", "KDR", "CLDN5", "EMCN", "ESAM", "RAMP2", "FLT1", "CDH5"],
    )
    blueprint.add_gene_module(
        "mesenchymal_stromal",
        ["CXCL12", "LEPR", "PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "DCN", "LUM", "PRRX1", "THY1"],
    )
    blueprint.add_gene_module(
        "osteolineage",
        ["RUNX2", "SP7", "ALPL", "BGLAP", "IBSP", "COL1A1", "COL1A2", "POSTN", "SPARC"],
    )
    blueprint.add_gene_module(
        "adipocyte",
        ["ADIPOQ", "PLIN1", "FABP4", "LPL", "CEBPA", "PPARG", "LEP", "APIPO", "CIDEC"],
    )
    blueprint.add_gene_module(
        "interferon_response",
        ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "IFI44L", "RSAD2", "IRF7"],
    )
    blueprint.add_gene_module(
        "inflammatory_activation",
        ["TNF", "NFKBIA", "JUN", "FOS", "DUSP1", "CXCL8", "CCL3", "CCL4", "IER3"],
    )

    blueprint.add_axis(
        "primitive_hspc",
        module="hspc_primitive",
    )
    blueprint.add_gated_axis(
        "cycling_hspc_progenitor",
        gate_module="hspc_primitive",
        positive_module="cell_cycle",
        negative_modules=[],
        negative_weight=0.5,
    )
    blueprint.add_contrast_axis(
        "lymphoid_progenitor_vs_primitive_hspc",
        positive_module="lymphoid_progenitor",
        negative_module="hspc_primitive",
    )
    blueprint.add_contrast_axis(
        "erythroid_vs_leukocyte",
        positive_module="erythroid_early",
        negative_module="leukocyte_core",
    )
    blueprint.add_contrast_axis(
        "erythroid_maturation",
        positive_module="erythroid_late",
        negative_module="erythroid_early",
    )
    blueprint.add_contrast_axis(
        "megakaryocyte_platelet_vs_leukocyte",
        positive_module="megakaryocyte_platelet",
        negative_module="leukocyte_core",
    )
    blueprint.add_contrast_axis(
        "lymphoid_vs_myeloid",
        positive_module="lymphocyte_core",
        negative_module="myeloid_core",
    )
    blueprint.add_contrast_axis(
        "granulocyte_vs_monocyte_macrophage",
        positive_module="granulocyte_immature",
        negative_module="monocyte_macrophage_core",
    )
    blueprint.add_gated_axis(
        "neutrophil_maturation",
        gate_module="myeloid_core",
        positive_module="neutrophil_mature",
        negative_modules=["granulocyte_immature"],
        negative_weight=0.7,
    )
    blueprint.add_gated_axis(
        "eosinophil_granulocyte_branch",
        gate_module="myeloid_core",
        positive_module="eosinophil",
        negative_modules=["neutrophil_mature", "basophil_mast", "classical_monocyte"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "basophil_mast_branch",
        gate_module="myeloid_core",
        positive_module="basophil_mast",
        negative_modules=["neutrophil_mature", "eosinophil", "classical_monocyte"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "classical_monocyte_branch",
        gate_module="monocyte_macrophage_core",
        positive_module="classical_monocyte",
        negative_modules=["nonclassical_monocyte", "marrow_macrophage"],
        negative_weight=0.6,
    )
    blueprint.add_gated_axis(
        "nonclassical_monocyte_branch",
        gate_module="monocyte_macrophage_core",
        positive_module="nonclassical_monocyte",
        negative_modules=["classical_monocyte", "marrow_macrophage"],
        negative_weight=0.6,
    )
    blueprint.add_gated_axis(
        "marrow_macrophage_branch",
        gate_module="monocyte_macrophage_core",
        positive_module="marrow_macrophage",
        negative_modules=["classical_monocyte", "nonclassical_monocyte"],
        negative_weight=0.6,
    )
    blueprint.add_gated_axis(
        "cdc2_branch",
        gate_module="dendritic_core",
        positive_module="cdc2",
        negative_modules=["cdc1", "pdc", "monocyte_macrophage_core"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "cdc1_branch",
        gate_module="dendritic_core",
        positive_module="cdc1",
        negative_modules=["cdc2", "pdc", "monocyte_macrophage_core"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "pdc_branch",
        gate_module="dendritic_core",
        positive_module="pdc",
        negative_modules=["cdc1", "cdc2", "b_lineage"],
        negative_weight=0.5,
    )
    blueprint.add_contrast_axis(
        "b_lineage_vs_t_nk",
        positive_module="b_lineage",
        negative_module="t_cell_core",
    )
    blueprint.add_gated_axis(
        "developing_b_branch",
        gate_module="b_lineage",
        positive_module="developing_b",
        negative_modules=["plasma_cell"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "plasma_cell_branch",
        gate_module="b_lineage",
        positive_module="plasma_cell",
        negative_modules=["developing_b"],
        negative_weight=0.5,
    )
    blueprint.add_contrast_axis(
        "t_cell_vs_nk",
        positive_module="t_cell_core",
        negative_module="nk_cell",
    )
    blueprint.add_gated_axis(
        "cd4_t_branch",
        gate_module="t_cell_core",
        positive_module="cd4_t",
        negative_modules=["cd8_cytotoxic_t", "treg", "innate_like_t"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "cd8_cytotoxic_t_branch",
        gate_module="t_cell_core",
        positive_module="cd8_cytotoxic_t",
        negative_modules=["cd4_t", "treg", "innate_like_t"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "treg_branch",
        gate_module="t_cell_core",
        positive_module="treg",
        negative_modules=["cd4_t", "cd8_cytotoxic_t", "innate_like_t"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "innate_like_t_branch",
        gate_module="t_cell_core",
        positive_module="innate_like_t",
        negative_modules=["cd4_t", "cd8_cytotoxic_t", "treg"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "nk_branch",
        gate_module="lymphocyte_core",
        positive_module="nk_cell",
        negative_modules=["t_cell_core", "b_lineage"],
        negative_weight=0.6,
    )
    blueprint.add_contrast_axis(
        "endothelial_vs_hematopoietic",
        positive_module="endothelial",
        negative_module="leukocyte_core",
    )
    blueprint.add_contrast_axis(
        "mesenchymal_stromal_vs_hematopoietic",
        positive_module="mesenchymal_stromal",
        negative_module="leukocyte_core",
    )
    blueprint.add_gated_axis(
        "osteolineage_branch",
        gate_module="mesenchymal_stromal",
        positive_module="osteolineage",
        negative_modules=["adipocyte"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "adipocyte_branch",
        gate_module="mesenchymal_stromal",
        positive_module="adipocyte",
        negative_modules=["osteolineage"],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "leukocyte_interferon_state",
        gate_module="leukocyte_core",
        positive_module="interferon_response",
        negative_modules=[],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "myeloid_inflammatory_state",
        gate_module="myeloid_core",
        positive_module="inflammatory_activation",
        negative_modules=[],
        negative_weight=0.5,
    )
    blueprint.add_gated_axis(
        "lymphocyte_activation_state",
        gate_module="lymphocyte_core",
        positive_module="inflammatory_activation",
        negative_modules=[],
        negative_weight=0.5,
    )

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
            "Build a count-matrix cell embedding from a declarative blueprint "
            "(no-ref variant: no literature-resource citations)."
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
