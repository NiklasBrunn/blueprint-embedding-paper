"""Assemble the client-side explore bundle shipped to the browser once.

Mirrors the standalone interactive-UMAP DATA shape so the Explore tab can
recolor entirely client-side. Reuses webapp.explorer / webapp.coloring rather
than re-porting file-driven builders.
"""
from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd

from webapp import explorer
from webapp.coloring import discrete_colors, scheme_for_obs

UNMATCHED = "unmatched"
DEFAULT_CELL_CAP = 75_000


def _round(arr: Any) -> list[float]:
    """Round a 1-D numeric array to 3 decimals and return a plain float list."""
    return [round(float(v), 3) for v in np.asarray(arr, dtype=float)]


def _downsample_index(total: int, cell_cap: int) -> np.ndarray:
    """Deterministic, uniform row indices into [0, total) capped at cell_cap."""
    if total <= cell_cap:
        return np.arange(total)
    idx = np.linspace(0, total - 1, cell_cap).round().astype(int)
    return np.unique(idx)


def displayed_cell_ids(parquet_path, cell_cap: int = DEFAULT_CELL_CAP) -> list[str]:
    """Cell ids shown in the explorer, in plot order (= embedding row order,
    uniformly downsampled to ``cell_cap``).

    Mirrors ``build()``'s ordering exactly so per-gene vectors fetched on demand
    can be aligned to the displayed points. ``build()`` merges the embedding with
    the UMAP coords ``how="left"`` on ``cell_id`` (preserving embedding row order
    and count), then applies the same uniform downsample; reading the embedding's
    ``cell_id`` column here reproduces that order without the coords/value work.
    """
    emb = pd.read_parquet(parquet_path, columns=["cell_id"])
    sel = _downsample_index(len(emb), cell_cap)
    return [str(c) for c in emb["cell_id"].to_numpy()[sel]]


def build(map_result, adata, summary, *, cell_cap: int = DEFAULT_CELL_CAP) -> dict:
    """Assemble the client-side explore bundle (see module docstring)."""
    bp = summary["blueprint"]
    axis_order = [ax["name"] for ax in bp["axes"]]

    emb = pd.read_parquet(map_result.parquet)
    coords = pd.read_parquet(map_result.umap_coords)
    merged = emb.merge(coords[["cell_id", "umap_0", "umap_1"]],
                       on="cell_id", how="left")
    if merged[["umap_0", "umap_1"]].isna().any().any():
        raise ValueError("UMAP coords contain NaN after merge — coords and "
                         "embedding cell_ids do not align")

    dim_cols = [c for c in emb.columns if c.startswith("dim_")]
    dim_cols.sort(key=lambda c: int(c.split("_")[1]))
    if len(dim_cols) != len(axis_order):
        raise ValueError(f"dim columns ({len(dim_cols)}) != axes "
                         f"({len(axis_order)}); embedding/summary mismatch")

    total = len(merged)
    sel = _downsample_index(total, cell_cap)
    shown = len(sel)

    coords_out = {
        "x": _round(merged["umap_0"].to_numpy()[sel]),
        "y": _round(merged["umap_1"].to_numpy()[sel]),
    }
    values = {
        name: _round(merged[dim].to_numpy()[sel])
        for name, dim in zip(axis_order, dim_cols)
    }

    axis_meta = {
        ax["name"]: {
            "kind": ax.get("kind"),
            "positive": ax.get("positive"),
            "gate": ax.get("gate"),
            "negatives": list(ax.get("negatives", [])),
        }
        for ax in bp["axes"]
    }

    branches = explorer.branches_from_summary(summary)

    modules = {m["name"]: [str(g) for g in m.get("genes", [])]
               for m in bp["modules"]}

    gene_list: list[str] = []
    seen: set[str] = set()
    for genes in modules.values():
        for g in genes:
            if g not in seen:
                seen.add(g)
                gene_list.append(g)
    for v in adata.var_names:
        g = str(v)
        if g not in seen:
            seen.add(g)
            gene_list.append(g)

    metadata = _build_metadata(adata, merged["cell_id"].to_numpy(), sel)

    # Tree-plot mode auto-derived from the blueprint's scenario notes (the
    # Generate-mode toggle's marker); the Explore download selector pre-selects
    # this and the user can override it. Local import keeps matplotlib out of
    # this module's import path.
    from webapp.tree_plot import detect_tree_kind
    tree_kind = detect_tree_kind((summary.get("scenario") or {}).get("notes"))

    return {
        "coords": coords_out,
        "axisOrder": axis_order,
        "values": values,
        "branches": branches,
        "axisMeta": axis_meta,
        "metadata": metadata,
        "modules": modules,
        "geneList": gene_list,
        "shown": shown,
        "total": total,
        "treeKind": tree_kind,
    }


def _build_metadata(adata, cell_ids: np.ndarray, sel: np.ndarray) -> dict:
    """Categorical obs columns aligned to embedding cell_id order, downsampled."""
    obs = adata.obs
    # Align obs to embedding cell_id order (obs index == cell barcodes).
    aligned = obs.reindex([str(c) for c in cell_ids])
    out: dict = {}
    for col in obs.columns:
        series = aligned[col]
        # Cast to object FIRST: pandas validates the `other` sentinel against a
        # Categorical column's existing categories and raises even when nothing
        # is missing ("Cannot setitem on a Categorical with a new category").
        # Real h5ads store obs (e.g. cell_type) as Categorical, so this path is
        # the common case, not the exception.
        raw = series.astype(object).where(series.notna(), UNMATCHED).tolist()
        if scheme_for_obs([v for v in raw if v != UNMATCHED]).get("type") != "categorical":
            continue
        labels = list(dict.fromkeys(raw))
        nunique = len([l for l in labels if l != UNMATCHED])
        if not (2 <= nunique <= 50):
            continue
        cmap = discrete_colors([l for l in labels if l != UNMATCHED])
        colors = [cmap.get(l, "#999999") for l in labels]
        label_idx = {l: i for i, l in enumerate(labels)}
        cell_idx = [int(label_idx[v]) for v in raw]
        cell_idx = [cell_idx[i] for i in sel]
        out[str(col)] = {
            "kind": "label",
            "labels": [str(l) for l in labels],
            "colors": colors,
            "cellIdx": cell_idx,
        }
    return out
