"""
Sweep UMAP plots across every embedding parquet under outs/embeddings/<dataset>/.

For each parquet, a 2D UMAP is computed once (cached as <stem>_umap_coords.parquet)
and PNG + PDF figures are rendered, one per dataset-specific obs column,
into outs/plots/umap/<dataset>/<model>/ — mirroring the embeddings layout.

Inputs:
    outs/embeddings/<dataset>/*.parquet         (cell_id + dim_* columns)
    data/preprocessed/<dataset>_pp.h5ad         (obs columns for coloring)

Outputs (per embedding parquet):
    outs/plots/umap/<dataset>/<model>/<stem>_umap_coords.parquet
    outs/plots/umap/<dataset>/<model>/<stem>_<column>.png
    outs/plots/umap/<dataset>/<model>/<stem>_<column>.pdf

Usage:
    python scripts/plots/umap_plots.py --dataset mouse_brain
    python scripts/plots/umap_plots.py --dataset human_pancreas
    python scripts/plots/umap_plots.py --dataset human_pancreas \
        --input outs/embeddings/human_pancreas/pca_human_pancreas_seed0.parquet
"""
from __future__ import annotations

import argparse
import colorsys
import math
import os
import re
import tempfile
import time
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap


REPO_ROOT = Path(__file__).resolve().parents[2]


# ---------------------------------------------------------------- categorical palettes
UNMATCHED_COLOR_HEX = "#d0d0d0"

_FAMILY_ORDER = [
    "Astro", "Endo", "Igtp", "Micro", "OPC", "Oligo",
    "L2", "L2/3", "L4", "L5", "L5a", "L5b", "L6a", "L6b",
    "Ndnf", "Vip", "Sncg", "Smad3", "Sst", "Pvalb",
    "Unclassified", "unmatched",
]
_FAMILY_RANK = {family: index for index, family in enumerate(_FAMILY_ORDER)}

# Tasic 2016 cell-type families. Cool hues (blue/teal/green) cover excitatory L2-L6,
# warm hues (pink/purple/orange) cover GABA interneurons; the remaining glia/vascular
# live in distinct neutral hues.
_FAMILY_BASE_COLORS = {
    "Astro": "#4c78a8", "Endo": "#72b7b2", "Igtp": "#54a24b",
    "Micro": "#9d6fb8", "OPC": "#9c755f", "Oligo": "#bab0ab",
    "L2": "#3454d1", "L2/3": "#4c78e0", "L4": "#48b6ff",
    "L5": "#00a6a6", "L5a": "#12b886", "L5b": "#52b788",
    "L6a": "#84cc16", "L6b": "#b8de29",
    "Ndnf": "#f2c6de", "Vip": "#ff4f9a", "Sncg": "#ef476f",
    "Smad3": "#cc5de8", "Sst": "#7b2cbf", "Pvalb": "#f77f00",
    "Unclassified": "#8d99ae", "unmatched": UNMATCHED_COLOR_HEX,
}

# Coarse Tasic broad types — each picks a representative subclass color so the
# hierarchy stays consistent (e.g. Pvalb-orange ↔ "GABA-ergic Neuron").
_BROAD_TYPE_COLORS = {
    "Glutamatergic Neuron": "#4c78e0",
    "GABA-ergic Neuron": "#f77f00",
    "Astrocyte": "#4c78a8",
    "Microglia": "#9d6fb8",
    "Oligodendrocyte": "#bab0ab",
    "Oligodendrocyte Precursor Cell": "#9c755f",
    "Endothelial Cell": "#72b7b2",
    "Unclassified": "#8d99ae",
}
_BROAD_TYPE_ORDER = list(_BROAD_TYPE_COLORS.keys())
_BROAD_TYPE_RANK = {name: i for i, name in enumerate(_BROAD_TYPE_ORDER)}

# Human pancreas celltypes grouped by function: endocrine (warm), exocrine (greens),
# stromal (browns), immune (purples), vascular (blue), glial (cyan).
_PANCREAS_CELLTYPE_COLORS = {
    "alpha":              "#f77f00",
    "beta":               "#d62828",
    "delta":              "#fcbf49",
    "gamma":              "#e85d04",
    "epsilon":            "#ffd60a",
    "acinar":             "#2a9d8f",
    "ductal":             "#52b788",
    "activated_stellate": "#7f4f24",
    "quiescent_stellate": "#c08552",
    "macrophage":         "#7209b7",
    "mast":               "#c77dff",
    "t_cell":             "#5a189a",
    "endothelial":        "#023e8a",
    "schwann":            "#06aed5",
}
_PANCREAS_CELLTYPE_ORDER = list(_PANCREAS_CELLTYPE_COLORS.keys())
_PANCREAS_CELLTYPE_RANK = {name: i for i, name in enumerate(_PANCREAS_CELLTYPE_ORDER)}

def _generic_palette_cycle():
    palette = []
    for name in ("tab20", "tab20b", "tab20c"):
        cmap = matplotlib.colormaps[name]
        for i in range(cmap.N):
            palette.append(cmap(i))
    return palette


def _batch_palette_cycle():
    # Set2/Set3 pastels — visually distinct from the saturated tab20-style cell-type
    # palette, so a glance at the figure tells you whether you're looking at biology
    # or a technical/batch axis.
    palette = []
    for name in ("Set2", "Set3"):
        cmap = matplotlib.colormaps[name]
        for i in range(cmap.N):
            palette.append(cmap(i))
    return palette


def _vary_family_color(base_hex: str, index: int, total: int):
    r, g, b = mcolors.to_rgb(base_hex)
    h, l, s = colorsys.rgb_to_hls(r, g, b)
    if total <= 1:
        return (r, g, b, 1.0)
    span = 0.22
    frac = index / max(total - 1, 1)
    lightness = min(0.82, max(0.32, l - span / 2 + frac * span))
    saturation = min(0.95, max(0.45, s + 0.08 * (0.5 - frac)))
    rr, gg, bb = colorsys.hls_to_rgb(h, lightness, saturation)
    return (rr, gg, bb, 1.0)


def _build_family_palette(labels):
    """Family-ordered palette for fine-grained labels like 'Pvalb Cpne5'."""
    unique_labels = sorted(str(value) for value in labels.dropna().unique())
    if not unique_labels:
        return {}, [], []

    def infer_family(label: str) -> str:
        token = label.split()[0] if label else label
        return token if token in _FAMILY_BASE_COLORS else label

    grouped: dict[str, list[str]] = {}
    for label in unique_labels:
        grouped.setdefault(infer_family(label), []).append(label)

    family_sort = sorted(
        grouped, key=lambda f: (_FAMILY_RANK.get(f, len(_FAMILY_RANK)), f),
    )

    ordered_labels: list[str] = []
    color_map: dict[str, tuple] = {}
    for family in family_sort:
        family_labels = sorted(grouped[family])
        ordered_labels.extend(family_labels)
        base_color = _FAMILY_BASE_COLORS.get(family, "#666666")
        for index, label in enumerate(family_labels):
            color_map[label] = _vary_family_color(base_color, index, len(family_labels))

    point_colors = [color_map.get(str(value), "#bdbdbd") for value in labels]
    return color_map, point_colors, ordered_labels


def _build_subclass_palette(labels):
    """Subclass labels (e.g. 'Pvalb', 'Vip') map directly to family base colors."""
    unique_labels = sorted(str(value) for value in labels.dropna().unique())
    if not unique_labels:
        return {}, [], []

    fallback = _generic_palette_cycle()
    fallback_idx = 0
    color_map: dict[str, tuple] = {}
    for label in unique_labels:
        if label in _FAMILY_BASE_COLORS:
            color_map[label] = mcolors.to_rgba(_FAMILY_BASE_COLORS[label])
        else:
            color_map[label] = fallback[fallback_idx % len(fallback)]
            fallback_idx += 1

    ordered_labels = sorted(
        unique_labels, key=lambda lbl: (_FAMILY_RANK.get(lbl, len(_FAMILY_RANK)), lbl),
    )
    point_colors = [color_map.get(str(value), UNMATCHED_COLOR_HEX) for value in labels]
    return color_map, point_colors, ordered_labels


def _build_broad_type_palette(labels):
    """Coarse broad_type → representative subclass color."""
    unique_labels = sorted(str(value) for value in labels.dropna().unique())
    if not unique_labels:
        return {}, [], []

    fallback = _generic_palette_cycle()
    fallback_idx = 0
    color_map: dict[str, tuple] = {}
    for label in unique_labels:
        if label in _BROAD_TYPE_COLORS:
            color_map[label] = mcolors.to_rgba(_BROAD_TYPE_COLORS[label])
        else:
            color_map[label] = fallback[fallback_idx % len(fallback)]
            fallback_idx += 1

    ordered_labels = sorted(
        unique_labels, key=lambda lbl: (_BROAD_TYPE_RANK.get(lbl, len(_BROAD_TYPE_RANK)), lbl),
    )
    point_colors = [color_map.get(str(value), UNMATCHED_COLOR_HEX) for value in labels]
    return color_map, point_colors, ordered_labels


def _build_pancreas_celltype_palette(labels):
    """human_pancreas celltype: functionally grouped warm/green/brown/purple."""
    unique_labels = sorted(str(value) for value in labels.dropna().unique())
    if not unique_labels:
        return {}, [], []

    fallback = _generic_palette_cycle()
    fallback_idx = 0
    color_map: dict[str, tuple] = {}
    for label in unique_labels:
        if label in _PANCREAS_CELLTYPE_COLORS:
            color_map[label] = mcolors.to_rgba(_PANCREAS_CELLTYPE_COLORS[label])
        else:
            color_map[label] = fallback[fallback_idx % len(fallback)]
            fallback_idx += 1

    ordered_labels = sorted(
        unique_labels,
        key=lambda lbl: (_PANCREAS_CELLTYPE_RANK.get(lbl, len(_PANCREAS_CELLTYPE_RANK)), lbl),
    )
    point_colors = [color_map.get(str(value), UNMATCHED_COLOR_HEX) for value in labels]
    return color_map, point_colors, ordered_labels


def _build_batch_palette(labels):
    """Qualitative pastel palette (Set2/Set3) for technical/batch columns."""
    unique_labels = sorted(str(value) for value in labels.dropna().unique())
    if not unique_labels:
        return {}, [], []

    palette = _batch_palette_cycle()
    color_map = {lbl: palette[idx % len(palette)] for idx, lbl in enumerate(unique_labels)}
    point_colors = [color_map.get(str(value), UNMATCHED_COLOR_HEX) for value in labels]
    return color_map, point_colors, unique_labels


def _build_generic_palette(labels):
    unique_labels = sorted(str(value) for value in labels.dropna().unique())
    if not unique_labels:
        return {}, [], []

    palette = _generic_palette_cycle()
    color_map = {label: palette[idx % len(palette)] for idx, label in enumerate(unique_labels)}
    point_colors = [color_map.get(str(value), UNMATCHED_COLOR_HEX) for value in labels]
    return color_map, point_colors, unique_labels


def _build_ordinal_palette(labels, *, cmap_name: str = "viridis"):
    """Sequential palette for ordinal labels (e.g. 'E6.5', 'E7.0', ...).

    Labels are sorted by their leading numeric component so colours track
    progression along the axis (early = dark purple, late = yellow under
    viridis). Falls back to alphabetical order if no number can be parsed.
    """
    unique_labels = [str(value) for value in labels.dropna().unique()]
    if not unique_labels:
        return {}, [], []

    def _to_number(s: str) -> float:
        m = re.search(r"E(\d+(?:\.\d+)?)", s)
        if m is None:
            m = re.search(r"(\d+(?:\.\d+)?)", s)
        return float(m.group(1)) if m else float("nan")

    pairs = [(s, _to_number(s)) for s in unique_labels]
    if any(math.isnan(v) for _, v in pairs):
        ordered_labels = sorted(unique_labels)
    else:
        ordered_labels = [s for s, _ in sorted(pairs, key=lambda x: x[1])]

    cmap = matplotlib.colormaps[cmap_name]
    n = len(ordered_labels)
    positions = [0.5] if n == 1 else list(np.linspace(0.05, 0.95, n))
    color_map = {label: cmap(pos) for label, pos in zip(ordered_labels, positions)}
    point_colors = [color_map.get(str(value), UNMATCHED_COLOR_HEX) for value in labels]
    return color_map, point_colors, ordered_labels


def build_categorical_colors(labels, dataset=None, column=None, kind=None):
    """Return (color_map, point_colors, ordered_labels) for a categorical Series.

    Dispatches by (dataset, column, kind):
      - kind == "batch" → Set2/Set3 pastel cycle (distinct from cell-type palettes)
      - mouse_brain primary_type/secondary_type/cluster_id → family-ordered, varied shades
      - mouse_brain subclass → family base colors looked up directly
      - mouse_brain broad_type → coarse-level map mirroring the family hierarchy
      - human_pancreas celltype → endocrine-warm / exocrine-green / stromal-brown / immune-purple
      - everything else → tab20/tab20b/tab20c cycle
    """
    if kind == "batch":
        return _build_batch_palette(labels)
    if kind == "ordinal":
        return _build_ordinal_palette(labels)
    if dataset == "mouse_brain":
        if column == "subclass":
            return _build_subclass_palette(labels)
        if column == "broad_type":
            return _build_broad_type_palette(labels)
        if column in {"primary_type", "secondary_type", "cluster_id"}:
            return _build_family_palette(labels)
    if dataset == "human_pancreas" and column == "celltype":
        return _build_pancreas_celltype_palette(labels)
    return _build_generic_palette(labels)


def style_axes(ax, title: str) -> None:
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1", labelpad=8)
    ax.set_ylabel("UMAP 2", labelpad=8)
    if title:
        ax.set_title(title, pad=8)


def build_legend_handles(color_map, ordered_labels):
    return [
        plt.Line2D(
            [0], [0], marker="o", linestyle="",
            markerfacecolor=color_map[label], markeredgecolor="none",
            markersize=5.5, label=label,
        )
        for label in ordered_labels
    ]


def add_inline_legend(fig, ax, handles, title: str) -> None:
    legend_cols = 1 if len(handles) <= 16 else 2 if len(handles) <= 32 else 3
    ax.legend(
        handles=handles, title=title,
        loc="upper left", bbox_to_anchor=(1.02, 1.0),
        frameon=False, ncol=legend_cols, borderaxespad=0.0,
        fontsize=7.5, title_fontsize=8,
        handletextpad=0.4, columnspacing=1.0, labelspacing=0.45,
    )
    fig.subplots_adjust(right=0.72)


def separate_legend_layout(n_items: int):
    ncol = 2 if n_items <= 20 else 3 if n_items <= 36 else 4 if n_items <= 56 else 5
    nrow = math.ceil(n_items / ncol)
    width = max(4.6, 1.6 * ncol + 0.8)
    height = 0.32 * nrow + 0.6
    return ncol, nrow, width, height


def add_bottom_legend(ax, handles, title: str) -> None:
    ncol, _, _, _ = separate_legend_layout(len(handles))
    ax.axis("off")
    ax.legend(
        handles=handles, title=title,
        loc="center", frameon=False, ncol=ncol,
        fontsize=7.3, title_fontsize=8,
        handletextpad=0.4, columnspacing=1.1, labelspacing=0.45,
        borderaxespad=0.0,
    )
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")

DEFAULT_PARAMS = {
    "mouse_brain": {
        "color_keys": ["primary_type", "subclass", "broad_type"],
        "column_kinds": {},  # all defaults to "label"
        "point_size": 14.0,  # 1.8k cells: prominent dots
    },
    "human_pancreas": {
        "color_keys": ["celltype", "tech"],
        "column_kinds": {"tech": "batch"},
        "point_size": 4.0,   # 16k cells: small dots to reduce overdraw
    },
    "mouse_gastrulation": {
        "color_keys": ["celltype", "stage", "sample_id"],
        "column_kinds": {"sample_id": "batch", "stage": "ordinal"},
        "point_size": 1.5,   # ~116k cells: tiny dots
    },
    "mouse_cortex": {
        # three-source integration: celltype_harmonized (~24 Tasic-2018 subclasses)
        # and class_harmonized (3 lineage roots) are the biology; study is the
        # study-of-origin batch; protocol is the technology batch (SMART-Seq /
        # SMART-Seq v4 / 10x) that the embedding should mix; area (VISp/ALM/cortex)
        # is a within-cortex biological axis.
        "color_keys": ["celltype_harmonized", "class_harmonized",
                       "study", "protocol", "area"],
        "column_kinds": {"study": "batch", "protocol": "batch"},
        "point_size": 1.5,   # ~45k cells: small dots to reduce overdraw
    },
    "human_immune": {
        # scIB "Immune (human)" benchmark: final_annotation (16 immune cell
        # types) is the biology; batch (10 donor/study batches) is what the
        # embedding should mix; study (5 source studies) and chemistry
        # (10x v2/v3, Smart-seq2) are coarser technical splits; tissue
        # (bone marrow vs PBMC) is a within-immune biological axis.
        "color_keys": ["final_annotation", "batch", "study", "chemistry", "tissue"],
        "column_kinds": {"batch": "batch", "study": "batch", "chemistry": "batch"},
        "point_size": 2.0,   # ~33k cells: small dots to reduce overdraw
    },
}

PLOT_RC = {
    "font.family": "DejaVu Sans",
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
    "savefig.facecolor": "white", "savefig.transparent": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
}


def parse_stem(stem: str, dataset: str) -> tuple[str, int | None]:
    """Recover (model, seed) from '<model>_<dataset>[_seed<N>]'."""
    if f"_{dataset}" not in stem:
        raise ValueError(f"Stem {stem!r} does not contain dataset {dataset!r}.")
    head, _, tail = stem.partition(f"_{dataset}")
    if not head:
        raise ValueError(f"Empty model name parsed from stem {stem!r}.")
    if not tail:
        return head, None
    if tail.startswith("_seed"):
        return head, int(tail[len("_seed"):])
    return head + tail, None


def discover_embeddings(embeddings_dir: Path, override: Path | None) -> list[Path]:
    if override is not None:
        if not override.exists():
            raise FileNotFoundError(f"--input not found: {override}")
        return [override.resolve()]
    if not embeddings_dir.exists():
        raise FileNotFoundError(f"Embeddings directory not found: {embeddings_dir}")
    parquets = [
        p for p in sorted(embeddings_dir.glob("*.parquet"))
        if not p.name.endswith("_umap_coords.parquet")
    ]
    if not parquets:
        raise FileNotFoundError(f"No embedding parquets in {embeddings_dir}")
    return parquets


def load_obs_columns(h5ad_path: Path, columns: list[str]) -> dict[str, pd.Series]:
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Labels h5ad not found: {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path, backed="r")
    missing = [c for c in columns if c not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"{h5ad_path} is missing obs columns {missing}. "
            f"Available: {list(adata.obs.columns)}"
        )
    obs_index = pd.Index(adata.obs_names.astype(str), name="cell_id")
    return {
        col: pd.Series(adata.obs[col].astype(str).to_numpy(), index=obs_index, name=col)
        for col in columns
    }


def compute_or_load_umap(parquet_path: Path, coords_path: Path, params: dict, overwrite: bool):
    df = pd.read_parquet(parquet_path)
    if "cell_id" not in df.columns:
        raise ValueError(f"{parquet_path} has no 'cell_id' column.")
    dim_cols = sorted(
        (c for c in df.columns if c.startswith("dim_")),
        key=lambda n: int(n.split("_")[1]),
    )
    if not dim_cols:
        raise ValueError(f"{parquet_path} has no dim_* columns.")
    df = df.sort_values("cell_id", kind="mergesort").reset_index(drop=True)
    cell_ids = df["cell_id"].astype(str)

    if coords_path.exists() and not overwrite:
        cached = pd.read_parquet(coords_path)
        if cached["cell_id"].astype(str).tolist() != cell_ids.tolist():
            raise ValueError(f"Cached coords {coords_path} disagree with {parquet_path}.")
        return cell_ids, cached[["umap_0", "umap_1"]].to_numpy(), len(dim_cols), True

    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=params["n_neighbors"], min_dist=params["min_dist"], metric=params["metric"],
        init="spectral",
        random_state=params["random_state"], transform_seed=params["random_state"],
        low_memory=True,
    )
    coords = reducer.fit_transform(df[dim_cols].to_numpy(dtype=float))

    coords_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"cell_id": cell_ids, "umap_0": coords[:, 0], "umap_1": coords[:, 1]}).to_parquet(
        coords_path, index=False
    )
    return cell_ids, coords, len(dim_cols), False


def render_panel(
    cell_ids: pd.Series, coords, label_series: pd.Series, column: str, dataset: str,
    out_png: Path, out_pdf: Path, params: dict, kind: str | None = None,
) -> None:
    label_values = pd.Series(
        cell_ids.map(label_series).fillna("unmatched").astype(str).to_numpy(),
        index=cell_ids.index, name=column,
    )
    color_map, point_colors, ordered_labels = build_categorical_colors(
        label_values, dataset=dataset, column=column, kind=kind,
    )
    legend_labels = [lbl for lbl in ordered_labels if lbl != "unmatched"]
    legend_color_map = {k: v for k, v in color_map.items() if k != "unmatched"}
    handles = build_legend_handles(legend_color_map, legend_labels)
    legend_title = column.replace("_", " ").capitalize()

    legend_mode = params["legend_mode"]
    if legend_mode == "auto":
        legend_mode = "inline" if len(legend_labels) <= 16 else "separate"

    panel = params["panel_size"]
    if legend_mode == "separate" and handles:
        _, _, lw, lh = separate_legend_layout(len(handles))
        fig = plt.figure(figsize=(max(panel, lw), panel + lh))
        gs = fig.add_gridspec(2, 1, height_ratios=[panel, lh], hspace=0.04)
        ax = fig.add_subplot(gs[0])
        legend_ax = fig.add_subplot(gs[1])
    else:
        fig, ax = plt.subplots(figsize=(panel, panel))
        legend_ax = None

    # Draw largest groups first so the smallest categories sit on top — keeps
    # minority clusters visible regardless of cell_id ordering. Cached coords
    # are untouched (this only reorders the scatter call).
    counts = label_values.value_counts()
    priority = np.array([counts[lbl] for lbl in label_values], dtype=int)
    order = np.argsort(-priority, kind="mergesort")
    ax.scatter(
        coords[order, 0], coords[order, 1],
        s=params["point_size"], alpha=params["alpha"], linewidths=0,
        c=[point_colors[i] for i in order], rasterized=True,
    )
    style_axes(ax, title="")
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0.03)

    if legend_mode == "inline" and handles:
        add_inline_legend(fig, ax, handles, legend_title)
    elif legend_mode == "separate" and handles:
        add_bottom_legend(legend_ax, handles, legend_title)

    if legend_ax is not None:
        fig.subplots_adjust(left=0.08, right=0.98, top=0.98, bottom=0.05, hspace=0.12)
    else:
        fig.tight_layout(pad=0.4)

    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=params["png_dpi"], bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, dpi=params["pdf_dpi"], bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _safe_col(col: str) -> str:
    """Composite color keys use '+' (e.g. 'typeSample+species'); replace it with
    '_x_' for the on-disk filename so the stem stays glob/URL-safe."""
    return col.replace("+", "_x_")


def existing_outputs(out_dir: Path, stem: str, color_keys: list[str], coords_path: Path) -> list[Path]:
    targets = [coords_path] if coords_path.exists() else []
    for col in color_keys:
        for ext in ("png", "pdf"):
            p = out_dir / f"{stem}_{_safe_col(col)}.{ext}"
            if p.exists():
                targets.append(p)
    return targets


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Batch UMAP plotting across outs/embeddings/.")
    parser.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    parser.add_argument("--input", type=Path, help="Single parquet to plot (overrides discovery).")
    parser.add_argument("--color-keys", type=str,
                        help="Comma-separated obs columns; overrides DEFAULT_PARAMS. "
                             "Composite keys via '+' (e.g. 'typeSample+species') colour by the joined "
                             "string of those obs columns; the on-disk filename uses '_x_' in place of '+'.")
    parser.add_argument("--n-neighbors", type=int, default=30)
    parser.add_argument("--min-dist", type=float, default=0.15)
    parser.add_argument("--metric", type=str, default="euclidean")
    parser.add_argument("--random-state", type=int, default=42)
    parser.add_argument("--point-size", type=float, default=None,
                        help="Override per-dataset DEFAULT_PARAMS point_size.")
    parser.add_argument("--alpha", type=float, default=0.9)
    parser.add_argument("--panel-size", type=float, default=4.6)
    parser.add_argument("--png-dpi", type=int, default=600)
    parser.add_argument("--pdf-dpi", type=int, default=600)
    parser.add_argument("--legend-mode", choices=("auto", "inline", "separate", "none"), default="separate",
                        help="Default 'separate' keeps a consistent below-panel layout across color keys.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    matplotlib.rcParams.update(PLOT_RC)

    color_keys = (
        [c.strip() for c in args.color_keys.split(",") if c.strip()]
        if args.color_keys else list(DEFAULT_PARAMS[args.dataset]["color_keys"])
    )
    if not color_keys:
        raise ValueError("At least one color key is required.")

    embeddings_dir = REPO_ROOT / "outs" / "embeddings" / args.dataset
    output_root    = REPO_ROOT / "outs" / "plots" / "umap" / args.dataset
    labels_h5ad    = REPO_ROOT / "data" / "preprocessed" / f"{args.dataset}_pp.h5ad"

    parquets = discover_embeddings(embeddings_dir, args.input)
    # Resolve composite color keys (e.g. "typeSample+species") into the underlying
    # obs columns to load, then build a joined string Series for each
    # composite. Atomic keys ("celltype", "species", …) flow through unchanged.
    needed_cols: list[str] = []
    for key in color_keys:
        for part in key.split("+"):
            if part not in needed_cols:
                needed_cols.append(part)
    obs_by_column = load_obs_columns(labels_h5ad, needed_cols)
    for key in color_keys:
        if "+" not in key:
            continue
        parts = key.split("+")
        joined = obs_by_column[parts[0]].astype(str)
        for p in parts[1:]:
            joined = joined.str.cat(obs_by_column[p].astype(str), sep=" | ")
        joined.name = key
        obs_by_column[key] = joined

    umap_params = {k: getattr(args, k) for k in ("n_neighbors", "min_dist", "metric", "random_state")}
    point_size = (
        args.point_size if args.point_size is not None
        else DEFAULT_PARAMS[args.dataset]["point_size"]
    )
    plot_params = {
        "point_size": point_size,
        "alpha": args.alpha, "panel_size": args.panel_size,
        "png_dpi": args.png_dpi, "pdf_dpi": args.pdf_dpi,
        "legend_mode": args.legend_mode,
    }
    column_kinds = DEFAULT_PARAMS[args.dataset].get("column_kinds", {})

    print(f"[ok] dataset={args.dataset} files={len(parquets)} color_keys={color_keys}")

    for parquet_path in parquets:
        stem = parquet_path.stem
        model, seed = parse_stem(stem, args.dataset)
        out_dir = output_root / model
        coords_path = out_dir / f"{stem}_umap_coords.parquet"

        existing = existing_outputs(out_dir, stem, color_keys, coords_path)
        if existing and not args.overwrite:
            raise FileExistsError(
                "Refusing to overwrite (pass --overwrite):\n  "
                + "\n  ".join(str(p) for p in existing)
            )

        t0 = time.perf_counter()
        cell_ids, coords, n_dims, cached = compute_or_load_umap(
            parquet_path, coords_path, umap_params, args.overwrite
        )
        umap_msg = "cached" if cached else f"{time.perf_counter() - t0:.1f}s"
        seed_repr = f"seed={seed}" if seed is not None else "seed=-"
        print(
            f"[ok] {parquet_path.name}: model={model} {seed_repr} "
            f"n_cells={len(cell_ids)} n_dims={n_dims} umap={umap_msg}"
        )

        for column in color_keys:
            safe = _safe_col(column)
            out_png = out_dir / f"{stem}_{safe}.png"
            out_pdf = out_dir / f"{stem}_{safe}.pdf"
            t1 = time.perf_counter()
            render_panel(
                cell_ids, coords, obs_by_column[column], column, args.dataset,
                out_png, out_pdf, plot_params, kind=column_kinds.get(column),
            )
            print(
                f"[ok]   {column}: render={time.perf_counter() - t1:.1f}s "
                f"-> {out_png.relative_to(REPO_ROOT)}"
            )


if __name__ == "__main__":
    main()
