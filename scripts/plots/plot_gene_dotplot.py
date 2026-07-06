"""
Marker-gene grounding figure for ONE blueprint run.

Visualises that the agent's *literature-cited* marker genes are real markers in
the data: each axis's positive-module genes (rows, grouped by the axis that
declares them) against the annotated cell types (columns). Two styles:

  --style dotplot  (default): one dot per (gene, cell type) over
      ALL cells -- dot COLOR = mean z-scored expression, dot SIZE = fraction of
      cells of that type expressing the gene. Axis row-blocks are ordered by the
      per-axis AUROC pattern (argmax cell type, reusing the axis x cell-type AUROC
      matrix) and the cell-type columns by the dataset's canonical biological
      order, so the grounded markers run down a clean diagonal.
  --style heatmap  (alternative): the original single-cell heatmap --
      cells (columns, subsampled) grouped by type, expression z-scored per gene.

The blueprint structure (axes in dim order, each module's gene list) is read from
the run's `<stem>_summary.json`; expression from the preprocessed AnnData
(`adata.X`, log1p, unless --layer overrides). Genes are matched to `var_names`
case-insensitively, mirroring the embedder.

Inputs:
    outs/embeddings/<dataset>/<stem>.parquet                 (cell_id order + axis activations)
    outs/embeddings/<dataset>/<stem>_summary.json            (axes + module genes)
    data/preprocessed/<dataset>_pp.h5ad                      (expression + labels)
    configs/axis_type_mappings/<dataset>/axis_types.csv      (curated axis->type, for boxing)

Outputs (file/dir tagged by --style: 'dotplot' default, or 'heatmap'):
    outs/plots/gene_<style>/<dataset>/<stem>_gene_<style>.{png,pdf}

Usage:
    python scripts/plots/plot_gene_dotplot.py --dataset mouse_cortex
    python scripts/plots/plot_gene_dotplot.py --dataset mouse_cortex --style heatmap \
        --cohort blueprint_gpt_5_5_unrefined_high --seed 2 --max-cells 4000
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]

# Default representative run per dataset: the literature-grounded Stage-A cohort and the
# composite-best seed used in the model comparison (umap_grid / frontier).
DATASET_DEFAULTS = {
    "mouse_cortex":   ("blueprint_gpt_5_5_unrefined_high", 2, "celltype_harmonized"),
    "human_pancreas": ("blueprint_gpt_5_5_unrefined_high", 1, "celltype"),
    "human_immune":   ("blueprint_gpt_5_5_unrefined_high", 1, "final_annotation"),
}

# Curated biological cell-type grouping for the dotplot rows, the source of truth
# for both the row order and the left-hand class band. Each dataset maps to an
# ordered list of (broad class, [cell types in within-class order]): neural first
# -- glutamatergic by cortical layer (superficial -> deep), then GABAergic by
# developmental origin (CGE -> MGE -> other) -- then non-neural by compartment
# (macroglia, immune, vascular/mural). Types absent from a dataset's list fall back
# to the palette order; nothing is dropped.
CELLTYPE_CLASSES = {
    "mouse_cortex": [
        ("Glutamatergic", ["CR", "L2/3 IT", "L4", "L5 IT", "L5 PT", "L5 NP",
                            "L6 IT", "L6 CT", "L6b"]),
        ("GABAergic",     ["Lamp5", "Sncg", "Vip", "Sst", "Sst Chodl", "Pvalb", "Meis2"]),
        ("Non-neural",    ["Astro", "OPC", "Oligo", "Microglia", "Macrophage",
                           "Endo", "Peri", "SMC", "VLMC"]),
    ],
}

# Colors for the broad-class band -- a hue family kept distinct from the axis-kind
# key (blue/orange/green) and the per-type palette so it reads as annotation.
CLASS_COLORS = {
    "Glutamatergic": "#756bb1",  # purple
    "GABAergic":     "#2ca7a0",  # teal
    "Non-neural":    "#969696",  # grey
}

# Cell outline boxes: red = curated named target, blue = the cell type this axis
# tops by AUROC, orange = both agree (named target is also the top-AUROC axis).
BOX_COLOR = {"named": "#d62728", "auroc": "#1f77b4", "both": "#ff7f0e"}


def _load(name: str):
    path = Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ----------------------------------------------------------------------
# Blueprint structure
# ----------------------------------------------------------------------
def axis_gene_blocks(summary: dict) -> list[tuple[str, list[str]]]:
    """Return [(axis_name, positive-module genes)] in dim order.

    Each axis's defining genes are its *positive* module's genes (the program
    the axis is named for). Genes already used by an earlier axis are dropped so
    every gene appears once, in the block of the first axis that defines it.
    """
    bp = summary["blueprint"]
    modules = {m["name"]: list(m["genes"]) for m in bp["modules"]}
    seen: set[str] = set()
    blocks: list[tuple[str, list[str]]] = []
    for axis in bp["axes"]:
        genes = []
        for g in modules.get(axis["positive"], []):
            if g not in seen:
                genes.append(g)
                seen.add(g)
        if genes:
            blocks.append((axis["name"], genes))
    return blocks


# ----------------------------------------------------------------------
# Expression
# ----------------------------------------------------------------------
def subsample_by_label(labels: pd.Series, max_cells: int, floor: int, rng) -> np.ndarray:
    """Stratified positional subsample: a floor per label so rare types stay
    visible, otherwise proportional. Returns positional indices into `labels`."""
    n = len(labels)
    if n <= max_cells:
        return np.arange(n)
    codes = labels.to_numpy()
    keep: list[np.ndarray] = []
    for val in pd.unique(codes):
        idx = np.flatnonzero(codes == val)
        take = min(len(idx), max(floor, round(max_cells * len(idx) / n)))
        keep.append(rng.choice(idx, size=take, replace=False))
    return np.sort(np.concatenate(keep))


def fetch_expression(dataset: str, genes: list[str], cell_ids: pd.Index,
                     layer: str) -> tuple[np.ndarray, list[str]]:
    """Return (cells x genes log-expression, genes_present) for `cell_ids`,
    matching gene symbols to var_names case-insensitively."""
    import anndata as ad
    import scipy.sparse as sp

    pp = REPO_ROOT / "data" / "preprocessed" / f"{dataset}_pp.h5ad"
    adata = ad.read_h5ad(pp, backed="r")
    lc_to_idx = {g.lower(): i for i, g in enumerate(adata.var_names.astype(str))}
    gidx, present = [], []
    for g in genes:
        j = lc_to_idx.get(g.lower())
        if j is not None:
            gidx.append(j)
            present.append(g)
    obs_names = pd.Index(adata.obs_names.astype(str))
    pos = obs_names.get_indexer(cell_ids)
    if (pos < 0).any():
        raise ValueError("Some embedding cell_ids are absent from the AnnData.")
    sub = adata[pos, :].to_memory()
    X = sub.layers[layer] if layer != "X" else sub.X
    X = X[:, gidx]
    if sp.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=float), present


def zscore_clip(X: np.ndarray, clip: float) -> np.ndarray:
    """Column-wise (per-gene) z-score with constant columns mapped to 0."""
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    sd[sd == 0] = 1.0
    return np.clip((X - mu) / sd, -clip, clip)


def aggregate_by_type(X: np.ndarray, labels: np.ndarray, types: list[str],
                      clip: float) -> tuple[np.ndarray, np.ndarray]:
    """Per (gene, cell-type) summaries over ALL cells, for the dotplot.

    Returns (mean_z, frac) each [n_genes x n_types]:
      mean_z[g,t] = mean log1p expression of gene g over cells of type t,
                    then z-scored per gene across the n_types group means and
                    clipped to +-clip (scanpy `standard_scale='var'` convention);
      frac[g,t]   = fraction of type-t cells with detectable expression (>0).
    """
    n_genes = X.shape[1]
    mean = np.zeros((n_genes, len(types)))
    frac = np.zeros((n_genes, len(types)))
    for j, t in enumerate(types):
        m = labels == t
        if not m.any():
            continue
        Xt = X[m]
        mean[:, j] = Xt.mean(axis=0)
        frac[:, j] = (Xt > 0).mean(axis=0)
    mu = mean.mean(axis=1, keepdims=True)
    sd = mean.std(axis=1, keepdims=True)
    sd[sd == 0] = 1.0
    mean_z = np.clip((mean - mu) / sd, -clip, clip)
    return mean_z, frac


# ----------------------------------------------------------------------
# Plot
# ----------------------------------------------------------------------
def plot(Z: np.ndarray, gene_names: list[str], gene_axis: list[str],
         cell_labels: pd.Series, color_map: dict, ordered_types: list[str],
         clip: float, title: str, out_stem: Path) -> None:
    n_cells, n_genes = Z.shape
    # cells (rows of Z) are already ordered by cell type; gene_names/columns by axis.
    type_codes = cell_labels.to_numpy()
    # boundaries between cell-type blocks (columns of the heatmap = cells)
    col_bounds = np.flatnonzero(type_codes[1:] != type_codes[:-1]) + 1
    # boundaries between axis blocks (rows of the heatmap = genes)
    row_bounds = np.flatnonzero(np.array(gene_axis[1:]) != np.array(gene_axis[:-1])) + 1

    fig_h = max(4.0, 0.105 * n_genes + 1.4)
    fig = plt.figure(figsize=(11.0, fig_h))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.4, fig_h - 0.4],
                          width_ratios=[1, 0.012], hspace=0.02, wspace=0.02)
    strip_ax = fig.add_subplot(gs[0, 0])
    ax = fig.add_subplot(gs[1, 0], sharex=strip_ax)
    cax = fig.add_subplot(gs[1, 1])

    # Heatmap: genes (rows) x cells (cols).
    im = ax.imshow(Z.T, aspect="auto", cmap="RdBu_r", vmin=-clip, vmax=clip,
                   interpolation="nearest", rasterized=True)
    ax.set_yticks(range(n_genes))
    # Gene symbols are italicised per standard nomenclature.
    ax.set_yticklabels(gene_names, fontsize=4.5, style="italic")
    ax.set_xticks([])
    for b in col_bounds:
        ax.axvline(b - 0.5, color="white", lw=0.4)
    for b in row_bounds:
        ax.axhline(b - 0.5, color="black", lw=0.5)
    ax.set_xlabel(f"cells (n={n_cells:,}, grouped by cell type)", fontsize=8)

    # Axis-block brackets + labels down the RIGHT edge (kept clear of gene names).
    starts = np.concatenate([[0], row_bounds])
    ends = np.concatenate([row_bounds, [n_genes]])
    x_lab = n_cells - 0.5 + 0.012 * n_cells
    for s, e in zip(starts, ends):
        ax.plot([n_cells - 0.5 + 0.004 * n_cells] * 2, [s - 0.4, e - 0.6],
                color="0.3", lw=0.8, clip_on=False)
        ax.text(x_lab, (s + e - 1) / 2, gene_axis[s], fontsize=5.5,
                ha="left", va="center", clip_on=False)

    # Top color strip = cell type.
    strip = np.array([[*color_map[t]] for t in type_codes])[None, :, :]
    strip_ax.imshow(strip, aspect="auto", interpolation="nearest")
    strip_ax.set_yticks([])
    strip_ax.set_xticks([])
    for b in col_bounds:
        strip_ax.axvline(b - 0.5, color="white", lw=0.4)
    # Cell-type names centered over their block.
    cstarts = np.concatenate([[0], col_bounds])
    cends = np.concatenate([col_bounds, [n_cells]])
    for s, e in zip(cstarts, cends):
        strip_ax.text((s + e - 1) / 2, -0.6, type_codes[s], rotation=90,
                      fontsize=4.8, ha="center", va="bottom")
    strip_ax.set_title(title, fontsize=10, pad=14)

    fig.colorbar(im, cax=cax).set_label(f"z-scored log1p (clip ±{clip:g})", fontsize=7)
    cax.tick_params(labelsize=6)

    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_stem.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"[ok] wrote {out_stem}.{{png,pdf}}")


def plot_dotplot(mean_z: np.ndarray, frac: np.ndarray, gene_names: list[str],
                 gene_axis: list[str], axis_kind_of: dict, kind_color: dict,
                 types: list[str], type_class: list[str], class_colors: dict,
                 color_map: dict, claimed: dict, best_axis_per_type: dict, clip: float,
                 title: str, out_stem: Path) -> None:
    """Genes (rows, grouped by the axis that declares them) x cell types (columns)
    dot plot: dot color = mean z-scored expression, dot size = fraction of cells of
    that type expressing the gene (portrait orientation; cell types share the x-axis
    with the AUROC heatmap). A broad-class band annotates the columns on top; axis
    brackets label the row blocks on the right; legends run along the bottom."""
    n_genes, n_types = mean_z.shape
    row_bounds = np.flatnonzero(np.array(gene_axis[1:]) != np.array(gene_axis[:-1])) + 1
    starts = np.concatenate([[0], row_bounds])
    ends = np.concatenate([row_bounds, [n_genes]])
    type_ix = {t: i for i, t in enumerate(types)}

    # Typography hierarchy (pt): class band > type labels > axis labels > genes.
    # NOTE: FS_GENE and row_h drive the figure HEIGHT (155 gene rows); raising them
    # pushes the whole composite past a page and floats it away from its references.
    # The other sizes (ticks/axis/class/legend) are free to grow -- they do not.
    FS_GENE, FS_AXIS, FS_TYPE, FS_CLASS = 7.5, 10.5, 11.5, 12.5
    FS_LEG, FS_LEG_TITLE, FS_CBAR, FS_CBAR_TICK = 10.0, 10.5, 10.0, 8.5

    s_max = 32.0  # points^2 for a fully-expressed (frac=1) dot
    row_h = 0.122                 # inches per gene row (sized to the gene-label font)
    fig_w = max(7.0, 0.36 * n_types + 5.0)
    fig_h = max(6.5, row_h * n_genes + 4.4)
    fig = plt.figure(figsize=(fig_w, fig_h))
    # rows: class+type strip | dots | spacer (type tick labels) | legends
    gs = fig.add_gridspec(4, 2, height_ratios=[0.55, row_h * n_genes, 0.9, 1.5],
                          width_ratios=[1.0, 0.30], hspace=0.02, wspace=0.015)
    top_ax = fig.add_subplot(gs[0, 0])               # class band + per-type colour strip
    ax = fig.add_subplot(gs[1, 0])                    # main dots
    brk_ax = fig.add_subplot(gs[1, 1]); brk_ax.axis("off")   # axis brackets/labels
    leg_ax = fig.add_subplot(gs[3, :]); leg_ax.axis("off")   # legends (below the spacer row)

    # --- main dots: x = cell types, y = genes (row 0 at the top) ---
    xx, yy = np.meshgrid(np.arange(n_types), np.arange(n_genes))   # both (n_genes, n_types)
    sc = ax.scatter(xx.ravel(), yy.ravel(), s=(frac.ravel() * s_max),
                    c=mean_z.ravel(), cmap="RdBu_r", vmin=-clip, vmax=clip,
                    edgecolors="0.55", linewidths=0.2, rasterized=True)
    ax.set_xlim(-0.5, n_types - 0.5)
    ax.set_ylim(n_genes - 0.5, -0.5)
    ax.set_xticks(range(n_types)); ax.set_xticklabels(types, rotation=90, fontsize=FS_TYPE)
    # Gene symbols (rows) are italicised per standard nomenclature.
    ax.set_yticks(range(n_genes)); ax.set_yticklabels(gene_names, fontsize=FS_GENE, style="italic")
    ax.tick_params(length=0)
    for b in row_bounds:
        ax.axhline(b - 0.5, color="0.9", lw=0.3)
    for s in ax.spines.values():
        s.set_visible(False)

    # Per axis-block, outline cells: RED = the axis's curated named target; BLUE = a
    # cell type for which this axis is the top-AUROC axis; ORANGE = both.
    def _box(c, s, e, color):
        ax.add_patch(mpatches.Rectangle((c - 0.5, s - 0.5), 1, e - s, fill=False,
                                        edgecolor=color, lw=2.6, zorder=4))

    for s, e in zip(starts, ends):
        ax_name = gene_axis[s]
        cc = type_ix.get(claimed.get(ax_name, ""))                  # named-target column
        best_cols = {type_ix[t] for t in types                      # columns this axis tops
                     if best_axis_per_type.get(t) == ax_name}
        if cc is not None and cc in best_cols:                      # overlap -> orange
            _box(cc, s, e, BOX_COLOR["both"]); best_cols.discard(cc)
        elif cc is not None:                                        # named only -> red
            _box(cc, s, e, BOX_COLOR["named"])
        for cb_ in best_cols:                                       # top-AUROC only -> blue
            _box(cb_, s, e, BOX_COLOR["auroc"])

    # --- top: per-type colour strip (lower) + broad-class band with names (upper) ---
    top_ax.set_xlim(-0.5, n_types - 0.5); top_ax.set_ylim(0, 1)
    top_ax.set_xticks([]); top_ax.set_yticks([])
    for s in top_ax.spines.values():
        s.set_visible(False)
    strip = np.array([[*color_map[t]] for t in types])[None, :, :]
    top_ax.imshow(strip, aspect="auto", interpolation="nearest", extent=(-0.5, n_types - 0.5, 0.0, 0.42))
    i = 0
    while i < n_types:
        c = type_class[i]; j = i
        while j + 1 < n_types and type_class[j + 1] == c:
            j += 1
        if c:
            top_ax.add_patch(mpatches.Rectangle((i - 0.5, 0.55), j - i + 1, 0.45,
                             facecolor=class_colors.get(c, "0.7"), edgecolor="white", lw=0.8))
            top_ax.text((i + j) / 2, 0.775, c, ha="center", va="center",
                        fontsize=FS_CLASS, color="white", fontweight="bold")
        i = j + 1

    # --- right: axis brackets + horizontal axis names, coloured by kind ---
    brk_ax.set_ylim(n_genes - 0.5, -0.5); brk_ax.set_xlim(0, 1)
    for s, e in zip(starts, ends):
        ax_name = gene_axis[s]
        kc = kind_color.get(axis_kind_of.get(ax_name, ""), "0.3")
        brk_ax.plot([0.06, 0.06], [s - 0.4, e - 0.6], color=kc, lw=1.2)
        brk_ax.text(0.16, (s + e - 1) / 2, ax_name, fontsize=FS_AXIS,
                    ha="left", va="center", color=kc)

    # --- bottom: colorbar + dot-size + axis-kind + cell-type-box legends ---
    cax = leg_ax.inset_axes([0.02, 0.50, 0.15, 0.14])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_label("mean log1p, z-scored per gene", fontsize=FS_CBAR)
    cb.ax.tick_params(labelsize=FS_CBAR_TICK, length=2)
    cb.outline.set_linewidth(0.5)
    size_handles = [plt.scatter([], [], s=f * s_max, color="0.45", edgecolors="0.55",
                                linewidths=0.2, label=f"{int(f * 100)}%")
                    for f in (0.1, 0.25, 0.5, 0.75, 1.0)]
    leg1 = leg_ax.legend(handles=size_handles, title="dot size = % cells expressing",
                         loc="upper left", bbox_to_anchor=(0.27, 0.98), frameon=False,
                         ncol=5, columnspacing=1.0, fontsize=FS_LEG,
                         title_fontsize=FS_LEG_TITLE, handletextpad=0.3)
    leg_ax.add_artist(leg1)
    kind_handles = [plt.Line2D([0], [0], marker="s", ls="", mfc=c, mec="none",
                               markersize=7, label=k) for k, c in kind_color.items()]
    leg2 = leg_ax.legend(handles=kind_handles, title="axis kind", loc="upper left",
                         bbox_to_anchor=(0.27, 0.45), frameon=False, ncol=3,
                         fontsize=FS_LEG, title_fontsize=FS_LEG_TITLE, handletextpad=0.3)
    leg_ax.add_artist(leg2)
    box_handles = [mpatches.Patch(facecolor="none", edgecolor=BOX_COLOR["named"],
                                  lw=2.6, label="named target"),
                   mpatches.Patch(facecolor="none", edgecolor=BOX_COLOR["auroc"],
                                  lw=2.6, label="top-AUROC type"),
                   mpatches.Patch(facecolor="none", edgecolor=BOX_COLOR["both"],
                                  lw=2.6, label="named = top-AUROC")]
    leg_ax.legend(handles=box_handles, title="cell-type box", loc="upper left",
                  bbox_to_anchor=(0.62, 0.98), frameon=False, ncol=3,
                  fontsize=FS_LEG, title_fontsize=FS_LEG_TITLE, handletextpad=0.4)

    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_stem.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"[ok] wrote {out_stem}.{{png,pdf}}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Marker-gene grounding figure for one blueprint run.")
    p.add_argument("--dataset", choices=tuple(DATASET_DEFAULTS), default="mouse_cortex")
    p.add_argument("--style", choices=("dotplot", "heatmap"), default="dotplot",
                   help="dotplot (aggregated, default) or heatmap (single-cell).")
    p.add_argument("--cohort", default=None, help="Blueprint cohort (model) stem prefix.")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--label-col", default=None)
    p.add_argument("--layer", default="X", help="AnnData layer ('X' = adata.X, log1p).")
    p.add_argument("--max-cells", type=int, default=4000)
    p.add_argument("--floor", type=int, default=25, help="Min cells kept per cell type.")
    p.add_argument("--clip", type=float, default=2.5)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    up = _load("umap_plots")
    matplotlib.rcParams.update(up.PLOT_RC)

    cohort, seed, label_col = DATASET_DEFAULTS[args.dataset]
    cohort = args.cohort or cohort
    seed = args.seed if args.seed is not None else seed
    label_col = args.label_col or label_col

    stem = f"{cohort}_{args.dataset}_seed{seed}"
    emb_dir = REPO_ROOT / "outs" / "embeddings" / args.dataset
    summary = json.loads((emb_dir / f"{stem}_summary.json").read_text())
    parquet = pd.read_parquet(emb_dir / f"{stem}.parquet")
    cell_ids_all = pd.Index(parquet["cell_id"].astype(str))

    kind = "dotplot" if args.style == "dotplot" else "heatmap"   # tag outputs by style
    out_stem = REPO_ROOT / "outs" / "plots" / f"gene_{kind}" / args.dataset / f"{stem}_gene_{kind}"
    if out_stem.with_suffix(".png").exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite {out_stem}.png (pass --overwrite).")

    # Labels (aligned to the embedding's cells).
    labels = up.load_obs_columns(
        REPO_ROOT / "data" / "preprocessed" / f"{args.dataset}_pp.h5ad", [label_col]
    )[label_col]
    cell_labels_all = cell_ids_all.map(labels).astype(str)
    keep = cell_labels_all.to_numpy() != "nan"
    cell_ids_all, cell_labels_all = cell_ids_all[keep], cell_labels_all[keep]

    blocks = axis_gene_blocks(summary)
    kind_of = {a["name"]: a["kind"] for a in summary["blueprint"]["axes"]}

    if args.style == "heatmap":
        # --- single-cell heatmap: subsample + order cells by type ---
        rng = np.random.default_rng(0)
        pos = subsample_by_label(pd.Series(cell_labels_all), args.max_cells, args.floor, rng)
        cell_ids = cell_ids_all[pos]
        cell_labels = pd.Series(cell_labels_all[pos])
        order = np.argsort(cell_labels.to_numpy(), kind="mergesort")
        cell_ids, cell_labels = cell_ids[order], cell_labels.iloc[order].reset_index(drop=True)

        X, present = fetch_expression(args.dataset, [g for _, gs in blocks for g in gs],
                                      cell_ids, args.layer)
        present_set = set(present)
        gene_names = [g for _, gs in blocks for g in gs if g in present_set]
        gene_axis = [ax for ax, gs in blocks for g in gs if g in present_set]
        Z = zscore_clip(X, args.clip)

        color_map, _, ordered = up.build_categorical_colors(
            cell_labels, dataset=args.dataset, column=label_col, kind=None)
        color_map = {k: matplotlib.colors.to_rgb(v) for k, v in color_map.items()}
        title = (f"{args.dataset} · blueprint seed{seed} · "
                 f"{len(gene_names)} cited marker genes / {len(blocks)} axes")
        plot(Z, gene_names, gene_axis, cell_labels, color_map, ordered, args.clip, title, out_stem)

    else:
        # --- aggregated dotplot (default): all cells, AUROC-sorted ---
        pf = _load("plot_axis_faithfulness")
        M, auroc_names, _kinds, type_counts = pf.axis_type_auroc(
            args.dataset, stem, label_col, up)
        claimed = pf.load_axis_type_map(args.dataset)   # curated axis -> named type
        auroc_types = list(type_counts.index)
        best_auroc = {}                                 # per axis: its peak AUROC (sort tiebreak)
        for i, name in enumerate(auroc_names):
            row = M[i]
            if not np.all(np.isnan(row)):
                best_auroc[name] = float(np.nanmax(row))
        best_axis_per_type = {}                         # per type: the top-AUROC axis (blue box)
        for j, t in enumerate(auroc_types):
            col = M[:, j]
            if not np.all(np.isnan(col)):
                best_axis_per_type[t] = auroc_names[int(np.nanargmax(col))]

        # Cell-type row order: curated biological grouping where available, else the
        # dataset palette order; any remaining present types appended (nothing dropped).
        present_set = set(pd.unique(cell_labels_all))
        color_map, _, ordered = up.build_categorical_colors(
            pd.Series(cell_labels_all), dataset=args.dataset, column=label_col, kind=None)
        classes_def = CELLTYPE_CLASSES.get(args.dataset, [])
        type_to_class = {t: cname for cname, ts in classes_def for t in ts}
        curated_order = [t for cname, ts in classes_def for t in ts]
        types = [t for t in curated_order if t in present_set]
        types += [t for t in ordered if t in present_set and t not in set(types)]
        types += [t for t in present_set if t not in set(types)]
        type_class = [type_to_class.get(t, "") for t in types]
        color_map = {t: matplotlib.colors.to_rgb(color_map[t]) for t in types}
        col_ix = {t: i for i, t in enumerate(types)}

        # Axis blocks ordered by each axis's curated named-target type (identical to
        # the faithfulness matrix); un-named program/contrast axes fall to the end,
        # ties broken by AUROC strength.
        blocks_sorted = sorted(
            blocks, key=lambda b: (col_ix.get(claimed.get(b[0], ""), len(types)),
                                   -best_auroc.get(b[0], 0.0)))

        X, present = fetch_expression(
            args.dataset, [g for _, gs in blocks_sorted for g in gs], cell_ids_all, args.layer)
        present_set = set(present)
        gene_names = [g for _, gs in blocks_sorted for g in gs if g in present_set]
        gene_axis = [ax for ax, gs in blocks_sorted for g in gs if g in present_set]

        mean_z, frac = aggregate_by_type(X, cell_labels_all.to_numpy(), types, args.clip)
        title = (f"{args.dataset} · blueprint seed{seed} · "
                 f"{len(gene_names)} cited marker genes / {len(blocks)} axes")
        plot_dotplot(mean_z, frac, gene_names, gene_axis, kind_of, pf.KIND_COLOR,
                     types, type_class, CLASS_COLORS, color_map, claimed,
                     best_axis_per_type, args.clip, title, out_stem)



if __name__ == "__main__":
    main()
