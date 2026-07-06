"""
UMAP grid for the held-out reference-to-query transfer experiment.

For one representative seed, computes a 2D UMAP of each method-variant's joint
ref+query embedding and renders a grid:

  row 1 : colored by cell type (does structure hold across the joint set?)
  row 2 : reference cells in grey, the held-out query cells colored by cell type
          and drawn on top -- the visual test of transfer (do unseen query cells
          land on the matching reference structure without retraining?).

UMAP coords are cached next to the figure so re-plots are instant. Reuses the
palettes / obs loading / UMAP helper from scripts/plots/umap_plots.py so colors
match the other UMAP figures.

Inputs:
    outs/embeddings/<dataset>/transfer/<variant>_<dataset>_joint_seed<N>.parquet
    data/preprocessed/<dataset>_pp.h5ad                       (obs labels)

Outputs:
    outs/plots/transfer/<dataset>/umap_coords/<variant>_seed<N>.parquet  (cache)
    outs/plots/transfer/<dataset>/transfer_umap.{svg,png}

Usage:
    python scripts/plots/plot_transfer_umap.py --dataset mouse_cortex --seed 1
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.parquet as pqt

REPO_ROOT = Path(__file__).resolve().parents[2]

VARIANT_ORDER = ["blueprint", "blueprint_noref", "scvi_frozen", "scvi_scarches",
                 "expimap_frozen", "expimap_scarches"]
VARIANT_LABEL = {"blueprint": "blueprint\n(literature-grounded)", "blueprint_noref": "blueprint (noref)",
                 "scvi_frozen": "scVI frozen", "scvi_scarches": "scVI scArches",
                 "expimap_frozen": "expiMap frozen", "expimap_scarches": "expiMap scArches"}
# label col + held-out study/tech per dataset (mirrors scripts/transfer/measure.py)
DATASET_CFG = {
    "mouse_cortex":   {"label_col": "celltype_harmonized", "query_key": "study", "query_values": ["tasic2016"]},
    "human_pancreas": {"label_col": "celltype", "query_key": "tech", "query_values": ["smartseq2"]},
    "human_immune":   {"label_col": "final_annotation", "query_key": "batch",
                       "query_values": ["Oetjen_U", "Sun_sample4_TC"]},
}
# Match umap_plots.py (the full-dataset UMAP-grid source) so recomputed panels are
# stylistically consistent with the reused blueprint coords.
UMAP_PARAMS = {"n_neighbors": 30, "min_dist": 0.15, "metric": "euclidean", "random_state": 42}
REF_GREY = (0.82, 0.82, 0.82)
# Blueprint arms reuse the exact cached full-dataset UMAP coords from umap_plots.py
# (their transfer-joint embedding == the full-dataset embedding, byte-identical).
BLUEPRINT_COHORT = {"blueprint": "blueprint_gpt_5_5_unrefined_high",
                    "blueprint_noref": "blueprint_gpt_5_5_noref_unrefined_high"}


def _load_umap_plots():
    path = REPO_ROOT / "scripts" / "plots" / "umap_plots.py"
    sys.path.insert(0, str(path.parent))
    try:
        spec = importlib.util.spec_from_file_location("_transfer_umap_plots", path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules["_transfer_umap_plots"] = mod
        spec.loader.exec_module(mod)
        return mod
    finally:
        if str(path.parent) in sys.path:
            sys.path.remove(str(path.parent))


def main() -> None:
    ap = argparse.ArgumentParser(description="UMAP grid for the held-out transfer variants.")
    ap.add_argument("--dataset", choices=list(DATASET_CFG), required=True)
    ap.add_argument("--seed", type=int, default=None,
                    help="Force a single shared seed for all variants. Default: pick each "
                         "variant's best seed by --metric (matches umap_grid.py).")
    ap.add_argument("--metric", default="composite_scib",
                    help="Per-variant best-seed selection metric from transfer_metrics.csv "
                         "(default composite_scib, matching the main frontier UMAP grid).")
    ap.add_argument("--overwrite", action="store_true", help="Recompute cached UMAP coords.")
    args = ap.parse_args()

    up = _load_umap_plots()
    cfg = DATASET_CFG[args.dataset]
    tdir = REPO_ROOT / "outs" / "embeddings" / args.dataset / "transfer"
    coord_dir = REPO_ROOT / "outs" / "plots" / "transfer" / args.dataset / "umap_coords"

    obs = up.load_obs_columns(REPO_ROOT / "data" / "preprocessed" / f"{args.dataset}_pp.h5ad",
                              [cfg["label_col"], cfg["query_key"]])
    label_series = obs[cfg["label_col"]]
    qv = set(map(str, cfg["query_values"]))
    is_query_by_id = obs[cfg["query_key"]].map(lambda v: str(v) in qv)

    # Resolve per-variant (seed, cached-coords-path). The blueprint arms REUSE the exact
    # cached full-dataset coords from the umap-grid (same byte-identical embedding, same
    # layout), at the seed those figures use (main metrics.csv best composite_scib).
    # scVI/expiMap have transfer-specific joint embeddings, so they are
    # recomputed on the joint at the model's transfer-best seed, shared across the
    # frozen/scArches columns so the surgery is the only difference between them.
    tdf = pd.read_csv(REPO_ROOT / "outs" / "metrics" / args.dataset / "transfer_metrics.csv")
    tdf = tdf[tdf.dataset == args.dataset]
    main_csv = REPO_ROOT / "outs" / "metrics" / args.dataset / "metrics.csv"
    mdf = pd.read_csv(main_csv) if main_csv.exists() else None
    umap_dir = REPO_ROOT / "outs" / "plots" / "umap" / args.dataset

    def _argmax_seed(df, col, val):
        if df is None:
            return None
        g = df[df[col] == val]
        return int(g.loc[g[args.metric].idxmax(), "seed"]) if len(g) else None

    def resolve(v):  # -> (seed, cached_coords_path|None)
        if args.seed is not None:
            return args.seed, None
        if v in BLUEPRINT_COHORT:                 # reuse the full-dataset umap-grid layout
            cohort = BLUEPRINT_COHORT[v]
            seed = _argmax_seed(mdf, "model", cohort) or _argmax_seed(tdf, "variant", v)
            cached = umap_dir / cohort / f"{cohort}_{args.dataset}_seed{seed}_umap_coords.parquet"
            return seed, (cached if (seed is not None and cached.exists()) else None)
        model = v.split("_")[0]                   # share seed across frozen/scArches
        return (_argmax_seed(tdf, "variant", f"{model}_scarches")
                or _argmax_seed(tdf, "variant", v)), None

    plan = {v: resolve(v) for v in VARIANT_ORDER}
    variants = [v for v in VARIANT_ORDER if plan[v][0] is not None and (
        plan[v][1] is not None or (tdir / f"{v}_{args.dataset}_joint_seed{plan[v][0]}.parquet").exists())]
    if not variants:
        raise FileNotFoundError(f"No coords/joint parquets found for {args.dataset}")
    print(f"[ok] {args.dataset} seeds -> " + ", ".join(
        f"{v}=seed{plan[v][0]}{' (reused coords)' if plan[v][1] else ''}" for v in variants))

    # Per-variant panel width kept near the umap_grid panel size so the figure is not
    # downscaled as hard when laid out full-width and its text renders at a comparable
    # size; height holds two square-ish UMAP rows plus the title and legend bands.
    fig, axes = plt.subplots(2, len(variants), figsize=(2.8 * len(variants), 6.5), squeeze=False)
    color_map = None
    for j, v in enumerate(variants):
        sd, cached_path = plan[v]
        if cached_path is not None:               # reuse published full-dataset coords
            cdf = pd.read_parquet(cached_path)
            cell_ids = cdf["cell_id"].astype(str).reset_index(drop=True)
            coords = cdf[["umap_0", "umap_1"]].to_numpy()
            full_pq = (REPO_ROOT / "outs" / "embeddings" / args.dataset
                       / f"{BLUEPRINT_COHORT[v]}_{args.dataset}_seed{sd}.parquet")
            ndim = (sum(c.startswith("dim_") for c in pqt.ParquetFile(full_pq).schema.names)
                    if full_pq.exists() else 0)
        else:                                     # recompute on the transfer joint
            cell_ids, coords, ndim, _ = up.compute_or_load_umap(
                tdir / f"{v}_{args.dataset}_joint_seed{sd}.parquet",
                coord_dir / f"{v}_seed{sd}.parquet", UMAP_PARAMS, args.overwrite)
        labels = pd.Series(cell_ids.map(label_series).fillna("unmatched").to_numpy(), name=cfg["label_col"])
        cmap, point_colors, ordered = up.build_categorical_colors(
            labels, dataset=args.dataset, column=cfg["label_col"], kind=None)
        if color_map is None:
            color_map, ordered_labels = cmap, [l for l in ordered if l != "unmatched"]
        is_q = cell_ids.map(is_query_by_id).fillna(False).to_numpy()

        # row 0: all cells by cell type
        ax0 = axes[0][j]
        ax0.scatter(coords[:, 0], coords[:, 1], s=3.5, alpha=0.6, linewidths=0,
                    c=point_colors, rasterized=True)
        ax0.set_title(f"{VARIANT_LABEL[v]}\nseed {sd} · {ndim}d", fontsize=11)
        # row 1: reference grey, query on top colored by cell type
        ax1 = axes[1][j]
        ax1.scatter(coords[~is_q, 0], coords[~is_q, 1], s=2.5, alpha=0.35, linewidths=0,
                    c=[REF_GREY], rasterized=True)
        ax1.scatter(coords[is_q, 0], coords[is_q, 1], s=7.0, alpha=0.9, linewidths=0,
                    c=[point_colors[i] for i in np.where(is_q)[0]], rasterized=True)
        for ax in (ax0, ax1):
            ax.set_xticks([]); ax.set_yticks([]); ax.set_aspect("equal", adjustable="datalim")
            for s in ax.spines.values():
                s.set_visible(False)
    axes[0][0].set_ylabel("by cell type", fontsize=12)
    axes[1][0].set_ylabel(f"reference (gray) vs\nquery ({','.join(cfg['query_values'])})", fontsize=12)

    handles = up.build_legend_handles(color_map, ordered_labels)
    ncol = min(7, len(handles))
    fig.legend(handles=handles, loc="lower center", ncol=ncol, fontsize=9,
               frameon=False, markerscale=2.2, columnspacing=1.2, handletextpad=0.4,
               bbox_to_anchor=(0.5, -0.01))
    # Block title is just the dataset name; the seed criterion and coordinate-reuse note
    # live in the caption, and the held-out query is already labeled on the bottom-row axis.
    fig.suptitle(args.dataset.replace('_', ' '), y=1.0, fontsize=12)
    # reserve bottom space for the (multi-row) cell-type legend
    n_legend_rows = int(np.ceil(len(handles) / ncol))
    fig.tight_layout(rect=(0, 0.03 + 0.022 * n_legend_rows, 1, 0.99))

    out_dir = REPO_ROOT / "outs" / "plots" / "transfer" / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(out_dir / f"transfer_umap.{ext}", dpi=200, bbox_inches="tight", facecolor="white")
    print(f"[ok] wrote {out_dir / 'transfer_umap.svg'} and .png ({len(variants)} variants)")


if __name__ == "__main__":
    main()
