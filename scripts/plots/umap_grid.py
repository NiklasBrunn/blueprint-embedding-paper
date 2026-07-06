"""
UMAP grid: 2 rows x N models, shared per-row legend.

Row 1 colors cells by the dataset's cell-type label (biology), row 2 by the batch
the embedding should mix (mouse_cortex: study; human_pancreas: tech; human_immune:
batch). One model per column, each at its best seed (highest composite_scib at the
cell-type label, from metrics.csv). The panel makes "batch-robust by construction"
visual: blueprint mixes batches (no batch channel) while PCA / foundation models
separate by batch.

Re-plots from the cached UMAP coordinates written by umap_plots.py
(<stem>_umap_coords.parquet) — no UMAP is recomputed — reusing that module's
categorical palettes and obs loading so colors match the standalone UMAP figures.

Inputs:
    data/preprocessed/mouse_cortex_pp.h5ad                       (obs columns)
    outs/plots/umap/mouse_cortex/<model>/<stem>_umap_coords.parquet
    outs/metrics/mouse_cortex/metrics.csv                        (best-seed pick)

Outputs:
    outs/plots/umap/mouse_cortex/_grid.{png,pdf}

Usage:
    # literature-grounded blueprint + baselines, each at its composite-best seed:
    python scripts/plots/umap_grid.py --dataset mouse_cortex
    python scripts/plots/umap_grid.py --dataset human_pancreas
    python scripts/plots/umap_grid.py --dataset human_immune
    # optional per-model seed override (e.g. force the blueprint to a chosen seed):
    python scripts/plots/umap_grid.py --seeds blueprint_constr=7
"""
from __future__ import annotations

import argparse
import importlib.util
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]

# Per-dataset rows: (obs column, palette kind, row label). Row 0 = biology
# (cell type), row 1 = the batch the embedding should mix.
DATASET_ROWS = {
    "mouse_cortex":   [("celltype_harmonized", None, "cell type"),
                       ("study", "batch", "study (batch)")],
    "human_pancreas": [("celltype", None, "cell type"),
                       ("tech", "batch", "tech (batch)")],
    "human_immune":   [("final_annotation", None, "cell type"),
                       ("batch", "batch", "batch")],
}

# Per-dataset legend relabeling, keyed by (dataset, column): append the sequencing
# technology to each study so the batch legend reads e.g. "tasic2016 (SMART-Seq)".
# Only the displayed label changes; colors stay keyed by the raw obs value.
LEGEND_RELABEL = {
    "mouse_cortex": {
        "study": {
            "tasic2016":  "tasic2016 (SMART-Seq)",
            "tasic2018":  "tasic2018 (SMART-Seq v4)",
            "zeisel2018": "zeisel2018 (10x v2)",
        },
    },
}


def _load(name: str):
    path = Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


up = _load("umap_plots")
pm = _load("plot_metrics")

# Friendly CLI aliases -> the `model` value used in metrics.csv / umap dir names.
ALIASES = {
    "blueprint_noref":  "blueprint_gpt_5_5_noref_unrefined_high",
    "blueprint_constr": "blueprint_gpt_5_5_unrefined_high",
}
# Grounded ('with references') blueprint, to match the frontier figure. Mixers
# (blueprint, batch-keyed scVI, expiMap) then prior-free separators (UCE, PCA);
# scVI-census is dropped here for legibility (it lives in the frontier + tables).
DEFAULT_MODELS = ["blueprint_constr", "scvi", "expimap", "uce", "pca"]


def resolve_model(name: str) -> str:
    return ALIASES.get(name, name)


def col_title(model: str) -> str:
    """Compact column title (shorter than prettify_model for blueprint variants)."""
    if model.startswith("blueprint_"):
        return "Blueprint"
    return pm.prettify_model(model)


def best_seed(metrics: pd.DataFrame, model: str, label_col: str) -> float | None:
    """Seed of the highest-composite run at the cell-type label; None if seedless."""
    sub = metrics[(metrics["model"] == model)
                  & (metrics["label_col"] == label_col)]
    if sub.empty:
        raise ValueError(f"{model}: no metrics rows at {label_col}.")
    row = sub.loc[sub["composite_scib"].idxmax()]
    seed = row["seed"]
    return None if pd.isna(seed) else int(seed)


def coords_path(dataset: str, model: str, seed: int | None) -> Path:
    """Locate the cached coords parquet (best seed, else seedless, else any)."""
    d = REPO_ROOT / "outs" / "plots" / "umap" / dataset / model
    if seed is not None:
        p = d / f"{model}_{dataset}_seed{seed}_umap_coords.parquet"
        if p.exists():
            return p
    p = d / f"{model}_{dataset}_umap_coords.parquet"
    if p.exists():
        return p
    found = sorted(d.glob("*_umap_coords.parquet"))
    if not found:
        raise FileNotFoundError(f"No coords parquet under {d}.")
    print(f"[warn] {model}: best-seed coords missing, using {found[0].name}")
    return found[0]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="UMAP grid (cell type / batch) per dataset.")
    p.add_argument("--dataset", choices=tuple(DATASET_ROWS), default="mouse_cortex")
    p.add_argument("--models", nargs="+", default=DEFAULT_MODELS,
                   help="Models (friendly aliases ok). Best seed chosen automatically.")
    p.add_argument("--seeds", nargs="*", default=None,
                   help="Per-model seed overrides as alias=seed (e.g. "
                        "blueprint_constr=7); otherwise the composite-best seed is used.")
    p.add_argument("--rows", type=str, default=None,
                   help="Comma-separated obs columns to use as grid rows, overriding "
                        "the dataset default (e.g. 'tissue'). The dataset's batch key "
                        "is colored as a batch axis; other columns as biology.")
    p.add_argument("--suffix", type=str, default=None,
                   help="Output filename suffix (e.g. 'tissue' -> umap_grid_<ds>_tissue).")
    p.add_argument("--point-size", type=float, default=1.8)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    matplotlib.rcParams.update(up.PLOT_RC)
    dataset = args.dataset
    rows = DATASET_ROWS[dataset]
    celltype_col = rows[0][0]  # always pick the seed by the cell-type composite
    if args.rows:
        batch_default = {r[0] for r in DATASET_ROWS[dataset] if r[1] == "batch"}
        rows = [(c.strip(), "batch" if c.strip() in batch_default else None, c.strip())
                for c in args.rows.split(",") if c.strip()]

    grid_name = "_grid" + (f"_{args.suffix}" if args.suffix else "")
    out_stem = REPO_ROOT / "outs" / "plots" / "umap" / dataset / grid_name
    if out_stem.with_suffix(".png").exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite {out_stem}.png (pass --overwrite).")

    raw = pd.read_csv(REPO_ROOT / "outs" / "metrics" / dataset / "metrics.csv")

    overrides = {}
    for item in (args.seeds or []):
        alias, _, val = item.partition("=")
        overrides[resolve_model(alias.strip())] = int(val)

    models = [resolve_model(m) for m in args.models]
    seeds = {m: (overrides[m] if m in overrides else best_seed(raw, m, celltype_col))
             for m in models}

    obs = up.load_obs_columns(
        REPO_ROOT / "data" / "preprocessed" / f"{dataset}_pp.h5ad",
        [col for col, _, _ in rows],
    )

    # Pre-load coords once per model (shared across the two rows).
    coords_by_model = {}
    for m in models:
        cdf = pd.read_parquet(coords_path(dataset, m, seeds[m]))
        coords_by_model[m] = (cdf["cell_id"].astype(str),
                              cdf[["umap_0", "umap_1"]].to_numpy())

    n = len(models)
    nr = len(rows)
    # gridspec: per data row a panel row + a legend row beneath it (so each row
    # keeps its own shared legend).
    panel = 2.0
    fig = plt.figure(figsize=(panel * n + 0.6, (panel + 1.1) * nr + 0.2))
    gs = fig.add_gridspec(2 * nr, n, height_ratios=[panel, 0.9] * nr,
                          hspace=0.18, wspace=0.06)

    for ri, (col, kind, row_label) in enumerate(rows):
        panel_r = ri * 2
        legend_r = ri * 2 + 1
        labels = obs[col]
        legend_handles = None
        for ci, m in enumerate(models):
            ax = fig.add_subplot(gs[panel_r, ci])
            cell_ids, coords = coords_by_model[m]
            label_values = pd.Series(
                cell_ids.map(labels).fillna("unmatched").astype(str).to_numpy(),
                index=cell_ids.index, name=col,
            )
            color_map, point_colors, ordered = up.build_categorical_colors(
                label_values, dataset=dataset, column=col, kind=kind)
            # Largest groups first so minority clusters sit on top.
            counts = label_values.value_counts()
            prio = np.array([counts[v] for v in label_values], dtype=int)
            order = np.argsort(-prio, kind="mergesort")
            ax.scatter(coords[order, 0], coords[order, 1], s=args.point_size,
                       alpha=0.8, linewidths=0,
                       c=[point_colors[i] for i in order], rasterized=True)
            ax.set_xticks([]); ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            ax.set_aspect("equal", adjustable="datalim")
            ax.margins(0.03)
            if ri == 0:
                ax.set_title(col_title(m), fontsize=9)
            if ci == 0:
                ax.set_ylabel(row_label, fontsize=9.5, labelpad=6)
            if legend_handles is None:
                leg_labels = [l for l in ordered if l != "unmatched"]
                leg_cmap = {k: v for k, v in color_map.items() if k != "unmatched"}
                legend_handles = up.build_legend_handles(leg_cmap, leg_labels)
                relabel = LEGEND_RELABEL.get(dataset, {}).get(col)
                if relabel:
                    for h in legend_handles:
                        h.set_label(relabel.get(h.get_label(), h.get_label()))

        # One shared legend for the row, spanning all columns.
        legend_ax = fig.add_subplot(gs[legend_r, :])
        legend_ax.axis("off")
        ncol = 2 if len(legend_handles) <= 8 else 4 if len(legend_handles) <= 24 else 6
        # No legend title: the row's identity (cell type / study) is already given by
        # the per-row y-axis label on the leftmost panel.
        legend_ax.legend(handles=legend_handles, loc="center", frameon=False,
                         ncol=ncol, fontsize=7.0,
                         handletextpad=0.4, columnspacing=1.0, labelspacing=0.4)

    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_stem.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"[ok] wrote {out_stem}.{{png,pdf}}  (best seeds: "
          + ", ".join(f"{m}={seeds[m]}" for m in models) + ")")



if __name__ == "__main__":
    main()
