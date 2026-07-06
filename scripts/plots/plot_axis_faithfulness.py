"""
Per-axis faithfulness figures for ONE blueprint run: an axis x cell-type AUROC
matrix plus a minimal across-seed reproducibility panel.

  (a) FULL axis x cell-type AUROC matrix -- every embedding axis scored
      (one-vs-rest AUROC of its activation) against every annotated cell type,
      computed from the run's parquet + labels. Each axis's *claimed* type
      (the label its name resolves to, read from
      configs/axis_type_mappings/<dataset>/axis_types.csv) is boxed; a faithful
      axis is brightest in its boxed cell and dim elsewhere.

and, as a separate figure, the minimal cherry-pick guard:

  across-seed AUROC distribution -- matched-axis AUROCs for every seed of the
      cohort, so the representative seed is shown to be typical, not lucky.

Inputs:
    outs/embeddings/<dataset>/<stem>.parquet          (axis activations)
    outs/embeddings/<dataset>/<stem>_summary.json     (axis names + kinds, dim order)
    configs/axis_type_mappings/<dataset>/axis_types.csv  (claimed type per axis)
    data/preprocessed/<dataset>_pp.h5ad               (labels)

Outputs:
    outs/plots/axis_metrics/<dataset>/axis_faithfulness_<stem>.{png,pdf}
    outs/plots/axis_metrics/<dataset>/axis_auroc_seeds_<cohort>.{png,pdf}

Usage:
    python scripts/plots/plot_axis_faithfulness.py --dataset mouse_cortex
    python scripts/plots/plot_axis_faithfulness.py --dataset mouse_cortex \
        --cohort blueprint_gpt_5_5_unrefined_high --seed 2 --label-col celltype_harmonized
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
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]

DATASET_DEFAULTS = {
    "mouse_cortex":   ("blueprint_gpt_5_5_unrefined_high", 2, "celltype_harmonized"),
    "human_pancreas": ("blueprint_gpt_5_5_unrefined_high", 1, "celltype"),
    "human_immune":   ("blueprint_gpt_5_5_unrefined_high", 1, "final_annotation"),
}
KIND_COLOR = {"plain": "#4c78a8", "contrast": "#f58518", "gated": "#54a24b"}


def _load(name: str):
    path = Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


# ----------------------------------------------------------------------
# Full axis x cell-type AUROC matrix (computed from the embedding)
# ----------------------------------------------------------------------
def axis_type_auroc(dataset: str, stem: str, label_col: str, up) -> tuple:
    """Return (M[n_axes x n_types], axis_names, axis_kinds, type_counts: Series).

    M[i, j] = one-vs-rest AUROC of axis i's activation for cell type j.
    """
    from sklearn.metrics import roc_auc_score

    emb_dir = REPO_ROOT / "outs" / "embeddings" / dataset
    summary = json.loads((emb_dir / f"{stem}_summary.json").read_text())
    axes = summary["blueprint"]["axes"]
    axis_names = [a["name"] for a in axes]
    axis_kinds = [a["kind"] for a in axes]

    df = pd.read_parquet(emb_dir / f"{stem}.parquet").set_index("cell_id")
    df.index = df.index.astype(str)
    dim_cols = [c for c in df.columns if c.startswith("dim_")]
    axis_names = axis_names[:len(dim_cols)]
    axis_kinds = axis_kinds[:len(dim_cols)]

    labels = up.load_obs_columns(
        REPO_ROOT / "data" / "preprocessed" / f"{dataset}_pp.h5ad", [label_col])[label_col]
    lab = df.index.map(labels).astype(str)
    keep = lab.to_numpy() != "nan"
    A = df[dim_cols].to_numpy()[keep]
    lab = pd.Series(lab.to_numpy()[keep])

    counts = lab.value_counts()
    types = list(counts.index)
    M = np.full((len(dim_cols), len(types)), np.nan)
    for j, t in enumerate(types):
        y = (lab.to_numpy() == t).astype(int)
        if y.sum() == 0 or y.sum() == len(y):
            continue
        for i in range(len(dim_cols)):
            M[i, j] = roc_auc_score(y, A[:, i])
    return M, axis_names, axis_kinds, counts


def load_axis_type_map(dataset: str) -> dict:
    """axis_name -> claimed cell type, from the hand-curated mapping CSV.

    The claimed (named-target) type is NOT auto-resolved; it is read from
    configs/axis_type_mappings/<dataset>/axis_types.csv (columns axis_name,
    cell_type, ...), which is edited by hand. Axes with a blank cell_type
    (e.g. pan-class or program axes with no single type) are simply not boxed.
    """
    path = REPO_ROOT / "configs" / "axis_type_mappings" / dataset / "axis_types.csv"
    if not path.exists():
        raise SystemExit(
            f"Missing axis->type mapping {path}. Create it (columns "
            f"axis_name,cell_type) to declare each axis's named cell type.")
    df = pd.read_csv(path).fillna("")
    return {str(r["axis_name"]): str(r["cell_type"]).strip()
            for _, r in df.iterrows() if str(r["cell_type"]).strip()}


# ----------------------------------------------------------------------
# (a) full matrix + (c) coverage bar
# ----------------------------------------------------------------------
def plot_main(M, axis_names, axis_kinds, type_counts, claimed,
              dataset: str, title: str, out_stem: Path) -> None:
    types_all = list(type_counts.index)
    # Cell-type (row) order: reuse the marker-gene dot plot's curated biological
    # grouping (glutamatergic -> GABAergic -> non-neural) so the two figures align;
    # fall back to claimed-types-first-then-orphans for datasets without a curated order.
    curated = [t for _c, ts in _load("plot_gene_dotplot").CELLTYPE_CLASSES.get(dataset, [])
               for t in ts]
    if curated:
        type_order = [t for t in curated if t in set(types_all)]
        type_order += [t for t in types_all if t not in set(type_order)]
    else:
        type_order, seen = [], set()
        for ax in axis_names:
            t = claimed.get(ax)
            if t in type_counts.index and t not in seen:
                type_order.append(t); seen.add(t)
        for t in types_all:
            if t not in seen:
                type_order.append(t); seen.add(t)
    col_ix = {t: i for i, t in enumerate(type_order)}
    Mc = M[:, [types_all.index(t) for t in type_order]]

    # Named-target column per axis, taken verbatim from the hand-curated mapping
    # (no auto-matching). Axes with no declared type get no box (len(type_order)).
    claimed_col = [col_ix.get(claimed.get(ax), len(type_order)) for ax in axis_names]
    row_best = np.nan_to_num(Mc, nan=0.0).max(axis=1)
    row_order = sorted(range(len(axis_names)),
                       key=lambda i: (claimed_col[i], -row_best[i]))
    Mc = Mc[row_order]
    names = [axis_names[i] for i in row_order]
    kinds = [axis_kinds[i] for i in row_order]
    claim_col = [claimed_col[i] for i in row_order]

    n_ax, n_t = Mc.shape
    # Axes on rows, cell types on columns -- cell types share the x-axis with the
    # marker-gene dot plot; colorbar + axis-kind legend stacked at the right.
    # Height factor sized so the panel fills down to the dot plot's bottom when the
    # two are stacked. Font sizes match the dot plot (FS_TYPE 11.5, FS_AXIS 10.5).
    fig, axm = plt.subplots(figsize=(max(7.0, 0.42 * n_t + 3.0), max(4.8, 0.32 * n_ax + 2.2)))

    im = axm.imshow(Mc, aspect="auto", cmap="viridis", vmin=0.0, vmax=1.0)
    axm.set_xticks(range(n_t)); axm.set_xticklabels(type_order, rotation=90, fontsize=11.5)
    axm.set_yticks(range(n_ax)); axm.set_yticklabels(names, fontsize=10.5)
    for tick, k in zip(axm.get_yticklabels(), kinds):
        tick.set_color(KIND_COLOR.get(k, "black"))
    # Box each axis's claimed (named-target) cell (claimed-type column, axis row).
    for r, c in enumerate(claim_col):
        if c < n_t:
            axm.add_patch(mpatches.Rectangle((c - 0.5, r - 0.5), 1, 1, fill=False,
                                             edgecolor="#d62728", lw=2.8))
    # Mark the single best axis per cell type (argmax AUROC over axes).
    for j in range(n_t):
        col = Mc[:, j]
        if np.all(np.isnan(col)):
            continue
        r = int(np.nanargmax(col))
        axm.scatter(j, r, marker="*", s=100, color="white", edgecolor="black",
                    linewidths=1.0, zorder=5)
    # Colorbar (upper right) + axis-kind legend stacked just below it.
    cb = fig.colorbar(im, ax=axm, fraction=0.030, pad=0.015, shrink=0.5, anchor=(0.0, 1.0))
    cb.set_label("AUROC", fontsize=10); cb.ax.tick_params(labelsize=8.5)
    handles = [plt.Line2D([0], [0], marker="s", ls="", mfc=c, mec="none", label=k)
               for k, c in KIND_COLOR.items()]
    axm.legend(handles=handles, title="axis kind", fontsize=10, title_fontsize=10.5,
               loc="upper left", bbox_to_anchor=(1.02, 0.46), frameon=False,
               labelspacing=0.5, handletextpad=0.4)

    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_stem.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"[ok] wrote {out_stem}.{{png,pdf}}")


# ----------------------------------------------------------------------
# across-seed faithfulness distribution
# ----------------------------------------------------------------------
def plot_seeds(dataset, cohort, label_col, seeds, seed_star, up, title, out_stem) -> None:
    # Per seed, the per-cell-type best AUROC = column-max of that run's axis x type
    # matrix (the best any axis discriminates each type). Purely empirical, no
    # name-matching -- so the representative seed is shown to be typical.
    data = []
    for s in seeds:
        M, _, _, type_counts = axis_type_auroc(
            dataset, f"{cohort}_{dataset}_seed{s}", label_col, up)
        col_best = np.nanmax(M, axis=0)
        data.append(col_best[~np.isnan(col_best)])
    fig, ax = plt.subplots(figsize=(max(5.0, 0.5 * len(seeds) + 1.5), 3.4))
    bp = ax.boxplot(data, positions=range(len(seeds)), widths=0.6,
                    patch_artist=True, showfliers=False)
    for i, (box, s) in enumerate(zip(bp["boxes"], seeds)):
        box.set(facecolor="#d62728" if int(s) == seed_star else "#9ecae1",
                edgecolor="0.3", lw=0.6)
        jit = np.random.default_rng(i).uniform(-0.12, 0.12, len(data[i]))
        ax.scatter(np.full(len(data[i]), i) + jit, data[i], s=5, color="0.25",
                   alpha=0.5, zorder=3)
    for m in bp["medians"]:
        m.set(color="black", lw=1.0)
    ax.axhline(0.5, color="k", ls=":", lw=0.7, label="chance")
    ax.set_xticks(range(len(seeds)))
    ax.set_xticklabels([f"seed{int(s)}" for s in seeds], rotation=90, fontsize=7)
    ax.set_ylabel("per-type best-axis AUROC", fontsize=8); ax.set_ylim(0.45, 1.02)
    ax.set_title(title, fontsize=9)
    ax.legend(fontsize=6, frameon=False, loc="lower right")
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_stem.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"[ok] wrote {out_stem}.{{png,pdf}}  (* = representative seed{seed_star})")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Per-axis faithfulness figures for one run.")
    p.add_argument("--dataset", choices=tuple(DATASET_DEFAULTS), default="mouse_cortex")
    p.add_argument("--cohort", default=None)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--label-col", default=None)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    up = _load("umap_plots")
    cohort, seed, label_col = DATASET_DEFAULTS[args.dataset]
    cohort = args.cohort or cohort
    seed = args.seed if args.seed is not None else seed
    label_col = args.label_col or label_col
    stem = f"{cohort}_{args.dataset}_seed{seed}"
    emb_dir = REPO_ROOT / "outs" / "embeddings" / args.dataset

    out_dir = REPO_ROOT / "outs" / "plots" / "axis_metrics" / args.dataset
    main_stem = out_dir / f"axis_faithfulness_{stem}"
    seeds_stem = out_dir / f"axis_auroc_seeds_{cohort}"
    if main_stem.with_suffix(".png").exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite {main_stem}.png (pass --overwrite).")

    claimed = load_axis_type_map(args.dataset)            # hand-curated; no auto-match
    M, axis_names, axis_kinds, type_counts = axis_type_auroc(args.dataset, stem, label_col, up)
    plot_main(M, axis_names, axis_kinds, type_counts, claimed, args.dataset,
              f"{args.dataset} · blueprint seed{seed} · axis faithfulness ({label_col})",
              main_stem)

    # All seeds of this cohort on disk -> across-seed reproducibility panel.
    import re
    seeds = sorted(int(re.search(r"_seed(\d+)\.parquet$", p.name).group(1))
                   for p in emb_dir.glob(f"{cohort}_{args.dataset}_seed*.parquet")
                   if re.search(r"_seed(\d+)\.parquet$", p.name))
    plot_seeds(args.dataset, cohort, label_col, seeds, seed, up,
               f"{args.dataset} · per-type best-axis AUROC across {len(seeds)} seeds", seeds_stem)



if __name__ == "__main__":
    main()
