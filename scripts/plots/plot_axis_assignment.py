"""
Argmax axis-assignment view for ONE blueprint run.

A secondary, intuitive complement to the AUROC faithfulness panels: assign each
cell to the axis on which it scores highest (argmax over the z-scored axes, the
same convention used for the developmental-tree gated blueprint) and ask where
each annotated cell type's cells land.

  (a) cell-type x axis assignment composition -- row-normalized fraction of each
      cell type's cells whose argmax is each axis. A near-block-diagonal means
      each type is captured by its own named axis.
  (b) assignment margin -- per cell, the gap between the top and second axis
      z-score; larger margins mean less ambiguous assignments.

Caveat: argmax over independently z-scored axes is a lossy, scale-sensitive
reduction (a cell can argmax onto a broadly-expressed program rather than its
identity axis); it is reported as a descriptive complement, not the faithfulness
instrument -- that is the threshold-free per-axis AUROC (plot_axis_faithfulness.py).

Inputs:
    outs/embeddings/<dataset>/<stem>.parquet
    outs/embeddings/<dataset>/<stem>_summary.json
    data/preprocessed/<dataset>_pp.h5ad

Outputs:
    outs/plots/axis_metrics/<dataset>/axis_assignment_<stem>.{png,pdf}

Usage:
    python scripts/plots/plot_axis_assignment.py --dataset mouse_cortex
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
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]

DATASET_DEFAULTS = {
    "mouse_cortex":   ("blueprint_gpt_5_5_unrefined_high", 2, "celltype_harmonized"),
    "human_pancreas": ("blueprint_gpt_5_5_unrefined_high", 1, "celltype"),
    "human_immune":   ("blueprint_gpt_5_5_unrefined_high", 1, "final_annotation"),
}


def _load(name: str):
    path = Path(__file__).resolve().parent / f"{name}.py"
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_run(dataset: str, stem: str, label_col: str, up):
    emb_dir = REPO_ROOT / "outs" / "embeddings" / dataset
    summary = json.loads((emb_dir / f"{stem}_summary.json").read_text())
    axis_names = [a["name"] for a in summary["blueprint"]["axes"]]
    df = pd.read_parquet(emb_dir / f"{stem}.parquet").set_index("cell_id")
    df.index = df.index.astype(str)
    dim_cols = [c for c in df.columns if c.startswith("dim_")]
    axis_names = axis_names[:len(dim_cols)]
    labels = up.load_obs_columns(
        REPO_ROOT / "data" / "preprocessed" / f"{dataset}_pp.h5ad", [label_col])[label_col]
    lab = df.index.map(labels).astype(str)
    keep = lab.to_numpy() != "nan"
    return df[dim_cols].to_numpy()[keep], pd.Series(lab.to_numpy()[keep]), axis_names


def plot(A: np.ndarray, lab: pd.Series, axis_names: list[str], title: str,
         out_stem: Path, curated: list[str], claimed: dict) -> None:
    # argmax assignment + margin (top1 - top2).
    order = np.argsort(-A, axis=1)
    assigned = order[:, 0]
    top1 = np.take_along_axis(A, order[:, :1], axis=1).ravel()
    top2 = np.take_along_axis(A, order[:, 1:2], axis=1).ravel()
    margin = top1 - top2

    # Cell-type (row) order: reuse the curated biological grouping that orders the
    # axis x cell-type AUROC heatmap's columns, so the two figures align; any type
    # not in the curated list falls in afterwards.
    types_all = list(lab.value_counts().index)
    type_order = [t for t in curated if t in set(types_all)]
    type_order += [t for t in types_all if t not in type_order]
    type_ix = {t: i for i, t in enumerate(type_order)}
    n_t, n_a = len(type_order), len(axis_names)
    C = np.zeros((n_t, n_a))
    codes = lab.map(type_ix).to_numpy()
    for ti, ai in zip(codes, assigned):
        C[ti, ai] += 1
    C = C / C.sum(axis=1, keepdims=True).clip(min=1)
    # Axis (column) order: by each axis's named cell type's position in type_order,
    # matching the AUROC heatmap's row order; axes that name no type fall last
    # (in declaration order).
    col_order = sorted(range(n_a),
                       key=lambda i: (type_ix.get(claimed.get(axis_names[i]), n_t), i))
    C = C[:, col_order]
    axis_lab = [axis_names[i] for i in col_order]

    fig = plt.figure(figsize=(max(9.5, 0.36 * n_a + 3.5), max(5.0, 0.22 * n_t + 1.8)))
    gs = fig.add_gridspec(1, 2, width_ratios=[n_a, max(12, 0.5 * n_a)], wspace=0.30)
    axm = fig.add_subplot(gs[0, 0])
    axh = fig.add_subplot(gs[0, 1])

    im = axm.imshow(C, aspect="auto", cmap="magma_r", vmin=0, vmax=1)
    axm.set_xticks(range(n_a)); axm.set_xticklabels(axis_lab, rotation=90, fontsize=8)
    axm.set_yticks(range(n_t)); axm.set_yticklabels(type_order, fontsize=8.5)
    axm.set_title("(a) cell-type → argmax-axis composition (row-normalized)", fontsize=11)
    cb = fig.colorbar(im, ax=axm, fraction=0.025, pad=0.01)
    cb.set_label("fraction of type's cells", fontsize=9.5); cb.ax.tick_params(labelsize=8)

    axh.hist(margin, bins=60, color="#6a51a3", alpha=0.85)
    axh.axvline(float(np.median(margin)), color="k", ls="--", lw=1.0,
                label=f"median = {np.median(margin):.2f}")
    axh.set_xlabel("assignment margin\n(top1 − top2 axis z-score)", fontsize=10.5)
    axh.set_ylabel("cells", fontsize=10.5)
    axh.set_title("(b) assignment margin", fontsize=11)
    axh.legend(fontsize=9.5, frameon=False)

    fig.suptitle(title, fontsize=12.5, y=1.005)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_stem.with_suffix(f".{ext}"), dpi=300, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)
    print(f"[ok] wrote {out_stem}.{{png,pdf}}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Argmax axis-assignment view.")
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

    out_stem = REPO_ROOT / "outs" / "plots" / "axis_metrics" / args.dataset / f"axis_assignment_{stem}"
    if out_stem.with_suffix(".png").exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite {out_stem}.png (pass --overwrite).")

    A, lab, axis_names = load_run(args.dataset, stem, label_col, up)
    curated = [t for _c, ts in _load("plot_gene_dotplot").CELLTYPE_CLASSES.get(args.dataset, [])
               for t in ts]
    claimed = _load("plot_axis_faithfulness").load_axis_type_map(args.dataset)
    plot(A, lab, axis_names,
         f"{args.dataset} · blueprint seed{seed} · argmax axis assignment ({label_col})",
         out_stem, curated, claimed)



if __name__ == "__main__":
    main()
