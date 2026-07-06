#!/usr/bin/env python3
"""Temporal-fidelity of a gated blueprint tree: does tree depth encode time?

The gastrulation blueprint is authored (via the scenario `notes`) so that its
gated axes form a time-forward tree --- depth = developmental time, root =
earliest state, leaves = latest. This script tests that structural property
directly against the atlas's physical stage label (E6.5--E8.5), independent of
any trajectory algorithm: the time ordering is *imposed by construction*, not
discovered from expression similarity.

For each blueprint embedding parquet under outs/embeddings/<dataset>/ we
  1. recover the gate tree from the run's `_report.md` (reusing
     blueprint_tree.parse_axis_table / build_edges), giving every axis a depth;
  2. assign each cell to its argmax axis over dim_*, and take the mean numeric
     embryonic stage of each axis's owned cells (parse_numeric_stage; cells with
     no resolvable numeric stage, e.g. `mixed_gastrulation`, are dropped);
  3. score the ordering: a pooled Spearman(axis depth, axis mean stage) over all
     axes, the mean Spearman along each root->leaf branch, and the fraction of
     parent->child steps along which mean stage is non-decreasing (monotonicity).

Outputs (swept over all blueprint runs of the dataset):
    outs/temporal_metrics/<dataset>/temporal_metrics.csv        per-run summary
    outs/temporal_metrics/<dataset>/temporal_metrics_axes.csv   per-axis detail

Usage:
    python scripts/temporal_metrics/metrics.py --dataset mouse_gastrulation
    python scripts/temporal_metrics/metrics.py --dataset mouse_gastrulation \
        --input outs/embeddings/mouse_gastrulation/blueprint_<...>_seed3.parquet
"""
from __future__ import annotations

import argparse
import importlib.util
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_gastrulation",)

# obs columns: the physical-time label whose ordering the tree should recover,
# and the cell-type label used only for the per-axis branch-point sanity check.
STAGE_COL = {"mouse_gastrulation": "stage"}
CELLTYPE_COL = {"mouse_gastrulation": "celltype"}


def _load_module(name: str, relpath: str):
    """Import a sibling repo script by file path (avoids the metrics.py name
    clash and keeps the heavy deps of those modules lazy)."""
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relpath)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


_BT = _load_module("bp_tree", "scripts/plots/blueprint_tree.py")
_EM = _load_module("emb_metrics", "scripts/embedding_metrics/metrics.py")
parse_axis_table = _BT.parse_axis_table
build_edges = _BT.build_edges
roots_and_children = _BT.roots_and_children
parse_numeric_stage = _EM.parse_numeric_stage


# ============================================================ tree geometry


def axis_depths(axes: list[dict], parent_of: dict[str, str]) -> dict[str, int]:
    """Depth of every axis (root = 0) over the gate/stitch tree."""
    roots, children = roots_and_children(axes, parent_of)
    depth: dict[str, int] = {}
    stack = [(r, 0) for r in roots]
    while stack:
        name, d = stack.pop()
        depth[name] = d
        for kid in children.get(name, []):
            stack.append((kid, d + 1))
    # disconnected safety net (should not happen): unseen axes -> depth 0
    for ax in axes:
        depth.setdefault(ax["name"], 0)
    return depth


def root_to_leaf_paths(axes: list[dict], parent_of: dict[str, str]) -> list[list[str]]:
    """Every root->leaf path as an ordered list of axis names."""
    roots, children = roots_and_children(axes, parent_of)
    paths: list[list[str]] = []

    def walk(name: str, acc: list[str]) -> None:
        acc = acc + [name]
        kids = children.get(name, [])
        if not kids:
            paths.append(acc)
            return
        for kid in kids:
            walk(kid, acc)

    for r in roots:
        walk(r, [])
    return paths


# ============================================================ per-run scoring


def monotonic_fraction(values: list[float]) -> float:
    """Fraction of consecutive steps that are non-decreasing (NaNs dropped)."""
    v = [x for x in values if not np.isnan(x)]
    if len(v) < 2:
        return float("nan")
    steps = np.diff(v)
    return float(np.mean(steps >= -1e-9))


def score_run(parquet: Path, stage_col: str, celltype_col: str,
              top_frac: float) -> tuple[dict, pd.DataFrame]:
    """Return (summary_row, per_axis_df) for one blueprint parquet."""
    stem = parquet.stem
    report = parquet.with_name(stem + "_report.md")
    axes = parse_axis_table(report)
    _, parent_of = build_edges(axes, stitch=True)
    depth = axis_depths(axes, parent_of)

    tbl = pq.read_table(parquet)
    df = tbl.to_pandas()
    dim_cols = [c for c in df.columns if c.startswith("dim_")]
    if len(dim_cols) != len(axes):
        raise SystemExit(
            f"{stem}: {len(dim_cols)} dim cols but {len(axes)} axis rows in report")
    X = df[dim_cols].to_numpy(dtype=np.float32)  # axis i == dim_i == axes[i]
    cell_ids = df["cell_id"].to_numpy()

    # align the physical-stage + celltype labels to the parquet's cells
    obs = ad.read_h5ad(parquet_dataset_h5ad(parquet), backed="r").obs
    obs = obs.reindex(cell_ids)
    stage_num = np.array([parse_numeric_stage(v) for v in obs[stage_col].to_numpy()],
                         dtype=float)
    celltype = obs[celltype_col].astype("object").to_numpy()

    owner = X.argmax(axis=1)  # axis each cell scores most strongly (z-scored axes)
    valid = ~np.isnan(stage_num)

    rows = []
    for i, ax in enumerate(axes):
        owned = owner == i
        owned_valid = owned & valid
        mean_stage = float(np.mean(stage_num[owned_valid])) if owned_valid.any() else float("nan")
        # robustness variant: mean stage over this axis's top-frac scoring cells
        k = max(1, int(top_frac * X.shape[0]))
        top_idx = np.argpartition(-X[:, i], k - 1)[:k]
        top_valid = top_idx[valid[top_idx]]
        mean_stage_top = float(np.mean(stage_num[top_valid])) if top_valid.size else float("nan")
        ct = celltype[owned]
        dominant = pd.Series(ct).dropna().mode()
        rows.append({
            "axis": ax["name"], "kind": ax["kind"], "depth": depth[ax["name"]],
            "parent": parent_of.get(ax["name"], ""),
            "n_owned": int(owned.sum()), "mean_stage": mean_stage,
            "mean_stage_top": mean_stage_top,
            "dominant_celltype": dominant.iloc[0] if len(dominant) else "",
        })
    axis_df = pd.DataFrame(rows)

    # pooled depth<->stage ordering over all axes that own valid cells
    pool = axis_df.dropna(subset=["mean_stage"])
    pooled_rho = (float(spearmanr(pool["depth"], pool["mean_stage"]).statistic)
                  if pool["depth"].nunique() > 1 and len(pool) > 2 else float("nan"))

    # per branch: Spearman(depth-order, mean stage) + monotonicity
    mean_by_axis = dict(zip(axis_df["axis"], axis_df["mean_stage"]))
    branch_rhos, mono_steps_ok, mono_steps_tot = [], 0, 0
    paths = root_to_leaf_paths(axes, parent_of)
    for path in paths:
        seq = [mean_by_axis[a] for a in path]
        valid_seq = [s for s in seq if not np.isnan(s)]
        if len(valid_seq) >= 3 and len(set(valid_seq)) > 1:
            order = list(range(len(valid_seq)))
            branch_rhos.append(float(spearmanr(order, valid_seq).statistic))
        if len(valid_seq) >= 2:
            steps = np.diff(valid_seq)
            mono_steps_ok += int(np.sum(steps >= -1e-9))
            mono_steps_tot += len(steps)

    cohort, seed = parse_cohort_seed(stem, parquet.parent.name)
    summary = {
        "model": cohort, "run": stem, "seed": seed,
        "n_axes": len(axes), "n_roots": len(roots_and_children(axes, parent_of)[0]),
        "n_branches": len(paths), "max_depth": int(axis_df["depth"].max()),
        "n_cells_used": int(valid.sum()),
        "pooled_depth_stage_spearman": pooled_rho,
        "mean_branch_spearman": float(np.nanmean(branch_rhos)) if branch_rhos else float("nan"),
        "monotonic_frac": (mono_steps_ok / mono_steps_tot) if mono_steps_tot else float("nan"),
    }
    axis_df.insert(0, "seed", seed)
    axis_df.insert(0, "model", cohort)
    return summary, axis_df


# ============================================================ small helpers


def parquet_dataset_h5ad(parquet: Path) -> Path:
    """The preprocessed AnnData the embedding was built from."""
    # outs/embeddings/<dataset>/<stem>.parquet -> data/preprocessed/<dataset>_pp.h5ad
    dataset = parquet.parent.name
    return REPO_ROOT / "data" / "preprocessed" / f"{dataset}_pp.h5ad"


def parse_cohort_seed(stem: str, dataset: str) -> tuple[str, int | None]:
    """`blueprint_<...>_<dataset>_seed<N>` -> ('blueprint_<...>', N)."""
    m = re.search(r"_seed(\d+)$", stem)
    seed = int(m.group(1)) if m else None
    cohort = stem[:m.start()] if m else stem
    suffix = f"_{dataset}"
    if cohort.endswith(suffix):
        cohort = cohort[:-len(suffix)]
    return cohort, seed


def discover(dataset: str) -> list[Path]:
    emb_dir = REPO_ROOT / "outs" / "embeddings" / dataset
    return sorted(emb_dir.glob("blueprint_*.parquet"))


# ============================================================ main


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", choices=DATASET_CHOICES, default="mouse_gastrulation")
    p.add_argument("--input", type=Path, default=None,
                   help="Single blueprint parquet (overrides discovery).")
    p.add_argument("--top-frac", type=float, default=0.05,
                   help="Robustness variant: mean stage over each axis's top fraction.")
    args = p.parse_args()

    stage_col = STAGE_COL[args.dataset]
    celltype_col = CELLTYPE_COL[args.dataset]
    parquets = [args.input] if args.input else discover(args.dataset)
    if not parquets:
        raise SystemExit(f"No blueprint parquets found for {args.dataset}")

    summaries, axis_frames = [], []
    for pq_path in parquets:
        print(f"[run] {pq_path.name}", flush=True)
        summary, axis_df = score_run(pq_path, stage_col, celltype_col, args.top_frac)
        summaries.append(summary)
        axis_frames.append(axis_df)
        print(f"      depth<->stage rho={summary['pooled_depth_stage_spearman']:.3f} "
              f"branch rho={summary['mean_branch_spearman']:.3f} "
              f"mono={summary['monotonic_frac']:.3f} "
              f"(axes={summary['n_axes']}, branches={summary['n_branches']})", flush=True)

    out_dir = REPO_ROOT / "outs" / "temporal_metrics" / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    sum_df = pd.DataFrame(summaries).sort_values(["model", "seed"])
    sum_df.to_csv(out_dir / "temporal_metrics.csv", index=False)
    pd.concat(axis_frames, ignore_index=True).to_csv(
        out_dir / "temporal_metrics_axes.csv", index=False)
    print(f"[ok] wrote {out_dir/'temporal_metrics.csv'} ({len(sum_df)} runs)")

    # quick cohort summary (mean +/- sd across seeds)
    print("\n[cohort mean +/- sd]")
    for model, g in sum_df.groupby("model"):
        for col in ["pooled_depth_stage_spearman", "monotonic_frac"]:
            print(f"  {model:42s} {col:28s} "
                  f"{g[col].mean():.3f} +/- {g[col].std(ddof=0):.3f} (n={len(g)})")


if __name__ == "__main__":
    main()
