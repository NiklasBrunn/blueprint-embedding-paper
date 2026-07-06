"""
Common preprocessing for QC'd single-cell data.

Standard workflow (no scaling, no PCA):
  1) normalize_total(target_sum)              -> layers["norm"]
  2) log1p                                    -> layers["log1p"]; also sets X
  3) highly_variable_genes(seurat_v3, n=2000, layer="counts", batch_key=...)
     -> var["highly_variable"]                (flag only; no subsetting)

Inputs:
    data/qc/raw_<dataset>_qc.h5ad             (produced by scripts/data_qc/qc.py)

Outputs (per dataset):
    data/preprocessed/<dataset>/preprocessing_report.md
    data/preprocessed/<dataset>/plots/{hvg_dispersion,highest_expr_genes}.png
    data/preprocessed/<dataset>_pp.h5ad       (counts + norm + log1p layers, HVGs flagged)

Usage:
    python scripts/data_preprocessing/preprocess.py --dataset mouse_brain
    python scripts/data_preprocessing/preprocess.py --dataset human_pancreas \
        --n-top-genes 3000 --target-sum 1e4
"""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from datetime import datetime
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")

# Per-dataset defaults. batch_key drives per-batch HVG ranking so no single
# study dominates the HVG set: mouse_brain has no natural batch; pancreas is
# multi-tech; mouse_gastrulation batches per embryo (obs["sample_id"]);
# mouse_cortex per study (obs["study"]); human_immune (scIB "Immune (human)"
# benchmark) per batch (obs["batch"]).
DEFAULT_PARAMS = {
    "mouse_brain":         {"n_top_genes": 2000, "target_sum": 1e4, "flavor": "seurat_v3", "batch_key": None},
    "human_pancreas":      {"n_top_genes": 2000, "target_sum": 1e4, "flavor": "seurat_v3", "batch_key": "tech"},
    "mouse_gastrulation":  {"n_top_genes": 2000, "target_sum": 1e4, "flavor": "seurat_v3", "batch_key": "sample_id"},
    "mouse_cortex":        {"n_top_genes": 2000, "target_sum": 1e4, "flavor": "seurat_v3", "batch_key": "study"},
    "human_immune":        {"n_top_genes": 2000, "target_sum": 1e4, "flavor": "seurat_v3", "batch_key": "batch"},
}


def normalize_and_log1p(adata: ad.AnnData, target_sum: float) -> None:
    if "counts" not in adata.layers:
        raise RuntimeError("Expected raw counts in adata.layers['counts'].")
    adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=target_sum)
    adata.layers["norm"] = adata.X.copy()
    sc.pp.log1p(adata)
    adata.layers["log1p"] = adata.X.copy()


def select_hvgs(adata: ad.AnnData, n_top_genes: int, flavor: str, batch_key: str | None) -> None:
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
        flavor=flavor,
        layer="counts",
        batch_key=batch_key,
        subset=False,
    )


def save_plots(adata: ad.AnnData, plots_dir: Path) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
    hv = adata.var["highly_variable"].to_numpy()
    score_col = "variances_norm" if "variances_norm" in adata.var else "dispersions_norm"
    means = adata.var["means"].to_numpy()
    score = adata.var[score_col].to_numpy()

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(means[~hv], score[~hv], s=4, alpha=0.4, c="grey", label="other")
    ax.scatter(means[hv], score[hv], s=6, alpha=0.7, c="C0", label="HVG")
    ax.set_xscale("log")
    ax.set_xlabel("mean expression")
    ax.set_ylabel(score_col.replace("_", " "))
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(plots_dir / "hvg_dispersion.png", dpi=120, bbox_inches="tight")
    plt.close(fig)

    sc.pl.highest_expr_genes(adata, n_top=20, show=False)
    plt.savefig(plots_dir / "highest_expr_genes.png", dpi=120, bbox_inches="tight")
    plt.close()


def hvg_summary(adata: ad.AnnData, top_n: int = 20) -> dict:
    score_col = "variances_norm" if "variances_norm" in adata.var else "dispersions_norm"
    rank_col = "highly_variable_rank" if "highly_variable_rank" in adata.var else None
    hvg = adata.var.loc[adata.var["highly_variable"]]
    ranked = hvg.sort_values(rank_col) if rank_col else hvg.sort_values(score_col, ascending=False)
    return {
        "score_col": score_col,
        "score_min": float(ranked[score_col].min()),
        "score_max": float(ranked[score_col].max()),
        "mean_min": float(ranked["means"].min()),
        "mean_max": float(ranked["means"].max()),
        "top_genes": ranked.head(top_n).index.astype(str).tolist(),
    }


def stats(adata: ad.AnnData) -> dict:
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_hvgs": int(adata.var["highly_variable"].sum()),
    }


def write_report(report: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    p, s, h = report["params"], report["stats"], report["hvg"]
    lines = [
        f"# Preprocessing report — {report['dataset']}",
        "",
        f"- Timestamp: `{report['timestamp']}`",
        f"- Input:  `{report['input_path']}`",
        f"- Output: `{report['output_path']}`",
        "",
        "## Parameters",
        *(f"- `{k}`: {v}" for k, v in p.items()),
        "",
        "## Steps",
        f"- normalize_total(target_sum={p['target_sum']}) → `layers['norm']`",
        "- log1p → `layers['log1p']`; also set `X`",
        f"- highly_variable_genes(flavor={p['flavor']}, n_top_genes={p['n_top_genes']}, "
        f"layer='counts', batch_key={p['batch_key']!r}, subset=False) → `var['highly_variable']`",
        "",
        "## Counts",
        f"- cells: {s['n_cells']}",
        f"- genes: {s['n_genes']}",
        f"- HVGs flagged: {s['n_hvgs']}",
        "",
        "## HVG summary",
        f"- score column: `{h['score_col']}`",
        f"- score range (HVGs): {h['score_min']:.4g} .. {h['score_max']:.4g}",
        f"- mean expression range (HVGs): {h['mean_min']:.4g} .. {h['mean_max']:.4g}",
        "",
        f"### Top {len(h['top_genes'])} HVGs",
        *(f"{i}. `{g}`" for i, g in enumerate(h["top_genes"], 1)),
    ]
    path.write_text("\n".join(lines) + "\n")


def resolve_params(args: argparse.Namespace) -> dict:
    p = DEFAULT_PARAMS[args.dataset].copy()
    if args.n_top_genes is not None:
        p["n_top_genes"] = args.n_top_genes
    if args.target_sum is not None:
        p["target_sum"] = args.target_sum
    if args.flavor is not None:
        p["flavor"] = args.flavor
    if args.batch_key is not None:
        p["batch_key"] = None if args.batch_key.lower() == "none" else args.batch_key
    return p


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Common preprocessing for QC'd single-cell data.")
    parser.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    parser.add_argument("--n-top-genes", type=int)
    parser.add_argument("--target-sum", type=float)
    parser.add_argument("--flavor", choices=("seurat_v3", "seurat", "cell_ranger"))
    parser.add_argument("--batch-key", type=str, help="Set to 'none' to disable.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = resolve_params(args)

    input_path  = REPO_ROOT / "data" / "qc" / f"raw_{args.dataset}_qc.h5ad"
    pp_dir      = REPO_ROOT / "data" / "preprocessed" / args.dataset
    plots_dir   = pp_dir / "plots"
    output_path = REPO_ROOT / "data" / "preprocessed" / f"{args.dataset}_pp.h5ad"
    report_path = pp_dir / "preprocessing_report.md"

    adata = ad.read_h5ad(input_path)
    normalize_and_log1p(adata, target_sum=params["target_sum"])
    select_hvgs(adata, params["n_top_genes"], params["flavor"], params["batch_key"])
    save_plots(adata, plots_dir)

    report = {
        "dataset": args.dataset,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "input_path": str(input_path),
        "output_path": str(output_path),
        "params": params,
        "stats": stats(adata),
        "hvg": hvg_summary(adata),
    }
    write_report(report, report_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.uns["preprocessing"] = json.dumps(report, default=str)
    adata.write_h5ad(output_path)

    print(f"[ok] {input_path.name}: {report['stats']['n_cells']}x{report['stats']['n_genes']}, "
          f"{report['stats']['n_hvgs']} HVGs flagged")
    print(f"[ok] wrote {output_path}")
    print(f"[ok] wrote {report_path}")
    print(f"[ok] wrote plots in {plots_dir}")


if __name__ == "__main__":
    main()
