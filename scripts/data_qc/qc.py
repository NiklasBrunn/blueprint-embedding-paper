"""
Common quality control for raw single-cell data.

Pure QC: compute QC metrics, filter cells/genes, generate before/after plots
and a markdown report. **No preprocessing** (no normalization, log1p, HVG
selection, dimensionality reduction).

Inputs:
    data/raw/raw_<dataset>.h5ad           (produced by scripts/IO/load_<dataset>.py)

Outputs (per dataset):
    data/qc/<dataset>/qc_report.md        (thresholds, before/after counts, summary stats)
    data/qc/<dataset>/plots/before_*.png  (pre-QC violin + scatter)
    data/qc/<dataset>/plots/after_*.png   (post-QC violin + scatter)
    data/qc/raw_<dataset>_qc.h5ad         (filtered raw counts; X and layers["counts"] unchanged)

Usage:
    python scripts/data_qc/qc.py --dataset mouse_brain
    python scripts/data_qc/qc.py --dataset human_pancreas \
        --min-genes 200 --max-pct-mt 20 --min-cells 3
"""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from datetime import datetime
from pathlib import Path

os.environ.setdefault("NUMBA_DISABLE_JIT", "1")
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")

# Sensible defaults per dataset (overridable via config or CLI). For
# mixed-protocol datasets (mouse_cortex, human_immune) thresholds are biased
# toward the droplet floor so QC does not become a covert per-protocol
# preprocessing step. The mt mask is case-insensitive (var_names lower-cased
# before the "mt-" prefix test), covering both mouse "mt-Co1" and human
# "MT-CO1" symbol forms.
DEFAULT_THRESHOLDS = {
    "mouse_brain":         {"min_genes": 1000, "max_pct_mt": 10, "min_cells": 3},
    "human_pancreas":      {"min_genes": 200,  "max_pct_mt": 20, "min_cells": 3},
    "mouse_gastrulation":  {"min_genes": 500,  "max_pct_mt": 5,  "min_cells": 3},
    "mouse_cortex":        {"min_genes": 500,  "max_pct_mt": 10, "min_cells": 3},
    "human_immune":        {"min_genes": 200,  "max_pct_mt": 20, "min_cells": 3},
}

QC_KEYS = ("n_genes_by_counts", "total_counts", "pct_counts_mt")


def calculate_qc_metrics(adata: ad.AnnData) -> None:
    var_lower = adata.var_names.str.lower()
    adata.var["mt"] = var_lower.str.startswith("mt-")
    adata.var["ribo"] = var_lower.str.startswith(("rps", "rpl"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True
    )


def save_qc_plots(adata: ad.AnnData, plots_dir: Path, label: str) -> None:
    plots_dir.mkdir(parents=True, exist_ok=True)
    sc.pl.violin(adata, keys=list(QC_KEYS), jitter=0.4, multi_panel=True, show=False)
    plt.savefig(plots_dir / f"{label}_violin.png", dpi=120, bbox_inches="tight")
    plt.close()

    sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", color="pct_counts_mt", show=False)
    plt.savefig(plots_dir / f"{label}_scatter_counts_vs_genes.png", dpi=120, bbox_inches="tight")
    plt.close()

    sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=False)
    plt.savefig(plots_dir / f"{label}_scatter_counts_vs_mt.png", dpi=120, bbox_inches="tight")
    plt.close()


def filter_data(
    adata: ad.AnnData, min_genes: int, max_pct_mt: float, min_cells: int
) -> tuple[ad.AnnData, list[dict]]:
    steps: list[dict] = []
    n = adata.n_obs
    sc.pp.filter_cells(adata, min_genes=min_genes)
    steps.append({"step": f"min_genes_per_cell >= {min_genes}", "removed_cells": int(n - adata.n_obs)})

    n = adata.n_obs
    adata = adata[adata.obs["pct_counts_mt"] <= max_pct_mt].copy()
    steps.append({"step": f"pct_counts_mt <= {max_pct_mt}", "removed_cells": int(n - adata.n_obs)})

    n = adata.n_vars
    sc.pp.filter_genes(adata, min_cells=min_cells)
    steps.append({"step": f"min_cells_per_gene >= {min_cells}", "removed_genes": int(n - adata.n_vars)})
    return adata, steps


def summary_stats(adata: ad.AnnData) -> dict:
    return {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "median_genes_per_cell": float(adata.obs["n_genes_by_counts"].median()),
        "median_total_counts": float(adata.obs["total_counts"].median()),
        "median_pct_mt": float(adata.obs["pct_counts_mt"].median()),
    }


def write_report(report: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        f"# QC report — {report['dataset']}",
        "",
        f"- Timestamp: `{report['timestamp']}`",
        f"- Input:  `{report['input_path']}`",
        f"- Output: `{report['output_path']}`",
        "",
        "## Thresholds",
    ]
    for k, v in report["thresholds"].items():
        lines.append(f"- `{k}`: {v}")
    lines += [
        "",
        "## Cell / gene counts",
        f"- before: {report['before']['n_cells']} cells × {report['before']['n_genes']} genes",
        f"- after:  {report['after']['n_cells']} cells × {report['after']['n_genes']} genes",
        "",
        "## Filtering steps",
    ]
    for step in report["steps"]:
        kind, n = ("removed_cells", step.get("removed_cells")) if "removed_cells" in step else ("removed_genes", step.get("removed_genes"))
        lines.append(f"- {step['step']} → {kind}={n}")
    lines += ["", "## Summary stats (before → after)"]
    for key in ("n_cells", "n_genes", "median_genes_per_cell", "median_total_counts", "median_pct_mt"):
        lines.append(f"- `{key}`: {report['before'][key]} → {report['after'][key]}")
    path.write_text("\n".join(lines) + "\n")


def resolve_thresholds(args: argparse.Namespace) -> dict:
    th = DEFAULT_THRESHOLDS[args.dataset].copy()
    if args.min_genes is not None:
        th["min_genes"] = args.min_genes
    if args.max_pct_mt is not None:
        th["max_pct_mt"] = args.max_pct_mt
    if args.min_cells is not None:
        th["min_cells"] = args.min_cells
    return th


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Common QC for raw single-cell data.")
    parser.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    parser.add_argument("--min-genes", type=int)
    parser.add_argument("--max-pct-mt", type=float)
    parser.add_argument("--min-cells", type=int)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    thresholds = resolve_thresholds(args)

    input_path = REPO_ROOT / "data" / "raw" / f"raw_{args.dataset}.h5ad"
    qc_dir = REPO_ROOT / "data" / "qc" / args.dataset
    plots_dir = qc_dir / "plots"
    output_path = REPO_ROOT / "data" / "qc" / f"raw_{args.dataset}_qc.h5ad"
    report_path = qc_dir / "qc_report.md"

    adata = ad.read_h5ad(input_path)
    calculate_qc_metrics(adata)
    save_qc_plots(adata, plots_dir, "before")
    before = summary_stats(adata)

    adata, steps = filter_data(
        adata,
        min_genes=thresholds["min_genes"],
        max_pct_mt=thresholds["max_pct_mt"],
        min_cells=thresholds["min_cells"],
    )
    save_qc_plots(adata, plots_dir, "after")
    after = summary_stats(adata)

    report = {
        "dataset": args.dataset,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "input_path": str(input_path),
        "output_path": str(output_path),
        "thresholds": thresholds,
        "steps": steps,
        "before": before,
        "after": after,
    }
    write_report(report, report_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.uns["qc"] = json.dumps(report, default=str)
    adata.write_h5ad(output_path)

    print(f"[ok] {input_path.name}: {before['n_cells']}x{before['n_genes']} -> "
          f"{after['n_cells']}x{after['n_genes']} cells x genes")
    print(f"[ok] wrote {output_path}")
    print(f"[ok] wrote {report_path}")
    print(f"[ok] wrote plots in {plots_dir}")


if __name__ == "__main__":
    main()
