"""
Canonical scIB results table across all models, from the saved Benchmarker output.

scripts/embedding_metrics/metrics.py writes outs/metrics/<dataset>/scib_table_<label>.csv
(the raw `Benchmarker.get_results()` frame: rows = embedding parquets + a
'Metric Type' row; columns = per-metric scores + Bio/Batch/Total aggregates). This
script reads that table, **aggregates the per-seed parquet rows into one row per
model** (mean across seeds), and renders the package's own `plot_results_table()`
styling — so no re-benchmarking is needed (the metrics run already did the work).

Absolute scores by default (matching metrics.csv); pass --min-max for the
panel-relative (within-figure) leaderboard view.

Inputs:
    outs/metrics/<dataset>/scib_table_<label_col>.csv   (written by metrics.py)

Outputs:
    outs/plots/metrics/<dataset>/scib_results_table_<label_col>.{svg,png}

Usage:
    python scripts/plots/plot_scib_table.py --dataset human_pancreas
    python scripts/plots/plot_scib_table.py --dataset mouse_brain --min-max
"""
from __future__ import annotations

import argparse
import importlib.util
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_metrics_module():
    """Load scripts/embedding_metrics/metrics.py (for parse_stem + prettify_model)."""
    path = REPO_ROOT / "scripts" / "embedding_metrics" / "metrics.py"
    spec = importlib.util.spec_from_file_location("embedding_metrics", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build_shim_benchmarker(table_csv: Path, dataset: str, parse_stem, prettify,
                           exclude: list[str] | None = None):
    """Reconstruct a minimal scib-metrics Benchmarker whose results are aggregated
    to one row per model, so its native plot_results_table() can render it.

    The saved CSV is get_results() output: index = parquet stems + 'Metric Type';
    columns = per-metric scores + Bio/Batch/Total aggregates. Inverted back to
    the Benchmarker's internal `_results` layout (metrics x models + a 'Metric Type'
    column), averaging the per-metric scores across seeds per model. Aggregate
    columns are dropped — get_results() recomputes them from the per-metric means.

    `exclude` is a list of tokens matched against the underscore-delimited words
    of the parsed model name (so 'refine' drops '..._refine_...' but NOT
    '..._unrefined_...'); any stem with a matching model is dropped before
    aggregation (e.g. restrict to the benchmark cohort set via 'freeform,refine').
    """
    from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
    from scib_metrics.benchmark._core import _METRIC_TYPE, _AGGREGATE_SCORE

    df = pd.read_csv(table_csv, index_col=0)
    mtype = df.loc[_METRIC_TYPE]                       # per-column metric type
    per_metric = [c for c in df.columns if mtype[c] != _AGGREGATE_SCORE]
    emb = df.drop(index=_METRIC_TYPE)[per_metric].astype(float)  # rows=stems

    model_of = {s: parse_stem(s, dataset)[0] for s in emb.index}
    if exclude:
        keep = [s for s in emb.index
                if not (set(exclude) & set(model_of[s].split("_")))]
        dropped = sorted({model_of[s] for s in emb.index if s not in keep})
        if dropped:
            print(f"[info] excluded models: {dropped}")
        emb = emb.loc[keep]
        model_of = {s: model_of[s] for s in keep}
    by_model = (emb.assign(__model__=[model_of[s] for s in emb.index])
                .groupby("__model__").mean())          # rows=models, cols=metrics
    # Display labels; keep raw on collision so rows stay distinct.
    pretty = {m: prettify(m) for m in by_model.index}
    dup = {k for k, c in Counter(pretty.values()).items() if c > 1}
    by_model.index = [m if pretty[m] in dup else pretty[m] for m in by_model.index]

    results = by_model.T                               # rows=metrics, cols=models
    results[_METRIC_TYPE] = mtype.reindex(results.index)
    has_batch = (mtype == "Batch correction").any()

    bm = Benchmarker.__new__(Benchmarker)              # bypass __init__ (no adata)
    bm._results = results
    bm._embedding_obsm_keys = [c for c in results.columns if c != _METRIC_TYPE]
    bm._bio_conservation_metrics = BioConservation()
    bm._batch_correction_metrics = BatchCorrection() if has_batch else None
    return bm, len(bm._embedding_obsm_keys)


def parse_args() -> argparse.Namespace:
    M = _load_metrics_module()
    p = argparse.ArgumentParser(description="scIB results-table figure from saved scib_table CSVs.")
    p.add_argument("--dataset", choices=M.DATASET_CHOICES, required=True)
    p.add_argument("--label-cols", type=str, help="Comma-separated subset (default: all saved tables).")
    p.add_argument("--min-max", action="store_true",
                   help="Min-max scale within this figure (panel-relative leaderboard view).")
    p.add_argument("--exclude", type=str, default=None,
                   help="Comma-separated model-name substrings to drop "
                        "(e.g. 'freeform,refine' for a benchmark-cohort-only table).")
    p.add_argument("--out-dir", type=Path, default=None)
    p.add_argument("--overwrite", action="store_true")
    args = p.parse_args()
    args._M = M
    return args


def main() -> None:
    args = parse_args()
    M = args._M
    dataset = args.dataset

    metrics_dir = REPO_ROOT / "outs" / "metrics" / dataset
    tables = sorted(metrics_dir.glob("scib_table_*.csv"))
    if not tables:
        raise FileNotFoundError(
            f"No scib_table_*.csv in {metrics_dir}. Run "
            f"scripts/embedding_metrics/metrics.py --dataset {dataset} first.")
    if args.label_cols:
        wanted = {c.strip() for c in args.label_cols.split(",") if c.strip()}
        tables = [t for t in tables if t.stem[len("scib_table_"):] in wanted]

    out_dir = args.out_dir or (REPO_ROOT / "outs" / "plots" / "metrics" / dataset)
    out_dir.mkdir(parents=True, exist_ok=True)
    suffix = "_minmax" if args.min_max else ""
    exclude = [t.strip() for t in args.exclude.split(",")] if args.exclude else None

    for table in tables:
        label_col = table.stem[len("scib_table_"):]
        bm, n_models = build_shim_benchmarker(table, dataset, M.parse_stem,
                                              M.prettify_model, exclude=exclude)
        dst_svg = out_dir / f"scib_results_table_{label_col}{suffix}.svg"
        if dst_svg.exists() and not args.overwrite:
            raise FileExistsError(f"{dst_svg} exists. Pass --overwrite.")
        bm.plot_results_table(min_max_scale=args.min_max, show=False, save_dir=str(out_dir))
        src_svg = out_dir / "scib_results.svg"
        if src_svg.exists():
            src_svg.replace(dst_svg)
        try:
            plt.gcf().savefig(out_dir / f"scib_results_table_{label_col}{suffix}.png",
                              dpi=200, bbox_inches="tight")
        except Exception as e:
            print(f"[warn] png render failed for {label_col}: {e}")
        plt.close("all")
        print(f"[ok] {label_col}: {n_models} model rows -> {dst_svg.name}")


if __name__ == "__main__":
    main()
