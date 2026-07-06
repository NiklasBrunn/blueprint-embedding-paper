"""
Embedding evaluation with the official scIB integration metric suite.

Scores every embedding parquet under outs/embeddings/<dataset>/ with the
`scib-metrics` Benchmarker (the scverse JAX reimplementation of Luecken et al.
2022), using the package defaults — *no deviation from the package*:

  Bio-conservation : isolated-label score, KMeans NMI, KMeans ARI, label
                     silhouette (ASW-label), cLISI.
  Batch-correction : batch-ASW (BRAS), iLISI, kBET, graph connectivity, PCR.
  Composite        : Total = 0.6*bio + 0.4*batch (the package aggregation),
                     read absolute via get_results(min_max_scale=False).

Neighbors are built by the Benchmarker via pynndescent; NMI/ARI use the
resolution-optimized Leiden (scib-metrics sweeps 10 resolutions),
not the Benchmarker's single-fixed-k kmeans default. Per the scib-metrics docs, values are internally
consistent but NOT directly comparable to the original `theislab/scib`
leaderboard — every method here is scored identically on our own datasets, so
the internal comparison is what matters.

PCR needs an *unintegrated* reference: we reuse the dataset's PCA baseline
(outs/embeddings/<dataset>/pca_<dataset>.parquet — PCA on HVG log1p, exactly the
scIB "Unintegrated" baseline) as `adata.obsm["Unintegrated"]`. If it is missing,
PCR is reported as NaN.

All scored embeddings of a dataset are placed in ONE AnnData (obsm keyed by
parquet stem) and benchmarked together per label column, so a single Benchmarker
covers every model — and `plot_results_table()` (see plot_scib_table.py) draws
the canonical scIB table across them.

Long-format CSV — one row per (parquet, label_col) — plus a per-embedding
markdown report. Aggregation across seeds is a downstream concern.

Inputs:
    outs/embeddings/<dataset>/*.parquet           (cell_id + dim_* columns)
    outs/embeddings/<dataset>/pca_<dataset>.parquet   (Unintegrated ref for PCR)
    data/preprocessed/<dataset>_pp.h5ad           (obs columns for labels/batch)

Outputs:
    outs/metrics/<dataset>/metrics.csv            (long format; rewritten on sweep)
    outs/metrics/<dataset>/scib_table_<label>.csv (raw Benchmarker get_results)
    outs/metrics/<dataset>/<stem>_metrics.md      (per-embedding markdown report)

Usage:
    python scripts/embedding_metrics/metrics.py --dataset mouse_brain
    python scripts/embedding_metrics/metrics.py --dataset human_pancreas --workers 8
    python scripts/embedding_metrics/metrics.py --dataset human_pancreas \
        --input outs/embeddings/human_pancreas/scvi_human_pancreas_seed0.parquet

References: Luecken et al. (2022) Nat. Methods 19:41-50,
doi:10.1038/s41592-021-01336-8 (scIB); implementation: scib-metrics,
https://github.com/YosefLab/scib-metrics. iLISI/cLISI: Korsunsky et al. (2019)
Harmony, doi:10.1038/s41592-019-0619-0. kBET: Büttner et al. (2019),
doi:10.1038/s41592-018-0254-1.
"""
from __future__ import annotations

import argparse
import re
import warnings
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq

# scib_metrics is imported lazily inside run_benchmarker() so the pure-Python
# helpers (parse_stem, participation_ratio, merge_csv, ...) and their unit tests
# don't require the heavy dep just to load this module.


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")

LABEL_COLS = {
    "mouse_brain":         ("primary_type", "subclass", "broad_type"),
    "human_pancreas":      ("celltype",),
    "mouse_gastrulation":  ("celltype", "stage"),
    "mouse_cortex":        ("celltype_harmonized", "class_harmonized"),
    "human_immune":        ("final_annotation",),
}
# mouse_cortex is cross-study/cross-protocol: `study` is the batch the embedding
# should mix (Tasic 2016 / Tasic 2018 / Zeisel 2018). human_immune is the scIB
# "Immune (human)" multi-study benchmark: `batch` is what the embedding should mix.
BATCH_KEY = {"mouse_brain": None, "human_pancreas": "tech", "mouse_gastrulation": "sample_id",
             "mouse_cortex": "study", "human_immune": "batch"}

DEFAULT_PARAMS = {
    "mouse_brain":         {"n_neighbors": 15},
    "human_pancreas":      {"n_neighbors": 15},
    "mouse_gastrulation":  {"n_neighbors": 15},
    "mouse_cortex":        {"n_neighbors": 15},
    "human_immune":        {"n_neighbors": 15},
}

CSV_COLUMNS = [
    "model", "dataset", "seed", "label_col",
    "n_cells", "n_cells_expected", "n_cells_dropped",
    "embedding_dim", "walltime_s", "participation_ratio",
    # scIB bio-conservation
    "ari", "nmi", "asw_label", "clisi", "isolated_label_f1",
    # scIB batch-correction (graph_connectivity is a BATCH metric in scIB)
    "asw_batch", "ilisi", "kbet_acceptance", "graph_connectivity", "pcr",
    # supplementary (not part of the scIB composite)
    "trajectory_corr",
    # derived scIB scores (straight from the Benchmarker aggregates)
    "bio_score", "batch_score_scib", "composite_scib",
]

# scib-metrics `get_results(clean_names=True)` column -> our CSV column.
RESULT_COL_MAP = {
    "Leiden NMI":         "nmi",   # resolution-optimized Leiden (our configured default)
    "Leiden ARI":         "ari",
    "KMeans NMI":         "nmi",   # kept so a kmeans-configured run still maps
    "KMeans ARI":         "ari",
    "Silhouette label":   "asw_label",
    "cLISI":              "clisi",
    "Isolated labels":    "isolated_label_f1",  # scIB silhouette isolated-label score
    "Batch ASW":          "asw_batch",          # alias if not clean-named to BRAS
    "BRAS":               "asw_batch",
    "iLISI":              "ilisi",
    "KBET":               "kbet_acceptance",
    "Graph connectivity": "graph_connectivity",
    "PCR comparison":     "pcr",
    "Bio conservation":   "bio_score",
    "Batch correction":   "batch_score_scib",
    "Total":              "composite_scib",
}
# Per-metric (non-aggregate) scIB columns we expect in a row, for NaN-filling.
SCIB_METRIC_COLS = ["ari", "nmi", "asw_label", "clisi", "isolated_label_f1",
                    "asw_batch", "ilisi", "kbet_acceptance", "graph_connectivity", "pcr"]

# Datasets with an ordinal-time obs column the embedding can be scored against.
# trajectory_corr = max_i |Spearman(embedding_axis_i, numeric_stage)|; NaN
# elsewhere. Numeric stage is parsed from the value via parse_numeric_stage.
TRAJECTORY_OBS_COL = {
    "mouse_brain":         None,
    "human_pancreas":      None,
    "mouse_gastrulation":  "stage",
}


# --------------------------------------------------------------------
# Pure-Python helpers (no heavy deps)
# --------------------------------------------------------------------
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


# Display names for the `model` column. Baselines map exactly; the
# agentic cohorts (blueprint_<engine>_[noref_]{unrefined|refine}_<effort>,
# freeform_<variant>_<engine>_<effort>) are compacted to e.g. "BP GPT-5.5 N U".
# Grounding/refinement flags: noref->N, unrefined->U, refine->R.
_BASELINE_PRETTY = {
    "pca": "PCA", "scvi": "scVI", "scvi_census": "scVI-census",
    "scvi_nobatch": "scVI (no batch)", "uce": "UCE", "ucell": "UCell", "expimap": "expiMap",
}
_ENGINE_PRETTY = {"gpt_5_5": "GPT-5.5", "gpt_5": "GPT-5"}
_VARIANT_PRETTY = {"constrained": "constrained", "naive": "naïve", "freeform": "freeform"}
_REFSTATE_PRETTY = {"unrefined": "U", "refine": "R"}
_EFFORT = ("low", "medium", "high", "xhigh")


def prettify_model(name: str) -> str:
    """Compact display label for a `model` value. Unknown names pass through."""
    if name in _BASELINE_PRETTY:
        return _BASELINE_PRETTY[name]
    if name.startswith("blueprint_"):
        family, s = "BP", name[len("blueprint_"):]
    elif name.startswith("freeform_"):
        family, s = "FF", name[len("freeform_"):]
    else:
        return name
    eng = ""
    for k, v in _ENGINE_PRETTY.items():
        if k in s:
            eng = v
            s = s.replace(k, "_").strip("_")
            break
    parts = [p for p in s.split("_") if p]
    effort = next((p for p in parts if p in _EFFORT), "")
    flags, variant = [], []
    for p in parts:
        if p in _EFFORT:
            continue
        if p in _REFSTATE_PRETTY:
            flags.append(_REFSTATE_PRETTY[p])
        elif p == "noref":
            flags.append("N")
        elif p == "userctx":
            flags.append("userctx")
        else:
            variant.append(_VARIANT_PRETTY.get(p, p))
    label = [family] + variant + ([eng] if eng else []) + flags
    if effort and effort != "high":  # all current runs are "high"; show others
        label.append(effort)
    return " ".join(label)


def discover_embeddings(emb_dir: Path, override: Path | None) -> list[Path]:
    if override is not None:
        if not override.exists():
            raise FileNotFoundError(f"--input not found: {override}")
        return [override.resolve()]
    if not emb_dir.exists():
        raise FileNotFoundError(f"Embeddings directory not found: {emb_dir}")
    parquets = [
        p for p in sorted(emb_dir.glob("*.parquet"))
        if not p.name.endswith("_umap_coords.parquet")
        # Exploratory local-model runs (OpenWebUI/Ollama/UFR) are excluded from
        # the aggregation for now; drop the guard to fold them back in.
        and "openwebui" not in p.name
    ]
    if not parquets:
        raise FileNotFoundError(f"No embedding parquets in {emb_dir}")
    return parquets


def load_embedding(parquet_path: Path) -> tuple[np.ndarray, list[str]]:
    table = pq.read_table(parquet_path)
    df = table.to_pandas()
    if "cell_id" not in df.columns:
        raise ValueError(f"{parquet_path} has no 'cell_id' column.")
    dim_cols = sorted(
        (c for c in df.columns if c.startswith("dim_")),
        key=lambda n: int(n.split("_")[1]),
    )
    if not dim_cols:
        raise ValueError(f"{parquet_path} has no dim_* columns.")
    X = df[dim_cols].to_numpy(dtype=np.float32)
    cell_ids = df["cell_id"].astype(str).tolist()
    return X, cell_ids


def load_obs(h5ad_path: Path, label_cols: list[str], batch_key: str | None,
             extra_cols: list[str] | None = None) -> pd.DataFrame:
    if not h5ad_path.exists():
        raise FileNotFoundError(f"Labels h5ad not found: {h5ad_path}")
    adata = ad.read_h5ad(h5ad_path, backed="r")
    needed = list(label_cols) + ([batch_key] if batch_key else []) + list(extra_cols or [])
    needed = list(dict.fromkeys(needed))  # dedupe preserving order
    missing = [c for c in needed if c not in adata.obs.columns]
    if missing:
        raise ValueError(
            f"{h5ad_path} is missing obs columns {missing}. "
            f"Available: {list(adata.obs.columns)}"
        )
    return pd.DataFrame(
        {c: adata.obs[c].astype(str).to_numpy() for c in needed},
        index=pd.Index(adata.obs_names.astype(str), name="cell_id"),
    )


def participation_ratio(X: np.ndarray) -> float:
    """Effective dimensionality: (sum sigma^2)^2 / sum sigma^4 of centered X."""
    X = X - X.mean(axis=0, keepdims=True)
    s = np.linalg.svd(X, compute_uv=False)
    s2 = (s ** 2).astype(np.float64)
    denom = float((s2 ** 2).sum())
    return float((s2.sum() ** 2) / denom) if denom > 0 else float("nan")


def parse_numeric_stage(value: object) -> float:
    """Map an ordinal-stage label (e.g. 'E6.5', 'embryonic day 7.25') to a float.

    Returns NaN if no numeric component can be parsed."""
    s = str(value)
    m = re.search(r"E(\d+(?:\.\d+)?)", s)
    if m is None:
        m = re.search(r"embryonic day\s+(\d+(?:\.\d+)?)", s, re.IGNORECASE)
    if m is None:
        m = re.search(r"(\d+(?:\.\d+)?)", s)
    return float(m.group(1)) if m else float("nan")


def trajectory_axis_spearman(X: np.ndarray, stage_labels: np.ndarray) -> float:
    """Maximum |Spearman corr| between any single embedding axis and the
    numeric ordinal-stage label. Returns NaN when fewer than 80% of stage
    labels parse to a number, or when the parsed series has zero variance."""
    parsed = np.array([parse_numeric_stage(v) for v in stage_labels], dtype=float)
    valid = ~np.isnan(parsed)
    if valid.mean() < 0.8 or np.unique(parsed[valid]).size < 2:
        return float("nan")
    from scipy.stats import spearmanr
    best = 0.0
    for i in range(X.shape[1]):
        rho, _ = spearmanr(X[valid, i], parsed[valid])
        if not np.isnan(rho):
            best = max(best, abs(float(rho)))
    return best


def read_sibling_report(parquet_path: Path) -> dict:
    """Best-effort regex extract of `walltime: <float> s` from the embedding's report."""
    md = parquet_path.with_name(parquet_path.stem + "_report.md")
    out: dict[str, float | str | None] = {"walltime_s": float("nan"), "report_path": None}
    if not md.exists():
        return out
    out["report_path"] = str(md)
    try:
        match = re.search(r"walltime:\s*([\d.]+)\s*s", md.read_text())
        if match:
            out["walltime_s"] = float(match.group(1))
    except Exception as e:
        print(f"[warn] could not parse walltime from {md}: {e}")
    return out


def merge_csv(csv_path: Path, new_rows: list[dict]) -> pd.DataFrame:
    """Replace existing rows that share (model, seed, label_col) with new_rows; append the rest."""
    new_df = pd.DataFrame(new_rows, columns=CSV_COLUMNS)
    if not csv_path.exists() or new_df.empty:
        return new_df
    existing = pd.read_csv(csv_path)
    mask = pd.Series(False, index=existing.index)
    for model, seed, label in {(r["model"], r.get("seed"), r["label_col"]) for r in new_rows}:
        seed_match = existing["seed"].isna() if pd.isna(seed) else existing["seed"] == seed
        mask |= (existing["model"] == model) & seed_match & (existing["label_col"] == label)
    return pd.concat([existing[~mask], new_df], ignore_index=True)


def write_csv(rows: list[dict], path: Path, overwrite: bool, full_sweep: bool) -> None:
    """Full-sweep mode rewrites the table; single-input mode merges in place."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if full_sweep:
        if path.exists() and not overwrite:
            raise FileExistsError(f"{path} exists. Pass --overwrite to replace it.")
        df = pd.DataFrame(rows, columns=CSV_COLUMNS)
    else:
        df = merge_csv(path, rows)
    df = df.sort_values(["model", "seed", "label_col"], na_position="last").reset_index(drop=True)
    df.to_csv(path, index=False)


# --------------------------------------------------------------------
# Dataset assembly
# --------------------------------------------------------------------
def align_to(ids: list[str], common: list[str]) -> np.ndarray:
    """Index array reordering a run's rows (given its `ids`) to `common` order."""
    pos = {c: i for i, c in enumerate(ids)}
    return np.array([pos[c] for c in common], dtype=int)


def build_obsm(parquets: list[Path], obs_df: pd.DataFrame, unint_path: Path | None,
               trajectory_col: str | None) -> tuple[dict, dict, list[str], np.ndarray | None, list[str]]:
    """Load every parquet's embedding, intersect cell_ids with obs (and the
    Unintegrated reference), and return:

      obsm    : {stem: matrix aligned to `common`} for the scored embeddings
      info    : {stem: {model, seed, n_total, n_cells, n_dropped, embedding_dim,
                        walltime_s, participation_ratio, trajectory_corr}}
      common  : the shared cell_id order (str)
      unint   : Unintegrated reference matrix aligned to `common`, or None
      order   : stem order (sorted) for stable obsm key iteration
    """
    raw: dict[str, tuple[np.ndarray, list[str]]] = {}
    for p in parquets:
        X, ids = load_embedding(p)
        raw[p.stem] = (X, ids)

    unint_ids = None
    if unint_path is not None and unint_path.exists():
        uX, unint_ids = load_embedding(unint_path)

    # Common cell order = obs ∩ every embedding (∩ Unintegrated when present).
    sets = [set(ids) for _, ids in raw.values()]
    common_set = set(obs_df.index).intersection(*sets) if sets else set()
    if unint_ids is not None:
        common_set &= set(unint_ids)
    common = sorted(common_set)
    if not common:
        raise ValueError("No cell_id overlap between embeddings and obs.")

    obsm: dict[str, np.ndarray] = {}
    info: dict[str, dict] = {}
    for p in parquets:
        X, ids = raw[p.stem]
        idx = align_to(ids, common)
        Xc = X[idx]
        obsm[p.stem] = Xc
        sib = read_sibling_report(p)
        info[p.stem] = {
            "n_total": len(ids),
            "n_cells": len(common),
            "n_dropped": len(ids) - len(common),
            "embedding_dim": int(X.shape[1]),
            "walltime_s": sib["walltime_s"],
            "participation_ratio": participation_ratio(Xc),
        }
        if trajectory_col and trajectory_col in obs_df.columns:
            info[p.stem]["trajectory_corr"] = trajectory_axis_spearman(
                Xc, obs_df.loc[common, trajectory_col].to_numpy())
        else:
            info[p.stem]["trajectory_corr"] = float("nan")

    unint = None
    if unint_ids is not None:
        unint = uX[align_to(unint_ids, common)]

    order = sorted(obsm)
    return obsm, info, common, unint, order


# --------------------------------------------------------------------
# Benchmarker driver
# --------------------------------------------------------------------
def run_benchmarker(obsm: dict, common: list[str], obs_sub: pd.DataFrame,
                    label_col: str, batch_key: str | None, unint: np.ndarray | None,
                    n_neighbors: int, n_jobs: int):
    """Build, prepare, and run one scib-metrics Benchmarker over all embeddings
    for a single label column; return the (benchmarked) Benchmarker. Callers get
    the table via `bm.get_results(...)` or `bm.plot_results_table(...)`."""
    from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

    n = len(common)
    a = ad.AnnData(X=np.zeros((n, 1), dtype=np.float32),
                   obs=obs_sub.loc[common].copy())
    for k, M in obsm.items():
        a.obsm[k] = M
    pre_key = None
    if unint is not None:
        a.obsm["Unintegrated"] = unint
        pre_key = "Unintegrated"

    has_batch = batch_key is not None
    batch_metrics = BatchCorrection() if has_batch else None
    # When there is no batch we cannot supply a batch_key; reuse the label as a
    # dummy key (unused because batch metrics are disabled).
    bkey = batch_key if has_batch else label_col

    # Leiden NMI/ARI (kmeans default under-clusters and deflates ARI/NMI); rationale in module docstring.
    bio_metrics = BioConservation(nmi_ari_cluster_labels_leiden=True,
                                  nmi_ari_cluster_labels_kmeans=False)
    bm = Benchmarker(
        a, batch_key=bkey, label_key=label_col,
        embedding_obsm_keys=list(obsm.keys()),
        bio_conservation_metrics=bio_metrics,
        batch_correction_metrics=batch_metrics,
        pre_integrated_embedding_obsm_key=pre_key,
        n_jobs=n_jobs, progress_bar=False,
    )
    bm.prepare()
    bm.benchmark()
    return bm


def results_to_rows(results: pd.DataFrame, dataset: str, label_col: str,
                    info: dict, has_batch: bool) -> list[dict]:
    """Map a Benchmarker results frame to long-format CSV rows (one per stem)."""
    rows: list[dict] = []
    scored = [ix for ix in results.index if ix != "Metric Type"]
    for stem in scored:
        model, seed = parse_stem(stem, dataset)
        meta = info[stem]
        row = {
            "_stem": stem,
            "model": model, "dataset": dataset, "seed": seed, "label_col": label_col,
            "n_cells": meta["n_cells"],
            "n_cells_expected": meta["n_total"],
            "n_cells_dropped": meta["n_dropped"],
            "embedding_dim": meta["embedding_dim"],
            "walltime_s": meta["walltime_s"],
            "participation_ratio": meta["participation_ratio"],
            "trajectory_corr": meta["trajectory_corr"],
        }
        for col in SCIB_METRIC_COLS + ["bio_score", "batch_score_scib", "composite_scib"]:
            row[col] = float("nan")
        for src, dst in RESULT_COL_MAP.items():
            if src in results.columns:
                val = results.at[stem, src]
                try:
                    row[dst] = float(val)
                except (TypeError, ValueError):
                    row[dst] = float("nan")
        # Bio-only (batchless) runs have no Total: composite falls back to bio.
        if not has_batch or np.isnan(row.get("composite_scib", float("nan"))):
            if not has_batch:
                row["composite_scib"] = row["bio_score"]
        rows.append(row)
    return rows


# --------------------------------------------------------------------
# Reporting
# --------------------------------------------------------------------
def _fmt(v) -> str:
    if isinstance(v, str):
        return v
    if v is None or (isinstance(v, float) and np.isnan(v)):
        return "—"
    return f"{v:.3f}"


def write_report(rows: list[dict], stem: str, dataset: str, out_dir: Path,
                 overwrite: bool) -> None:
    path = out_dir / f"{stem}_metrics.md"
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists. Pass --overwrite to replace it.")
    r0 = rows[0]
    lines = [
        f"# Embedding metrics — {dataset} / {r0['model']}"
        + (f" (seed={r0['seed']})" if r0["seed"] is not None else ""),
        "",
        f"- Timestamp:   `{datetime.now().isoformat(timespec='seconds')}`",
        f"- Engine:      `scib-metrics Benchmarker` (Leiden NMI/ARI [res-optimized], cLISI, PCR; pynndescent neighbors)",
        f"- cells:       {r0['n_cells']} (dropped {r0['n_cells_dropped']})",
        f"- embedding_dim: {r0['embedding_dim']}",
        f"- participation_ratio: {_fmt(r0['participation_ratio'])}",
        "",
        "## Metrics (scIB, scib-metrics implementation)",
        "",
        "| label | ari | nmi | asw_label | cLISI | iso | asw_batch | iLISI | kBET | graph_conn | PCR | traj | bio | batch | composite |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        cells = [r["label_col"], r["ari"], r["nmi"], r["asw_label"], r["clisi"],
                 r["isolated_label_f1"], r["asw_batch"], r["ilisi"], r["kbet_acceptance"],
                 r["graph_connectivity"], r["pcr"], r["trajectory_corr"],
                 r["bio_score"], r["batch_score_scib"], r["composite_scib"]]
        lines.append("| " + " | ".join(_fmt(c) for c in cells) + " |")
    lines += ["", "## Notes",
              "- scib-metrics reimplementation of Luecken et al. (2022); internally "
              "consistent but not directly comparable to the original scib leaderboard.",
              f"- batch metrics: {'computed' if not np.isnan(r0['batch_score_scib']) else 'NaN (no batch in dataset)'}"]
    path.write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------
# Main
# --------------------------------------------------------------------
def resolve_params(args: argparse.Namespace) -> dict:
    p = DEFAULT_PARAMS[args.dataset].copy()
    if args.n_neighbors is not None:
        p["n_neighbors"] = args.n_neighbors
    if args.label_cols is not None:
        p["label_cols"] = tuple(c.strip() for c in args.label_cols.split(",") if c.strip())
    else:
        p["label_cols"] = LABEL_COLS[args.dataset]
    p["batch_key"] = BATCH_KEY[args.dataset]
    return p


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="scIB metrics (scib-metrics Benchmarker) for embedding parquets.")
    parser.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    parser.add_argument("--input", type=Path, help="Single parquet (overrides discovery).")
    parser.add_argument("--label-cols", type=str, help="Comma-separated obs columns.")
    parser.add_argument("--n-neighbors", type=int)
    parser.add_argument("--seed", type=int, default=0,
                        help="Accepted for CLI compatibility; the scoring engine uses a fixed internal state.")
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--skip-existing", action="store_true",
                        help="Skip parquets whose <stem>_metrics.md already exists.")
    parser.add_argument("--workers", type=int, default=1,
                        help="n_jobs for the Benchmarker (default: 1).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = resolve_params(args)

    embeddings_dir = REPO_ROOT / "outs" / "embeddings" / args.dataset
    out_dir        = REPO_ROOT / "outs" / "metrics" / args.dataset
    csv_path       = out_dir / "metrics.csv"
    obs_h5ad       = REPO_ROOT / "data" / "preprocessed" / f"{args.dataset}_pp.h5ad"
    unint_path     = embeddings_dir / f"pca_{args.dataset}.parquet"

    parquets = discover_embeddings(embeddings_dir, args.input)
    if args.skip_existing:
        before = len(parquets)
        parquets = [p for p in parquets if not (out_dir / f"{p.stem}_metrics.md").exists()]
        print(f"[ok] skip-existing: {before - len(parquets)}/{before} already have reports; "
              f"processing {len(parquets)}")
        if not parquets:
            print("[ok] nothing to do.")
            return

    trajectory_col = TRAJECTORY_OBS_COL.get(args.dataset)
    obs_df = load_obs(obs_h5ad, list(params["label_cols"]), params["batch_key"],
                      extra_cols=[trajectory_col] if trajectory_col else None)
    obs_df.attrs["dataset"] = args.dataset

    has_unint = unint_path.exists()
    print(f"[ok] dataset={args.dataset} files={len(parquets)} labels={list(params['label_cols'])} "
          f"batch={params['batch_key']} trajectory={trajectory_col or '—'} "
          f"unintegrated_ref={'pca' if has_unint else 'MISSING (PCR=NaN)'} n_jobs={args.workers}")
    if not has_unint:
        print(f"[warn] {unint_path.name} not found — PCR will be NaN for all embeddings.")

    obsm, info, common, unint, _order = build_obsm(
        parquets, obs_df, unint_path if has_unint else None, trajectory_col)
    # Backfill (model, seed) now that we know the dataset (build_obsm placeholder).
    for stem in obsm:
        model, seed = parse_stem(stem, args.dataset)
        info[stem]["model"], info[stem]["seed"] = model, seed
    print(f"[ok] common cells: {len(common)}; embeddings: {len(obsm)}")

    has_batch = params["batch_key"] is not None
    all_rows: list[dict] = []
    per_stem_rows: dict[str, list[dict]] = {stem: [] for stem in obsm}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")  # silence package FutureWarning / anndata noise
        for label_col in params["label_cols"]:
            obs_sub = obs_df[[c for c in ([label_col] + ([params["batch_key"]] if has_batch else []))]]
            print(f"[run] benchmarking label_col={label_col} ...", flush=True)
            bm = run_benchmarker(obsm, common, obs_sub, label_col, params["batch_key"],
                                 unint, params["n_neighbors"], args.workers)
            results = bm.get_results(min_max_scale=False, clean_names=True)
            results.to_csv(out_dir / f"scib_table_{label_col}.csv")
            rows = results_to_rows(results, args.dataset, label_col, info, has_batch)
            for r in rows:
                per_stem_rows.setdefault(r["_stem"], []).append(r)
            all_rows.extend(rows)
            print(f"[ok] label_col={label_col}: {len(rows)} rows")

    out_dir.mkdir(parents=True, exist_ok=True)
    for stem, rows in per_stem_rows.items():
        if rows:
            write_report(rows, stem, args.dataset, out_dir, overwrite=True)

    full_sweep = args.input is None and not args.skip_existing
    write_csv(all_rows, csv_path, overwrite=args.overwrite, full_sweep=full_sweep)
    print(f"[ok] wrote {csv_path} ({len(all_rows)} rows)")


if __name__ == "__main__":
    main()
