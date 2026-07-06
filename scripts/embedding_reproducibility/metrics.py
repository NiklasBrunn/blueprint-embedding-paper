"""
Embedding-level reproducibility across replicate seeds.

The blueprint-variability metrics (scripts/blueprint_variability_metrics/) measure
how much the *specification* drifts across runs. This stage measures reproducibility at the
representation level: do replicate runs produce the *same embedding* and the *same
conclusions*? Both are computed cell-by-cell (all seeds of a cohort share identical
cell_ids), so they work uniformly for blueprint, freeform, and conventional
baselines — enabling the key comparison "is re-authoring a blueprint any less
reproducible than re-training a stochastic neural method (scVI)?".

For every embedding cohort with >= 2 seeds under outs/embeddings/<dataset>/:

  cka           (Level 2, representation — HEADLINE): mean pairwise linear CKA
                (Kornblith et al. 2019) between runs' embeddings — agreement of the
                representational geometry up to rotation/isotropic scale, the right
                notion of "same embedding". Run-level (cluster) bootstrap CI.
  knn_overlap   (Level 2, representation — fine-grained secondary): for each seed
                pair, the mean over cells of the Jaccard between the two runs'
                k-nearest-neighbor sets. A deliberately harsh local metric (even
                scVI-vs-scVI scores low); reported for honesty, not as the headline.
  cluster_ari   (Level 3, conclusion — cross-run agreement): KMeans(k = #cell-type
                labels) on each run's embedding, then the mean pairwise adjusted Rand
                index between runs' cluster assignments — do runs group the same cells
                together. A pure between-run agreement metric.
  scib_ari/nmi  (Level 3, conclusion — quality consistency): the per-seed scIB
                ARI/NMI-vs-ground-truth (reused from outs/metrics/<dataset>/metrics.csv)
                summarised as mean ± SD across seeds. The SD is the reproducibility
                signal: do replicate runs achieve consistent biological quality. Uses
                the same scIB metrics as the rest of the pipeline, so it is directly
                comparable to the quality figures. (Note: low SD is necessary but not
                sufficient for reproducibility — two runs can score equally yet
                partition cells differently; cluster_ari is the agreement complement.)

Anchors:
  knn_chance_floor           = k / (n_cells - 1): kNN overlap expected by chance.
  spec_jaccard_random_floor  = mean Jaccard of two random gene draws of the median
                               blueprint gene-set size from the gene universe — the
                               floor for the spec-Jaccard bar in the companion figure.

Inputs:
    outs/embeddings/<dataset>/*.parquet              (cell_id + dim_* columns)
    data/preprocessed/<dataset>_pp.h5ad              (label column for #clusters; var)
    outs/embeddings/<dataset>/blueprint_*_seed*_summary.json  (gene-set sizes for the null)

Outputs:
    outs/embedding_reproducibility/<dataset>/reproducibility.csv   (one row per cohort)
    outs/embedding_reproducibility/<dataset>/pairs_long.csv        (per seed-pair)

Usage:
    python scripts/embedding_reproducibility/metrics.py --dataset human_pancreas
    python scripts/embedding_reproducibility/metrics.py --dataset human_pancreas --k 15
"""
from __future__ import annotations

import argparse
import json
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")
# Label column whose cardinality sets the cluster count k for cross-run ARI.
LABEL_COL = {"mouse_brain": "subclass", "mouse_cortex": "celltype_harmonized",
             "human_immune": "final_annotation"}  # default below is "celltype"
DEFAULT_K = 15
DEFAULT_BOOTSTRAP = 1000
NULL_REPS = 500
RNG_SEED = 0


# ============================================================
# Loading / discovery
# ============================================================
def parse_stem(stem: str, dataset: str) -> tuple[str, int | None]:
    """Recover (model, seed) from '<model>_<dataset>[_seed<N>]'. Mirrors
    scripts/embedding_metrics/metrics.py:parse_stem."""
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


def load_embedding(parquet_path: Path) -> tuple[np.ndarray, list[str]]:
    df = pq.read_table(parquet_path).to_pandas()
    dim_cols = sorted((c for c in df.columns if c.startswith("dim_")),
                      key=lambda n: int(n.split("_")[1]))
    if "cell_id" not in df.columns or not dim_cols:
        raise ValueError(f"{parquet_path} missing cell_id or dim_* columns.")
    X = df[dim_cols].to_numpy(dtype=np.float32)
    return X, df["cell_id"].astype(str).tolist()


def discover_cohorts(emb_dir: Path, dataset: str) -> dict[str, list[tuple[int, Path]]]:
    """model stem -> sorted list of (seed, parquet_path), only stems with >= 2 seeds."""
    by_stem: dict[str, list[tuple[int, Path]]] = defaultdict(list)
    for p in sorted(emb_dir.glob("*.parquet")):
        if p.name.endswith("_umap_coords.parquet"):
            continue
        # Exploratory local-model runs (OpenWebUI/Ollama/UFR) are excluded from
        # the aggregation for now; drop the guard to fold them back in.
        if "openwebui" in p.name:
            continue
        model, seed = parse_stem(p.stem, dataset)
        if seed is not None:
            by_stem[model].append((seed, p))
    return {m: sorted(v) for m, v in by_stem.items() if len(v) >= 2}


def cohort_kind(stem: str) -> str:
    if stem.startswith("blueprint_"):
        return "blueprint"
    if stem.startswith("freeform_"):
        return "freeform"
    return "baseline"


def align_to_reference(mats: list[np.ndarray], ids: list[list[str]]
                       ) -> list[np.ndarray]:
    """Reorder every run's rows to a common cell order (the sorted intersection of
    all runs' cell_ids). Fast-paths the common case where all runs already share an
    identical cell_id ordering."""
    ref = ids[0]
    if all(other == ref for other in ids[1:]):
        return mats
    common = sorted(set(ref).intersection(*[set(i) for i in ids[1:]]))
    out = []
    for X, cid in zip(mats, ids):
        pos = {c: i for i, c in enumerate(cid)}
        out.append(X[[pos[c] for c in common]])
    return out


# ============================================================
# Reproducibility statistics
# ============================================================
def linear_cka(X: np.ndarray, Y: np.ndarray) -> float:
    """Linear CKA (Kornblith et al. 2019) between two column-centered embeddings on
    the same cells: ||YᵀX||_F² / (||XᵀX||_F · ||YᵀY||_F). Invariant to orthogonal
    rotation and isotropic scale, so it measures whether two runs represent the same
    geometry regardless of arbitrary axis orientation. In [0, 1], 1 = identical."""
    X = X.astype(np.float64); Y = Y.astype(np.float64)
    X = X - X.mean(axis=0, keepdims=True)
    Y = Y - Y.mean(axis=0, keepdims=True)
    xtx = X.T @ X; yty = Y.T @ Y; ytx = Y.T @ X
    denom = np.sqrt((xtx ** 2).sum()) * np.sqrt((yty ** 2).sum())
    return float((ytx ** 2).sum() / denom) if denom > 0 else float("nan")


def knn_indices(X: np.ndarray, k: int) -> np.ndarray:
    """Self-excluded k-nearest-neighbor index matrix, shape (n_cells, k)."""
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="auto").fit(X)
    return nn.kneighbors(X, return_distance=False)[:, 1:]


def pair_knn_overlap(idx_a: np.ndarray, idx_b: np.ndarray) -> float:
    """Mean over cells of Jaccard(neighbors_a, neighbors_b). With fixed k,
    |A∪B| = 2k - |A∩B|, so only the per-cell intersection count is needed."""
    k = idx_a.shape[1]
    inter = (idx_a[:, :, None] == idx_b[:, None, :]).any(axis=2).sum(axis=1)
    return float((inter / (2 * k - inter)).mean())


def bootstrap_pairmean_ci(M: np.ndarray, rng: np.random.Generator,
                          n_boot: int = DEFAULT_BOOTSTRAP) -> tuple[float, float]:
    """Run-level (cluster) bootstrap of a mean-over-pairs statistic: resample the
    seeds, average M over the distinct-seed pairs in the resample. The unit is the
    run, not the (dependent) pair, so the CI is honest. M is a symmetric seed×seed
    matrix of a pairwise statistic (here CKA)."""
    n = M.shape[0]
    means = []
    for _ in range(n_boot):
        s = rng.integers(0, n, size=n)
        vals = [M[i, j] for a, i in enumerate(s) for j in s[a + 1:] if i != j]
        if vals:
            means.append(float(np.mean(vals)))
    if not means:
        return float("nan"), float("nan")
    return float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def random_gene_jaccard_floor(gene_universe: int, m: int,
                              rng: np.random.Generator, reps: int = NULL_REPS) -> float:
    """Expected Jaccard of two independent uniform random m-gene draws from a
    universe of `gene_universe` genes. The floor for the spec gene-Jaccard."""
    if gene_universe <= 0 or m <= 0 or m > gene_universe:
        return float("nan")
    vals = []
    for _ in range(reps):
        a = set(rng.choice(gene_universe, size=m, replace=False).tolist())
        b = set(rng.choice(gene_universe, size=m, replace=False).tolist())
        vals.append(len(a & b) / len(a | b))
    return float(np.mean(vals))


def scib_quality_dispersion(dataset: str, label_col: str) -> dict[str, dict]:
    """Per-cohort cross-seed mean ± SD of the scIB ARI/NMI-vs-truth already computed
    by scripts/embedding_metrics/metrics.py. Reused (not recomputed) so the Level-3
    quality-consistency numbers are identical to the quality figures. Keyed by the
    `model` stem, which matches our cohort stem."""
    path = REPO_ROOT / "outs" / "metrics" / dataset / "metrics.csv"
    if not path.is_file():
        return {}
    df = pd.read_csv(path)
    df = df[df["label_col"] == label_col]
    out: dict[str, dict] = {}
    for stem, sub in df.groupby("model"):
        a = sub["ari"].dropna().to_numpy()
        n = sub["nmi"].dropna().to_numpy()
        out[str(stem)] = {
            "scib_ari_mean": float(a.mean()) if a.size else float("nan"),
            "scib_ari_sd":   float(a.std(ddof=1)) if a.size > 1 else 0.0,
            "scib_nmi_mean": float(n.mean()) if n.size else float("nan"),
            "scib_nmi_sd":   float(n.std(ddof=1)) if n.size > 1 else 0.0,
        }
    return out


def blueprint_gene_set_sizes(emb_dir: Path) -> list[int]:
    """Union-of-module-genes size per blueprint run, from the Embed summary
    mirrors. Used to pick a representative m for the random-gene null."""
    sizes = []
    for sj in sorted(emb_dir.glob("blueprint_*_seed*_summary.json")):
        try:
            data = json.loads(sj.read_text())
        except Exception:
            continue
        gs: set[str] = set()
        for mod in (data.get("blueprint") or {}).get("modules") or []:
            gs.update(g for g in (mod.get("genes") or []) if isinstance(g, str))
        if gs:
            sizes.append(len(gs))
    return sizes


# ============================================================
# Driver
# ============================================================
def compute_cohort(stem: str, runs: list[tuple[int, Path]], n_clusters: int,
                   k: int, rng: np.random.Generator) -> tuple[dict, list[dict]]:
    seeds = [s for s, _ in runs]
    mats, ids = zip(*[load_embedding(p) for _, p in runs])
    mats = align_to_reference(list(mats), list(ids))
    n_cells = mats[0].shape[0]

    # Per-seed precompute: neighbor indices + KMeans cluster labels.
    print(f"  [{stem}] {len(seeds)} seeds, {n_cells} cells: kNN + KMeans...",
          flush=True)
    knn = [knn_indices(X, k) for X in mats]
    labels = [KMeans(n_clusters=n_clusters, random_state=RNG_SEED, n_init=10)
              .fit_predict(X) for X in mats]

    n = len(seeds)
    cka_M = np.full((n, n), np.nan)  # seed×seed CKA, for the run-level bootstrap
    pair_rows: list[dict] = []
    ckas, overlaps, aris = [], [], []
    for a, b in combinations(range(n), 2):
        cka = linear_cka(mats[a], mats[b])
        ov = pair_knn_overlap(knn[a], knn[b])
        ari = float(adjusted_rand_score(labels[a], labels[b]))
        cka_M[a, b] = cka_M[b, a] = cka
        ckas.append(cka)
        overlaps.append(ov)
        aris.append(ari)
        pair_rows.append({"cohort": stem, "seed_i": seeds[a], "seed_j": seeds[b],
                          "cka": cka, "knn_overlap": ov, "cluster_ari": ari})

    ci_lo, ci_hi = bootstrap_pairmean_ci(cka_M, rng)
    row = {
        "cohort": stem,
        "kind": cohort_kind(stem),
        "n_seeds": n,
        "n_cells": n_cells,
        "k": k,
        "cka_mean": float(np.mean(ckas)),
        "cka_ci_lo": ci_lo,
        "cka_ci_hi": ci_hi,
        "knn_overlap_mean": float(np.mean(overlaps)),
        "knn_chance_floor": k / (n_cells - 1),
        "cluster_ari_mean": float(np.mean(aris)),
    }
    return row, pair_rows


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Embedding-level reproducibility across replicate seeds.")
    ap.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    ap.add_argument("--k", type=int, default=DEFAULT_K,
                    help=f"kNN neighbors for the overlap metric (default {DEFAULT_K}).")
    ap.add_argument("--overwrite", action="store_true")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    emb_dir = REPO_ROOT / "outs" / "embeddings" / args.dataset
    h5ad = REPO_ROOT / "data" / "preprocessed" / f"{args.dataset}_pp.h5ad"
    out_dir = REPO_ROOT / "outs" / "embedding_reproducibility" / args.dataset
    out_csv = out_dir / "reproducibility.csv"
    if out_csv.exists() and not args.overwrite:
        raise SystemExit(f"{out_csv} exists. Pass --overwrite to replace it.")

    cohorts = discover_cohorts(emb_dir, args.dataset)
    if not cohorts:
        raise SystemExit(f"No multi-seed cohorts under {emb_dir}.")

    # Cluster count = #cell-type labels; gene universe for the null = #genes.
    label_col = LABEL_COL.get(args.dataset, "celltype")
    adata = ad.read_h5ad(h5ad, backed="r")
    if label_col not in adata.obs.columns:
        raise SystemExit(f"{h5ad} has no obs column {label_col!r}.")
    n_clusters = int(pd.Series(adata.obs[label_col].astype(str)).nunique())
    gene_universe = int(adata.n_vars)

    rng = np.random.default_rng(RNG_SEED)
    qual = scib_quality_dispersion(args.dataset, label_col)
    sizes = blueprint_gene_set_sizes(emb_dir)
    m_repr = int(np.median(sizes)) if sizes else 0
    spec_floor = random_gene_jaccard_floor(gene_universe, m_repr, rng)
    print(f"[info] {args.dataset}: {len(cohorts)} cohorts, k_clusters={n_clusters}, "
          f"gene_universe={gene_universe}, spec null m={m_repr} -> "
          f"floor={spec_floor:.4f}", flush=True)

    rows, pair_rows = [], []
    for stem in sorted(cohorts):
        row, pairs = compute_cohort(stem, cohorts[stem], n_clusters, args.k, rng)
        row["spec_jaccard_random_floor"] = spec_floor
        row["spec_null_gene_set_size"] = m_repr
        row.update(qual.get(stem, {"scib_ari_mean": float("nan"), "scib_ari_sd": float("nan"),
                                   "scib_nmi_mean": float("nan"), "scib_nmi_sd": float("nan")}))
        rows.append(row)
        pair_rows.extend(pairs)
        print(f"  -> CKA={row['cka_mean']:.3f} "
              f"[{row['cka_ci_lo']:.3f}, {row['cka_ci_hi']:.3f}] "
              f"knn_overlap={row['knn_overlap_mean']:.3f} "
              f"cluster_ari={row['cluster_ari_mean']:.3f}", flush=True)

    out_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).sort_values("cohort").to_csv(out_csv, index=False)
    pd.DataFrame(pair_rows).to_csv(out_dir / "pairs_long.csv", index=False)
    print(f"[ok] wrote {out_csv} ({len(rows)} cohorts)")


if __name__ == "__main__":
    main()
