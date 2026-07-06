#!/usr/bin/env python3
"""Pseudotime-substrate comparison: how well does each embedding support a
trajectory algorithm at recovering developmental time?

The blueprint encodes depth = time *directly in its axis layout* (see
metrics.py). A black-box embedding (PCA, scVI) only yields time when a separate
trajectory algorithm is run on top of it. For a fair comparison, the *same*
diffusion pseudotime (scanpy DPT) is run on every embedding's kNN graph, rooted
at the earliest state (an Epiblast cell at E6.5), and scored on how well the
resulting pseudotime recovers (i) the physical embryonic stage globally and (ii) a
published reference gut-endoderm DPT (`obs['endo_gutDPT']`, Nowotschin) on the
endoderm cells it was derived for.

For each embedding parquet under outs/embeddings/<dataset>/:
  neighbors(use_rep=embedding) -> diffmap -> iroot=Epiblast@E6.5 medoid -> dpt,
  then Spearman(dpt, numeric stage) [excl. mixed_gastrulation] and
  Spearman(dpt, endo_gutDPT) on endoderm-trajectory cells.

Output:
    outs/temporal_metrics/<dataset>/pseudotime_substrate.csv

Usage:
    python scripts/temporal_metrics/pseudotime_substrate.py --dataset mouse_gastrulation
"""
from __future__ import annotations

import argparse
import importlib.util
import re
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow.parquet as pq
import scanpy as sc
from scipy.stats import spearmanr

REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_gastrulation",)
RANDOM_STATE = 0
N_NEIGHBORS = 15
N_COMPS = 15

STAGE_COL = {"mouse_gastrulation": "stage"}
CELLTYPE_COL = {"mouse_gastrulation": "celltype"}
ROOT_CELLTYPE = {"mouse_gastrulation": "Epiblast"}
ROOT_STAGE = {"mouse_gastrulation": "E6.5"}
REF_DPT_COL = {"mouse_gastrulation": "endo_gutDPT"}  # published reference pseudotime


def _load_module(name: str, relpath: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relpath)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


parse_numeric_stage = _load_module("emb_metrics",
                                   "scripts/embedding_metrics/metrics.py").parse_numeric_stage


def dataset_h5ad(dataset: str) -> Path:
    return REPO_ROOT / "data" / "preprocessed" / f"{dataset}_pp.h5ad"


def cohort_of(stem: str, dataset: str) -> str:
    s = re.sub(r"_seed\d+$", "", stem)
    suffix = f"_{dataset}"
    return s[:-len(suffix)] if s.endswith(suffix) else s


def choose_iroot(X: np.ndarray, is_root_pop: np.ndarray) -> int:
    """Index of the root-population cell closest to that population's centroid
    in the embedding (a stable medoid-approx; identical rule for every method)."""
    idx = np.flatnonzero(is_root_pop)
    centroid = X[idx].mean(axis=0)
    d = np.linalg.norm(X[idx] - centroid, axis=1)
    return int(idx[np.argmin(d)])


def run_dpt(X: np.ndarray, is_root_pop: np.ndarray) -> np.ndarray:
    adata = ad.AnnData(X=X.astype(np.float32))
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, use_rep="X",
                    random_state=RANDOM_STATE)
    sc.tl.diffmap(adata, n_comps=N_COMPS + 1, random_state=RANDOM_STATE)
    adata.uns["iroot"] = choose_iroot(X, is_root_pop)
    sc.tl.dpt(adata)
    return np.asarray(adata.obs["dpt_pseudotime"], dtype=float)


def score_embedding(parquet: Path, obs: pd.DataFrame, dataset: str) -> dict:
    df = pq.read_table(parquet).to_pandas()
    dim_cols = [c for c in df.columns if c.startswith("dim_")]
    X = df[dim_cols].to_numpy(dtype=np.float32)
    o = obs.reindex(df["cell_id"].to_numpy())

    stage_num = np.array([parse_numeric_stage(v) for v in o[STAGE_COL[dataset]].to_numpy()],
                         dtype=float)
    is_root = ((o[CELLTYPE_COL[dataset]].to_numpy() == ROOT_CELLTYPE[dataset])
               & (o[STAGE_COL[dataset]].astype(str).to_numpy() == ROOT_STAGE[dataset]))
    if not is_root.any():
        raise SystemExit(f"No {ROOT_CELLTYPE[dataset]}@{ROOT_STAGE[dataset]} root cells")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        pt = run_dpt(X, is_root)

    finite = np.isfinite(pt)
    stage_ok = finite & ~np.isnan(stage_num)
    rho_stage = float(spearmanr(pt[stage_ok], stage_num[stage_ok]).statistic)

    ref = pd.to_numeric(o[REF_DPT_COL[dataset]], errors="coerce").to_numpy()
    endo = finite & ~np.isnan(ref)
    rho_ref = (float(spearmanr(pt[endo], ref[endo]).statistic)
               if endo.sum() > 2 else float("nan"))

    return {
        "model": cohort_of(parquet.stem, dataset), "run": parquet.stem,
        "embedding_dim": len(dim_cols),
        "spearman_stage": rho_stage,        # signed: + means root->late ordering recovered
        "spearman_stage_abs": abs(rho_stage),
        "spearman_ref_endo": rho_ref,       # vs published gut-endoderm DPT
        "spearman_ref_endo_abs": abs(rho_ref) if not np.isnan(rho_ref) else float("nan"),
        "n_endo_cells": int(endo.sum()),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", choices=DATASET_CHOICES, default="mouse_gastrulation")
    p.add_argument("--input", type=Path, default=None,
                   help="Single parquet (overrides discovery).")
    args = p.parse_args()

    emb_dir = REPO_ROOT / "outs" / "embeddings" / args.dataset
    parquets = [args.input] if args.input else sorted(emb_dir.glob("*.parquet"))
    if not parquets:
        raise SystemExit(f"No parquets for {args.dataset}")

    obs = ad.read_h5ad(dataset_h5ad(args.dataset), backed="r").obs

    rows = []
    for pqf in parquets:
        print(f"[run] {pqf.name}", flush=True)
        r = score_embedding(pqf, obs, args.dataset)
        rows.append(r)
        print(f"      DPT-vs-stage rho={r['spearman_stage']:+.3f}  "
              f"DPT-vs-refEndo rho={r['spearman_ref_endo']:+.3f} "
              f"(endo n={r['n_endo_cells']})", flush=True)

    out_dir = REPO_ROOT / "outs" / "temporal_metrics" / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    out = out_dir / "pseudotime_substrate.csv"
    pd.DataFrame(rows).sort_values("model").to_csv(out, index=False)
    print(f"[ok] wrote {out} ({len(rows)} embeddings)")


if __name__ == "__main__":
    main()
