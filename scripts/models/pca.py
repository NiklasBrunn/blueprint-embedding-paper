"""
PCA embedding for preprocessed single-cell data.

Linear baseline: PCA on the log1p-normalized matrix, restricted to HVGs
(var["highly_variable"]). Uses scikit-learn's PCA (n_comps=50, zero-centered)
so the *fit* (computing the principal-component loadings — the "model") and the
*transform* (projecting the cells onto them — constructing the embedding) are
two separate, separately-timed calls; the trained loadings are persisted under
models/<dataset>/pca/pca_model.npz (components + feature mean + variance ratios)
so the model footprint is non-trivial like the other parametric methods.

PCA is effectively deterministic on a fixed input (truncated SVD up to sign
flips on individual components, which do not affect any downstream distance-
or angle-based metric we compute). We therefore do not expose --seed; the
random_state below is hard-coded for reproducibility of the underlying
randomized SVD path. Embeddings differ from a previous run only by per-component
sign (and column-permutation in the pathological eigenvalue-degenerate case);
all metrics in scripts/embedding_metrics/metrics.py are invariant to that.

Inputs:
    data/preprocessed/<dataset>_pp.h5ad   (produced by scripts/data_preprocessing/preprocess.py)

Outputs:
    outs/embeddings/<dataset>/pca_<dataset>.parquet
    outs/embeddings/<dataset>/pca_<dataset>_report.md
    models/<dataset>/pca/pca_model.npz                   (trained PC loadings)
(Re-running just hits FileExistsError; delete the existing parquet to overwrite,
e.g. after a preprocessing change.)

Usage:
    python scripts/models/pca.py --dataset mouse_brain
    python scripts/models/pca.py --dataset human_pancreas
"""
from __future__ import annotations

import argparse
import platform
import time
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp
from sklearn.decomposition import PCA as SKLearnPCA

from _resource_meta import path_bytes as _path_bytes, resource_lines as _resource_lines
from _resource_meta import rss_mb as _rss_mb, time_split as _time_split


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")
MODEL_NAME = "pca"

# Hard-coded random_state for the (randomized) SVD path. PCA is effectively
# deterministic; this just pins the per-component sign convention.
PCA_RANDOM_STATE = 0

N_COMPS = 50
DEFAULT_PARAMS = {
    "mouse_brain":         {"n_comps": N_COMPS, "use_hvgs": True},
    "human_pancreas":      {"n_comps": N_COMPS, "use_hvgs": True},
    "mouse_gastrulation":  {"n_comps": N_COMPS, "use_hvgs": True},
    "mouse_cortex":        {"n_comps": N_COMPS, "use_hvgs": True},
    "human_immune":        {"n_comps": N_COMPS, "use_hvgs": True},
}


def fit(adata: ad.AnnData, params: dict, model_dir: Path) -> tuple[np.ndarray, list[float], dict]:
    """Fit the PCA model (compute loadings), persist it, then project the data.

    Returns (embedding, variance_ratio, meta) where meta carries obtain_model_s
    (time to compute + save the loadings — the "obtain the model" step) and
    build_embedding_s (time to project the cells onto them)."""
    if params["use_hvgs"]:
        if "highly_variable" not in adata.var:
            raise RuntimeError("var['highly_variable'] missing; rerun preprocessing.")
        X = adata[:, adata.var["highly_variable"].to_numpy()].X
    else:
        X = adata.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)  # sklearn PCA needs dense, zero-centered

    _t_fit0 = time.perf_counter()
    pca = SKLearnPCA(n_components=params["n_comps"], random_state=PCA_RANDOM_STATE)
    pca.fit(X)
    model_dir.mkdir(parents=True, exist_ok=True)
    np.savez(model_dir / "pca_model.npz",
             components=pca.components_.astype(np.float32),       # (n_comps, n_features) loadings
             mean=pca.mean_.astype(np.float32),                   # (n_features,) feature means
             explained_variance_ratio=np.asarray(pca.explained_variance_ratio_, dtype=np.float64))
    _t_fit1 = time.perf_counter()
    rss_after_obtain_mb = _rss_mb()

    emb = np.asarray(pca.transform(X), dtype=np.float32)
    _t_tr1 = time.perf_counter()

    meta = {
        "model_source":        "trained",
        "obtain_model_s":      float(_t_fit1 - _t_fit0),
        "build_embedding_s":   float(_t_tr1 - _t_fit1),
        "rss_after_obtain_mb": rss_after_obtain_mb,
        "model_dir":           str(model_dir),
    }
    return emb, np.asarray(pca.explained_variance_ratio_, dtype=float).tolist(), meta


def write_parquet(emb: np.ndarray, cell_ids: list[str], path: Path, kv: dict) -> None:
    cols = {f"dim_{i}": emb[:, i].astype(np.float32) for i in range(emb.shape[1])}
    table = pa.Table.from_pydict({"cell_id": cell_ids, **cols})
    md = {**(table.schema.metadata or {}),
          **{k.encode(): str(v).encode() for k, v in kv.items()}}
    table = table.replace_schema_metadata(md)
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, path, compression="zstd")


def write_report(report: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    p, s = report["params"], report["stats"]
    var_ratio = report["variance_ratio"]
    lines = [
        f"# PCA embedding report — {report['dataset']}",
        "",
        f"- Timestamp:  `{report['timestamp']}`",
        f"- Input:      `{report['input_path']}`",
        f"- Output:     `{report['output_path']}`",
        f"- Model:      `{report['model_dir']}/pca_model.npz`",
        "",
        "## Parameters",
        f"- `random_state`: {PCA_RANDOM_STATE} (hard-coded; PCA is deterministic)",
        *(f"- `{k}`: {v}" for k, v in p.items()),
        "",
        "## Stats",
        f"- cells: {s['n_cells']}",
        f"- input genes (after HVG mask): {s['n_input_genes']}",
        f"- embedding dim: {s['embedding_dim']}",
        *_resource_lines(s, storage_what="pca_model.npz",
                         obtain_what="compute PC loadings + save model",
                         build_what="project cells onto the loadings"),
        "",
        "## Notes",
        "- The persisted model is the PC loadings (genes × n_comps), the per-feature",
        "  means used for zero-centering, and the explained-variance ratios — enough to",
        "  project new cells onto the same axes. It is cheap to recompute, so this is a",
        "  convenience artifact, not a hard dependency.",
        "",
        "## Variance explained",
        f"- first component:        {var_ratio[0]:.4g}",
        f"- top 10 cumulative:      {sum(var_ratio[:10]):.4g}",
        f"- top {len(var_ratio)} cumulative:    {sum(var_ratio):.4g}",
    ]
    path.write_text("\n".join(lines) + "\n")


def resolve_params(args: argparse.Namespace) -> dict:
    p = DEFAULT_PARAMS[args.dataset].copy()
    if args.n_comps is not None:
        p["n_comps"] = args.n_comps
    return p


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="PCA embedding for preprocessed AnnData.")
    parser.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    parser.add_argument("--n-comps", type=int, help=f"Override n_comps (default {N_COMPS}).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = resolve_params(args)

    input_path  = REPO_ROOT / "data" / "preprocessed" / f"{args.dataset}_pp.h5ad"
    out_dir     = REPO_ROOT / "outs" / "embeddings" / args.dataset
    stem        = f"{MODEL_NAME}_{args.dataset}"
    output_path = out_dir / f"{stem}.parquet"
    report_path = out_dir / f"{stem}_report.md"
    model_dir   = REPO_ROOT / "models" / args.dataset / "pca"

    if output_path.exists():
        raise FileExistsError(
            f"{output_path} exists. Delete or rename it to overwrite."
        )

    # Start the clock before the disk read; walltime_s is the full
    # "file-on-disk → embedding-on-disk" cost.
    t0 = time.perf_counter()
    adata = ad.read_h5ad(input_path)
    cell_ids = adata.obs_names.astype(str).tolist()
    t_after_load = time.perf_counter()
    rss_after_load_mb = _rss_mb()

    emb, var_ratio, fit_meta = fit(adata, params, model_dir)
    t_after_build = time.perf_counter()
    rss_after_build_mb = _rss_mb()

    n_input_genes = (
        int(adata.var["highly_variable"].sum()) if params["use_hvgs"] else int(adata.n_vars)
    )
    model_storage_bytes = _path_bytes(model_dir)
    load_data_s = t_after_load - t0
    obtain_model_s = fit_meta.get("obtain_model_s")
    build_embedding_s = fit_meta.get("build_embedding_s")

    stats = {
        "n_cells":             int(adata.n_obs),
        "n_input_genes":       n_input_genes,
        "embedding_dim":       int(emb.shape[1]),
        "model_source":        fit_meta.get("model_source", "trained"),
        "load_data_s":         float(load_data_s),
        "obtain_model_s":      obtain_model_s,
        "build_embedding_s":   build_embedding_s,
        "model_storage_bytes": model_storage_bytes,
        "model_download_s":    0.0,    # PCA fits locally; nothing downloaded
        "rss_after_load_mb":   rss_after_load_mb,
        "rss_after_obtain_mb": fit_meta.get("rss_after_obtain_mb", rss_after_load_mb),
        "rss_after_build_mb":  rss_after_build_mb,
        "device":              "cpu",   # scikit-learn — no GPU path
        "host":                platform.node(),
    }
    kv = {
        "model": MODEL_NAME, "dataset": args.dataset,
        "embedding_dim": stats["embedding_dim"], "n_input_genes": stats["n_input_genes"],
        "random_state": PCA_RANDOM_STATE,
        "model_source": stats["model_source"], "model_storage_bytes": model_storage_bytes,
        "device": "cpu", "host": stats["host"],
    }
    write_parquet(emb, cell_ids, output_path, kv)
    walltime_s = time.perf_counter() - t0
    stats["walltime_s"] = float(walltime_s)
    stats["other_s"] = _time_split(walltime_s, load_data_s, obtain_model_s, build_embedding_s)
    stats["peak_rss_mb"] = _rss_mb()
    stats["embedding_bytes"] = _path_bytes(output_path)
    kv["walltime_s"] = round(walltime_s, 3)
    kv["peak_rss_mb"] = round(stats["peak_rss_mb"], 1)

    report = {
        "dataset": args.dataset,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "input_path": str(input_path),
        "output_path": str(output_path),
        "model_dir": str(model_dir),
        "params": params,
        "stats": stats,
        "variance_ratio": var_ratio,
    }
    write_report(report, report_path)

    print(f"[ok] {input_path.name}: ({stats['n_cells']}, {stats['embedding_dim']}) "
          f"PCA embedding in {walltime_s:.2f}s "
          f"(load {load_data_s:.2f}s + fit {obtain_model_s:.2f}s + transform {build_embedding_s:.2f}s)")
    print(f"[ok] wrote {output_path}")
    print(f"[ok] wrote {report_path}")
    print(f"[ok] saved model in {model_dir}")


if __name__ == "__main__":
    main()
