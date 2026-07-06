"""
Pretrained scVI embedding from CZ CELLxGENE Census.

Encoder-only inference against the Census-trained scVI checkpoint for the
appropriate organism. Symbol -> Ensembl mapping uses the Census `var` table
(canonical). The encoder forward pass is deterministic — no --seed is exposed
because the output is fully determined by (dataset, CENSUS_VERSION, organism).

Inputs:
    data/preprocessed/<dataset>_pp.h5ad

Outputs:
    outs/embeddings/<dataset>/scvi_census_<dataset>.parquet
    outs/embeddings/<dataset>/scvi_census_<dataset>_report.md
(Re-running just hits FileExistsError; delete the existing file to overwrite,
e.g. after a dependency bump or a CENSUS_VERSION change.)

Caches (shared across datasets/runs):
    models/cellxgene_census/<organism>/scvi_<census_version>/   (downloaded model)

Usage:
    python scripts/models/scvi_census.py --dataset mouse_brain
    python scripts/models/scvi_census.py --dataset human_pancreas
"""
from __future__ import annotations

# Import the shared resource-cost helper *before* the sys.path surgery below.
from _resource_meta import path_bytes as _path_bytes, resource_lines as _resource_lines
from _resource_meta import rss_mb as _rss_mb, time_split as _time_split

# IMPORTANT: scripts/models/scvi.py shadows the installed `scvi` package because
# Python prepends the script's own directory to sys.path. Strip that entry before
# any third-party imports so `import scvi` (in fit()) resolves to scvi-tools.
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:] = [p for p in _sys.path if _os.path.abspath(p) != _HERE]
del _os, _sys, _HERE

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


def _resolve_device(requested: str) -> str:
    """Resolve the requested compute device for the encoder forward pass,
    downgrading to 'cpu' (with a warning) if the asked-for backend isn't
    available. 'auto' -> 'cuda' if present else 'cpu' (Apple MPS is opt-in via
    --accelerator mps; for a single inference pass it tends to be slower than
    CPU on these small networks)."""
    req = (requested or "auto").lower()
    try:
        import torch
        if req in ("auto", "gpu", "cuda") and torch.cuda.is_available():
            return "cuda"
        if req == "mps" and torch.backends.mps.is_available():
            return "mps"
        if req in ("gpu", "cuda"):
            print("[warn] --accelerator cuda requested but no CUDA device; using 'cpu'.")
        elif req == "mps":
            print("[warn] --accelerator mps requested but MPS not available; using 'cpu'.")
    except Exception:
        if req not in ("auto", "cpu"):
            print(f"[warn] could not query torch for device {req!r}; using 'cpu'.")
    return "cpu"


def _round_to_integer_counts(X):
    """Round nonzero entries to nearest integer in place.

    No-op for true UMI data (integer-valued floats); only matters for SMART-Seq
    RSEM expected counts and similar fractional quantifications. Matches the
    preprocessing CZI applies before training Census scVI checkpoints.
    """
    if sp.issparse(X):
        X.data = np.rint(X.data).astype(X.data.dtype)
        return X
    return np.rint(X).astype(X.dtype)


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")
MODEL_NAME = "scvi_census"
CENSUS_VERSION = "2025-11-08"

DEFAULT_PARAMS = {
    "mouse_brain":         {"organism": "mus_musculus"},
    "human_pancreas":      {"organism": "homo_sapiens"},
    "mouse_gastrulation":  {"organism": "mus_musculus"},
    "mouse_cortex":        {"organism": "mus_musculus"},
    "human_immune":        {"organism": "homo_sapiens"},
}


def _ensure_model_cache(organism: str, cache: Path) -> None:
    """Download the Census scVI `model.pt` into `cache` if not already present.

    For releases >= 2024-07-01, CZI ships only the bare model.pt at
    `s3://cellxgene-contrib-public/models/scvi/<date>/<organism>/model.pt`
    (no `modelhub/` wrapper directory). scvi-tools 1.4+ reads everything it
    needs (state, registry, var_names) from model.pt itself.
    """
    import botocore
    import botocore.session
    import botocore.config

    target = cache / "model.pt"
    if target.exists() and target.stat().st_size > 0:
        return
    if target.exists():
        target.unlink()  # 0-byte leftover from a previous failed run

    cache.mkdir(parents=True, exist_ok=True)
    bucket = "cellxgene-contrib-public"
    key = f"models/scvi/{CENSUS_VERSION}/{organism}/model.pt"
    s3 = botocore.session.Session().create_client(
        "s3", config=botocore.config.Config(signature_version=botocore.UNSIGNED)
    )
    resp = s3.get_object(Bucket=bucket, Key=key)
    with target.open("wb") as f:
        for chunk in resp["Body"].iter_chunks(8 * 1024 * 1024):
            f.write(chunk)
    if not target.exists() or target.stat().st_size == 0:
        raise RuntimeError(
            f"Failed to download s3://{bucket}/{key}. "
            f"Check that CENSUS_VERSION='{CENSUS_VERSION}' has a published scVI model "
            "for this organism."
        )


def _symbol_to_ensembl(organism: str):
    """Build symbol → Ensembl ID map from the Census var table.

    Empirically the Census duplicates are dominated by ncRNA gene families
    (Y_RNA, Metazoa_SRP, U3/U6, SNORA*) where there's no canonical "right"
    Ensembl ID anyway. We keep the first occurrence — for proper protein-coding
    genes that have <5 paralogs the first row in Census order is typically the
    canonical one; for RNA-family symbols the choice is arbitrary by nature.
    """
    import cellxgene_census
    with cellxgene_census.open_soma(census_version=CENSUS_VERSION) as census:
        var_tbl = (
            census["census_data"][organism].ms["RNA"].var
            .read(column_names=["feature_id", "feature_name"]).concat().to_pandas()
        )
    return (var_tbl.drop_duplicates("feature_name", keep="first")
                   .set_index("feature_name")["feature_id"])


def fit(adata: ad.AnnData, params: dict, cache: Path) -> tuple[np.ndarray, dict]:
    import scvi
    import scvi.model  # noqa: F401  (force submodule registration)

    if "counts" not in adata.layers:
        raise RuntimeError("layers['counts'] missing; rerun preprocessing.")

    organism = params["organism"]
    _t_fit0 = time.perf_counter()
    _ensure_model_cache(organism, cache)
    sym2ens = _symbol_to_ensembl(organism)

    keep = adata.var_names.isin(sym2ens.index)
    n_input = int(adata.n_vars)
    n_dropped = int((~keep).sum())

    sub = adata[:, keep].copy()
    sub.X = sub.layers["counts"].copy()
    sub.X = _round_to_integer_counts(sub.X)
    sub.var["ensembl_id"] = sub.var_names.map(sym2ens).values
    sub.var.index = sub.var["ensembl_id"]
    sub.var.index.name = None
    sub.obs["batch"] = "unassigned"
    sub.obs["n_counts"] = np.asarray(sub.X.sum(axis=1)).ravel()

    scvi.model.SCVI.prepare_query_anndata(sub, str(cache))
    vae = scvi.model.SCVI.load_query_data(sub, str(cache))
    vae.is_trained = True
    device = _resolve_device(params.get("accelerator", "auto"))
    if device != "cpu":
        try:
            vae.to_device(device)
        except Exception as exc:
            print(f"[warn] vae.to_device({device!r}) failed ({exc}); staying on cpu.")
            device = "cpu"
    _t_fit1 = time.perf_counter()
    rss_after_obtain_mb = _rss_mb()
    emb = np.asarray(vae.get_latent_representation(), dtype=np.float32)
    _t_tr1 = time.perf_counter()
    try:
        device = str(next(vae.module.parameters()).device)
    except Exception:
        pass

    meta = {
        "census_version":  CENSUS_VERSION,
        "organism":        organism,
        "n_input_genes":   n_input,
        "n_mapped_genes":  int(keep.sum()),
        "n_dropped":       n_dropped,
        # The model is pretrained by others; "obtain model" here = (cold-cache
        # download, if any) + Census var-table query + query-anndata adaptation +
        # pretrained-checkpoint load. "build embedding" = the encoder forward pass.
        "model_source":        "loaded",
        "obtain_model_s":      float(_t_fit1 - _t_fit0),
        "build_embedding_s":   float(_t_tr1 - _t_fit1),
        "rss_after_obtain_mb": rss_after_obtain_mb,
        "accelerator_requested": params.get("accelerator", "auto"),
        "device":          device,
    }
    return emb, meta


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
    p, s, m = report["params"], report["stats"], report["model"]
    overlap_pct = 100.0 * m["n_mapped_genes"] / max(m["n_input_genes"], 1)
    lines = [
        f"# scVI Census embedding report — {report['dataset']}",
        "",
        f"- Timestamp: `{report['timestamp']}`",
        f"- Input:     `{report['input_path']}`",
        f"- Output:    `{report['output_path']}`",
        f"- Cache:     `{report['cache']}`",
        "",
        "## Parameters",
        *(f"- `{k}`: {v}" for k, v in p.items()),
        "",
        "## Model",
        f"- census_version: `{m['census_version']}`",
        f"- organism:       `{m['organism']}`",
        "",
        "## Stats",
        f"- cells: {s['n_cells']}",
        f"- input genes (preprocessed):           {m['n_input_genes']}",
        f"- mapped to model gene set (Ensembl):   {m['n_mapped_genes']} ({overlap_pct:.1f}%)",
        f"- dropped (unmapped or ambiguous symbol): {m['n_dropped']}",
        f"- embedding dim: {s['embedding_dim']}",
        *_resource_lines(s, storage_what="pretrained model.pt",
                         obtain_what="cold-cache download (if any) + Census var query + query-anndata adaptation + checkpoint load",
                         build_what="encoder forward pass"),
    ]
    path.write_text("\n".join(lines) + "\n")


def resolve_params(args: argparse.Namespace) -> dict:
    p = DEFAULT_PARAMS[args.dataset].copy()
    p["accelerator"] = args.accelerator or "auto"
    return p


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Pretrained scVI embedding from CZ CELLxGENE Census."
    )
    parser.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    parser.add_argument("--accelerator", choices=["auto", "cpu", "mps", "cuda", "gpu"],
                        default=None,
                        help="Device for the encoder forward pass (default: auto = CUDA if "
                             "present else CPU; 'mps' to try Apple Metal — usually slower than "
                             "CPU for a single inference pass on these small networks).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = resolve_params(args)

    input_path  = REPO_ROOT / "data" / "preprocessed" / f"{args.dataset}_pp.h5ad"
    out_dir     = REPO_ROOT / "outs" / "embeddings" / args.dataset
    stem        = f"{MODEL_NAME}_{args.dataset}"
    output_path = out_dir / f"{stem}.parquet"
    report_path = out_dir / f"{stem}_report.md"
    cache       = (REPO_ROOT / "models" / "cellxgene_census" /
                   params["organism"] / f"scvi_{CENSUS_VERSION}")

    if output_path.exists():
        raise FileExistsError(
            f"{output_path} exists. Delete or rename it to overwrite."
        )

    # Start the clock before the disk read; walltime_s = load_data + obtain_model
    # (= checkpoint load + gene-space adaptation; this is "loaded", not trained —
    # the model was pretrained by CZI) + build_embedding (= encoder forward) + other.
    t0 = time.perf_counter()
    adata = ad.read_h5ad(input_path)
    cell_ids = adata.obs_names.astype(str).tolist()
    t_after_load = time.perf_counter()
    rss_after_load_mb = _rss_mb()

    emb, model_meta = fit(adata, params, cache=cache)
    t_after_build = time.perf_counter()
    rss_after_build_mb = _rss_mb()

    model_storage_bytes = _path_bytes(cache / "model.pt")
    device = model_meta.get("device", "cpu")
    host = platform.node()
    load_data_s = t_after_load - t0
    obtain_model_s = model_meta.get("obtain_model_s")
    build_embedding_s = model_meta.get("build_embedding_s")
    stats = {
        "n_cells":             int(adata.n_obs),
        "embedding_dim":       int(emb.shape[1]),
        "model_source":        model_meta.get("model_source", "loaded"),
        "load_data_s":         float(load_data_s),
        "obtain_model_s":      obtain_model_s,
        "build_embedding_s":   build_embedding_s,
        "model_storage_bytes": model_storage_bytes,
        "model_download_s":    0.0,    # model.pt cached under models/cellxgene_census/...
        "rss_after_load_mb":   rss_after_load_mb,
        "rss_after_obtain_mb": model_meta.get("rss_after_obtain_mb", rss_after_load_mb),
        "rss_after_build_mb":  rss_after_build_mb,
        "device":              device,
        "host":                host,
    }
    kv = {
        "model": MODEL_NAME, "dataset": args.dataset, "embedding_dim": stats["embedding_dim"],
        "census_version": CENSUS_VERSION, "organism": params["organism"],
        "model_source": stats["model_source"], "model_storage_bytes": model_storage_bytes,
        "device": device, "host": host,
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
        "cache": str(cache),
        "params": params,
        "stats": stats,
        "model": model_meta,
    }
    write_report(report, report_path)

    print(f"[ok] {input_path.name}: ({stats['n_cells']}, {stats['embedding_dim']}) "
          f"Census scVI ({params['organism']} {CENSUS_VERSION}) in {walltime_s:.1f}s "
          f"(load {load_data_s:.1f}s + acquire {obtain_model_s:.1f}s + forward {build_embedding_s:.2f}s)")
    print(f"[ok] wrote {output_path}")
    print(f"[ok] wrote {report_path}")


if __name__ == "__main__":
    main()
