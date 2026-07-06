"""
scVI embedding (de novo) for preprocessed single-cell data.

Trains scVI from scratch on raw counts (layers["counts"]) restricted to HVGs.
scvi.model.SCVI defaults (n_latent=10, n_layers=1, n_hidden=128); model.train()
is called with no args so scvi-tools picks max_epochs heuristically. The seed
controls weight init + minibatch order (different seed -> different embedding).

Inputs:
    data/preprocessed/<dataset>_pp.h5ad

Outputs (per run):
    outs/embeddings/<dataset>/<model>_<dataset>_<run_tag>.parquet
    outs/embeddings/<dataset>/<model>_<dataset>_<run_tag>_report.md
    models/<dataset>/<model>/<run_tag>/                  (trained model checkpoint)

<model> is "scvi" by default, or "scvi_nobatch" when `--batch-key none` overrides
a dataset whose default batch_key is non-null (i.e., human_pancreas). This keeps
the no-batch run as a distinct cohort that `scripts/embedding_metrics/metrics.py`
parses as a separate method (via parse_stem on `_<dataset>`).

Different seeds produce distinct files (...seed{N}.parquet). Same-seed re-runs
will fail with FileExistsError so historical results aren't silently clobbered
— delete or rename the existing file if you really mean to overwrite.

Usage:
    python scripts/models/scvi.py --dataset mouse_brain
    python scripts/models/scvi.py --dataset human_pancreas --seed 1
    python scripts/models/scvi.py --dataset human_pancreas --batch-key none --seed 0
"""
from __future__ import annotations

# Import the shared resource-cost helper *before* the sys.path surgery below
# (it lives in this dir, so it needs the script dir on sys.path).
from _resource_meta import path_bytes as _path_bytes, resource_lines as _resource_lines
from _resource_meta import rss_mb as _rss_mb, time_split as _time_split

# IMPORTANT: this script's own filename (scvi.py) shadows the installed `scvi`
# package because Python prepends the script's directory to sys.path. Strip that
# entry before any third-party imports so `import scvi` (in fit()) resolves to
# scvi-tools.
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


def _resolve_accelerator(requested: str) -> str:
    """Validate the requested Lightning accelerator, downgrading to 'auto' (and
    warning) if the asked-for backend isn't available on this machine.

    'auto' lets scvi-tools / Lightning pick: CUDA when present, else CPU
    (scvi-tools 1.x does NOT auto-select Apple's MPS — pass --accelerator mps to
    try it; on these small networks it's a marginal-at-best speedup, see README/
    repo notes). Returns the accelerator string to hand to model.train()."""
    req = (requested or "auto").lower()
    if req in ("gpu", "cuda"):
        try:
            import torch
            if torch.cuda.is_available():
                return "cuda"
        except Exception:
            pass
        print("[warn] --accelerator cuda requested but no CUDA device; using 'auto'.")
        return "auto"
    if req == "mps":
        try:
            import torch
            if torch.backends.mps.is_available():
                return "mps"
        except Exception:
            pass
        print("[warn] --accelerator mps requested but MPS not available; using 'auto'.")
        return "auto"
    return req if req in ("auto", "cpu") else "auto"


def _round_to_integer_counts(X):
    """Round nonzero entries to nearest integer in place.

    No-op for true UMI data (integer-valued floats); only matters for SMART-Seq
    RSEM expected counts and similar fractional quantifications. Aligns the
    input with scVI's Negative Binomial likelihood, which assumes count data.
    """
    if sp.issparse(X):
        X.data = np.rint(X.data).astype(X.data.dtype)
        return X
    return np.rint(X).astype(X.dtype)


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")
MODEL_NAME = "scvi"

N_LATENT = 10
DEFAULT_PARAMS = {
    "mouse_brain":         {"batch_key": None,        "n_latent": N_LATENT},
    "human_pancreas":      {"batch_key": "tech",      "n_latent": N_LATENT},
    "mouse_gastrulation":  {"batch_key": "sample_id", "n_latent": N_LATENT},
    # mouse_cortex is the three-source isocortex integration — `study`
    # (tasic2016 / tasic2018 / zeisel2018) captures both protocol AND
    # cohort-of-origin batch effects (Tasic 2016 and 2018 share lab but
    # differ in SMART-Seq / SMART-Seq v4 protocol versions). Matches the
    # HVG batch_key used at preprocess time and is the natural eval batch
    # for the cross-protocol integration claim.
    "mouse_cortex":        {"batch_key": "study",     "n_latent": N_LATENT},
    # human_immune is the scIB Immune (human) multi-study integration benchmark —
    # `batch` (10 studies: 10X / Freytag / Oetjen_* / Sun_* / Villani) is the
    # cross-study batch the embedding should mix, matching the eval BATCH_KEY in
    # embedding_metrics/metrics.py and the HVG batch_key used at preprocess time.
    # --batch-key none yields the `scvi_nobatch` cohort.
    "human_immune":        {"batch_key": "batch",     "n_latent": N_LATENT},
}


def fit(adata: ad.AnnData, params: dict, seed: int, ckpt_dir: Path) -> tuple[np.ndarray, dict]:
    import scvi  # heavy import; defer until needed
    scvi.settings.seed = seed

    if "highly_variable" not in adata.var:
        raise RuntimeError("var['highly_variable'] missing; rerun preprocessing.")
    if "counts" not in adata.layers:
        raise RuntimeError("layers['counts'] missing; rerun preprocessing.")

    sub = adata[:, adata.var["highly_variable"].to_numpy()].copy()
    sub.X = sub.layers["counts"].copy()  # scVI requires raw integer counts in X
    sub.X = _round_to_integer_counts(sub.X)

    _t_fit0 = time.perf_counter()
    scvi.model.SCVI.setup_anndata(sub, batch_key=params["batch_key"])
    model = scvi.model.SCVI(sub)
    accelerator = _resolve_accelerator(params.get("accelerator", "auto"))
    train_kw = {"accelerator": accelerator}
    if accelerator != "auto":
        train_kw["devices"] = 1
    model.train(**train_kw)  # heuristic max_epochs

    ckpt_dir.parent.mkdir(parents=True, exist_ok=True)
    model.save(str(ckpt_dir), overwrite=False)
    _t_fit1 = time.perf_counter()
    rss_after_obtain_mb = _rss_mb()

    emb = np.asarray(model.get_latent_representation(), dtype=np.float32)
    _t_tr1 = time.perf_counter()

    history_last: dict[str, float] = {}
    for k, v in (model.history or {}).items():
        if hasattr(v, "iloc") and len(v):
            history_last[k] = float(v.iloc[-1, 0])
    max_epochs = None
    if model.trainer is not None and hasattr(model.trainer, "max_epochs"):
        max_epochs = int(model.trainer.max_epochs)

    try:
        device = str(next(model.module.parameters()).device)
    except Exception:
        device = accelerator
    train_meta = {
        "max_epochs":   max_epochs,
        "history_last": history_last,
        "ckpt_dir":     str(ckpt_dir),
        "model_source":        "trained",
        "obtain_model_s":      float(_t_fit1 - _t_fit0),
        "build_embedding_s":   float(_t_tr1 - _t_fit1),
        "rss_after_obtain_mb": rss_after_obtain_mb,
        "accelerator_requested": params.get("accelerator", "auto"),
        "device":       device,
    }
    return emb, train_meta


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
    p, s, t = report["params"], report["stats"], report["train"]
    lines = [
        f"# scVI embedding report — {report['dataset']}",
        "",
        f"- Timestamp:  `{report['timestamp']}`",
        f"- Input:      `{report['input_path']}`",
        f"- Output:     `{report['output_path']}`",
        f"- Checkpoint: `{t['ckpt_dir']}`",
        "",
        "## Parameters",
        f"- `seed`: {report['seed']}",
        *(f"- `{k}`: {v}" for k, v in p.items()),
        "",
        "## Stats",
        f"- cells: {s['n_cells']}",
        f"- input genes (HVGs): {s['n_input_genes']}",
        f"- embedding dim: {s['embedding_dim']}",
        *_resource_lines(s, storage_what="scVI checkpoint dir",
                         obtain_what="setup_anndata + construct + train + save",
                         build_what="encoder forward pass"),
        "",
        "## Training",
        f"- max_epochs: {t['max_epochs']}",
        *(f"- `{k}`: {v:.4g}" for k, v in t["history_last"].items()),
    ]
    path.write_text("\n".join(lines) + "\n")


def resolve_params(args: argparse.Namespace) -> dict:
    p = DEFAULT_PARAMS[args.dataset].copy()
    if args.batch_key is not None:
        p["batch_key"] = None if args.batch_key.lower() == "none" else args.batch_key
    p["accelerator"] = args.accelerator or "auto"
    return p


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="scVI (de novo) embedding for preprocessed AnnData.")
    parser.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--batch-key", type=str,
                        help="Override the dataset default. Pass 'none' to disable.")
    parser.add_argument("--accelerator", choices=["auto", "cpu", "mps", "cuda", "gpu"],
                        default=None,
                        help="Lightning accelerator (default: auto = CUDA if present else CPU; "
                             "scvi-tools does not auto-pick Apple MPS — pass 'mps' to try it).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = resolve_params(args)
    seed_tag = f"seed{args.seed}"

    is_nobatch_override = (
        params["batch_key"] is None
        and DEFAULT_PARAMS[args.dataset]["batch_key"] is not None
    )
    model_name = "scvi_nobatch" if is_nobatch_override else MODEL_NAME

    input_path  = REPO_ROOT / "data" / "preprocessed" / f"{args.dataset}_pp.h5ad"
    out_dir     = REPO_ROOT / "outs" / "embeddings" / args.dataset
    stem        = f"{model_name}_{args.dataset}_{seed_tag}"
    output_path = out_dir / f"{stem}.parquet"
    report_path = out_dir / f"{stem}_report.md"
    ckpt_dir    = REPO_ROOT / "models" / args.dataset / model_name / seed_tag

    if output_path.exists():
        raise FileExistsError(
            f"{output_path} exists. Delete or rename it to overwrite, "
            "or pass --seed N to write a new file."
        )

    # Start the clock before the disk read; walltime_s = load_data + obtain_model
    # (= train) + build_embedding (= encoder forward) + other (the four sum to it).
    t0 = time.perf_counter()
    adata = ad.read_h5ad(input_path)
    cell_ids = adata.obs_names.astype(str).tolist()
    t_after_load = time.perf_counter()
    rss_after_load_mb = _rss_mb()

    emb, train_meta = fit(adata, params, seed=args.seed, ckpt_dir=ckpt_dir)
    t_after_build = time.perf_counter()
    rss_after_build_mb = _rss_mb()

    device = train_meta.get("device", "cpu")
    host = platform.node()
    model_storage_bytes = _path_bytes(ckpt_dir)
    load_data_s = t_after_load - t0
    obtain_model_s = train_meta.get("obtain_model_s")
    build_embedding_s = train_meta.get("build_embedding_s")
    stats = {
        "n_cells":             int(adata.n_obs),
        "n_input_genes":       int(adata.var["highly_variable"].sum()),
        "embedding_dim":       int(emb.shape[1]),
        "model_source":        train_meta.get("model_source", "trained"),
        "load_data_s":         float(load_data_s),
        "obtain_model_s":      obtain_model_s,
        "build_embedding_s":   build_embedding_s,
        "model_storage_bytes": model_storage_bytes,
        "model_download_s":    0.0,    # trained de novo; nothing downloaded
        "rss_after_load_mb":   rss_after_load_mb,
        "rss_after_obtain_mb": train_meta.get("rss_after_obtain_mb", rss_after_load_mb),
        "rss_after_build_mb":  rss_after_build_mb,
        "device":              device,
        "host":                host,
    }
    kv = {
        "model": model_name, "dataset": args.dataset, "seed": args.seed,
        "embedding_dim": stats["embedding_dim"], "n_latent": params["n_latent"],
        "batch_key": params["batch_key"] if params["batch_key"] is not None else "none",
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
        "seed": args.seed,
        "input_path": str(input_path),
        "output_path": str(output_path),
        "params": params,
        "stats": stats,
        "train": train_meta,
    }
    write_report(report, report_path)

    print(f"[ok] {input_path.name}: ({stats['n_cells']}, {stats['embedding_dim']}) "
          f"scVI embedding in {walltime_s:.1f}s "
          f"(load {load_data_s:.1f}s + train {obtain_model_s:.1f}s + forward {build_embedding_s:.2f}s)")
    print(f"[ok] wrote {output_path}")
    print(f"[ok] wrote {report_path}")
    print(f"[ok] saved model in {ckpt_dir}")


if __name__ == "__main__":
    main()
