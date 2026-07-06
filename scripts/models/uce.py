"""
UCE (Universal Cell Embedding) inference for preprocessed single-cell data.

Subprocess wrapper around snap-stanford/UCE's `eval_single_anndata.py` using
the 33-layer checkpoint (the largest released). The same checkpoint embeds both species
— only `--species` differs. Output dim = 1280. The encoder is in eval mode
(deterministic forward pass), but UCE's input pipeline draws a random gene
subset per cell (`np.random.choice` weighted by expression) and shuffles the
chromosome blocks (`np.random.shuffle`), and the upstream eval script never
seeds the RNGs — so by default UCE is non-deterministic across runs. We make
it reproducible by launching the eval script through a tiny seed-bootstrap that
sets `random` / `numpy` / `torch` seeds in the child process first; the seed is
exposed via `--seed` (default 0) and recorded in the report/metadata. (Output
stem stays `uce_<dataset>` — a 30-seed UCE cohort is impractical on CPU; if a
multi-seed cohort is ever wanted on GPU, switch to a `_seed{N}` stem then.)

Inputs:
    data/preprocessed/<dataset>_pp.h5ad

Outputs:
    outs/embeddings/<dataset>/uce_<dataset>.parquet
    outs/embeddings/<dataset>/uce_<dataset>_report.md
(Re-running just hits FileExistsError; delete the existing file to overwrite,
e.g. after a dependency bump or a checkpoint swap.)

Caches (shared across datasets/runs):
    models/uce/repo/                          (cloned from snap-stanford/UCE)
    models/uce/model_files/                   (33-layer checkpoint + support files)

Required one-time setup:
    1. The repo is cloned automatically on first run (or clone manually if offline).
    2. Download the 33-layer checkpoint from FigShare article 24320806
       ("Universal Cell Embedding Model Files") and place it at
       models/uce/model_files/33l_8ep_1024t_1280.torch

Usage:
    python scripts/models/uce.py --dataset mouse_brain
    python scripts/models/uce.py --dataset human_pancreas
"""
from __future__ import annotations

import argparse
import platform
import subprocess
import sys
import tempfile
import time
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from _resource_meta import path_bytes as _path_bytes, resource_lines as _resource_lines
from _resource_meta import rss_mb as _rss_mb, time_split as _time_split


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")
MODEL_NAME = "uce"


UCE_REPO_URL = "https://github.com/snap-stanford/UCE.git"
UCE_DIR     = REPO_ROOT / "models" / "uce" / "repo"
MODEL_FILES = REPO_ROOT / "models" / "uce" / "model_files"
CKPT_NAME   = "33l_8ep_1024t_1280.torch"
UCE_NLAYERS = 33
OUTPUT_DIM  = 1280

DEFAULT_PARAMS = {
    "mouse_brain":         {"species": "mouse"},
    "human_pancreas":      {"species": "human"},
    "mouse_gastrulation":  {"species": "mouse"},
    "mouse_cortex":        {"species": "mouse"},
    "human_immune":        {"species": "human"},
}


def _ensure_repo() -> None:
    if UCE_DIR.exists():
        return
    UCE_DIR.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(
        ["git", "clone", "--depth", "1", UCE_REPO_URL, str(UCE_DIR)],
        check=True,
    )


def _ensure_checkpoint() -> Path:
    ckpt = MODEL_FILES / CKPT_NAME
    if not ckpt.exists():
        raise FileNotFoundError(
            f"33-layer UCE checkpoint missing at {ckpt}. "
            "Download it manually from FigShare article 24320806 "
            "('Universal Cell Embedding Model Files') and place it there. "
            "We do not silently fall back to the 4-layer model — its embedding "
            "space is not comparable to the 33-layer one."
        )
    return ckpt


def _detect_device() -> str:
    try:
        import torch
        return "cuda" if torch.cuda.is_available() else "cpu"
    except Exception:
        return "cpu"


def _seeded_eval_argv(in_h5ad: Path, tmpdir: Path, ckpt: Path, species: str,
                      device: str, seed: int) -> list[str]:
    """Build the `python -c <bootstrap>` invocation that runs UCE's
    `eval_single_anndata.py` with `random`/`numpy`/`torch` seeded first.

    `runpy.run_path(..., run_name='__main__')` triggers the script's
    `if __name__ == '__main__'` block; we prepend the cwd (= UCE_DIR) to
    sys.path so its sibling-module imports (`utils`, `evaluate`, …) resolve,
    exactly as `python eval_single_anndata.py` would."""
    argv = [
        "eval_single_anndata.py",
        "--adata_path", str(in_h5ad),
        "--dir", str(tmpdir) + "/",
        "--species", species,
        "--model_loc", str(ckpt),
        "--nlayers", str(UCE_NLAYERS),
        "--batch_size", "25" if device == "cuda" else "4",
    ]
    bootstrap = (
        "import sys, random, numpy, runpy\n"
        "sys.path.insert(0, '.')\n"
        f"random.seed({seed}); numpy.random.seed({seed})\n"
        "try:\n"
        f"    import torch; torch.manual_seed({seed})\n"
        "except Exception:\n    pass\n"
        f"sys.argv = {argv!r}\n"
        "runpy.run_path('eval_single_anndata.py', run_name='__main__')\n"
    )
    return [sys.executable, "-c", bootstrap]


def fit(adata: ad.AnnData, params: dict, dataset: str) -> tuple[np.ndarray, dict]:
    if "counts" not in adata.layers:
        raise RuntimeError("layers['counts'] missing; rerun preprocessing.")

    _ensure_repo()
    ckpt = _ensure_checkpoint()
    device = _detect_device()
    seed = int(params.get("seed", 0))

    sub = adata.copy()
    sub.X = sub.layers["counts"].copy()  # raw counts in X; symbols already in var_names

    with tempfile.TemporaryDirectory() as tmp:
        tmpdir = Path(tmp)
        in_h5ad = tmpdir / f"{dataset}.h5ad"
        sub.write_h5ad(in_h5ad)

        subprocess.run(
            _seeded_eval_argv(in_h5ad, tmpdir, ckpt, params["species"], device, seed),
            cwd=UCE_DIR, check=True,
        )

        # eval_single_anndata writes "<input_stem>_uce_adata.h5ad" inside --dir.
        out_h5ad = tmpdir / f"{dataset}_uce_adata.h5ad"
        if not out_h5ad.exists():
            candidates = list(tmpdir.glob("*_uce_adata.h5ad"))
            if not candidates:
                raise RuntimeError(f"UCE produced no '*_uce_adata.h5ad' in {tmpdir}.")
            out_h5ad = candidates[0]
        out = ad.read_h5ad(out_h5ad)
        if "X_uce" not in out.obsm:
            raise RuntimeError(f"UCE output missing obsm['X_uce'] in {out_h5ad}.")
        emb = np.asarray(out.obsm["X_uce"], dtype=np.float32)

    meta = {
        "nlayers":       UCE_NLAYERS,
        "checkpoint":    CKPT_NAME,
        "species":       params["species"],
        "device":        device,
        "seed":          seed,
        "embedding_dim": int(emb.shape[1]),
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
    slow = "" if m["device"] == "cuda" else "  (CPU fallback — 33-layer UCE on CPU is slow)"
    lines = [
        f"# UCE embedding report — {report['dataset']}",
        "",
        f"- Timestamp: `{report['timestamp']}`",
        f"- Input:     `{report['input_path']}`",
        f"- Output:    `{report['output_path']}`",
        "",
        "## Parameters",
        *(f"- `{k}`: {v}" for k, v in p.items()),
        "",
        "## Model",
        f"- nlayers:    `{m['nlayers']}`",
        f"- checkpoint: `{m['checkpoint']}`",
        f"- species:    `{m['species']}`",
        f"- device:     `{m['device']}`{slow}",
        f"- seed:       `{m.get('seed', 0)}` (passed to random/numpy/torch in the eval subprocess)",
        "",
        "## Stats",
        f"- cells: {s['n_cells']}",
        f"- input genes (preprocessed): {s['n_input_genes']}",
        f"- embedding dim: {s['embedding_dim']}",
        *_resource_lines(s, storage_what="33-layer UCE checkpoint (.torch)",
                         obtain_what="loaded & run atomically inside UCE's eval subprocess",
                         build_what="write temp h5ad + UCE eval subprocess (checkpoint load + forward) + read output"),
    ]
    path.write_text("\n".join(lines) + "\n")


def resolve_params(args: argparse.Namespace) -> dict:
    p = DEFAULT_PARAMS[args.dataset].copy()
    p["seed"] = args.seed
    return p


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="UCE inference for preprocessed AnnData.")
    parser.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    parser.add_argument("--seed", type=int, default=0,
                        help="Seed for random/numpy/torch in UCE's eval subprocess "
                             "(makes the per-cell gene subsampling reproducible; default 0).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = resolve_params(args)

    input_path  = REPO_ROOT / "data" / "preprocessed" / f"{args.dataset}_pp.h5ad"
    out_dir     = REPO_ROOT / "outs" / "embeddings" / args.dataset
    stem        = f"{MODEL_NAME}_{args.dataset}"
    output_path = out_dir / f"{stem}.parquet"
    report_path = out_dir / f"{stem}_report.md"

    if output_path.exists():
        raise FileExistsError(
            f"{output_path} exists. Delete or rename it to overwrite."
        )

    # Start the clock before the disk read; walltime_s = load_data + obtain_model
    # + build_embedding + other. UCE's eval runs as one subprocess that loads the
    # pretrained checkpoint and forwards atomically, so obtain_model_s is n/a and
    # the whole subprocess (+ temp-h5ad I/O) counts as build_embedding. peak RSS
    # now includes the child process (self+children once it has been waited on).
    t0 = time.perf_counter()
    adata = ad.read_h5ad(input_path)
    cell_ids = adata.obs_names.astype(str).tolist()
    t_after_load = time.perf_counter()
    rss_after_load_mb = _rss_mb()

    emb, model_meta = fit(adata, params, dataset=args.dataset)
    t_after_build = time.perf_counter()
    rss_after_build_mb = _rss_mb()

    model_storage_bytes = _path_bytes(MODEL_FILES / CKPT_NAME)
    host = platform.node()
    load_data_s = t_after_load - t0
    build_embedding_s = t_after_build - t_after_load
    stats = {
        "n_cells":             int(adata.n_obs),
        "n_input_genes":       int(adata.n_vars),
        "embedding_dim":       int(emb.shape[1]),
        "model_source":        "loaded",   # pretrained checkpoint
        "load_data_s":         float(load_data_s),
        "obtain_model_s":      None,        # n/a — subprocess loads + forwards atomically
        "build_embedding_s":   float(build_embedding_s),
        "model_storage_bytes": model_storage_bytes,
        "model_download_s":    0.0,    # checkpoint is a one-time manual FigShare download; cached
        "rss_after_load_mb":   rss_after_load_mb,
        "rss_after_obtain_mb": rss_after_load_mb,   # no separable obtain step
        "rss_after_build_mb":  rss_after_build_mb,
        "device":              model_meta["device"],
        "host":                host,
    }
    kv = {
        "model": MODEL_NAME, "dataset": args.dataset, "embedding_dim": stats["embedding_dim"],
        "uce_nlayers": UCE_NLAYERS, "uce_checkpoint": CKPT_NAME, "species": params["species"],
        "seed": model_meta.get("seed", 0), "device": model_meta["device"],
        "model_source": "loaded", "model_storage_bytes": model_storage_bytes, "host": host,
    }
    write_parquet(emb, cell_ids, output_path, kv)
    walltime_s = time.perf_counter() - t0
    stats["walltime_s"] = float(walltime_s)
    stats["other_s"] = _time_split(walltime_s, load_data_s, 0.0, build_embedding_s)
    stats["peak_rss_mb"] = _rss_mb()
    stats["embedding_bytes"] = _path_bytes(output_path)
    kv["walltime_s"] = round(walltime_s, 3)
    kv["peak_rss_mb"] = round(stats["peak_rss_mb"], 1)

    report = {
        "dataset": args.dataset,
        "timestamp": datetime.now().isoformat(timespec="seconds"),
        "input_path": str(input_path),
        "output_path": str(output_path),
        "params": params,
        "stats": stats,
        "model": model_meta,
    }
    write_report(report, report_path)

    print(f"[ok] {input_path.name}: ({stats['n_cells']}, {stats['embedding_dim']}) "
          f"UCE-{UCE_NLAYERS} ({params['species']}) on {model_meta['device']} in {walltime_s:.1f}s")
    print(f"[ok] wrote {output_path}")
    print(f"[ok] wrote {report_path}")


if __name__ == "__main__":
    main()
