"""
expiMap embedding (de novo) for preprocessed single-cell data.

Conditional VAE from scArches whose latent dimensions are anchored to a
user-supplied collection of gene programs. Trained from scratch on raw counts
(layers["counts"]) restricted to the union of HVGs and program-member genes
(so program members are never dropped during HVG selection). Latent dim equals
the number of active programs (those not zeroed by expiMap's group-lasso
pruning), plus any extra de-novo dims when --n-ext > 0. The seed controls weight
init + minibatch order.

Default gene programs are MSigDB Reactome (organism-native: `c2.cp.reactome`
for human, `m2.cp.reactome` for mouse). The .gmt is downloaded directly from
the Broad MSigDB server and cached under models/genesets/<organism>/. Override
with --gene-sets <path-to-gmt>.

Inputs:
    data/preprocessed/<dataset>_pp.h5ad

Outputs (per run):
    outs/embeddings/<dataset>/expimap_<dataset>_<run_tag>.parquet
    outs/embeddings/<dataset>/expimap_<dataset>_<run_tag>_report.md
    models/<dataset>/expimap/<run_tag>/                  (trained model checkpoint)

Caches (shared across datasets/runs):
    models/genesets/<organism>/<category>.v<dbver>.symbols.gmt

Different seeds produce distinct files (...seed{N}.parquet). Same-seed re-runs
will fail with FileExistsError so historical results aren't silently clobbered
— delete or rename the existing file if you really mean to overwrite.

Usage:
    python scripts/models/expimap.py --dataset mouse_brain
    python scripts/models/expimap.py --dataset human_pancreas --seed 1
    python scripts/models/expimap.py --dataset human_pancreas --gene-sets path/to.gmt
"""
from __future__ import annotations

# Import the shared resource-cost helper *before* the sys.path surgery below.
from _resource_meta import path_bytes as _path_bytes, resource_lines as _resource_lines
from _resource_meta import rss_mb as _rss_mb, time_split as _time_split

# IMPORTANT: scripts/models/scvi.py shadows the installed `scvi` package because
# Python prepends the script's directory to sys.path. scarches imports scvi
# transitively, so strip that entry before any third-party imports.
import os as _os
import sys as _sys
_HERE = _os.path.dirname(_os.path.abspath(__file__))
_sys.path[:] = [p for p in _sys.path if _os.path.abspath(p) != _HERE]
del _os, _sys, _HERE

import argparse
import json
import platform
import time
import urllib.request
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp

# scArches 0.6.1 imports `from anndata import read`, which was removed in
# anndata 0.10. Restore the alias so `import scarches` succeeds.
if not hasattr(ad, "read"):
    ad.read = ad.read_h5ad


def _round_to_integer_counts(X):
    """Round nonzero entries to nearest integer in place.

    No-op for true UMI data; only matters for SMART-Seq RSEM expected counts
    and similar fractional quantifications. Aligns the input with expiMap's
    Negative Binomial likelihood, which assumes count data.
    """
    if sp.issparse(X):
        X.data = np.rint(X.data).astype(X.data.dtype)
        return X
    return np.rint(X).astype(X.dtype)


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation",
                   "mouse_cortex", "human_immune")
MODEL_NAME = "expimap"

MSIGDB_URL = ("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
              "{dbver}/{category}.v{dbver}.symbols.gmt")

# Standard scArches expiMap recipe (matches the public tutorial).
EARLY_STOPPING_KWARGS = {
    "early_stopping_metric": "val_unweighted_loss",
    "threshold":             0,
    "patience":              50,
    "reduce_lr":             True,
    "lr_patience":           13,
    "lr_factor":             0.1,
}

DEFAULT_PARAMS = {
    "mouse_brain": {
        "organism":     "mouse",
        "msigdb_dbver": "2024.1.Mm",
        "category":     "m2.cp.reactome",
        "batch_key":    None,
        "n_ext":        0,
        "min_genes":    12,
        "max_genes":    None,
        "n_epochs":     400,
        "alpha":        0.7,
        "alpha_kl":     0.5,
        "hidden_layer_sizes": [256, 256, 256],
    },
    # Benchmark datasets (mouse_cortex, human_pancreas, human_immune) are run
    # author-faithfully against Reactome ∪ PanglaoDB via --gene-sets (see
    # scripts/IO/build_panglao_gmt.py), using the scArches expiMap recipe:
    # restrict to ~2000 HVGs (in fit()) and keep programs with >= min_genes=12
    # members among them, which lands the latent at the published ~hundreds-of-
    # programs scale rather than >1000.
    "human_pancreas": {
        "organism":     "human",
        "msigdb_dbver": "2024.1.Hs",
        "category":     "c2.cp.reactome",
        "batch_key":    "tech",
        "n_ext":        0,
        "min_genes":    12,
        "max_genes":    None,
        "n_epochs":     400,
        "alpha":        0.7,
        "alpha_kl":     0.5,
        "hidden_layer_sizes": [256, 256, 256],
    },
    "mouse_cortex": {
        "organism":     "mouse",
        "msigdb_dbver": "2024.1.Mm",
        "category":     "m2.cp.reactome",
        "batch_key":    "study",
        "n_ext":        0,
        "min_genes":    12,
        "max_genes":    None,
        "n_epochs":     400,
        "alpha":        0.7,
        "alpha_kl":     0.5,
        "hidden_layer_sizes": [256, 256, 256],
    },
    "human_immune": {
        "organism":     "human",
        "msigdb_dbver": "2024.1.Hs",
        "category":     "c2.cp.reactome",
        "batch_key":    "batch",
        "n_ext":        0,
        "min_genes":    12,
        "max_genes":    None,
        "n_epochs":     400,
        "alpha":        0.7,
        "alpha_kl":     0.5,
        "hidden_layer_sizes": [256, 256, 256],
    },
    "mouse_gastrulation": {
        "organism":     "mouse",
        "msigdb_dbver": "2024.1.Mm",
        "category":     "m2.cp.reactome",
        "batch_key":    "sample_id",
        "n_ext":        0,
        "min_genes":    12,
        "max_genes":    None,
        "n_epochs":     400,
        "alpha":        0.7,
        "alpha_kl":     0.5,
        "hidden_layer_sizes": [256, 256, 256],
    },
}


def _download_msigdb_gmt(category: str, dbver: str, cache_dir: Path) -> Path:
    """Download <category>.v<dbver>.symbols.gmt from the Broad MSigDB server."""
    gmt_path = cache_dir / f"{category}.v{dbver}.symbols.gmt"
    if gmt_path.exists():
        return gmt_path
    cache_dir.mkdir(parents=True, exist_ok=True)
    url = MSIGDB_URL.format(dbver=dbver, category=category)
    try:
        urllib.request.urlretrieve(url, gmt_path)
    except Exception as e:
        if gmt_path.exists():
            gmt_path.unlink()
        raise RuntimeError(
            f"Failed to download MSigDB gmt from {url}: {e}. "
            "Check --category/--msigdb-dbver, or supply --gene-sets <local.gmt>."
        ) from e
    return gmt_path


def _load_gmt(path: Path) -> dict[str, list[str]]:
    """Read a .gmt; return {program_name: [genes]}. Skips the URL/description column."""
    out: dict[str, list[str]] = {}
    for line in path.read_text().splitlines():
        if not line.strip():
            continue
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        out[parts[0]] = parts[2:]
    return out


def fit(adata: ad.AnnData, params: dict, gmt_path: Path, seed: int,
        ckpt_dir: Path) -> tuple[np.ndarray, dict]:
    import scarches as sca  # heavy import; defer until needed

    if "highly_variable" not in adata.var:
        raise RuntimeError("var['highly_variable'] missing; rerun preprocessing.")
    if "counts" not in adata.layers:
        raise RuntimeError("layers['counts'] missing; rerun preprocessing.")

    # Author-faithful gene universe: restrict to the highly variable genes (the
    # scArches expiMap tutorial selects ~2000 HVGs), then keep programs with
    # >= min_genes members among them. Taking the union with all program-member
    # genes instead lets almost every Reactome pathway survive and inflates the
    # latent to >1000 dims, unlike the published expiMap (~hundreds of programs).
    sub = adata[:, adata.var["highly_variable"].to_numpy()].copy()
    sub.X = sub.layers["counts"].copy()  # expiMap requires raw integer counts in X
    sub.X = _round_to_integer_counts(sub.X)

    # expiMap's forward pass conditions theta on a one-hot batch vector and
    # cannot handle condition_key=None. When the dataset has no batch column,
    # inject a synthetic single-level "_expimap_batch" so the conditional
    # plumbing has something to one-hot. With one level, the conditioning
    # collapses to a constant — i.e. effectively no batch correction.
    if params["batch_key"] is None:
        sub.obs["_expimap_batch"] = "all"
        condition_key = "_expimap_batch"
    else:
        condition_key = params["batch_key"]

    _t_fit0 = time.perf_counter()
    # clean=False keeps full, unique program names. With clean=True scArches
    # "removes the word before the first underscore and cuts to 30 chars" (to
    # strip REACTOME_ prefixes), which collapses multi-word PanglaoDB cell-type
    # names onto a shared suffix (e.g. Acinar_cells / Alpha_cells -> CELLS),
    # colliding and silently dropping the cell-type programs.
    sca.utils.add_annotations(
        sub, str(gmt_path),
        min_genes=params["min_genes"],
        max_genes=params["max_genes"],
        clean=False,
    )
    n_programs = int(sub.varm["I"].shape[1]) if "I" in sub.varm else 0
    program_names_full = list(sub.uns.get("terms", []))
    if n_programs == 0:
        raise RuntimeError(
            f"No gene programs survived min_genes={params['min_genes']} on "
            f"{int(sub.n_vars)} retained genes. Lower --min-genes or supply a "
            "different --gene-sets file."
        )

    model = sca.models.EXPIMAP(
        adata=sub,
        condition_key=condition_key,
        hidden_layer_sizes=params["hidden_layer_sizes"],
        recon_loss="nb",
        n_ext=params["n_ext"],
    )
    model.train(
        n_epochs=params["n_epochs"],
        alpha=params["alpha"],
        alpha_kl=params["alpha_kl"],
        alpha_epoch_anneal=100,
        weight_decay=0.0,
        early_stopping_kwargs=EARLY_STOPPING_KWARGS,
        use_early_stopping=True,
        monitor_only_val=False,
        seed=seed,
    )

    ckpt_dir.parent.mkdir(parents=True, exist_ok=True)
    model.save(str(ckpt_dir), overwrite=False, save_anndata=False)
    _t_fit1 = time.perf_counter()
    rss_after_obtain_mb = _rss_mb()

    # scarches.BaseMixin.get_latent does sparse indexing with a torch tensor,
    # which scipy 1.12+ rejects. Pass x densely (and c explicitly) to bypass
    # the broken code path.
    x_dense = sub.X.toarray() if sp.issparse(sub.X) else np.asarray(sub.X)
    c_arr   = sub.obs[condition_key].to_numpy()
    emb = np.asarray(
        model.get_latent(x=x_dense, c=c_arr, mean=True, only_active=True),
        dtype=np.float32,
    )
    _t_tr1 = time.perf_counter()

    # Build program-name labels for each surviving latent dim. nonzero_terms()
    # is a bool array of length n_programs (+ n_ext_m + n_ext when used). The
    # first n_programs slots correspond to named Reactome programs; trailing
    # ones are unnamed de-novo dims (only present when n_ext / n_ext_m > 0).
    nz = np.asarray(model.nonzero_terms())
    n_extra = len(nz) - n_programs
    all_labels = list(program_names_full) + [f"ext_{i}" for i in range(n_extra)]
    active_indices = np.where(nz)[0].tolist()
    program_names_active = [all_labels[i] for i in active_indices]
    if len(program_names_active) != emb.shape[1]:
        raise RuntimeError(
            f"Program-name count {len(program_names_active)} does not match "
            f"latent dim count {emb.shape[1]}; scarches API drift?"
        )

    # scArches uses CUDA when available, else CPU (no MPS path); read the
    # actual device off the trained module, falling back to that heuristic.
    try:
        device = str(next(model.model.parameters()).device)
    except Exception:
        try:
            import torch
            device = "cuda" if torch.cuda.is_available() else "cpu"
        except Exception:
            device = "cpu"
    train_meta = {
        "n_programs":         n_programs,
        "n_active":           int(emb.shape[1]),
        "n_inactive":         int(n_programs + n_extra - emb.shape[1]),
        "program_names":      program_names_active,
        "n_ext":              params["n_ext"],
        "n_epochs":           params["n_epochs"],
        "ckpt_dir":           str(ckpt_dir),
        "n_input_genes":      int(sub.n_vars),
        "n_hvgs":             int(adata.var["highly_variable"].sum()),
        "model_source":        "trained",
        "obtain_model_s":      float(_t_fit1 - _t_fit0),
        "build_embedding_s":   float(_t_tr1 - _t_fit1),
        "rss_after_obtain_mb": rss_after_obtain_mb,
        "device":             device,
    }
    return emb, train_meta


def write_parquet(emb: np.ndarray, cell_ids: list[str], program_names: list[str],
                  path: Path, kv: dict) -> None:
    cols = {f"dim_{i}": emb[:, i].astype(np.float32) for i in range(emb.shape[1])}
    table = pa.Table.from_pydict({"cell_id": cell_ids, **cols})
    md = {**(table.schema.metadata or {}),
          **{k.encode(): str(v).encode() for k, v in kv.items()},
          b"programs_json": json.dumps(program_names).encode()}
    table = table.replace_schema_metadata(md)
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, path, compression="zstd")


def write_report(report: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    p, s, t = report["params"], report["stats"], report["train"]
    lines = [
        f"# expiMap embedding report — {report['dataset']}",
        "",
        f"- Timestamp:  `{report['timestamp']}`",
        f"- Input:      `{report['input_path']}`",
        f"- Output:     `{report['output_path']}`",
        f"- Gene sets:  `{report['gmt_path']}`",
        f"- Checkpoint: `{t['ckpt_dir']}`",
        "",
        "## Parameters",
        f"- `seed`: {report['seed']}",
        *(f"- `{k}`: {v}" for k, v in p.items()),
        "",
        "## Stats",
        f"- cells: {s['n_cells']}",
        f"- HVGs (preprocessed):                  {t['n_hvgs']}",
        f"- input genes (HVGs ∪ program members): {t['n_input_genes']}",
        f"- gene programs (after min_genes filter): {t['n_programs']}",
        f"- active programs (after group lasso):   {t['n_active']}",
        f"- inactive programs (zeroed by lasso):   {t['n_inactive']}",
        f"- embedding dim (= active programs):     {s['embedding_dim']}",
        *_resource_lines(s, storage_what="expiMap checkpoint dir",
                         obtain_what="annotate programs + construct + train + save",
                         build_what="get_latent forward pass"),
        "",
        "## Notes",
        "- Each `dim_i` corresponds to one Reactome program; the mapping is",
        "  stored as `programs_json` in the parquet schema metadata.",
        "- Inactive programs (zeroed out by group-lasso regularization) are",
        "  dropped from the output. To recover the full layout, retrain with",
        "  `only_active=False` in `model.get_latent`.",
    ]
    path.write_text("\n".join(lines) + "\n")


def resolve_params(args: argparse.Namespace) -> dict:
    p = DEFAULT_PARAMS[args.dataset].copy()
    if args.batch_key is not None:
        p["batch_key"] = None if args.batch_key.lower() == "none" else args.batch_key
    if args.category is not None:
        p["category"] = args.category
    if args.msigdb_dbver is not None:
        p["msigdb_dbver"] = args.msigdb_dbver
    if args.n_ext is not None:
        p["n_ext"] = args.n_ext
    if args.min_genes is not None:
        p["min_genes"] = args.min_genes
    if args.n_epochs is not None:
        p["n_epochs"] = args.n_epochs
    return p


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="expiMap (scArches) embedding for preprocessed AnnData."
    )
    parser.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--batch-key", type=str,
                        help="Override the dataset default. Pass 'none' to disable.")
    parser.add_argument("--category", type=str,
                        help="MSigDB category code (e.g. c2.cp.reactome / m2.cp.reactome / "
                             "h.all / mh.all). Ignored if --gene-sets is set.")
    parser.add_argument("--msigdb-dbver", type=str,
                        help="MSigDB release (e.g. 2024.1.Hs / 2024.1.Mm). "
                             "Ignored if --gene-sets is set.")
    parser.add_argument("--gene-sets", type=str,
                        help="Path to a .gmt file. Bypasses the MSigDB download.")
    parser.add_argument("--n-ext", type=int,
                        help="Extra de-novo latent dims (default 0).")
    parser.add_argument("--min-genes", type=int,
                        help="Drop programs with fewer than this many member genes "
                             "present in the input (default 12).")
    parser.add_argument("--n-epochs", type=int,
                        help="Max training epochs; early stopping may halt sooner "
                             "(default 400).")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    params = resolve_params(args)
    seed_tag = f"seed{args.seed}"

    input_path  = REPO_ROOT / "data" / "preprocessed" / f"{args.dataset}_pp.h5ad"
    out_dir     = REPO_ROOT / "outs" / "embeddings" / args.dataset
    stem        = f"{MODEL_NAME}_{args.dataset}_{seed_tag}"
    output_path = out_dir / f"{stem}.parquet"
    report_path = out_dir / f"{stem}_report.md"
    ckpt_dir    = REPO_ROOT / "models" / args.dataset / "expimap" / seed_tag

    if output_path.exists():
        raise FileExistsError(
            f"{output_path} exists. Delete or rename it to overwrite, "
            "or pass --seed N to write a new file."
        )

    if args.gene_sets is not None:
        gmt_path = Path(args.gene_sets).resolve()
        if not gmt_path.exists():
            raise FileNotFoundError(f"--gene-sets file not found: {gmt_path}")
    else:
        cache = REPO_ROOT / "models" / "genesets" / params["organism"]
        gmt_path = _download_msigdb_gmt(params["category"], params["msigdb_dbver"], cache)

    # Start the clock before the disk read; walltime_s = load_data + obtain_model
    # (= annotate + construct + train + save) + build_embedding (= get_latent) + other.
    t0 = time.perf_counter()
    adata = ad.read_h5ad(input_path)
    cell_ids = adata.obs_names.astype(str).tolist()
    t_after_load = time.perf_counter()
    rss_after_load_mb = _rss_mb()

    emb, train_meta = fit(adata, params, gmt_path=gmt_path, seed=args.seed, ckpt_dir=ckpt_dir)
    t_after_build = time.perf_counter()
    rss_after_build_mb = _rss_mb()

    model_storage_bytes = _path_bytes(ckpt_dir)
    device = train_meta.get("device", "cpu")
    host = platform.node()
    load_data_s = t_after_load - t0
    obtain_model_s = train_meta.get("obtain_model_s")
    build_embedding_s = train_meta.get("build_embedding_s")
    stats = {
        "n_cells":             int(adata.n_obs),
        "embedding_dim":       int(emb.shape[1]),
        "model_source":        train_meta.get("model_source", "trained"),
        "load_data_s":         float(load_data_s),
        "obtain_model_s":      obtain_model_s,
        "build_embedding_s":   build_embedding_s,
        "model_storage_bytes": model_storage_bytes,
        "model_download_s":    0.0,    # trained de novo; the .gmt gene set is tiny/cached
        "rss_after_load_mb":   rss_after_load_mb,
        "rss_after_obtain_mb": train_meta.get("rss_after_obtain_mb", rss_after_load_mb),
        "rss_after_build_mb":  rss_after_build_mb,
        "device":              device,
        "host":                host,
    }
    kv = {
        "model": MODEL_NAME, "dataset": args.dataset, "seed": args.seed,
        "embedding_dim": stats["embedding_dim"], "n_programs": train_meta["n_programs"],
        "n_active": train_meta["n_active"], "n_ext": params["n_ext"],
        "batch_key": params["batch_key"] if params["batch_key"] is not None else "none",
        "gene_sets": str(gmt_path), "msigdb_dbver": params["msigdb_dbver"], "category": params["category"],
        "model_source": stats["model_source"], "model_storage_bytes": model_storage_bytes,
        "device": device, "host": host,
    }
    write_parquet(emb, cell_ids, train_meta["program_names"], output_path, kv)
    walltime_s = time.perf_counter() - t0
    stats["walltime_s"] = float(walltime_s)
    stats["other_s"] = _time_split(walltime_s, load_data_s, obtain_model_s, build_embedding_s)
    stats["peak_rss_mb"] = _rss_mb()
    stats["embedding_bytes"] = _path_bytes(output_path)
    kv["walltime_s"] = round(walltime_s, 3)
    kv["peak_rss_mb"] = round(stats["peak_rss_mb"], 1)

    report = {
        "dataset":     args.dataset,
        "timestamp":   datetime.now().isoformat(timespec="seconds"),
        "seed":        args.seed,
        "input_path":  str(input_path),
        "output_path": str(output_path),
        "gmt_path":    str(gmt_path),
        "params":      params,
        "stats":       stats,
        "train":       train_meta,
    }
    write_report(report, report_path)

    print(f"[ok] {input_path.name}: ({stats['n_cells']}, {stats['embedding_dim']}) "
          f"expiMap embedding in {walltime_s:.1f}s "
          f"(load {load_data_s:.1f}s + train {obtain_model_s:.1f}s + forward {build_embedding_s:.2f}s)")
    print(f"[ok] wrote {output_path}")
    print(f"[ok] wrote {report_path}")
    print(f"[ok] saved model in {ckpt_dir}")


if __name__ == "__main__":
    main()
