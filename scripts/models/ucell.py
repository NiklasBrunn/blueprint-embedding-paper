"""
UCell embedding for preprocessed single-cell data.

Within-cell rank-based gene-set scoring (Mann-Whitney U statistic), via the
official Python port pyUCell (carmonalab/pyUCell). Non-learned and
deterministic — no --seed exposed because the output is fully determined by
(dataset, gene_set source). Embedding dim equals the number of supplied
signatures.

Default gene programs are MSigDB Reactome (organism-native: `c2.cp.reactome`
for human, `m2.cp.reactome` for mouse). The .gmt is downloaded directly from
the Broad MSigDB server and cached under models/genesets/<organism>/. Override
with --gene-sets <path-to-gmt>.

UCell's score is invariant to monotone transformations of expression (it ranks
genes within each cell), so log1p-normalized X gives the same result as raw
counts. We use adata.X (log1p) directly — no HVG restriction, since dropping
genes can break gene-set membership.

Inputs:
    data/preprocessed/<dataset>_pp.h5ad

Outputs:
    outs/embeddings/<dataset>/ucell_<dataset>.parquet
    outs/embeddings/<dataset>/ucell_<dataset>_report.md
(Re-running just hits FileExistsError; delete the existing file to overwrite,
e.g. after a gene-set source change.)

Caches (shared across datasets/runs):
    models/genesets/<organism>/<category>.v<dbver>.symbols.gmt

Usage:
    python scripts/models/ucell.py --dataset mouse_brain
    python scripts/models/ucell.py --dataset human_pancreas
    python scripts/models/ucell.py --dataset human_pancreas --gene-sets path/to.gmt
"""
from __future__ import annotations

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

from _resource_meta import path_bytes as _path_bytes, resource_lines as _resource_lines
from _resource_meta import rss_mb as _rss_mb, time_split as _time_split


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation",
                   "mouse_cortex", "human_immune")
MODEL_NAME = "ucell"

MSIGDB_URL = ("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
              "{dbver}/{category}.v{dbver}.symbols.gmt")

DEFAULT_PARAMS = {
    "mouse_brain":         {"organism": "mouse", "msigdb_dbver": "2024.1.Mm",
                            "category": "m2.cp.reactome"},
    "human_pancreas":      {"organism": "human", "msigdb_dbver": "2024.1.Hs",
                            "category": "c2.cp.reactome"},
    "mouse_gastrulation":  {"organism": "mouse", "msigdb_dbver": "2024.1.Mm",
                            "category": "m2.cp.reactome"},
    # Benchmark datasets are run author-faithfully against PanglaoDB cell-type
    # markers via --gene-sets (see scripts/IO/build_panglao_gmt.py); the MSigDB
    # fields below are the fallback only when --gene-sets is omitted.
    "mouse_cortex":        {"organism": "mouse", "msigdb_dbver": "2024.1.Mm",
                            "category": "m2.cp.reactome"},
    "human_immune":        {"organism": "human", "msigdb_dbver": "2024.1.Hs",
                            "category": "c2.cp.reactome"},
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


def fit(adata: ad.AnnData, gmt_path: Path) -> tuple[np.ndarray, list[str], dict]:
    import pyucell  # heavy import; defer until needed

    programs = _load_gmt(gmt_path)
    sig_names_full = list(programs.keys())

    program_genes = {g for genes in programs.values() for g in genes}
    n_input = int(adata.n_vars)
    n_mapped = int(adata.var_names.isin(program_genes).sum())
    if n_mapped == 0:
        raise RuntimeError(
            f"No gene from {gmt_path} maps to adata.var_names. Symbols/organism mismatch?"
        )

    work = adata.copy()  # pyucell appends columns to adata.obs in place
    pyucell.compute_ucell_scores(work, signatures=programs)

    score_cols = [f"{name}_UCell" for name in sig_names_full]
    missing = [c for c in score_cols if c not in work.obs.columns]
    if missing:
        raise RuntimeError(
            f"pyUCell did not write {len(missing)} expected score column(s); "
            f"first missing: {missing[0]!r}. Check pyucell version (>=0.6)."
        )
    emb_full = work.obs[score_cols].to_numpy(dtype=np.float32)

    # Drop signatures whose score is constant across cells (typically because
    # no member genes are present in the dataset, or all member genes are
    # imputed identically as max-rank). They carry no information and inflate
    # the embedding dim with dead axes.
    active_mask = emb_full.std(axis=0) > 0
    emb = emb_full[:, active_mask]
    sig_names = [n for n, keep in zip(sig_names_full, active_mask) if keep]

    meta = {
        "n_input_genes":     n_input,
        "n_mapped_genes":    n_mapped,
        "n_signatures_full": len(sig_names_full),
        "n_signatures":      int(emb.shape[1]),
        "n_dropped_zero":    int((~active_mask).sum()),
    }
    return emb, sig_names, meta


def write_parquet(emb: np.ndarray, cell_ids: list[str], sig_names: list[str],
                  path: Path, kv: dict) -> None:
    cols = {f"dim_{i}": emb[:, i].astype(np.float32) for i in range(emb.shape[1])}
    table = pa.Table.from_pydict({"cell_id": cell_ids, **cols})
    md = {**(table.schema.metadata or {}),
          **{k.encode(): str(v).encode() for k, v in kv.items()},
          b"programs_json": json.dumps(sig_names).encode()}
    table = table.replace_schema_metadata(md)
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, path, compression="zstd")


def write_report(report: dict, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    p, s, m = report["params"], report["stats"], report["model"]
    overlap_pct = 100.0 * m["n_mapped_genes"] / max(m["n_input_genes"], 1)
    lines = [
        f"# UCell embedding report — {report['dataset']}",
        "",
        f"- Timestamp: `{report['timestamp']}`",
        f"- Input:     `{report['input_path']}`",
        f"- Output:    `{report['output_path']}`",
        f"- Gene sets: `{report['gmt_path']}`",
        "",
        "## Parameters",
        *(f"- `{k}`: {v}" for k, v in p.items()),
        "",
        "## Stats",
        f"- cells: {s['n_cells']}",
        f"- input genes (preprocessed):              {m['n_input_genes']}",
        f"- mapped to any gene-set member symbol:    {m['n_mapped_genes']} ({overlap_pct:.1f}%)",
        f"- signatures (full GMT):                   {m['n_signatures_full']}",
        f"- signatures dropped (zero variance):      {m['n_dropped_zero']}",
        f"- signatures kept (= embedding dim):       {m['n_signatures']}",
        *_resource_lines(s, storage_what="gene-set .gmt — defines the axes UCell scores against",
                         obtain_what="none — UCell does no model fitting (within-cell rank scoring only)",
                         build_what="read .gmt + compute_ucell_scores (Mann-Whitney rank scoring) + drop dead axes"),
        "",
        "## Notes",
        "- Each `dim_i` corresponds to one gene program (signature) from the",
        "  supplied .gmt; the mapping is stored as `programs_json` in the parquet",
        "  schema metadata.",
        "- Coverage above is **reference → dataset**: what fraction of the",
        "  preprocessed gene set appears in any program. Not directly comparable",
        "  to scVI Census coverage, which is **dataset → model** (fraction of",
        "  preprocessed genes that map to the encoder's reference Ensembl space).",
    ]
    path.write_text("\n".join(lines) + "\n")


def resolve_params(args: argparse.Namespace) -> dict:
    p = DEFAULT_PARAMS[args.dataset].copy()
    if args.category is not None:
        p["category"] = args.category
    if args.msigdb_dbver is not None:
        p["msigdb_dbver"] = args.msigdb_dbver
    return p


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="UCell embedding for preprocessed AnnData.")
    parser.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    parser.add_argument("--category", type=str,
                        help="MSigDB category code (e.g. c2.cp.reactome / m2.cp.reactome / "
                             "h.all / mh.all). Ignored if --gene-sets is set.")
    parser.add_argument("--msigdb-dbver", type=str,
                        help="MSigDB release (e.g. 2024.1.Hs / 2024.1.Mm). "
                             "Ignored if --gene-sets is set.")
    parser.add_argument("--gene-sets", type=str,
                        help="Path to a .gmt file. Bypasses the MSigDB download.")
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

    if args.gene_sets is not None:
        gmt_path = Path(args.gene_sets).resolve()
        if not gmt_path.exists():
            raise FileNotFoundError(f"--gene-sets file not found: {gmt_path}")
    else:
        cache = REPO_ROOT / "models" / "genesets" / params["organism"]
        gmt_path = _download_msigdb_gmt(params["category"], params["msigdb_dbver"], cache)

    # Start the clock before the disk read; walltime_s = load_data + obtain_model
    # (= 0; UCell does no model fitting) + build_embedding (read .gmt + rank scoring)
    # + other. The .gmt that defines the axes is recorded as model_storage_bytes.
    t0 = time.perf_counter()
    adata = ad.read_h5ad(input_path)
    cell_ids = adata.obs_names.astype(str).tolist()
    t_after_load = time.perf_counter()
    rss_after_load_mb = _rss_mb()

    emb, sig_names, meta = fit(adata, gmt_path=gmt_path)
    t_after_build = time.perf_counter()
    rss_after_build_mb = _rss_mb()

    model_storage_bytes = _path_bytes(gmt_path)  # the gene-set .gmt defines the axes
    load_data_s = t_after_load - t0
    build_embedding_s = t_after_build - t_after_load
    stats = {
        "n_cells":             int(adata.n_obs),
        "embedding_dim":       int(emb.shape[1]),
        "model_source":        "none",
        "load_data_s":         float(load_data_s),
        "obtain_model_s":      0.0,
        "build_embedding_s":   float(build_embedding_s),
        "model_storage_bytes": model_storage_bytes,
        "model_download_s":    0.0,    # the .gmt is tiny and cached under models/genesets/
        "rss_after_load_mb":   rss_after_load_mb,
        "rss_after_obtain_mb": rss_after_load_mb,   # no obtain step
        "rss_after_build_mb":  rss_after_build_mb,
        "device":              "cpu",   # pyUCell rank scoring — no GPU path
        "host":                platform.node(),
    }
    kv = {
        "model": MODEL_NAME, "dataset": args.dataset, "embedding_dim": stats["embedding_dim"],
        "n_signatures": meta["n_signatures"], "n_signatures_full": meta["n_signatures_full"],
        "n_dropped_zero": meta["n_dropped_zero"], "gene_sets": str(gmt_path),
        "msigdb_dbver": params["msigdb_dbver"], "category": params["category"],
        "model_source": "none", "model_storage_bytes": model_storage_bytes,
        "device": "cpu", "host": stats["host"],
    }
    write_parquet(emb, cell_ids, sig_names, output_path, kv)
    walltime_s = time.perf_counter() - t0
    stats["walltime_s"] = float(walltime_s)
    stats["other_s"] = _time_split(walltime_s, load_data_s, 0.0, build_embedding_s)
    stats["peak_rss_mb"] = _rss_mb()
    stats["embedding_bytes"] = _path_bytes(output_path)
    kv["walltime_s"] = round(walltime_s, 3)
    kv["peak_rss_mb"] = round(stats["peak_rss_mb"], 1)

    report = {
        "dataset":     args.dataset,
        "timestamp":   datetime.now().isoformat(timespec="seconds"),
        "input_path":  str(input_path),
        "output_path": str(output_path),
        "gmt_path":    str(gmt_path),
        "params":      params,
        "stats":       stats,
        "model":       meta,
    }
    write_report(report, report_path)

    print(f"[ok] {input_path.name}: ({stats['n_cells']}, {stats['embedding_dim']}) "
          f"UCell embedding in {walltime_s:.1f}s")
    print(f"[ok] wrote {output_path}")
    print(f"[ok] wrote {report_path}")


if __name__ == "__main__":
    main()
