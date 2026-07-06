"""
Load the scIB human immune integration benchmark dataset.

Benchmark assembly:
    Luecken et al. (2022) "Benchmarking atlas-level data integration in
    single-cell genomics", Nature Methods, doi:10.1038/s41592-021-01336-8.

This is a composite of several published human immune-cell scRNA-seq studies
(peripheral blood and bone marrow), harmonized into a single AnnData. It is the
"Immune (human)" integration task from the scIB benchmark — a deliberately hard,
multi-study batch problem spanning 10x and Smart-seq2 chemistries. Batches live
in obs["batch"] (ten batches) and harmonized cell-type labels in
obs["final_annotation"] (T / NK / B / monocyte / DC / progenitor / erythroid /
megakaryocyte lineages, etc.).

Source: figshare ndownloader file 25717328
    - Immune_ALL_human.h5ad  (~2.06 GB)
from the same collection as the scIB pancreas file (24539828):
https://figshare.com/articles/dataset/12420968

The downloaded H5AD ships with normalized expression in adata.X and raw counts
in adata.layers["counts"]. We move counts into adata.X so the output file is a
clean raw-counts AnnData, with the original counts also retained in
layers["counts"]. Native obs columns (notably ``batch`` for batches and
``final_annotation`` for labels) are preserved unchanged.

Outputs:
    - data/raw/raw_human_immune.h5ad   (raw counts in adata.X and adata.layers["counts"])
    - data/raw/raw_human_immune.json   (dataset metadata and provenance, sibling of the h5ad)

Usage:
    python scripts/IO/load_human_immune.py [--force-download] [--overwrite]
"""
from __future__ import annotations

import argparse
import json
import urllib.request
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import scipy.sparse as sp


REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = REPO_ROOT / "data" / "raw"
RAW_OUTPUT = RAW_DIR / "raw_human_immune.h5ad"
METADATA_OUTPUT = RAW_DIR / "raw_human_immune.json"

DATASET_NAME = "human_immune"
CITATION = "Luecken et al. (2022) Nature Methods, doi:10.1038/s41592-021-01336-8"
SPECIES = "Homo sapiens"
TISSUE = "immune (peripheral blood + bone marrow)"
PROTOCOL = "mixed scRNA-seq technologies (10x, Smart-seq2)"

# Use the ndownloader.figshare.com subdomain (not figshare.com/ndownloader/...);
# the latter goes through a CDN that intermittently serves an AWS-WAF JS
# challenge, while the former 302-redirects straight to a presigned S3 URL.
SOURCE_URL = "https://ndownloader.figshare.com/files/25717328"
SOURCE_FILENAME = "Immune_ALL_human.h5ad"
BATCH_KEY = "batch"
LABEL_KEY = "final_annotation"


def download_source(force: bool = False) -> Path:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    out_path = RAW_DIR / SOURCE_FILENAME
    if out_path.exists() and out_path.stat().st_size > 0 and not force:
        print(f"[skip] {out_path}")
        return out_path
    print(f"[download] {SOURCE_URL} -> {out_path}")
    # figshare's CDN sometimes serves an AWS-WAF JS challenge (HTTP 202 with
    # an empty body, header x-amzn-waf-action=challenge) that a CLI client
    # cannot solve. We detect that and fail loudly instead of leaving a
    # zero-byte file behind.
    request = urllib.request.Request(
        SOURCE_URL, headers={"User-Agent": "Mozilla/5.0 (compatible; agentic-embeddings/1.0)"}
    )
    with urllib.request.urlopen(request) as response:
        if response.headers.get("x-amzn-waf-action") == "challenge":
            raise RuntimeError(
                f"figshare returned a WAF challenge for {SOURCE_URL}. "
                f"Download the file manually in a browser and save it to {out_path}, "
                "then re-run this script."
            )
        with out_path.open("wb") as out:
            while chunk := response.read(1 << 20):
                out.write(chunk)
    if out_path.stat().st_size == 0:
        out_path.unlink(missing_ok=True)
        raise RuntimeError(
            f"Downloaded file from {SOURCE_URL} is empty. Likely a transient WAF "
            f"challenge; download manually in a browser to {out_path}."
        )
    return out_path


def as_csr_float32(matrix) -> sp.csr_matrix:
    if sp.issparse(matrix):
        return matrix.tocsr().astype(np.float32)
    return sp.csr_matrix(np.asarray(matrix, dtype=np.float32))


def load_human_immune_adata(force_download: bool = False) -> ad.AnnData:
    source_path = download_source(force=force_download)
    adata = ad.read_h5ad(source_path)
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    if "counts" not in adata.layers:
        raise RuntimeError(
            f"Expected raw counts in adata.layers['counts'] in {source_path}; "
            "the layer is missing. Verify the downloaded file."
        )
    # Fail fast (and helpfully) if the expected obs columns are not present, so
    # downstream QC / preprocessing / metrics never silently run on a wrong key.
    for key in (BATCH_KEY, LABEL_KEY):
        if key not in adata.obs:
            raise RuntimeError(
                f"Expected obs column {key!r} in {source_path}; available columns: "
                f"{sorted(adata.obs.columns.astype(str))}"
            )

    counts = as_csr_float32(adata.layers["counts"])
    adata.X = counts.copy()
    adata.layers["counts"] = counts

    if "gene_symbol" not in adata.var:
        adata.var["gene_symbol"] = adata.var_names.astype(str)

    adata.uns["source"] = {
        "dataset": DATASET_NAME,
        "citation": CITATION,
        "source_url": SOURCE_URL,
        "source_files": [SOURCE_FILENAME],
        "download_date": datetime.today().strftime("%Y-%m-%d"),
        "batch_key": BATCH_KEY,
        "label_key": LABEL_KEY,
    }
    return adata


def build_metadata(adata: ad.AnnData) -> dict:
    batch_counts = (
        adata.obs[BATCH_KEY].astype(str).value_counts().to_dict()
        if BATCH_KEY in adata.obs
        else {}
    )
    celltype_counts = (
        adata.obs[LABEL_KEY].astype(str).value_counts().to_dict()
        if LABEL_KEY in adata.obs
        else {}
    )
    return {
        "dataset": DATASET_NAME,
        "citation": CITATION,
        "species": SPECIES,
        "tissue": TISSUE,
        "protocol": PROTOCOL,
        "n_cells_raw": int(adata.n_obs),
        "n_genes_raw": int(adata.n_vars),
        "source_url": SOURCE_URL,
        "source_files": [SOURCE_FILENAME],
        "download_date": datetime.today().strftime("%Y-%m-%d"),
        "counts_layer": "counts",
        "batch_key": BATCH_KEY,
        "label_key": LABEL_KEY,
        "batch_counts": batch_counts,
        "celltype_counts": celltype_counts,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=f"Load {DATASET_NAME} (scIB immune).")
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if RAW_OUTPUT.exists() and not args.overwrite:
        raise FileExistsError(f"{RAW_OUTPUT} exists. Pass --overwrite to replace it.")

    adata = load_human_immune_adata(force_download=args.force_download)

    RAW_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    METADATA_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(RAW_OUTPUT)
    with METADATA_OUTPUT.open("w") as handle:
        json.dump(build_metadata(adata), handle, indent=2)

    print(f"[ok] wrote {RAW_OUTPUT} ({adata.n_obs} cells x {adata.n_vars} genes)")
    print(f"[ok] wrote {METADATA_OUTPUT}")


if __name__ == "__main__":
    main()
