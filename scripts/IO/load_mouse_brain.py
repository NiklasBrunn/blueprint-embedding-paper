"""
Load the Tasic 2016 mouse primary visual cortex (V1) SMART-Seq dataset.

Source:
    Tasic et al. (2016) "Adult mouse cortical cell taxonomy revealed by single
    cell transcriptomics", Nature Neuroscience, doi:10.1038/nn.4216.

Source: GEO accession GSE71585 supplementary files
    - GSE71585_RefSeq_counts.csv.gz
    - GSE71585_Clustering_Results.csv.gz

Outputs:
    - data/raw/raw_mouse_brain.h5ad   (raw counts in adata.X and adata.layers["counts"])
    - data/raw/raw_mouse_brain.json   (dataset metadata and provenance, sibling of the h5ad)

Usage:
    python scripts/IO/load_mouse_brain.py [--force-download] [--overwrite]
"""
from __future__ import annotations

import argparse
import json
import urllib.request
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = REPO_ROOT / "data" / "raw"
RAW_OUTPUT = RAW_DIR / "raw_mouse_brain.h5ad"
METADATA_OUTPUT = RAW_DIR / "raw_mouse_brain.json"

DATASET_NAME = "mouse_brain"
CITATION = "Tasic et al. (2016) Nature Neuroscience, doi:10.1038/nn.4216"
SPECIES = "Mus musculus"
TISSUE = "primary visual cortex (V1)"
PROTOCOL = "SMART-Seq (single-cell, full-length)"

GEO_BASE = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE71585&format=file&file="
COUNTS_FILENAME = "GSE71585_RefSeq_counts.csv.gz"
CLUSTERING_FILENAME = "GSE71585_Clustering_Results.csv.gz"
COUNTS_URL = f"{GEO_BASE}{COUNTS_FILENAME}"
CLUSTERING_URL = f"{GEO_BASE}{CLUSTERING_FILENAME}"


def download_source(force: bool = False) -> tuple[Path, Path]:
    RAW_DIR.mkdir(parents=True, exist_ok=True)
    counts_path = RAW_DIR / COUNTS_FILENAME
    clustering_path = RAW_DIR / CLUSTERING_FILENAME
    for url, path in ((COUNTS_URL, counts_path), (CLUSTERING_URL, clustering_path)):
        if path.exists() and not force:
            print(f"[skip] {path}")
            continue
        print(f"[download] {url} -> {path}")
        urllib.request.urlretrieve(url, path)
    return counts_path, clustering_path


def _derive_subclass(primary_type: pd.Series) -> pd.Series:
    # The Tasic taxonomy "subclass" is the leading whitespace-separated token of
    # primary_type (e.g. "Vip Mybpc1" -> "Vip", "L5a Hsd11b1" -> "L5a"). It is
    # not in the source CSV, so we derive it here once at load time.
    subclass = primary_type.fillna("Unclassified").astype(str).str.split().str[0]
    return subclass.replace({"nan": "Unclassified"})


def parse_geo_to_anndata(counts_path: Path, clustering_path: Path) -> ad.AnnData:
    counts_df = pd.read_csv(counts_path, compression="gzip")
    gene_names = pd.Index(counts_df.pop("gene").astype(str), name="gene")
    sample_ids = counts_df.columns.astype(str).tolist()

    metadata = pd.read_csv(clustering_path, compression="gzip")
    metadata = metadata.rename(columns={"sample_title": "sample_id"}).set_index("sample_id")
    metadata = metadata.loc[sample_ids].copy()
    metadata.index = metadata.index.astype(str)
    metadata.index.name = "cell_id"
    metadata["sample_id"] = metadata.index
    metadata["secondary_type"] = metadata["secondary_type"].fillna("")
    metadata["core_intermediate"] = metadata["core_intermediate"].fillna("unclassified")
    metadata["cluster_id"] = metadata["primary_type"].fillna("Unclassified").astype(str)
    metadata["subclass"] = _derive_subclass(metadata["primary_type"])

    # Drop cells the Tasic 2016 taxonomy could not place ("Unclassified" in
    # broad_type / primary_type / subclass — same 82-cell set). They lack a
    # confident cell-type assignment, so they only add noise to downstream
    # cluster-quality metrics and UMAPs.
    keep = metadata["broad_type"].astype(str) != "Unclassified"
    n_dropped = int((~keep).sum())
    metadata = metadata.loc[keep].copy()
    counts_df = counts_df.loc[:, metadata.index]
    print(f"[filter] dropped {n_dropped} cells labelled 'Unclassified'; kept {metadata.shape[0]}")

    matrix = counts_df.to_numpy(dtype=np.float32, copy=False).T
    adata = ad.AnnData(
        X=sp.csr_matrix(matrix),
        obs=metadata,
        var=pd.DataFrame(index=gene_names),
        dtype=np.float32,
    )
    adata.var["gene_symbol"] = adata.var_names
    adata.layers["counts"] = adata.X.copy()
    return adata


def load_mouse_brain_adata(force_download: bool = False) -> ad.AnnData:
    counts_path, clustering_path = download_source(force=force_download)
    adata = parse_geo_to_anndata(counts_path, clustering_path)
    adata.uns["source"] = {
        "dataset": DATASET_NAME,
        "citation": CITATION,
        "source_url": [COUNTS_URL, CLUSTERING_URL],
        "source_files": [COUNTS_FILENAME, CLUSTERING_FILENAME],
        "download_date": datetime.today().strftime("%Y-%m-%d"),
    }
    return adata


def build_metadata(adata: ad.AnnData) -> dict:
    return {
        "dataset": DATASET_NAME,
        "citation": CITATION,
        "species": SPECIES,
        "tissue": TISSUE,
        "protocol": PROTOCOL,
        "n_cells_raw": int(adata.n_obs),
        "n_genes_raw": int(adata.n_vars),
        "source_url": [COUNTS_URL, CLUSTERING_URL],
        "source_files": [COUNTS_FILENAME, CLUSTERING_FILENAME],
        "download_date": datetime.today().strftime("%Y-%m-%d"),
        "counts_layer": "counts",
        "label_keys": ["broad_type", "subclass", "primary_type", "cluster_id"],
        "label_derivations": {
            "subclass": "first whitespace-separated token of primary_type "
                        "(e.g. 'Vip Mybpc1' -> 'Vip', 'L5a Hsd11b1' -> 'L5a'); "
                        "NaNs become 'Unclassified'.",
            "cluster_id": "primary_type with NaNs filled as 'Unclassified', cast to str.",
        },
        "cell_filter": "drop rows where broad_type == 'Unclassified' "
                       "(no confident Tasic taxonomy assignment).",
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=f"Load {DATASET_NAME} (Tasic 2016).")
    parser.add_argument("--force-download", action="store_true")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if RAW_OUTPUT.exists() and not args.overwrite:
        raise FileExistsError(f"{RAW_OUTPUT} exists. Pass --overwrite to replace it.")

    adata = load_mouse_brain_adata(force_download=args.force_download)

    RAW_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    METADATA_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(RAW_OUTPUT)
    with METADATA_OUTPUT.open("w") as handle:
        json.dump(build_metadata(adata), handle, indent=2)

    print(f"[ok] wrote {RAW_OUTPUT} ({adata.n_obs} cells x {adata.n_vars} genes)")
    print(f"[ok] wrote {METADATA_OUTPUT}")


if __name__ == "__main__":
    main()
