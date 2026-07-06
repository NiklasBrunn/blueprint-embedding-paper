"""
Load the Pijuan-Sala 2019 mouse gastrulation atlas as one AnnData object.

Source:
    Pijuan-Sala et al. (2019) "A single-cell molecular map of mouse gastrulation
    and early organogenesis", Nature 566:490-495,
    doi:10.1038/s41586-019-0933-9.

Source: Marioni Lab institutional storage (Cambridge / CRUK), the same processed
companion release that Bioconductor's `MouseGastrulationData::EmbryoAtlasData()`
wraps. Hosted at https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz; the
tarball contains raw_counts.mtx (genes x cells, MatrixMarket), genes.tsv,
meta.tab, and barcodes.tsv. The atlas is a single coherent experiment (multiple
embryos sampled across E6.5-E8.5 by 10x Genomics); we drop the doublet- and
stripped-nuclei calls the authors flagged plus any unannotated cells, leaving
the ~116k cells used in the published figures.

obs columns are harmonized to a stable triple:
    celltype   — manually curated annotation (37 classes; meta.tab "celltype")
    stage      — developmental day (E6.5, E6.75, ..., E8.5; canonical in source)
    sample_id  — per-embryo identifier (meta.tab "sample"), the natural batch key

Outputs:
    - data/raw/raw_mouse_gastrulation.h5ad   (raw counts in adata.X and adata.layers["counts"])
    - data/raw/raw_mouse_gastrulation.json   (provenance + per-stage / per-sample counts)

Caches (shared across reruns):
    - data/raw/_marionilab/atlas_data.tar.gz   (1.6 GB; downloaded on first use)
    - data/raw/_marionilab/atlas/{raw_counts.mtx,genes.tsv,meta.tab}

Usage:
    python scripts/IO/load_mouse_gastrulation.py [--overwrite]
"""
from __future__ import annotations

import argparse
import json
import re
import tarfile
import urllib.request
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = REPO_ROOT / "data" / "raw"
RAW_OUTPUT = RAW_DIR / "raw_mouse_gastrulation.h5ad"
METADATA_OUTPUT = RAW_DIR / "raw_mouse_gastrulation.json"

DATASET_NAME = "mouse_gastrulation"
CITATION = "Pijuan-Sala et al. (2019) Nature, doi:10.1038/s41586-019-0933-9"
SPECIES = "Mus musculus"
TISSUE = "gastrulating embryo (E6.5-E8.5)"
PROTOCOL = "10x Genomics Chromium 3'"
COLLECTION_DOI = "10.1038/s41586-019-0933-9"

ATLAS_URL = "https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz"
ATLAS_CACHE = RAW_DIR / "_marionilab"
TARBALL = ATLAS_CACHE / "atlas_data.tar.gz"
EXTRACT_FILES = ("atlas/raw_counts.mtx", "atlas/genes.tsv", "atlas/meta.tab")

BATCH_KEY = "sample_id"
LABEL_KEY = "celltype"


def _parse_stage(value: object) -> str:
    """Extract canonical 'E<x>' label from verbose stage strings (e.g. 'embryonic day 6.5')."""
    if value is None:
        return ""
    s = str(value)
    m = re.search(r"\bE\d+(?:\.\d+)?\b", s)
    if m:
        return m.group(0)
    m = re.search(r"embryonic day\s+(\d+(?:\.\d+)?)", s, re.IGNORECASE)
    if m:
        return f"E{m.group(1)}"
    return s


def _ensure_atlas_files(force_download: bool) -> None:
    ATLAS_CACHE.mkdir(parents=True, exist_ok=True)
    if force_download and TARBALL.exists():
        TARBALL.unlink()
    if not TARBALL.exists():
        print(f"[fetch] {ATLAS_URL} -> {TARBALL}", flush=True)
        urllib.request.urlretrieve(ATLAS_URL, TARBALL)
    missing = [name for name in EXTRACT_FILES if not (ATLAS_CACHE / name).exists()]
    if missing:
        print(f"[extract] {len(missing)} files from {TARBALL.name}", flush=True)
        with tarfile.open(TARBALL, "r:gz") as tf:
            for name in missing:
                tf.extract(name, path=ATLAS_CACHE)


def fetch_from_marionilab(force_download: bool) -> ad.AnnData:
    _ensure_atlas_files(force_download=force_download)
    print("[read] meta.tab", flush=True)
    meta = pd.read_csv(ATLAS_CACHE / "atlas/meta.tab", sep="\t",
                       low_memory=False, na_values=["NA"])
    print("[read] genes.tsv", flush=True)
    genes = pd.read_csv(ATLAS_CACHE / "atlas/genes.tsv", sep="\t",
                        header=None, names=["ensembl_id", "feature_name"])
    print("[read] raw_counts.mtx (genes x cells; transposing)", flush=True)
    counts_gxc = scipy.io.mmread(str(ATLAS_CACHE / "atlas/raw_counts.mtx"))
    counts = counts_gxc.tocsr().T.tocsr().astype(np.float32)

    if counts.shape != (len(meta), len(genes)):
        raise RuntimeError(
            f"shape mismatch: counts={counts.shape} vs "
            f"(cells={len(meta)}, genes={len(genes)})"
        )

    var = pd.DataFrame({
        "ensembl_id": genes["ensembl_id"].astype(str).values,
        "feature_name": genes["feature_name"].astype(str).values,
    })
    var.index = pd.Index(genes["ensembl_id"].astype(str).values, name="ensembl_id")

    obs = meta.copy()
    obs.index = pd.Index(obs["cell"].astype(str).values, name="cell_id")

    adata = ad.AnnData(X=counts, obs=obs, var=var)
    return adata


def _drop_authors_excluded(adata: ad.AnnData) -> ad.AnnData:
    flag = lambda col: adata.obs.get(col, pd.Series(False, index=adata.obs_names)) \
                                  .astype(str).str.upper() == "TRUE"
    keep = adata.obs["celltype"].notna() & ~flag("doublet") & ~flag("stripped")
    n_before = adata.n_obs
    adata = adata[keep.to_numpy(), :].copy()
    print(f"[filter] kept {adata.n_obs}/{n_before} cells "
          "(dropped doublets / stripped nuclei / unannotated)", flush=True)
    return adata


def harmonize_obs(adata: ad.AnnData) -> None:
    obs = adata.obs.copy()

    if "celltype" in obs.columns:
        obs["celltype"] = obs["celltype"].astype(str)
    elif "cell_type" in obs.columns:
        obs["celltype"] = obs["cell_type"].astype(str)
    else:
        raise RuntimeError(
            "Neither 'celltype' nor 'cell_type' present in obs. "
            f"Available columns: {list(obs.columns)}"
        )

    stage_src = next(
        (c for c in ("stage", "development_stage", "developmental_stage") if c in obs.columns),
        None,
    )
    if stage_src is None:
        raise RuntimeError(
            "No developmental-stage column found. "
            f"Available columns: {list(obs.columns)}"
        )
    obs["stage"] = obs[stage_src].astype(str).map(_parse_stage)

    if "sample_id" not in obs.columns:
        for cand in ("donor_id", "sample", "embryo_id", "library_id"):
            if cand in obs.columns:
                obs["sample_id"] = obs[cand].astype(str)
                break
        else:
            raise RuntimeError(
                "No batch column found in obs. Expected one of: "
                "sample_id, donor_id, sample, embryo_id, library_id. "
                f"Available columns: {list(obs.columns)}"
            )
    else:
        obs["sample_id"] = obs["sample_id"].astype(str)

    adata.obs = obs


def harmonize_var(adata: ad.AnnData) -> None:
    """Use gene symbols (not Ensembl IDs) as var_names for compatibility with
    the rest of the pipeline (HGNC/MGI symbol matching in QC, preprocessing,
    scoring functions). MarioniLab genes.tsv supplies symbols in `feature_name`."""
    if "feature_name" in adata.var.columns:
        symbols = adata.var["feature_name"].astype(str)
        adata.var["ensembl_id"] = adata.var_names.astype(str)
        adata.var_names = pd.Index(symbols.values, name=None)
        adata.var_names_make_unique()
    else:
        adata.var["ensembl_id"] = adata.var_names.astype(str)


def stash_counts(adata: ad.AnnData) -> None:
    X = adata.X
    if sp.issparse(X):
        counts = X.tocsr().astype(np.float32)
    else:
        counts = sp.csr_matrix(np.asarray(X, dtype=np.float32))
    adata.X = counts.copy()
    adata.layers["counts"] = counts


def load_adata(force_download: bool = False) -> ad.AnnData:
    adata = fetch_from_marionilab(force_download=force_download)
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    adata = _drop_authors_excluded(adata)
    harmonize_obs(adata)
    harmonize_var(adata)
    stash_counts(adata)

    adata.uns["source"] = {
        "dataset": DATASET_NAME,
        "citation": CITATION,
        "source": "marionilab/atlas_data.tar.gz",
        "collection_doi": COLLECTION_DOI,
        "atlas_url": ATLAS_URL,
        "download_date": datetime.today().strftime("%Y-%m-%d"),
        "batch_key": BATCH_KEY,
        "label_key": LABEL_KEY,
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
        "source": "marionilab/atlas_data.tar.gz",
        "collection_doi": COLLECTION_DOI,
        "atlas_url": ATLAS_URL,
        "download_date": datetime.today().strftime("%Y-%m-%d"),
        "counts_layer": "counts",
        "batch_key": BATCH_KEY,
        "label_key": LABEL_KEY,
        "stage_counts": adata.obs["stage"].astype(str).value_counts().to_dict(),
        "sample_counts": adata.obs[BATCH_KEY].astype(str).value_counts().to_dict(),
        "celltype_counts": adata.obs[LABEL_KEY].astype(str).value_counts().to_dict(),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=f"Load {DATASET_NAME} (Pijuan-Sala 2019) from MarioniLab atlas."
    )
    parser.add_argument("--force-download", action="store_true",
                        help="Re-download atlas_data.tar.gz even if cached.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if RAW_OUTPUT.exists() and not args.overwrite:
        raise FileExistsError(f"{RAW_OUTPUT} exists. Pass --overwrite to replace it.")

    adata = load_adata(force_download=args.force_download)

    RAW_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    METADATA_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(RAW_OUTPUT)
    with METADATA_OUTPUT.open("w") as handle:
        json.dump(build_metadata(adata), handle, indent=2)

    print(f"[ok] wrote {RAW_OUTPUT} ({adata.n_obs} cells x {adata.n_vars} genes)")
    print(f"[ok] wrote {METADATA_OUTPUT}")


if __name__ == "__main__":
    main()
