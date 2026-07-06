#!/usr/bin/env python3
"""Literature-grounded marker embedding for human pancreas scRNA-seq cells.

The scenario_config.jsonc for this run specifies human pancreas with no
restriction to a single cell type. This script therefore uses curated marker
modules for the major healthy human pancreatic epithelial, endocrine, stromal,
vascular, and immune compartments. It does not train a model and does not use
pretrained-model inference.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
from typing import Iterable

import numpy as np
import pandas as pd
from scipy import sparse

try:
    import anndata as ad
except ImportError as exc:  # pragma: no cover - depends on runtime env
    raise SystemExit("embed_freeform.py requires the 'anndata' package") from exc


PANCREAS_MARKER_MODULES: "OrderedDict[str, list[str]]" = OrderedDict(
    [
        (
            "acinar_digestive_enzyme",
            [
                "PRSS1",
                "PRSS2",
                "PRSS3",
                "CPA1",
                "CPA2",
                "CPB1",
                "CTRB1",
                "CTRB2",
                "CELA2A",
                "CELA3A",
                "CELA3B",
                "CEL",
                "PNLIP",
                "CLPS",
                "REG1A",
                "REG1B",
            ],
        ),
        (
            "ductal_epithelial",
            [
                "KRT7",
                "KRT8",
                "KRT18",
                "KRT19",
                "SOX9",
                "EPCAM",
                "MUC1",
                "CFTR",
                "CLDN4",
                "SLC4A4",
                "SPP1",
                "MMP7",
            ],
        ),
        (
            "centroacinar_ductal",
            ["ANXA4", "TSPAN8", "CLDN10", "KRT19", "SOX9", "SPP1", "LCN2"],
        ),
        (
            "endocrine_pan_islet",
            [
                "CHGA",
                "CHGB",
                "SCG5",
                "PCSK1N",
                "NEUROD1",
                "ISL1",
                "PAX6",
                "NKX2-2",
                "RFX6",
                "FEV",
            ],
        ),
        (
            "beta_cell",
            [
                "INS",
                "IAPP",
                "MAFA",
                "PDX1",
                "NKX6-1",
                "PCSK1",
                "SLC30A8",
                "G6PC2",
                "UCN3",
                "HADH",
            ],
        ),
        (
            "alpha_cell",
            ["GCG", "ARX", "IRX2", "MAFB", "PCSK2", "TTR", "TM4SF4", "LOXL4"],
        ),
        ("delta_cell", ["SST", "HHEX", "LEPR", "RBP4", "BCHE", "GHSR"]),
        ("pp_gamma_cell", ["PPY", "PYY", "FAP", "SERTM1", "SLITRK6"]),
        ("epsilon_cell", ["GHRL", "MLN", "KCNK16", "ISL1"]),
        (
            "stellate_fibroblast_matrix",
            [
                "COL1A1",
                "COL1A2",
                "COL3A1",
                "COL6A1",
                "DCN",
                "LUM",
                "MGP",
                "PDGFRA",
                "FBLN1",
                "THY1",
            ],
        ),
        (
            "activated_stellate_smooth_muscle",
            ["ACTA2", "TAGLN", "MYL9", "MYH11", "TPM2", "CNN1", "COL1A1", "CTGF"],
        ),
        (
            "endothelial",
            ["PECAM1", "VWF", "KDR", "FLT1", "ENG", "ESAM", "RAMP2", "PLVAP"],
        ),
        (
            "pericyte",
            ["RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3", "ABCC9", "KCNJ8"],
        ),
        (
            "macrophage_monocyte",
            ["LYZ", "C1QA", "C1QB", "C1QC", "CD68", "MSR1", "MRC1", "FCGR3A"],
        ),
        (
            "dendritic_antigen_presenting",
            ["HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "FCER1A", "CLEC10A", "LILRA4"],
        ),
        ("t_cell", ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "IL7R", "CCR7"]),
        ("nk_cell", ["NKG7", "GNLY", "KLRD1", "KLRF1", "PRF1", "GZMB", "GZMA"]),
        ("b_cell", ["MS4A1", "CD79A", "CD79B", "BANK1", "CD22", "PAX5"]),
        ("plasma_cell", ["MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1", "IGKC"]),
        ("mast_cell", ["TPSAB1", "TPSB2", "CPA3", "KIT", "HDC", "MS4A2"]),
        ("cycling", ["MKI67", "TOP2A", "UBE2C", "HMGB2", "CENPF", "STMN1"]),
        (
            "epithelial_identity",
            ["EPCAM", "KRT8", "KRT18", "KRT19", "KRT7", "CLDN4", "MUC1"],
        ),
        (
            "islet_maturity_secretory",
            ["SLC30A8", "UCN3", "MAFA", "PCSK1", "PCSK2", "CHGA", "CHGB", "SCGN"],
        ),
        (
            "inflammatory_interferon",
            ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "STAT1", "CXCL10"],
        ),
    ]
)


CONTRASTS: tuple[tuple[str, str, str], ...] = (
    ("endocrine_vs_exocrine", "endocrine_pan_islet", "acinar_digestive_enzyme"),
    ("beta_vs_alpha", "beta_cell", "alpha_cell"),
    ("beta_vs_delta", "beta_cell", "delta_cell"),
    ("alpha_vs_delta", "alpha_cell", "delta_cell"),
    ("minor_endocrine_vs_major", "pp_gamma_cell", "beta_cell"),
    ("acinar_vs_ductal", "acinar_digestive_enzyme", "ductal_epithelial"),
    ("ductal_vs_acinar", "ductal_epithelial", "acinar_digestive_enzyme"),
    ("epithelial_vs_matrix", "epithelial_identity", "stellate_fibroblast_matrix"),
    ("vascular_vs_matrix", "endothelial", "stellate_fibroblast_matrix"),
    ("pericyte_vs_endothelial", "pericyte", "endothelial"),
    ("immune_vs_epithelial", "macrophage_monocyte", "epithelial_identity"),
    ("lymphoid_vs_myeloid", "t_cell", "macrophage_monocyte"),
    ("nk_vs_t", "nk_cell", "t_cell"),
    ("plasma_vs_b", "plasma_cell", "b_cell"),
    ("antigen_presenting_vs_macrophage", "dendritic_antigen_presenting", "macrophage_monocyte"),
    ("cycling_vs_secretory_epithelium", "cycling", "epithelial_identity"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a fixed literature-marker embedding for human pancreas cells."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file")
    parser.add_argument("--output", required=True, help="Output parquet path")
    parser.add_argument("--seed", type=int, required=True, help="Run tag; accepted but unused")
    return parser.parse_args()


def as_count_matrix(adata: ad.AnnData):
    if "counts" not in adata.layers:
        raise KeyError('Expected raw counts in adata.layers["counts"]')
    counts = adata.layers["counts"]
    if sparse.issparse(counts):
        return counts.tocsr()
    return np.asarray(counts)


def gene_index(var_names: Iterable[str]) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for i, gene in enumerate(var_names):
        index.setdefault(str(gene).upper(), []).append(i)
    return index


def marker_columns(markers: Iterable[str], index: dict[str, list[int]]) -> list[int]:
    cols: list[int] = []
    seen: set[int] = set()
    for marker in markers:
        for col in index.get(marker.upper(), []):
            if col not in seen:
                cols.append(col)
                seen.add(col)
    return cols


def row_sum(matrix, cols: list[int]) -> np.ndarray:
    n_obs = matrix.shape[0]
    if not cols:
        return np.zeros(n_obs, dtype=np.float32)
    if sparse.issparse(matrix):
        return np.asarray(matrix[:, cols].sum(axis=1)).ravel().astype(np.float32)
    return np.asarray(matrix[:, cols].sum(axis=1)).ravel().astype(np.float32)


def module_score(counts, cols: list[int], library_size: np.ndarray) -> np.ndarray:
    marker_counts = row_sum(counts, cols)
    normalized = (marker_counts * 10000.0) / library_size
    return np.log1p(normalized).astype(np.float32)


def centered_ratio(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return ((left - right) / (left + right + 1.0)).astype(np.float32)


def build_embedding(adata: ad.AnnData) -> pd.DataFrame:
    counts = as_count_matrix(adata)
    if sparse.issparse(counts):
        library_size = np.asarray(counts.sum(axis=1)).ravel().astype(np.float32)
    else:
        library_size = np.asarray(counts.sum(axis=1)).ravel().astype(np.float32)
    library_size = np.maximum(library_size, 1.0)

    index = gene_index(adata.var_names)
    scores: "OrderedDict[str, np.ndarray]" = OrderedDict()
    for module_name, markers in PANCREAS_MARKER_MODULES.items():
        cols = marker_columns(markers, index)
        scores[module_name] = module_score(counts, cols, library_size)

    dims: list[np.ndarray] = []
    dim_names: list[str] = []

    for module_name, values in scores.items():
        dim_names.append(module_name)
        dims.append(values)

    for contrast_name, left_name, right_name in CONTRASTS:
        dim_names.append(contrast_name)
        dims.append(centered_ratio(scores[left_name], scores[right_name]))

    # A fixed aggregate endocrine-composition axis helps separate rare islet
    # hormone-producing populations without fitting a data-dependent basis.
    minor_endocrine = (
        scores["delta_cell"]
        + scores["pp_gamma_cell"]
        + scores["epsilon_cell"]
    ) / 3.0
    major_endocrine = (scores["beta_cell"] + scores["alpha_cell"]) / 2.0
    dim_names.append("minor_vs_major_endocrine")
    dims.append(centered_ratio(minor_endocrine, major_endocrine))

    embedding = np.vstack(dims).T.astype(np.float32)
    out = pd.DataFrame({"cell_id": adata.obs_names.astype(str)})
    for i in range(embedding.shape[1]):
        out[f"dim_{i}"] = embedding[:, i]
    return out


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.input)
    embedding = build_embedding(adata)
    embedding.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
