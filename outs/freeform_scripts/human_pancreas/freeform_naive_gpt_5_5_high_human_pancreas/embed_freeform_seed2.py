#!/usr/bin/env python3
"""Literature-marker embedding for human pancreas scRNA-seq cells.

The embedding is intentionally naive: each coordinate is a fixed,
literature-grounded marker/module score computed from observed counts.
It does not fit a model to the dataset and does not call any pretrained
model.  The scenario_config.jsonc for this run specifies human pancreas,
with all tissue cell types in scope.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


# Human pancreas markers from canonical scRNA-seq/cell-atlas literature:
# endocrine hormones and transcription factors, exocrine digestive enzymes
# and duct genes, plus stromal/vascular/immune/neural compartment markers.
MARKER_MODULES: "OrderedDict[str, list[str]]" = OrderedDict(
    [
        (
            "beta_cell_insulin",
            [
                "INS",
                "IAPP",
                "PCSK1",
                "PCSK2",
                "G6PC2",
                "SLC30A8",
                "MAFA",
                "PDX1",
                "NKX6-1",
                "NKX2-2",
                "PAX6",
                "ABCC8",
                "KCNJ11",
                "HADH",
            ],
        ),
        (
            "alpha_cell_glucagon",
            [
                "GCG",
                "TTR",
                "ARX",
                "IRX2",
                "MAFB",
                "PCSK2",
                "PAX6",
                "LOXL4",
                "TM4SF4",
                "GC",
            ],
        ),
        (
            "delta_cell_somatostatin",
            ["SST", "HHEX", "RBP4", "LEPR", "PCSK1", "PCSK2", "PAX6", "PRG4"],
        ),
        (
            "pp_cell_pancreatic_polypeptide",
            ["PPY", "PYY", "SLC18A2", "REG1A", "REG1B", "CHGA", "CHGB", "PAX6"],
        ),
        (
            "epsilon_cell_ghrelin",
            ["GHRL", "IAPP", "CHGA", "CHGB", "PCSK1", "PAX6"],
        ),
        (
            "endocrine_secretory_granule",
            [
                "CHGA",
                "CHGB",
                "SCG2",
                "SCG3",
                "SCG5",
                "VGF",
                "PCSK1N",
                "SNAP25",
                "SYT7",
                "SYP",
                "NEUROD1",
                "ISL1",
                "PAX6",
            ],
        ),
        (
            "beta_maturity_function",
            [
                "MAFA",
                "UCN3",
                "SLC2A2",
                "GCK",
                "G6PC2",
                "SLC30A8",
                "ABCC8",
                "KCNJ11",
                "PCSK1",
                "HADH",
                "GLP1R",
            ],
        ),
        (
            "immature_endocrine_progenitor",
            [
                "SOX9",
                "HES1",
                "NEUROG3",
                "RFX6",
                "NKX2-2",
                "PAX4",
                "MNX1",
                "ONECUT1",
                "SOX4",
                "INSM1",
            ],
        ),
        (
            "acinar_digestive_enzyme",
            [
                "PRSS1",
                "PRSS2",
                "PRSS3",
                "CTRB1",
                "CTRB2",
                "CPA1",
                "CPA2",
                "CPB1",
                "CELA3A",
                "CELA3B",
                "REG1A",
                "REG1B",
                "REG3A",
                "PNLIP",
                "CLPS",
                "AMY2A",
                "AMY2B",
            ],
        ),
        (
            "acinar_secretory_stress",
            [
                "XBP1",
                "PDIA2",
                "HSPA5",
                "HYOU1",
                "ERP27",
                "AGR2",
                "AGR3",
                "MIST1",
                "BHLHA15",
                "RBPJL",
            ],
        ),
        (
            "ductal_epithelium",
            [
                "KRT19",
                "KRT7",
                "KRT8",
                "KRT18",
                "EPCAM",
                "MUC1",
                "CFTR",
                "SLC4A4",
                "SLC26A3",
                "SLC26A6",
                "TFF1",
                "TFF2",
                "TFF3",
                "SOX9",
                "ANXA4",
            ],
        ),
        (
            "centroacinar_ductal",
            ["CFTR", "SOX9", "KRT19", "SPP1", "PROM1", "MUC1", "ANXA4", "CLDN10"],
        ),
        (
            "tuft_chemosensory",
            ["POU2F3", "TRPM5", "DCLK1", "GNAT3", "AVIL", "LRMP", "SH2D6", "IL25"],
        ),
        (
            "stellate_fibroblast_matrix",
            [
                "COL1A1",
                "COL1A2",
                "COL3A1",
                "COL5A1",
                "DCN",
                "LUM",
                "PDGFRA",
                "PDGFRB",
                "ACTA2",
                "TAGLN",
                "RGS5",
                "SPARC",
                "MMP2",
            ],
        ),
        (
            "endothelial",
            [
                "PECAM1",
                "VWF",
                "KDR",
                "FLT1",
                "EMCN",
                "CDH5",
                "ESAM",
                "ENG",
                "RAMP2",
                "PLVAP",
                "CLEC14A",
            ],
        ),
        (
            "pericyte_smooth_muscle",
            ["RGS5", "PDGFRB", "CSPG4", "MCAM", "ACTA2", "TAGLN", "MYH11", "NOTCH3"],
        ),
        (
            "immune_pan_leukocyte",
            ["PTPRC", "LST1", "TYROBP", "CORO1A", "CXCR4", "CD52", "HLA-DRA"],
        ),
        (
            "myeloid_macrophage_dc",
            [
                "LYZ",
                "LST1",
                "AIF1",
                "TYROBP",
                "C1QA",
                "C1QB",
                "C1QC",
                "FCGR3A",
                "MS4A7",
                "CD14",
                "FCER1A",
                "CLEC10A",
                "ITGAX",
            ],
        ),
        (
            "t_nk_cell",
            [
                "CD3D",
                "CD3E",
                "CD2",
                "TRAC",
                "IL7R",
                "CCR7",
                "NKG7",
                "GNLY",
                "PRF1",
                "GZMB",
                "KLRD1",
            ],
        ),
        (
            "b_plasma_cell",
            [
                "MS4A1",
                "CD79A",
                "CD79B",
                "BANK1",
                "MZB1",
                "JCHAIN",
                "SDC1",
                "IGHG1",
                "IGKC",
                "XBP1",
            ],
        ),
        (
            "mast_cell",
            ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC", "HPGDS"],
        ),
        (
            "schwann_neural",
            ["S100B", "SOX10", "PLP1", "MPZ", "GFRA3", "NRXN1", "NCAM1", "NGFR"],
        ),
        (
            "cycling",
            [
                "MKI67",
                "TOP2A",
                "UBE2C",
                "CENPF",
                "STMN1",
                "PCNA",
                "MCM2",
                "MCM5",
                "TYMS",
            ],
        ),
        (
            "interferon_inflammatory",
            [
                "ISG15",
                "IFIT1",
                "IFIT2",
                "IFIT3",
                "IFI6",
                "IFI27",
                "MX1",
                "OAS1",
                "OAS2",
                "STAT1",
                "CXCL10",
            ],
        ),
        (
            "hypoxia_stress",
            ["HIF1A", "VEGFA", "BNIP3", "NDRG1", "DDIT4", "SLC2A1", "ENO1", "LDHA"],
        ),
        (
            "apoptosis_damage",
            ["FOS", "JUN", "DUSP1", "ATF3", "BTG2", "IER2", "DNAJB1", "HSPA1A", "HSPA1B"],
        ),
        (
            "ribosomal_translation",
            [
                "RPL3",
                "RPL5",
                "RPL7",
                "RPL10",
                "RPL13A",
                "RPLP0",
                "RPS3",
                "RPS6",
                "RPS12",
                "RPS18",
            ],
        ),
        (
            "mitochondrial_activity",
            [
                "MT-CO1",
                "MT-CO2",
                "MT-CO3",
                "MT-CYB",
                "MT-ND1",
                "MT-ND2",
                "MT-ND3",
                "MT-ND4",
                "MT-ATP6",
            ],
        ),
    ]
)


CONTRASTS: "OrderedDict[str, tuple[str, str]]" = OrderedDict(
    [
        ("endocrine_over_exocrine", ("endocrine_secretory_granule", "acinar_digestive_enzyme")),
        ("beta_over_alpha", ("beta_cell_insulin", "alpha_cell_glucagon")),
        ("duct_over_acinar", ("ductal_epithelium", "acinar_digestive_enzyme")),
        ("stromal_over_epithelial", ("stellate_fibroblast_matrix", "ductal_epithelium")),
        ("vascular_over_stromal", ("endothelial", "stellate_fibroblast_matrix")),
        ("immune_over_parenchymal", ("immune_pan_leukocyte", "endocrine_secretory_granule")),
        ("maturity_over_progenitor", ("beta_maturity_function", "immature_endocrine_progenitor")),
        ("stress_over_secretory", ("apoptosis_damage", "endocrine_secretory_granule")),
    ]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Naive literature-marker embedding for human pancreas cells.")
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", required=True, type=int, help="Run tag; accepted for interface compatibility.")
    return parser.parse_args()


def _as_1d_array(x) -> np.ndarray:
    if sparse.issparse(x):
        x = x.toarray()
    return np.asarray(x).reshape(-1)


def _matrix_subset(matrix, idx: list[int]):
    if sparse.issparse(matrix):
        return matrix[:, idx]
    return np.asarray(matrix)[:, idx]


def _available_indices(var_names: Iterable[str], genes: Iterable[str]) -> list[int]:
    upper_to_index: dict[str, int] = {}
    for i, name in enumerate(var_names):
        key = str(name).upper()
        if key not in upper_to_index:
            upper_to_index[key] = i
    return [upper_to_index[g.upper()] for g in genes if g.upper() in upper_to_index]


def counts_to_log_cpm(counts) -> tuple[object, np.ndarray]:
    """Return log1p(CP10K) and per-cell library sizes."""
    lib_size = _as_1d_array(counts.sum(axis=1)).astype(np.float64)
    scale = np.divide(10000.0, lib_size, out=np.zeros_like(lib_size), where=lib_size > 0)
    if sparse.issparse(counts):
        norm = counts.astype(np.float64).multiply(scale[:, None]).tocsr(copy=False)
        norm.data = np.log1p(norm.data)
        return norm, lib_size
    norm = np.asarray(counts, dtype=np.float64) * scale[:, None]
    return np.log1p(norm), lib_size


def module_score(log_cpm, indices: list[int], n_cells: int) -> np.ndarray:
    """Mean log-normalized expression over available marker genes."""
    if not indices:
        return np.zeros(n_cells, dtype=np.float32)
    sub = _matrix_subset(log_cpm, indices)
    score = _as_1d_array(sub.mean(axis=1)).astype(np.float64)
    return np.log1p(score).astype(np.float32)


def robust_standardize(values: np.ndarray) -> np.ndarray:
    """Compress heavy-tailed module scores without learning dataset structure."""
    values = np.asarray(values, dtype=np.float64)
    median = np.median(values)
    q25, q75 = np.percentile(values, [25, 75])
    iqr = q75 - q25
    if not np.isfinite(iqr) or iqr < 1e-8:
        std = np.std(values)
        scale = std if std > 1e-8 else 1.0
    else:
        scale = iqr / 1.349
    z = (values - median) / scale
    return np.tanh(z / 3.0).astype(np.float32)


def build_embedding(adata: ad.AnnData) -> pd.DataFrame:
    if "counts" not in adata.layers:
        raise KeyError('Input AnnData must contain raw counts in layers["counts"].')

    counts = adata.layers["counts"]
    n_cells = adata.n_obs
    log_cpm, lib_size = counts_to_log_cpm(counts)

    raw_scores: dict[str, np.ndarray] = {}
    standardized: dict[str, np.ndarray] = {}
    for name, genes in MARKER_MODULES.items():
        idx = _available_indices(adata.var_names, genes)
        raw = module_score(log_cpm, idx, n_cells)
        raw_scores[name] = raw
        standardized[name] = robust_standardize(raw)

    features: list[np.ndarray] = []
    feature_names: list[str] = []

    for name in MARKER_MODULES:
        features.append(standardized[name])
        feature_names.append(name)

    for name, (positive, negative) in CONTRASTS.items():
        contrast = standardized[positive] - standardized[negative]
        features.append(np.clip(contrast, -2.0, 2.0).astype(np.float32))
        feature_names.append(name)

    detected_genes = _as_1d_array((counts > 0).sum(axis=1)).astype(np.float64)
    qc_features = OrderedDict(
        [
            ("count_depth", np.log1p(lib_size)),
            ("detected_gene_count", np.log1p(detected_genes)),
        ]
    )
    for name, values in qc_features.items():
        features.append(robust_standardize(values))
        feature_names.append(name)

    matrix = np.vstack(features).T.astype(np.float32, copy=False)
    data = {"cell_id": adata.obs_names.astype(str).to_numpy()}
    for i in range(matrix.shape[1]):
        data[f"dim_{i}"] = matrix[:, i]

    # Keep feature_names local for maintainability; the required output schema is
    # intentionally only cell_id plus dim_* columns.
    _ = feature_names
    return pd.DataFrame(data)


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.input)
    embedding = build_embedding(adata)
    embedding.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
