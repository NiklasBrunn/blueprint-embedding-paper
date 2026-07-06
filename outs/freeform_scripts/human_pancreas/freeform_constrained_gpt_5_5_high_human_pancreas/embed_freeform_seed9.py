#!/usr/bin/env python3
"""Literature-derived embedding for healthy human pancreas scRNA-seq."""

from __future__ import annotations

import argparse
from collections.abc import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


# Dimension map:
# dim_0 endocrine-vs-exocrine/stromal identity, dim_1 beta, dim_2 alpha,
# dim_3 delta, dim_4 pancreatic-polypeptide/gamma, dim_5 epsilon,
# dim_6 acinar, dim_7 ductal/centroacinar, dim_8 REG3A+ acinar stress,
# dim_9 stellate/fibroblast, dim_10 endothelial, dim_11 pericyte/smooth muscle,
# dim_12 myeloid, dim_13 T/NK lymphoid, dim_14 B/plasma, dim_15 cycling.


# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Muraro et al., A Single-Cell Transcriptome Atlas
# of the Human Pancreas, DOI 10.1016/j.cels.2016.09.002; Segerstolpe et al.,
# Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and
# Type 2 Diabetes, DOI 10.1016/j.cmet.2016.08.020.
ENDOCRINE_CORE = [
    "CHGA",
    "CHGB",
    "SCG2",
    "SCG3",
    "PCSK1",
    "PCSK2",
    "ISL1",
    "PAX6",
    "NEUROD1",
    "SLC18A2",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Muraro et al., A Single-Cell Transcriptome Atlas
# of the Human Pancreas, DOI 10.1016/j.cels.2016.09.002; Segerstolpe et al.,
# Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and
# Type 2 Diabetes, DOI 10.1016/j.cmet.2016.08.020.
BETA = [
    "INS",
    "IAPP",
    "MAFA",
    "PDX1",
    "NKX6-1",
    "SLC30A8",
    "PCSK1",
    "G6PC2",
    "DLK1",
    "HADH",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Muraro et al., A Single-Cell Transcriptome Atlas
# of the Human Pancreas, DOI 10.1016/j.cels.2016.09.002; Segerstolpe et al.,
# Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and
# Type 2 Diabetes, DOI 10.1016/j.cmet.2016.08.020.
ALPHA = [
    "GCG",
    "ARX",
    "IRX2",
    "MAFB",
    "TTR",
    "LOX",
    "PCSK2",
    "TM4SF4",
    "CRYBA2",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Muraro et al., A Single-Cell Transcriptome Atlas
# of the Human Pancreas, DOI 10.1016/j.cels.2016.09.002; Segerstolpe et al.,
# Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and
# Type 2 Diabetes, DOI 10.1016/j.cmet.2016.08.020.
DELTA = [
    "SST",
    "HHEX",
    "RBP4",
    "LEPR",
    "BCHE",
    "PCSK1",
    "PRG4",
    "GHSR",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Muraro et al., A Single-Cell Transcriptome Atlas
# of the Human Pancreas, DOI 10.1016/j.cels.2016.09.002; Segerstolpe et al.,
# Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and
# Type 2 Diabetes, DOI 10.1016/j.cmet.2016.08.020.
PANCREATIC_POLYPEPTIDE = [
    "PPY",
    "PNOC",
    "PYY",
    "MEIS2",
    "FAP",
    "TMEM176B",
    "SLC18A2",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Muraro et al., A Single-Cell Transcriptome Atlas
# of the Human Pancreas, DOI 10.1016/j.cels.2016.09.002; Segerstolpe et al.,
# Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and
# Type 2 Diabetes, DOI 10.1016/j.cmet.2016.08.020.
EPSILON = [
    "GHRL",
    "MLN",
    "RBP4",
    "ISL1",
    "SLC18A2",
    "PCSK1",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Muraro et al., A Single-Cell Transcriptome Atlas
# of the Human Pancreas, DOI 10.1016/j.cels.2016.09.002; Tabula Sapiens Consortium,
# The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans,
# DOI 10.1126/science.abl4896.
ACINAR = [
    "PRSS1",
    "PRSS2",
    "CPA1",
    "CPA2",
    "CPB1",
    "CTRB1",
    "CTRB2",
    "CTRC",
    "CELA3A",
    "CELA3B",
    "AMY2A",
    "AMY2B",
    "PNLIP",
    "REG1A",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Muraro et al., A Single-Cell Transcriptome Atlas
# of the Human Pancreas, DOI 10.1016/j.cels.2016.09.002; Tabula Sapiens Consortium,
# The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans,
# DOI 10.1126/science.abl4896.
DUCTAL_CENTROACINAR = [
    "KRT19",
    "KRT7",
    "KRT8",
    "KRT18",
    "EPCAM",
    "SOX9",
    "MUC1",
    "CFTR",
    "SLC4A4",
    "ANXA4",
    "SPP1",
    "TFF1",
]

# Muraro et al., A Single-Cell Transcriptome Atlas of the Human Pancreas,
# DOI 10.1016/j.cels.2016.09.002; Baron et al., A Single-Cell Transcriptomic Map
# of the Human and Mouse Pancreas, DOI 10.1016/j.cels.2016.08.011.
REG3A_ACINAR_STRESS = [
    "REG3A",
    "REG3G",
    "REG1A",
    "REG1B",
    "SPINK1",
    "CTRB1",
    "PRSS2",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Tabula Sapiens Consortium, The Tabula Sapiens:
# A multiple-organ, single-cell transcriptomic atlas of humans, DOI
# 10.1126/science.abl4896.
STELLATE_FIBROBLAST = [
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "PDGFRA",
    "PDGFRB",
    "RGS5",
    "ACTA2",
    "MMP2",
    "TAGLN",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Tabula Sapiens Consortium, The Tabula Sapiens:
# A multiple-organ, single-cell transcriptomic atlas of humans, DOI
# 10.1126/science.abl4896.
ENDOTHELIAL = [
    "PECAM1",
    "VWF",
    "KDR",
    "FLT1",
    "EMCN",
    "RAMP2",
    "PLVAP",
    "CDH5",
    "ESAM",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Tabula Sapiens Consortium, The Tabula Sapiens:
# A multiple-organ, single-cell transcriptomic atlas of humans, DOI
# 10.1126/science.abl4896.
PERICYTE_SMOOTH_MUSCLE = [
    "RGS5",
    "PDGFRB",
    "CSPG4",
    "MCAM",
    "NOTCH3",
    "ACTA2",
    "TAGLN",
    "MYH11",
    "MYLK",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Tabula Sapiens Consortium, The Tabula Sapiens:
# A multiple-organ, single-cell transcriptomic atlas of humans, DOI
# 10.1126/science.abl4896.
MYELOID = [
    "LYZ",
    "AIF1",
    "CSF1R",
    "C1QA",
    "C1QB",
    "C1QC",
    "CD68",
    "FCGR3A",
    "MS4A7",
    "HLA-DRA",
    "HLA-DRB1",
    "CD74",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Tabula Sapiens Consortium, The Tabula Sapiens:
# A multiple-organ, single-cell transcriptomic atlas of humans, DOI
# 10.1126/science.abl4896.
T_NK_LYMPHOID = [
    "CD3D",
    "CD3E",
    "CD2",
    "TRAC",
    "IL7R",
    "NKG7",
    "GNLY",
    "GZMB",
    "PRF1",
    "KLRD1",
]

# Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas,
# DOI 10.1016/j.cels.2016.08.011; Tabula Sapiens Consortium, The Tabula Sapiens:
# A multiple-organ, single-cell transcriptomic atlas of humans, DOI
# 10.1126/science.abl4896.
B_PLASMA = [
    "MS4A1",
    "CD79A",
    "CD79B",
    "BANK1",
    "MZB1",
    "JCHAIN",
    "XBP1",
    "IGHG1",
    "IGHM",
    "SDC1",
]

# MSigDB HALLMARK_G2M_CHECKPOINT and HALLMARK_E2F_TARGETS; Reactome R-HSA-69278
# Cell Cycle, Mitotic.
CYCLING = [
    "MKI67",
    "TOP2A",
    "UBE2C",
    "CENPF",
    "NUSAP1",
    "BIRC5",
    "HMGB2",
    "TYMS",
    "PCNA",
    "MCM2",
    "MCM5",
    "CDK1",
    "CCNB1",
    "CCNB2",
]


MODULES = {
    "endocrine": ENDOCRINE_CORE,
    "beta": BETA,
    "alpha": ALPHA,
    "delta": DELTA,
    "pp": PANCREATIC_POLYPEPTIDE,
    "epsilon": EPSILON,
    "acinar": ACINAR,
    "ductal": DUCTAL_CENTROACINAR,
    "reg3a_stress": REG3A_ACINAR_STRESS,
    "stellate_fibroblast": STELLATE_FIBROBLAST,
    "endothelial": ENDOTHELIAL,
    "pericyte_smooth_muscle": PERICYTE_SMOOTH_MUSCLE,
    "myeloid": MYELOID,
    "t_nk": T_NK_LYMPHOID,
    "b_plasma": B_PLASMA,
    "cycling": CYCLING,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Score literature-curated human pancreas signatures into an embedding."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, default=0, help="Run tag; unused.")
    return parser.parse_args()


def build_gene_index(var_names: Iterable[str]) -> dict[str, int]:
    index: dict[str, int] = {}
    for i, name in enumerate(var_names):
        key = str(name).upper()
        index.setdefault(key, i)
    return index


def as_dense_float32(matrix) -> np.ndarray:
    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32)


def module_score(matrix, gene_index: dict[str, int], genes: list[str]) -> np.ndarray:
    cols = [gene_index[g.upper()] for g in genes if g.upper() in gene_index]
    n_cells = matrix.shape[0]
    if not cols:
        return np.zeros(n_cells, dtype=np.float32)

    values = as_dense_float32(matrix[:, cols])
    means = np.asarray(values.mean(axis=0), dtype=np.float32)
    stds = np.asarray(values.std(axis=0), dtype=np.float32)
    valid = stds > 1e-6
    if not np.any(valid):
        return np.zeros(n_cells, dtype=np.float32)

    z = (values[:, valid] - means[valid]) / stds[valid]
    return np.asarray(z.mean(axis=1), dtype=np.float32)


def robust_scale(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float32)
    median = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - median))
    if not np.isfinite(mad) or mad < 1e-6:
        scale = np.nanstd(values)
    else:
        scale = 1.4826 * mad
    if not np.isfinite(scale) or scale < 1e-6:
        return np.zeros_like(values, dtype=np.float32)
    scaled = (values - median) / scale
    return np.clip(scaled, -8.0, 8.0).astype(np.float32)


def make_embedding(adata: ad.AnnData) -> pd.DataFrame:
    gene_index = build_gene_index(adata.var_names)
    matrix = adata.X

    scores = {
        name: robust_scale(module_score(matrix, gene_index, genes))
        for name, genes in MODULES.items()
    }

    exocrine = 0.5 * (scores["acinar"] + scores["ductal"])
    vascular_stromal = np.maximum.reduce(
        [
            scores["stellate_fibroblast"],
            scores["endothelial"],
            scores["pericyte_smooth_muscle"],
        ]
    )
    immune = np.maximum.reduce([scores["myeloid"], scores["t_nk"], scores["b_plasma"]])

    dims = [
        scores["endocrine"] - 0.35 * exocrine - 0.25 * vascular_stromal - 0.25 * immune,
        scores["beta"],
        scores["alpha"],
        scores["delta"],
        scores["pp"],
        scores["epsilon"],
        scores["acinar"],
        scores["ductal"],
        scores["reg3a_stress"] - 0.25 * scores["ductal"],
        scores["stellate_fibroblast"],
        scores["endothelial"],
        scores["pericyte_smooth_muscle"],
        scores["myeloid"],
        scores["t_nk"],
        scores["b_plasma"],
        scores["cycling"],
    ]

    data = {"cell_id": adata.obs_names.astype(str).to_numpy()}
    for i, dim in enumerate(dims):
        data[f"dim_{i}"] = robust_scale(dim).astype(np.float32)
    return pd.DataFrame(data)


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.input)
    embedding = make_embedding(adata)
    embedding.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
