#!/usr/bin/env python3
"""
Literature-derived embedding for baseline human pancreas single-cell data.

The script scores hand-curated marker modules from the human pancreas and
human cell-atlas literature. It does not train a model, cluster cells, or run
pretrained inference. Each output dimension is a robustly standardized module
score over the log1p-normalized expression matrix supplied in the input h5ad.
"""

from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from typing import Iterable

import numpy as np
import pandas as pd


# Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
ACINAR_DIGESTIVE_GENES = [
    "PRSS1",
    "PRSS2",
    "CPA1",
    "CPA2",
    "CPB1",
    "CELA3A",
    "CELA3B",
    "CTRB1",
    "CTRB2",
    "PNLIP",
    "CLPS",
    "AMY2A",
    "AMY2B",
    "REG1A",
]

# Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
ACINAR_REG_STRESS_GENES = [
    "REG1A",
    "REG1B",
    "REG3A",
    "REG3G",
    "SPINK1",
    "PAPPA2",
    "CLU",
    "GP2",
]

# Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
DUCTAL_CENTROACINAR_GENES = [
    "KRT19",
    "KRT7",
    "KRT8",
    "KRT18",
    "SOX9",
    "SLC4A4",
    "SLC26A6",
    "CFTR",
    "MUC1",
    "ANXA4",
    "TSPAN8",
    "SPP1",
]

# Segerstolpe et al., Cell Metabolism 2016, 10.1016/j.cmet.2016.08.020; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
BETA_CELL_GENES = [
    "INS",
    "IAPP",
    "MAFA",
    "PDX1",
    "NKX6-1",
    "PCSK1",
    "SLC30A8",
    "G6PC2",
    "ADCYAP1",
    "DLK1",
]

# Segerstolpe et al., Cell Metabolism 2016, 10.1016/j.cmet.2016.08.020; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
ALPHA_CELL_GENES = [
    "GCG",
    "ARX",
    "IRX2",
    "MAFB",
    "TTR",
    "TM4SF4",
    "GC",
    "LOXL4",
    "PCSK2",
]

# Segerstolpe et al., Cell Metabolism 2016, 10.1016/j.cmet.2016.08.020; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
DELTA_CELL_GENES = [
    "SST",
    "HHEX",
    "RBP4",
    "LEPR",
    "BCHE",
    "GHSR",
    "PCSK1N",
]

# Segerstolpe et al., Cell Metabolism 2016, 10.1016/j.cmet.2016.08.020; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
PP_GAMMA_CELL_GENES = [
    "PPY",
    "PYY",
    "SERTM1",
    "ID4",
    "MEIS2",
    "PCSK1",
    "CHGA",
]

# Segerstolpe et al., Cell Metabolism 2016, 10.1016/j.cmet.2016.08.020; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
EPSILON_CELL_GENES = [
    "GHRL",
    "MLN",
    "ISL1",
    "CHGA",
    "CHGB",
    "PCSK1",
]

# Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
ENDOCRINE_SECRETORY_CORE_GENES = [
    "CHGA",
    "CHGB",
    "SCG2",
    "SCG3",
    "SCG5",
    "PCSK1",
    "PCSK2",
    "PAX6",
    "ISL1",
    "NEUROD1",
]

# Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
ENDOTHELIAL_GENES = [
    "PECAM1",
    "VWF",
    "KDR",
    "FLT1",
    "EMCN",
    "ESAM",
    "CDH5",
    "RAMP2",
    "CLEC14A",
]

# Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
PERICYTE_SMOOTH_MUSCLE_GENES = [
    "RGS5",
    "PDGFRB",
    "MCAM",
    "CSPG4",
    "ACTA2",
    "TAGLN",
    "MYH11",
    "NOTCH3",
]

# Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
STELLATE_QUIESCENT_GENES = [
    "RGS5",
    "PDGFRB",
    "DCN",
    "LUM",
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "VIM",
    "MGP",
    "APOD",
]

# Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
STELLATE_ACTIVATED_GENES = [
    "ACTA2",
    "TAGLN",
    "COL1A1",
    "COL1A2",
    "FN1",
    "SPARC",
    "CTGF",
    "THY1",
    "POSTN",
]

# Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
MYELOID_MAST_GENES = [
    "PTPRC",
    "LST1",
    "TYROBP",
    "AIF1",
    "FCER1G",
    "C1QA",
    "C1QB",
    "CD68",
    "TPSAB1",
    "CPA3",
    "MS4A2",
]

# Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
LYMPHOID_GENES = [
    "PTPRC",
    "CD3D",
    "CD3E",
    "TRAC",
    "IL7R",
    "NKG7",
    "GNLY",
    "GZMB",
    "MS4A1",
    "CD79A",
    "MZB1",
]

# Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896; Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
SCHWANN_NEURAL_GENES = [
    "S100B",
    "SOX10",
    "PLP1",
    "MPZ",
    "MBP",
    "NRXN1",
    "NCAM1",
    "GFRA3",
]


@dataclass(frozen=True)
class Signature:
    name: str
    genes: tuple[str, ...]


SIGNATURES = [
    Signature("acinar_digestive", tuple(ACINAR_DIGESTIVE_GENES)),
    Signature("acinar_reg_stress", tuple(ACINAR_REG_STRESS_GENES)),
    Signature("ductal_centroacinar", tuple(DUCTAL_CENTROACINAR_GENES)),
    Signature("beta_cell", tuple(BETA_CELL_GENES)),
    Signature("alpha_cell", tuple(ALPHA_CELL_GENES)),
    Signature("delta_cell", tuple(DELTA_CELL_GENES)),
    Signature("pp_gamma_cell", tuple(PP_GAMMA_CELL_GENES)),
    Signature("epsilon_cell", tuple(EPSILON_CELL_GENES)),
    Signature("endocrine_secretory_core", tuple(ENDOCRINE_SECRETORY_CORE_GENES)),
    Signature("endothelial", tuple(ENDOTHELIAL_GENES)),
    Signature("pericyte_smooth_muscle", tuple(PERICYTE_SMOOTH_MUSCLE_GENES)),
    Signature("stellate_quiescent", tuple(STELLATE_QUIESCENT_GENES)),
    Signature("stellate_activated", tuple(STELLATE_ACTIVATED_GENES)),
    Signature("myeloid_mast", tuple(MYELOID_MAST_GENES)),
    Signature("lymphoid", tuple(LYMPHOID_GENES)),
    Signature("schwann_neural", tuple(SCHWANN_NEURAL_GENES)),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Score literature-derived human-pancreas modules into a parquet embedding."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, required=True, help="Run tag; accepted for contract compatibility.")
    return parser.parse_args()


def read_h5ad(path: str):
    try:
        import anndata as ad
    except ImportError as exc:
        raise SystemExit("anndata is required to read .h5ad input files") from exc
    return ad.read_h5ad(path)


def build_gene_index(var_names: Iterable[object]) -> dict[str, int]:
    gene_to_idx: dict[str, int] = {}
    for idx, gene in enumerate(var_names):
        key = str(gene).upper()
        if key not in gene_to_idx:
            gene_to_idx[key] = idx
    return gene_to_idx


def module_score(matrix, gene_indices: list[int], n_cells: int) -> np.ndarray:
    if not gene_indices:
        return np.zeros(n_cells, dtype=np.float32)

    try:
        from scipy import sparse
    except ImportError:
        sparse = None

    subset = matrix[:, gene_indices]
    if sparse is not None and sparse.issparse(subset):
        values = np.asarray(subset.mean(axis=1)).ravel()
    else:
        values = np.asarray(subset)
        if values.ndim == 1:
            values = values.reshape(-1, 1)
        values = values.mean(axis=1)
    return np.asarray(values, dtype=np.float32)


def robust_standardize(values: np.ndarray) -> np.ndarray:
    finite = np.isfinite(values)
    if not np.any(finite):
        return np.zeros(values.shape[0], dtype=np.float32)

    clean = values[finite].astype(np.float64, copy=False)
    center = float(np.median(clean))
    mad = float(np.median(np.abs(clean - center)))
    scale = 1.4826 * mad
    if scale < 1e-6:
        scale = float(np.std(clean))
    if scale < 1e-6:
        scale = 1.0

    standardized = (values.astype(np.float64, copy=False) - center) / scale
    standardized[~finite] = 0.0
    return np.clip(standardized, -8.0, 8.0).astype(np.float32)


def make_embedding(adata) -> pd.DataFrame:
    n_cells = int(adata.n_obs)
    gene_index = build_gene_index(adata.var_names)
    matrix = adata.X

    output = {"cell_id": np.asarray(adata.obs_names, dtype=str)}
    for dim, signature in enumerate(SIGNATURES):
        indices = [gene_index[g.upper()] for g in signature.genes if g.upper() in gene_index]
        raw = module_score(matrix, indices, n_cells)
        output[f"dim_{dim}"] = robust_standardize(raw)

    columns = ["cell_id"] + [f"dim_{dim}" for dim in range(len(SIGNATURES))]
    frame = pd.DataFrame(output, columns=columns)
    for column in columns[1:]:
        frame[column] = frame[column].astype(np.float32, copy=False)
    return frame


def main() -> None:
    args = parse_args()
    adata = read_h5ad(args.input)
    embedding = make_embedding(adata)

    parent = os.path.dirname(os.path.abspath(args.output))
    if parent:
        os.makedirs(parent, exist_ok=True)
    embedding.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
