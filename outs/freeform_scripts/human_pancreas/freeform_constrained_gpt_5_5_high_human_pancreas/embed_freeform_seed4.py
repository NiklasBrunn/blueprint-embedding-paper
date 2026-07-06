#!/usr/bin/env python3
"""
Literature-derived embedding for baseline human pancreas single-cell data.

The scenario is human pancreas with no age, sex, condition, or cell-type
restriction. Adult human pancreas single-cell atlases resolve major exocrine,
endocrine, stromal, vascular, immune, and neural compartments. Fifteen
dimensions are used because the expected landscape needs separate axes for the
dominant acinar/ductal epithelium, the endocrine core and principal endocrine
subtypes, pancreatic stellate/fibroblast states, vascular mural/endothelial
cells, immune lineages, and Schwann/neural cells. The script scores fixed gene
modules only; it does not fit, train, cluster, or run model inference.
"""

from __future__ import annotations

import argparse
from typing import Iterable


SIGNATURES = [
    {
        "name": "acinar_secretory",
        # Muraro et al., A Single-Cell Transcriptome Atlas of the Human Pancreas, 10.1016/j.cels.2016.09.002;
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011.
        "genes": [
            "PRSS1",
            "PRSS2",
            "CPA1",
            "CPA2",
            "CPB1",
            "CTRB1",
            "CTRB2",
            "CELA3A",
            "CELA3B",
            "PNLIP",
            "CLPS",
            "AMY2A",
            "AMY2B",
            "CEL",
            "GP2",
            "SYCN",
            "REG1A",
            "REG1B",
            "SPINK1",
        ],
    },
    {
        "name": "ductal_centroacinar",
        # Muraro et al., A Single-Cell Transcriptome Atlas of the Human Pancreas, 10.1016/j.cels.2016.09.002;
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011.
        "genes": [
            "KRT19",
            "KRT7",
            "KRT8",
            "KRT18",
            "EPCAM",
            "SOX9",
            "CFTR",
            "SLC4A4",
            "SPP1",
            "MMP7",
            "ANXA4",
            "TFF1",
            "TFF2",
            "MUC1",
            "MUC13",
        ],
    },
    {
        "name": "endocrine_core",
        # Segerstolpe et al., Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes, 10.1016/j.cmet.2016.08.020;
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011.
        "genes": [
            "CHGA",
            "CHGB",
            "SCG5",
            "SCG3",
            "PCSK1",
            "PCSK2",
            "PCSK1N",
            "NEUROD1",
            "RFX6",
            "ISL1",
            "PAX6",
            "INSM1",
            "NKX2-2",
            "UCHL1",
            "SLC30A8",
        ],
    },
    {
        "name": "beta_cell",
        # Segerstolpe et al., Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes, 10.1016/j.cmet.2016.08.020;
        # Muraro et al., A Single-Cell Transcriptome Atlas of the Human Pancreas, 10.1016/j.cels.2016.09.002.
        "genes": [
            "INS",
            "IAPP",
            "PDX1",
            "MAFA",
            "NKX6-1",
            "PAX4",
            "G6PC2",
            "SLC2A2",
            "HADH",
            "ADCYAP1",
            "FFAR1",
            "DLK1",
            "PCSK1",
        ],
    },
    {
        "name": "alpha_cell",
        # Segerstolpe et al., Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes, 10.1016/j.cmet.2016.08.020;
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011.
        "genes": [
            "GCG",
            "ARX",
            "IRX1",
            "IRX2",
            "MAFB",
            "PCSK2",
            "TTR",
            "GC",
            "LOXL4",
            "TM4SF4",
            "CRYBA2",
        ],
    },
    {
        "name": "delta_cell",
        # Segerstolpe et al., Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes, 10.1016/j.cmet.2016.08.020;
        # Muraro et al., A Single-Cell Transcriptome Atlas of the Human Pancreas, 10.1016/j.cels.2016.09.002.
        "genes": [
            "SST",
            "HHEX",
            "RBP4",
            "LEPR",
            "GHSR",
            "BCHE",
            "RGS2",
            "FEV",
            "PCSK1",
            "PCSK2",
        ],
    },
    {
        "name": "gamma_pp_cell",
        # Muraro et al., A Single-Cell Transcriptome Atlas of the Human Pancreas, 10.1016/j.cels.2016.09.002;
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011.
        "genes": [
            "PPY",
            "PYY",
            "FEV",
            "SLC18A2",
            "MEIS2",
            "PCSK1N",
            "CHGB",
            "TSPAN8",
        ],
    },
    {
        "name": "epsilon_ghrelin_cell",
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011;
        # Segerstolpe et al., Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes, 10.1016/j.cmet.2016.08.020.
        "genes": [
            "GHRL",
            "MLN",
            "PCSK1N",
            "CHGA",
            "CHGB",
            "ISL1",
            "PAX6",
            "RBP4",
        ],
    },
    {
        "name": "quiescent_stellate_fibroblast",
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011;
        # Tabula Sapiens Consortium, The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans, 10.1126/science.abl4896.
        "genes": [
            "COL1A1",
            "COL1A2",
            "COL3A1",
            "DCN",
            "LUM",
            "MGP",
            "CFD",
            "C7",
            "CXCL12",
            "PDGFRA",
            "VIM",
            "DES",
            "PDGFRB",
        ],
    },
    {
        "name": "activated_stellate_myofibroblast",
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011;
        # Tabula Sapiens Consortium, The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans, 10.1126/science.abl4896.
        "genes": [
            "ACTA2",
            "TAGLN",
            "MYL9",
            "TPM2",
            "FN1",
            "SPARC",
            "MMP2",
            "TIMP1",
            "SERPINE1",
            "POSTN",
            "CTGF",
            "COL1A1",
            "COL3A1",
        ],
    },
    {
        "name": "endothelial",
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011;
        # Tabula Sapiens Consortium, The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans, 10.1126/science.abl4896.
        "genes": [
            "PECAM1",
            "VWF",
            "CDH5",
            "KDR",
            "FLT1",
            "ENG",
            "ESAM",
            "CLDN5",
            "EMCN",
            "PLVAP",
            "RAMP2",
            "ACKR1",
        ],
    },
    {
        "name": "mural_pericyte_smooth_muscle",
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011;
        # Tabula Sapiens Consortium, The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans, 10.1126/science.abl4896.
        "genes": [
            "RGS5",
            "MCAM",
            "CSPG4",
            "PDGFRB",
            "NOTCH3",
            "ACTA2",
            "MYH11",
            "TAGLN",
            "MYL9",
            "PLN",
            "TPM2",
            "DES",
        ],
    },
    {
        "name": "lymphoid",
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011;
        # Tabula Sapiens Consortium, The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans, 10.1126/science.abl4896.
        "genes": [
            "PTPRC",
            "CD3D",
            "CD3E",
            "CD2",
            "TRAC",
            "IL7R",
            "NKG7",
            "GNLY",
            "KLRD1",
            "MS4A1",
            "CD79A",
            "BANK1",
            "MZB1",
            "IGKC",
        ],
    },
    {
        "name": "myeloid_antigen_presenting",
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011;
        # Tabula Sapiens Consortium, The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans, 10.1126/science.abl4896.
        "genes": [
            "PTPRC",
            "LST1",
            "AIF1",
            "LYZ",
            "CD14",
            "FCGR3A",
            "TYROBP",
            "FCER1G",
            "C1QA",
            "C1QB",
            "MS4A7",
            "HLA-DRA",
            "HLA-DPA1",
            "CD74",
            "ITGAM",
        ],
    },
    {
        "name": "schwann_neural",
        # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas Reveals Inter- and Intra-cell Population Structure, 10.1016/j.cels.2016.08.011;
        # Tabula Sapiens Consortium, The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans, 10.1126/science.abl4896.
        "genes": [
            "S100B",
            "SOX10",
            "PLP1",
            "PMP22",
            "MPZ",
            "NGFR",
            "GFRA3",
            "NRXN1",
            "CDH19",
            "CRYAB",
            "FOXD3",
        ],
    },
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Score fixed literature-derived human pancreas marker modules and "
            "write a parquet embedding."
        )
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, default=0, help="Accepted run tag; unused.")
    return parser.parse_args()


def first_matching_indices(var_names: Iterable[str]) -> dict[str, int]:
    mapping: dict[str, int] = {}
    for idx, name in enumerate(var_names):
        key = str(name).upper()
        mapping.setdefault(key, idx)
    return mapping


def dense_subset(matrix, indices: list[int]):
    import numpy as np
    from scipy import sparse

    subset = matrix[:, indices]
    if sparse.issparse(subset):
        subset = subset.toarray()
    return np.asarray(subset, dtype=np.float32)


def module_score(matrix, indices: list[int], n_cells: int):
    import numpy as np

    if not indices:
        return np.zeros(n_cells, dtype=np.float32)

    values = dense_subset(matrix, indices)
    means = values.mean(axis=0, dtype=np.float64)
    stds = values.std(axis=0, dtype=np.float64)
    usable = stds > 1e-6
    if not usable.any():
        return np.zeros(n_cells, dtype=np.float32)

    z = (values[:, usable] - means[usable]) / stds[usable]
    z = np.clip(z, -5.0, 5.0)
    raw = z.mean(axis=1, dtype=np.float64)

    median = np.median(raw)
    mad = np.median(np.abs(raw - median))
    if mad > 1e-6:
        scaled = (raw - median) / (1.4826 * mad)
    else:
        sd = raw.std()
        scaled = (raw - raw.mean()) / sd if sd > 1e-6 else raw * 0.0
    return np.clip(scaled, -8.0, 8.0).astype(np.float32)


def build_embedding(input_path: str):
    import anndata as ad
    import numpy as np
    import pandas as pd

    adata = ad.read_h5ad(input_path)
    matrix = adata.X
    if "counts" in adata.layers:
        # The contract provides counts for provenance, but the intended scoring
        # scale is X: normalize_total(1e4) followed by log1p.
        matrix = adata.X

    gene_to_idx = first_matching_indices(adata.var_names)
    output = {"cell_id": np.asarray(adata.obs_names, dtype=str)}
    n_cells = adata.n_obs

    for dim, signature in enumerate(SIGNATURES):
        indices = [gene_to_idx[g.upper()] for g in signature["genes"] if g.upper() in gene_to_idx]
        output[f"dim_{dim}"] = module_score(matrix, indices, n_cells)

    return pd.DataFrame(output, columns=["cell_id"] + [f"dim_{i}" for i in range(len(SIGNATURES))])


def main() -> None:
    args = parse_args()
    embedding = build_embedding(args.input)
    embedding.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
