#!/usr/bin/env python3
"""Knowledge-driven embedding for a human pancreas cell landscape.

Scenario read from scenario_config.jsonc before authoring:
  species: human
  tissue: pancreas
  celltype: all pancreatic cell types/states
  age/sex: unspecified
  condition: baseline / healthy

The coordinate system is intentionally interpretable rather than learned.
Dimensions are marker-module scores spanning the expected pancreas landscape:
endocrine islet subtypes, exocrine/ductal epithelium, developmental endocrine
and epithelial progenitor programs, stromal/vascular/immune compartments, and
general cell states that occur in baseline tissue single-cell profiles.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


# The number of dimensions is driven by the scenario breadth: all human pancreas
# cell types with no age restriction needs mature adult compartments plus
# progenitor/immature axes, while condition=null excludes disease-specific
# cancer/diabetes/fibrosis programs.
SIGNATURES: "OrderedDict[str, list[str]]" = OrderedDict(
    [
        (
            "pan_endocrine_islet",
            [
                "CHGA",
                "CHGB",
                "SCG2",
                "SCG3",
                "PCSK1N",
                "ISL1",
                "PAX6",
                "NEUROD1",
                "RFX6",
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
                "SLC30A8",
                "PCSK1",
                "PCSK2",
                "G6PC2",
                "UCN3",
                "ABCC8",
                "KCNJ11",
            ],
        ),
        (
            "alpha_cell",
            [
                "GCG",
                "ARX",
                "IRX2",
                "TTR",
                "GC",
                "TM4SF4",
                "LOXL4",
                "MAFB",
            ],
        ),
        (
            "delta_cell",
            [
                "SST",
                "HHEX",
                "RBP4",
                "LEPR",
                "BCHE",
                "GHSR",
                "PCSK1",
            ],
        ),
        (
            "pp_gamma_cell",
            [
                "PPY",
                "PYY",
                "SLC6A4",
                "ID4",
                "FAP",
                "CHGA",
            ],
        ),
        (
            "epsilon_cell",
            [
                "GHRL",
                "MLN",
                "CHGA",
                "CHGB",
                "ISL1",
            ],
        ),
        (
            "immature_endocrine_progenitor",
            [
                "NEUROG3",
                "NEUROD1",
                "NKX2-2",
                "PAX4",
                "PAX6",
                "FEV",
                "SOX4",
                "SOX11",
                "ISL1",
                "RFX6",
            ],
        ),
        (
            "acinar_exocrine",
            [
                "PRSS1",
                "PRSS2",
                "CPA1",
                "CPA2",
                "CTRB1",
                "CTRB2",
                "CELA3A",
                "CELA3B",
                "PNLIP",
                "CLPS",
                "AMY2A",
                "REG1A",
                "REG1B",
                "SPINK1",
            ],
        ),
        (
            "ductal_epithelial",
            [
                "KRT19",
                "KRT7",
                "KRT8",
                "KRT18",
                "EPCAM",
                "SOX9",
                "CFTR",
                "SLC4A4",
                "MUC1",
                "SPP1",
                "ANXA4",
                "CLDN10",
            ],
        ),
        (
            "centroacinar_ion_transport",
            [
                "CFTR",
                "SLC4A4",
                "ANXA4",
                "SOX9",
                "HNF1B",
                "ONECUT1",
                "ALDH1A1",
                "KRT19",
                "AQP1",
                "TFF1",
            ],
        ),
        (
            "pancreatic_epithelial_progenitor",
            [
                "PDX1",
                "SOX9",
                "HNF1B",
                "ONECUT1",
                "ONECUT2",
                "HES1",
                "PTF1A",
                "MNX1",
                "NKX6-1",
                "FOXA2",
            ],
        ),
        (
            "fibroblast_stellate_matrix",
            [
                "COL1A1",
                "COL1A2",
                "COL3A1",
                "COL6A1",
                "COL6A2",
                "DCN",
                "LUM",
                "THY1",
                "PDGFRA",
                "MMP2",
            ],
        ),
        (
            "activated_stellate_myofibroblast",
            [
                "ACTA2",
                "TAGLN",
                "MYL9",
                "TPM2",
                "POSTN",
                "CTHRC1",
                "COL1A1",
                "COL1A2",
                "TIMP1",
                "SERPINE1",
            ],
        ),
        (
            "blood_endothelial",
            [
                "PECAM1",
                "VWF",
                "KDR",
                "FLT1",
                "CDH5",
                "ESAM",
                "RAMP2",
                "PLVAP",
                "EMCN",
                "CLDN5",
            ],
        ),
        (
            "lymphatic_endothelial",
            [
                "PROX1",
                "PDPN",
                "LYVE1",
                "FLT4",
                "CCL21",
                "MMRN1",
                "RELN",
            ],
        ),
        (
            "mural_pericyte_smooth_muscle",
            [
                "RGS5",
                "PDGFRB",
                "MCAM",
                "CSPG4",
                "NOTCH3",
                "ACTA2",
                "MYH11",
                "TAGLN",
                "MYLK",
                "CNN1",
            ],
        ),
        (
            "myeloid_macrophage",
            [
                "LST1",
                "C1QA",
                "C1QB",
                "C1QC",
                "AIF1",
                "CD68",
                "CSF1R",
                "TYROBP",
                "FCGR3A",
                "MS4A7",
                "LYZ",
            ],
        ),
        (
            "dendritic_antigen_presentation",
            [
                "CLEC10A",
                "FCER1A",
                "CD1C",
                "ITGAX",
                "HLA-DRA",
                "HLA-DPA1",
                "HLA-DPB1",
                "LAMP3",
                "CCR7",
            ],
        ),
        (
            "t_cell_nk",
            [
                "CD3D",
                "CD3E",
                "TRAC",
                "IL7R",
                "NKG7",
                "GNLY",
                "KLRD1",
                "PRF1",
                "GZMB",
                "CCL5",
            ],
        ),
        (
            "b_cell_plasma",
            [
                "MS4A1",
                "CD79A",
                "CD79B",
                "BANK1",
                "MZB1",
                "JCHAIN",
                "IGKC",
                "IGHG1",
                "XBP1",
            ],
        ),
        (
            "mast_cell",
            [
                "TPSAB1",
                "TPSB2",
                "CPA3",
                "KIT",
                "MS4A2",
                "HDC",
                "GATA2",
            ],
        ),
        (
            "blood_erythroid_platelet",
            [
                "HBB",
                "HBA1",
                "HBA2",
                "ALAS2",
                "AHSP",
                "PF4",
                "PPBP",
                "GP9",
                "ITGA2B",
            ],
        ),
        (
            "neural_glial",
            [
                "S100B",
                "SOX10",
                "PLP1",
                "MPZ",
                "GFRA3",
                "NRXN1",
                "NCAM1",
                "ELAVL4",
                "PHOX2B",
            ],
        ),
        (
            "cell_cycle",
            [
                "MKI67",
                "TOP2A",
                "UBE2C",
                "PCNA",
                "MCM2",
                "MCM5",
                "STMN1",
                "HMGB2",
                "CENPF",
            ],
        ),
        (
            "interferon_antiviral",
            [
                "ISG15",
                "IFIT1",
                "IFIT2",
                "IFIT3",
                "MX1",
                "OAS1",
                "OAS2",
                "IFI6",
                "STAT1",
                "IRF7",
            ],
        ),
        (
            "stress_upr_immediate_early",
            [
                "HSPA1A",
                "HSPA1B",
                "DNAJB1",
                "FOS",
                "JUN",
                "DDIT3",
                "ATF4",
                "XBP1",
                "HSP90B1",
                "HERPUD1",
            ],
        ),
        (
            "hypoxia_glycolytic",
            [
                "VEGFA",
                "CA9",
                "BNIP3",
                "LDHA",
                "PDK1",
                "SLC2A1",
                "ENO1",
                "PGK1",
                "NDRG1",
            ],
        ),
    ]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a knowledge-based human pancreas embedding from an h5ad file."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, default=0, help="Accepted run tag; not used.")
    return parser.parse_args()


def make_gene_index(var_names: Iterable[str]) -> dict[str, int]:
    gene_to_index: dict[str, int] = {}
    for idx, gene in enumerate(var_names):
        key = str(gene).strip().upper()
        if key and key not in gene_to_index:
            gene_to_index[key] = idx
    return gene_to_index


def present_indices(genes: list[str], gene_to_index: dict[str, int]) -> list[int]:
    return [gene_to_index[g.upper()] for g in genes if g.upper() in gene_to_index]


def module_mean(matrix, indices: list[int]) -> np.ndarray:
    n_cells = matrix.shape[0]
    if not indices:
        return np.zeros(n_cells, dtype=np.float64)

    subset = matrix[:, indices]
    if sparse.issparse(subset):
        return np.asarray(subset.mean(axis=1)).ravel().astype(np.float64, copy=False)
    return np.asarray(subset, dtype=np.float64).mean(axis=1)


def robust_scale(values: np.ndarray) -> np.ndarray:
    x = np.asarray(values, dtype=np.float64)
    finite = np.isfinite(x)
    if not finite.any():
        return np.zeros_like(x, dtype=np.float32)

    clean = x[finite]
    center = float(np.median(clean))
    q05, q95 = np.percentile(clean, [5.0, 95.0])
    scale = float(q95 - q05)

    if scale < 1e-8:
        positive = clean[clean > 0]
        if positive.size >= 2:
            q10, q90 = np.percentile(positive, [10.0, 90.0])
            scale = float(q90 - q10)

    if scale < 1e-8:
        scale = float(np.std(clean))

    if scale < 1e-8:
        scale = 1.0

    z = (x - center) / scale
    z[~np.isfinite(z)] = 0.0
    return np.clip(z, -6.0, 6.0).astype(np.float32)


def build_embedding(adata: ad.AnnData) -> pd.DataFrame:
    matrix = adata.X
    gene_to_index = make_gene_index(adata.var_names)
    columns: dict[str, np.ndarray] = {
        "cell_id": np.asarray(adata.obs_names, dtype=str),
    }

    for dim, (_name, genes) in enumerate(SIGNATURES.items()):
        raw_score = module_mean(matrix, present_indices(genes, gene_to_index))
        columns[f"dim_{dim}"] = robust_scale(raw_score)

    return pd.DataFrame(columns)


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.input)
    embedding = build_embedding(adata)

    expected_columns = ["cell_id"] + [f"dim_{i}" for i in range(len(SIGNATURES))]
    embedding = embedding.loc[:, expected_columns]
    for col in expected_columns[1:]:
        embedding[col] = embedding[col].astype(np.float32, copy=False)

    output_path = Path(args.output)
    embedding.to_parquet(output_path, index=False)


if __name__ == "__main__":
    main()
