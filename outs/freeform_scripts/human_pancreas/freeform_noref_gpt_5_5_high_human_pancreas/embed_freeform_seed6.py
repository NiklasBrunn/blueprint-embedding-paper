#!/usr/bin/env python3
"""Knowledge-defined embedding for a human pancreas cell landscape.

Scenario interpreted from scenario_config.jsonc:
human pancreas, healthy/baseline, no age/sex restriction, all expected cell
types.  The axes below therefore cover endocrine, acinar, ductal/progenitor,
vascular, stromal, immune, neural, cycling, and stress/secretory-state programs.
They are intentionally literal gene signatures, not learned components.
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


SIGNATURES: "OrderedDict[str, list[str]]" = OrderedDict(
    [
        (
            "endocrine_core",
            [
                "CHGA",
                "CHGB",
                "SCG2",
                "SCG3",
                "SYP",
                "PCSK1",
                "PCSK2",
                "NEUROD1",
                "INSM1",
                "ISL1",
                "PAX6",
                "RFX6",
                "PAX4",
            ],
        ),
        (
            "acinar_secretory",
            [
                "PRSS1",
                "PRSS2",
                "PRSS3",
                "CTRB1",
                "CTRB2",
                "CPA1",
                "CPA2",
                "CELA2A",
                "CELA3A",
                "CELA3B",
                "CLPS",
                "PNLIP",
                "AMY2A",
                "AMY2B",
                "REG1A",
                "REG1B",
                "SYCN",
                "GP2",
            ],
        ),
        (
            "ductal_centroacinar",
            [
                "KRT19",
                "KRT7",
                "KRT8",
                "KRT18",
                "SOX9",
                "CFTR",
                "MUC1",
                "SLC4A4",
                "ANXA4",
                "EPCAM",
                "CLDN4",
                "CLDN10",
                "TFF1",
                "SPP1",
                "KLF4",
            ],
        ),
        (
            "beta_cell",
            [
                "INS",
                "IAPP",
                "SLC30A8",
                "MAFA",
                "UCN3",
                "G6PC2",
                "PDX1",
                "NKX6-1",
                "HADH",
                "PCSK1",
                "ADCYAP1",
            ],
        ),
        (
            "alpha_cell",
            [
                "GCG",
                "TTR",
                "ARX",
                "IRX2",
                "PCSK2",
                "TM4SF4",
                "GC",
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
                "CCKBR",
                "PCSK1",
                "PCSK2",
            ],
        ),
        (
            "pp_gamma_cell",
            [
                "PPY",
                "PYY",
                "FEV",
                "SLC18A2",
                "PCSK1",
                "CHGA",
            ],
        ),
        (
            "epsilon_ghrelin_cell",
            [
                "GHRL",
                "MLN",
                "GPX3",
                "GPR39",
                "CHGA",
                "SYP",
            ],
        ),
        (
            "enterochromaffin_serotonergic",
            [
                "TPH1",
                "DDC",
                "SLC18A1",
                "SLC6A4",
                "TAC1",
                "NTS",
                "LMX1A",
                "CHGA",
            ],
        ),
        (
            "developmental_pancreatic_progenitor",
            [
                "NEUROG3",
                "SOX9",
                "HNF1B",
                "ONECUT1",
                "ONECUT2",
                "HES1",
                "DLL1",
                "NOTCH2",
                "PROM1",
                "ALDH1A1",
                "PDX1",
                "PTF1A",
                "CPA1",
                "RFX6",
            ],
        ),
        (
            "blood_endothelial",
            [
                "PECAM1",
                "VWF",
                "KDR",
                "FLT1",
                "ESAM",
                "CDH5",
                "EMCN",
                "RAMP2",
                "CLDN5",
                "PLVAP",
                "ACKR1",
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
                "TBX1",
            ],
        ),
        (
            "fibroblast_stellate_matrix",
            [
                "COL1A1",
                "COL1A2",
                "COL3A1",
                "DCN",
                "LUM",
                "COL6A1",
                "COL6A2",
                "PDGFRA",
                "FBLN1",
                "MMP2",
                "SPARC",
                "CXCL12",
            ],
        ),
        (
            "mural_smooth_muscle_pericyte",
            [
                "ACTA2",
                "TAGLN",
                "MYH11",
                "MCAM",
                "RGS5",
                "PDGFRB",
                "CSPG4",
                "NOTCH3",
                "DES",
                "TPM2",
            ],
        ),
        (
            "myeloid_macrophage",
            [
                "PTPRC",
                "LST1",
                "AIF1",
                "TYROBP",
                "FCER1G",
                "CD68",
                "CSF1R",
                "C1QA",
                "C1QB",
                "C1QC",
                "CD14",
                "FCGR3A",
                "MS4A7",
            ],
        ),
        (
            "t_nk_lymphoid",
            [
                "PTPRC",
                "CD3D",
                "CD3E",
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
                "PTPRC",
                "MS4A1",
                "CD79A",
                "CD79B",
                "CD74",
                "MZB1",
                "JCHAIN",
                "IGKC",
                "IGHG1",
                "IGHA1",
                "XBP1",
            ],
        ),
        (
            "mast_cell",
            [
                "PTPRC",
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
            "schwann_neural",
            [
                "SOX10",
                "S100B",
                "MPZ",
                "PLP1",
                "GFRA3",
                "NRXN1",
                "NGFR",
                "GAP43",
                "FOXD3",
            ],
        ),
        (
            "cell_cycle",
            [
                "MKI67",
                "TOP2A",
                "PCNA",
                "MCM2",
                "MCM5",
                "TYMS",
                "CENPF",
                "UBE2C",
                "CCNB1",
                "CDC20",
            ],
        ),
        (
            "stress_interferon_secretory_load",
            [
                "HSPA1A",
                "HSPA1B",
                "DNAJB1",
                "JUN",
                "FOS",
                "DDIT3",
                "ATF3",
                "HSPA5",
                "XBP1",
                "PDIA4",
                "IFIT1",
                "IFIT3",
                "ISG15",
                "MX1",
                "OAS1",
                "CXCL8",
                "NFKBIA",
            ],
        ),
    ]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a biological-knowledge embedding for human pancreas cells "
            "from literal marker-gene signatures."
        )
    )
    parser.add_argument("--input", required=True, help="Input .h5ad file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, required=True, help="Run tag; not used.")
    return parser.parse_args()


def build_gene_index(var_names: Iterable[str]) -> dict[str, int]:
    """Map gene symbol variants to the first observed column index."""
    gene_to_index: dict[str, int] = {}
    for idx, gene in enumerate(var_names):
        text = str(gene)
        if text not in gene_to_index:
            gene_to_index[text] = idx
        upper = text.upper()
        if upper not in gene_to_index:
            gene_to_index[upper] = idx
    return gene_to_index


def present_indices(genes: Iterable[str], gene_to_index: dict[str, int]) -> list[int]:
    seen: set[int] = set()
    indices: list[int] = []
    for gene in genes:
        idx = gene_to_index.get(gene)
        if idx is None:
            idx = gene_to_index.get(gene.upper())
        if idx is not None and idx not in seen:
            indices.append(idx)
            seen.add(idx)
    return indices


def row_mean(matrix, indices: list[int], n_obs: int) -> np.ndarray:
    if not indices:
        return np.zeros(n_obs, dtype=np.float64)
    values = matrix[:, indices]
    if sparse.issparse(values):
        values = values.mean(axis=1)
        return np.asarray(values).ravel().astype(np.float64, copy=False)
    return np.asarray(values, dtype=np.float64).mean(axis=1)


def robust_standardize(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float64)
    if values.size == 0:
        return values

    median = np.nanmedian(values)
    q25, q75 = np.nanpercentile(values, [25, 75])
    scale = (q75 - q25) / 1.349

    if not np.isfinite(scale) or scale <= 1e-8:
        scale = np.nanstd(values)
    if not np.isfinite(scale) or scale <= 1e-8:
        return np.zeros_like(values, dtype=np.float32)

    standardized = (values - median) / scale
    standardized = np.nan_to_num(standardized, nan=0.0, posinf=8.0, neginf=-8.0)
    return np.clip(standardized, -8.0, 8.0).astype(np.float32)


def compute_embedding(adata: ad.AnnData) -> np.ndarray:
    matrix = adata.X
    gene_to_index = build_gene_index(adata.var_names)
    columns: list[np.ndarray] = []

    for genes in SIGNATURES.values():
        indices = present_indices(genes, gene_to_index)
        raw_score = row_mean(matrix, indices, adata.n_obs)
        columns.append(robust_standardize(raw_score))

    return np.column_stack(columns).astype(np.float32, copy=False)


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.input)
    embedding = compute_embedding(adata)

    result = pd.DataFrame(
        {"cell_id": np.asarray(adata.obs_names, dtype=str)},
    )
    for dim_idx in range(embedding.shape[1]):
        result[f"dim_{dim_idx}"] = embedding[:, dim_idx].astype(np.float32, copy=False)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output, index=False)


if __name__ == "__main__":
    main()
