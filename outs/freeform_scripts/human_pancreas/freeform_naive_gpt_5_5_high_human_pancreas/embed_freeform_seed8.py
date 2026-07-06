#!/usr/bin/env python3
"""Literature-marker embedding for human pancreas scRNA-seq cells.

The embedding is intentionally naive: each coordinate is a fixed score for a
published marker module or a fixed contrast between marker modules. No model is
trained on the input data, and no pretrained model inference is used.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


# Human pancreas marker modules drawn from standard scRNA-seq annotations and
# pancreatic biology literature: endocrine islet hormones/transcription factors,
# exocrine acinar and duct markers, stromal/vascular markers, immune lineage
# markers, and stress/proliferation programs commonly resolved in pancreas data.
SIGNATURES: "OrderedDict[str, list[str]]" = OrderedDict(
    [
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
                "G6PC2",
                "ADCYAP1",
                "HADH",
            ],
        ),
        (
            "alpha_cell",
            ["GCG", "ARX", "IRX2", "MAFB", "TTR", "PCSK2", "LOXL4", "GC"],
        ),
        (
            "delta_cell",
            ["SST", "HHEX", "LEPR", "RBP4", "BCHE", "PCSK1", "PCSK2"],
        ),
        (
            "gamma_pp_cell",
            ["PPY", "PYY", "SLC18A2", "PCSK1N", "FEV", "CHGA", "CHGB"],
        ),
        (
            "epsilon_ghrelin_cell",
            ["GHRL", "MLN", "RBP4", "ISL1", "CHGA", "CHGB"],
        ),
        (
            "serotonin_endocrine",
            ["TPH1", "DDC", "SLC18A1", "SLC18A2", "LMX1A", "CHGA"],
        ),
        (
            "pan_endocrine",
            [
                "CHGA",
                "CHGB",
                "SCG2",
                "SCG3",
                "PCSK1N",
                "NEUROD1",
                "ISL1",
                "PAX6",
                "RFX6",
            ],
        ),
        (
            "endocrine_progenitor",
            [
                "NEUROG3",
                "SOX4",
                "SOX11",
                "NKX2-2",
                "RFX6",
                "PAX4",
                "ONECUT1",
                "HES6",
            ],
        ),
        (
            "mature_beta",
            [
                "MAFA",
                "UCN3",
                "G6PC2",
                "SLC30A8",
                "HADH",
                "SLC2A2",
                "GLP1R",
                "ABCC8",
                "KCNJ11",
            ],
        ),
        (
            "immature_or_stressed_beta",
            ["ALDH1A3", "LDHA", "SOX9", "MYC", "HES1", "DUSP6", "HK2", "REST"],
        ),
        (
            "hormone_processing",
            ["PCSK1", "PCSK2", "CPE", "SCG2", "SCG3", "SLC18A2", "CHGA", "CHGB"],
        ),
        (
            "acinar",
            [
                "PRSS1",
                "PRSS2",
                "CPA1",
                "CPA2",
                "CPB1",
                "CLPS",
                "CTRB1",
                "CTRB2",
                "CELA3A",
                "CELA3B",
                "PNLIP",
                "CEL",
                "GP2",
                "AMY2A",
                "AMY2B",
            ],
        ),
        (
            "ductal",
            [
                "KRT19",
                "KRT7",
                "KRT8",
                "KRT18",
                "SLC4A4",
                "CFTR",
                "SOX9",
                "MMP7",
                "MUC1",
                "EPCAM",
                "TACSTD2",
                "CLDN4",
                "ANXA4",
            ],
        ),
        (
            "centroacinar",
            ["CLU", "REG3A", "ANXA4", "SPINK1", "KRT19", "SOX9", "PDPN", "MMP7"],
        ),
        (
            "acinar_ductal_metaplasia",
            ["REG3A", "REG1A", "REG1B", "KRT19", "MMP7", "TFF1", "TFF2", "ANXA2", "CLU"],
        ),
        (
            "epithelial",
            ["EPCAM", "KRT8", "KRT18", "KRT19", "TACSTD2", "CLDN4", "MUC1"],
        ),
        (
            "stellate_fibroblast",
            [
                "COL1A1",
                "COL1A2",
                "COL3A1",
                "DCN",
                "LUM",
                "COL6A1",
                "COL6A2",
                "PDGFRA",
                "PDGFRB",
            ],
        ),
        (
            "activated_stellate_myofibroblast",
            ["ACTA2", "TAGLN", "MYL9", "MYH11", "COL1A1", "COL1A2", "SPARC", "MMP2", "TIMP1"],
        ),
        (
            "endothelial",
            ["PECAM1", "VWF", "KDR", "ENG", "ESAM", "CDH5", "CLDN5", "EMCN", "RAMP2", "PLVAP"],
        ),
        (
            "pericyte_smooth_muscle",
            ["RGS5", "PDGFRB", "CSPG4", "MCAM", "ACTA2", "MYH11", "NOTCH3", "NDUFA4L2"],
        ),
        (
            "mesothelial",
            ["MSLN", "WT1", "UPK3B", "KRT19", "KRT8", "KRT18", "PDPN", "LUM"],
        ),
        (
            "macrophage_monocyte",
            [
                "LST1",
                "C1QA",
                "C1QB",
                "C1QC",
                "CD68",
                "CSF1R",
                "AIF1",
                "FCGR3A",
                "MS4A7",
                "TYROBP",
            ],
        ),
        (
            "dendritic_apc",
            ["HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "FCER1A", "CLEC10A", "ITGAX", "LILRA4"],
        ),
        (
            "t_nk",
            ["CD3D", "CD3E", "TRAC", "IL7R", "NKG7", "GNLY", "KLRD1", "GZMB", "PRF1"],
        ),
        (
            "b_plasma",
            ["MS4A1", "CD79A", "CD79B", "BANK1", "MZB1", "JCHAIN", "IGKC", "IGHG1", "XBP1"],
        ),
        (
            "mast_cell",
            ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC"],
        ),
        (
            "neutrophil_granulocyte",
            ["S100A8", "S100A9", "FCN1", "VCAN", "CSF3R", "CXCR2", "MPO", "ELANE"],
        ),
        (
            "proliferation",
            ["MKI67", "TOP2A", "CENPF", "PCNA", "HMGB2", "STMN1", "UBE2C", "TYMS"],
        ),
        (
            "interferon_inflammation",
            ["IFIT1", "IFIT2", "IFIT3", "ISG15", "MX1", "OAS1", "OAS2", "STAT1", "CXCL10"],
        ),
        (
            "heat_shock_er_stress",
            ["FOS", "JUN", "DDIT3", "HSPA1A", "HSPA1B", "ATF3", "HSP90AA1", "XBP1", "HERPUD1"],
        ),
    ]
)


CONTRASTS: "OrderedDict[str, tuple[dict[str, float], float]]" = OrderedDict(
    [
        (
            "endocrine_minus_exocrine",
            ({"pan_endocrine": 1.0, "acinar": -0.5, "ductal": -0.5}, 1.0),
        ),
        (
            "beta_minus_alpha",
            ({"beta_cell": 1.0, "mature_beta": 0.5, "alpha_cell": -1.0}, 1.5),
        ),
        (
            "alpha_delta_balance",
            ({"alpha_cell": 1.0, "delta_cell": -1.0}, 1.0),
        ),
        (
            "ductal_minus_acinar",
            ({"ductal": 1.0, "centroacinar": 0.5, "acinar": -1.0}, 1.5),
        ),
        (
            "adm_regeneration_axis",
            (
                {
                    "acinar_ductal_metaplasia": 1.0,
                    "centroacinar": 0.5,
                    "ductal": 0.25,
                    "acinar": -0.5,
                },
                1.75,
            ),
        ),
        (
            "stromal_vascular_axis",
            ({"stellate_fibroblast": 1.0, "activated_stellate_myofibroblast": 0.5, "endothelial": -1.0}, 1.5),
        ),
        (
            "immune_inflammatory_axis",
            (
                {
                    "macrophage_monocyte": 0.7,
                    "dendritic_apc": 0.5,
                    "t_nk": 0.5,
                    "b_plasma": 0.5,
                    "interferon_inflammation": 0.7,
                },
                2.4,
            ),
        ),
        (
            "beta_maturity_axis",
            ({"mature_beta": 1.0, "immature_or_stressed_beta": -1.0, "endocrine_progenitor": -0.5}, 1.5),
        ),
    ]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Embed human pancreas scRNA-seq cells with fixed literature marker count scores."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file")
    parser.add_argument("--output", required=True, help="Output parquet path")
    parser.add_argument("--seed", type=int, required=True, help="Run tag; accepted but not used")
    return parser.parse_args()


def normalize_symbol(symbol: object) -> str:
    return str(symbol).strip().upper()


def gene_indices(var_names: Iterable[object]) -> dict[str, list[int]]:
    lookup: dict[str, list[int]] = {}
    for idx, gene in enumerate(var_names):
        lookup.setdefault(normalize_symbol(gene), []).append(idx)
    return lookup


def marker_indices(markers: Iterable[str], lookup: dict[str, list[int]]) -> list[int]:
    seen: set[int] = set()
    indices: list[int] = []
    for marker in markers:
        for idx in lookup.get(normalize_symbol(marker), []):
            if idx not in seen:
                seen.add(idx)
                indices.append(idx)
    return indices


def row_sum(matrix) -> np.ndarray:
    return np.asarray(matrix.sum(axis=1)).ravel().astype(np.float64)


def score_signature(counts, total_counts: np.ndarray, indices: list[int]) -> np.ndarray:
    if not indices:
        return np.zeros(total_counts.shape[0], dtype=np.float32)

    module_counts = row_sum(counts[:, indices])
    library = np.maximum(total_counts, 1.0)
    per_marker_cp10k = 10000.0 * module_counts / (library * float(len(indices)))
    return np.log1p(per_marker_cp10k).astype(np.float32)


def build_embedding(adata: ad.AnnData) -> pd.DataFrame:
    if "counts" not in adata.layers:
        raise ValueError('Input AnnData must contain raw counts in layers["counts"].')

    counts = adata.layers["counts"]
    if sparse.issparse(counts):
        counts = counts.tocsr()
    else:
        counts = np.asarray(counts)

    total_counts = row_sum(counts)
    lookup = gene_indices(adata.var_names)

    scores: "OrderedDict[str, np.ndarray]" = OrderedDict()
    for name, markers in SIGNATURES.items():
        indices = marker_indices(markers, lookup)
        scores[name] = score_signature(counts, total_counts, indices)

    for name, (weights, scale) in CONTRASTS.items():
        contrast = np.zeros(adata.n_obs, dtype=np.float32)
        for signature_name, weight in weights.items():
            contrast += np.float32(weight) * scores[signature_name]
        scores[name] = (contrast / np.float32(scale)).astype(np.float32)

    data = {"cell_id": np.asarray(adata.obs_names, dtype=str)}
    for dim, values in enumerate(scores.values()):
        data[f"dim_{dim}"] = values.astype(np.float32, copy=False)
    return pd.DataFrame(data)


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.input)
    embedding = build_embedding(adata)
    embedding.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
