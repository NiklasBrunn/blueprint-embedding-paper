#!/usr/bin/env python3
"""
Literature-marker embedding for human pancreas scRNA-seq data.

This script does not train a model, infer clusters, select genes from the
dataset, or run pretrained inference.  It maps each cell to fixed coordinates
defined by pancreas marker modules from established scRNA-seq / histology
literature: endocrine islet hormones and transcription factors, exocrine
acinar/ductal genes, vascular, stromal, neural, immune, stress, and cell-cycle
programs.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import scipy.sparse as sp


PANCREAS_MODULES: "OrderedDict[str, list[str]]" = OrderedDict(
    [
        # Endocrine islet lineages and maturation programs.
        ("beta_insulin", ["INS", "IAPP", "PCSK1", "PCSK2", "SLC30A8", "G6PC2", "ADCYAP1"]),
        ("beta_maturity", ["MAFA", "PDX1", "NKX6-1", "NKX2-2", "PAX6", "ABCC8", "KCNJ11", "RBP4"]),
        ("alpha_glucagon", ["GCG", "ARX", "IRX2", "IRX1", "PCSK2", "TTR", "LOX"]),
        ("delta_somatostatin", ["SST", "HHEX", "RBP4", "LEPR", "BCHE", "RGS2"]),
        ("pp_gamma", ["PPY", "PAX6", "SLC30A8", "PCSK1", "REG1A"]),
        ("epsilon_ghrelin", ["GHRL", "RBP4", "ISL1", "PAX6", "PCSK1"]),
        ("endocrine_pan_islet", ["CHGA", "CHGB", "SCG2", "SCG3", "NEUROD1", "ISL1", "PAX6", "PCSK1N"]),
        ("islet_developmental", ["NEUROG3", "INSM1", "FEV", "MYT1", "PAX4", "RFX6", "NKX2-2"]),
        ("islet_secretory_granule", ["CHGA", "CHGB", "SCG2", "SCG5", "PCSK1N", "VGF", "SYT7", "RAB3A"]),
        ("islet_metabolic_sensing", ["GCK", "SLC2A2", "SLC30A8", "ABCC8", "KCNJ11", "PCSK1", "PCSK2"]),
        # Exocrine pancreas.
        ("acinar_digestive_enzyme", ["PRSS1", "PRSS2", "CPA1", "CPA2", "CTRB1", "CTRB2", "CELA3A", "CELA3B", "PNLIP", "CLPS"]),
        ("acinar_secretory_rer", ["REG1A", "REG1B", "REG3A", "SPINK1", "SYCN", "GP2", "AMY2A", "AMY2B"]),
        ("ductal_bicarbonate_epithelium", ["KRT19", "KRT7", "KRT8", "KRT18", "SOX9", "SLC4A4", "SLC26A6", "CFTR", "MUC1"]),
        ("centroacinar_ductal_progenitor", ["KRT19", "SOX9", "SPP1", "ANXA4", "TSPAN8", "MMP7", "LCN2"]),
        ("epithelial_general", ["EPCAM", "KRT8", "KRT18", "KRT19", "TACSTD2", "CDH1", "CLDN4"]),
        # Mesenchymal, vascular, neural, and immune compartments present in pancreas.
        ("endothelial_blood", ["PECAM1", "VWF", "KDR", "FLT1", "EMCN", "CDH5", "ESAM", "RAMP2"]),
        ("endothelial_arterial", ["EFNB2", "HEY1", "GJA5", "SOX17", "DLL4", "SEMA3G"]),
        ("endothelial_venous", ["ACKR1", "NR2F2", "VCAM1", "SELE", "PLVAP"]),
        ("lymphatic_endothelial", ["PROX1", "PDPN", "LYVE1", "FLT4", "CCL21", "MMRN1"]),
        ("pericyte_smooth_muscle", ["RGS5", "PDGFRB", "CSPG4", "MCAM", "ACTA2", "TAGLN", "MYH11", "NOTCH3"]),
        ("fibroblast_stellate", ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL6A1", "MMP2"]),
        ("activated_stellate_myofibroblast", ["ACTA2", "TAGLN", "COL1A1", "COL1A2", "POSTN", "THY1", "SPARC", "TIMP1"]),
        ("schwann_neural", ["SOX10", "S100B", "PLP1", "MPZ", "NRXN1", "GFRA3"]),
        ("immune_pan_leukocyte", ["PTPRC", "LST1", "TYROBP", "CXCR4", "CORO1A", "HCST"]),
        ("myeloid_macrophage", ["LYZ", "LST1", "AIF1", "C1QA", "C1QB", "C1QC", "CD68", "FCGR3A"]),
        ("dendritic_antigen_presenting", ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "CD74", "FCER1A", "CLEC10A"]),
        ("t_nk_cells", ["CD3D", "CD3E", "TRAC", "IL7R", "NKG7", "GNLY", "KLRD1", "PRF1"]),
        ("b_plasma_cells", ["MS4A1", "CD79A", "CD79B", "MZB1", "JCHAIN", "XBP1", "IGHG1", "IGKC"]),
        ("mast_cells", ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC"]),
        # Cross-cutting technical/biological states useful for preserving neighborhoods.
        ("erythroid_contamination", ["HBB", "HBA1", "HBA2", "ALAS2", "AHSP", "SLC4A1"]),
        ("cell_cycle_s", ["MKI67", "PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "TYMS"]),
        ("cell_cycle_g2m", ["TOP2A", "UBE2C", "BIRC5", "CENPF", "NUSAP1", "CDK1", "CCNB1", "CCNB2"]),
        ("immediate_stress", ["FOS", "JUN", "JUNB", "JUND", "EGR1", "IER2", "ATF3", "DDIT3"]),
        ("interferon_response", ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "STAT1", "IRF7"]),
        ("hypoxia_glycolysis", ["VEGFA", "SLC2A1", "LDHA", "ENO1", "PGK1", "CA9", "BNIP3", "NDRG1"]),
        ("ribosomal_translation", ["RPLP0", "RPL3", "RPL7", "RPL13A", "RPS3", "RPS6", "RPS18", "EEF1A1"]),
        ("mitochondrial_content", ["MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ATP6"]),
    ]
)


CONTRASTS: list[tuple[str, str]] = [
    ("beta_insulin", "alpha_glucagon"),
    ("beta_insulin", "delta_somatostatin"),
    ("alpha_glucagon", "delta_somatostatin"),
    ("pp_gamma", "epsilon_ghrelin"),
    ("endocrine_pan_islet", "acinar_digestive_enzyme"),
    ("endocrine_pan_islet", "ductal_bicarbonate_epithelium"),
    ("acinar_digestive_enzyme", "ductal_bicarbonate_epithelium"),
    ("ductal_bicarbonate_epithelium", "centroacinar_ductal_progenitor"),
    ("fibroblast_stellate", "endothelial_blood"),
    ("pericyte_smooth_muscle", "endothelial_blood"),
    ("immune_pan_leukocyte", "epithelial_general"),
    ("myeloid_macrophage", "t_nk_cells"),
    ("b_plasma_cells", "t_nk_cells"),
    ("activated_stellate_myofibroblast", "fibroblast_stellate"),
    ("cell_cycle_g2m", "cell_cycle_s"),
    ("immediate_stress", "endocrine_pan_islet"),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Embed cells with fixed human-pancreas literature marker modules.")
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", required=True, type=int, help="Run tag; accepted for interface compatibility.")
    return parser.parse_args()


def read_h5ad(path: str):
    try:
        import anndata as ad
    except ImportError as exc:  # pragma: no cover - environment dependent
        raise SystemExit("embed_freeform.py requires the 'anndata' package to read .h5ad input.") from exc
    return ad.read_h5ad(path)


def make_gene_index(var_names: Iterable[object]) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for i, name in enumerate(var_names):
        key = str(name).upper()
        index.setdefault(key, []).append(i)
    return index


def marker_indices(markers: Iterable[str], gene_index: dict[str, list[int]]) -> list[int]:
    idx: list[int] = []
    seen: set[int] = set()
    for marker in markers:
        for gene_idx in gene_index.get(marker.upper(), []):
            if gene_idx not in seen:
                idx.append(gene_idx)
                seen.add(gene_idx)
    return idx


def as_csr_counts(matrix):
    if sp.issparse(matrix):
        return matrix.tocsr().astype(np.float32)
    return sp.csr_matrix(np.asarray(matrix, dtype=np.float32))


def log1p_cp10k(counts_csr: sp.csr_matrix) -> sp.csr_matrix:
    libsize = np.asarray(counts_csr.sum(axis=1)).ravel().astype(np.float32)
    scale = np.zeros_like(libsize, dtype=np.float32)
    nonzero = libsize > 0
    scale[nonzero] = 10000.0 / libsize[nonzero]

    norm = counts_csr.copy().astype(np.float32)
    norm = norm.multiply(scale[:, None]).tocsr()
    norm.data = np.log1p(norm.data, dtype=np.float32)
    return norm


def module_mean(log_norm: sp.csr_matrix, idx: list[int]) -> np.ndarray:
    if not idx:
        return np.zeros(log_norm.shape[0], dtype=np.float32)
    values = log_norm[:, idx].mean(axis=1)
    return np.asarray(values).ravel().astype(np.float32)


def module_detected_fraction(counts_csr: sp.csr_matrix, idx: list[int]) -> np.ndarray:
    if not idx:
        return np.zeros(counts_csr.shape[0], dtype=np.float32)
    detected = counts_csr[:, idx].copy()
    detected.data = np.ones_like(detected.data, dtype=np.float32)
    frac = detected.mean(axis=1)
    return np.asarray(frac).ravel().astype(np.float32)


def signed_log_contrast(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return (left - right).astype(np.float32)


def bounded_score(score: np.ndarray) -> np.ndarray:
    # Log-normalized module means are non-negative; this keeps dimensions on a
    # compact, stable scale without using dataset-fitted normalization.
    return np.tanh(score / 2.5).astype(np.float32)


def build_embedding(adata) -> pd.DataFrame:
    if "counts" not in adata.layers:
        raise SystemExit('Input AnnData must contain raw counts in layers["counts"].')

    counts = as_csr_counts(adata.layers["counts"])
    log_norm = log1p_cp10k(counts)
    gene_index = make_gene_index(adata.var_names)

    module_scores: OrderedDict[str, np.ndarray] = OrderedDict()
    module_fracs: OrderedDict[str, np.ndarray] = OrderedDict()

    for name, markers in PANCREAS_MODULES.items():
        idx = marker_indices(markers, gene_index)
        module_scores[name] = bounded_score(module_mean(log_norm, idx))
        module_fracs[name] = module_detected_fraction(counts, idx)

    dims: list[np.ndarray] = []

    # Primary coordinates: fixed literature marker module intensities.
    dims.extend(module_scores.values())

    # Secondary coordinates: sparse count detection fractions for major modules.
    for name in [
        "beta_insulin",
        "alpha_glucagon",
        "delta_somatostatin",
        "pp_gamma",
        "epsilon_ghrelin",
        "acinar_digestive_enzyme",
        "ductal_bicarbonate_epithelium",
        "endothelial_blood",
        "fibroblast_stellate",
        "immune_pan_leukocyte",
    ]:
        dims.append(module_fracs[name].astype(np.float32))

    # Fixed literature-informed axes separating closely related pancreas states.
    for left, right in CONTRASTS:
        dims.append(signed_log_contrast(module_scores[left], module_scores[right]))

    embedding = np.column_stack(dims).astype(np.float32, copy=False)
    columns = {f"dim_{i}": embedding[:, i] for i in range(embedding.shape[1])}
    out = pd.DataFrame({"cell_id": adata.obs_names.astype(str).to_numpy()})
    for col, values in columns.items():
        out[col] = values
    return out


def main() -> None:
    args = parse_args()
    adata = read_h5ad(args.input)
    out = build_embedding(adata)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(output, index=False)


if __name__ == "__main__":
    main()
