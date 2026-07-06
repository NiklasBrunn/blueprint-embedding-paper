#!/usr/bin/env python3
"""Literature-marker embedding for human pancreas scRNA-seq cells.

This script uses fixed, hand-curated marker gene modules for major human
pancreas cell populations and states. It does not train a model, cluster cells,
or run pretrained-model inference.
"""

from __future__ import annotations

import argparse
import os
from collections import OrderedDict

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


PANCREAS_SIGNATURES = OrderedDict(
    [
        (
            "endocrine",
            [
                "CHGA",
                "CHGB",
                "SCG2",
                "SCG3",
                "PCSK1",
                "PCSK2",
                "ISL1",
                "NEUROD1",
                "PAX6",
            ],
        ),
        (
            "beta_cell",
            [
                "INS",
                "IAPP",
                "SLC30A8",
                "PDX1",
                "NKX6-1",
                "MAFA",
                "PCSK1",
                "G6PC2",
                "HADH",
                "UCN3",
                "RBP4",
            ],
        ),
        (
            "mature_beta",
            [
                "INS",
                "IAPP",
                "UCN3",
                "MAFA",
                "G6PC2",
                "SLC2A2",
                "SLC30A8",
                "SIX2",
                "SIX3",
            ],
        ),
        (
            "alpha_cell",
            [
                "GCG",
                "ARX",
                "IRX2",
                "MAFB",
                "PCSK2",
                "TTR",
                "DPP4",
            ],
        ),
        (
            "delta_cell",
            [
                "SST",
                "HHEX",
                "RBP4",
                "LEPR",
                "GHSR",
            ],
        ),
        (
            "pp_gamma_cell",
            [
                "PPY",
                "PYY",
                "ISL1",
                "REG1A",
            ],
        ),
        (
            "epsilon_cell",
            [
                "GHRL",
                "MLN",
                "IAPP",
            ],
        ),
        (
            "endocrine_progenitor_like",
            [
                "NEUROG3",
                "RFX6",
                "NKX2-2",
                "PAX4",
                "SOX4",
                "HES6",
                "INSM1",
            ],
        ),
        (
            "acinar",
            [
                "PRSS1",
                "PRSS2",
                "CPA1",
                "CPA2",
                "CTRB1",
                "CTRB2",
                "CELA3A",
                "CELA3B",
                "CLPS",
                "PNLIP",
                "AMY2A",
                "AMY2B",
                "REG1A",
                "REG1B",
            ],
        ),
        (
            "ductal",
            [
                "KRT19",
                "KRT7",
                "KRT8",
                "KRT18",
                "SOX9",
                "CFTR",
                "EPCAM",
                "MUC1",
                "MSLN",
                "CLDN10",
                "ANXA4",
            ],
        ),
        (
            "centroacinar_ductal",
            [
                "ALDH1A1",
                "CA2",
                "SOX9",
                "KRT19",
                "KRT7",
                "SPINK1",
            ],
        ),
        (
            "activated_ductal",
            [
                "SPP1",
                "KRT17",
                "MMP7",
                "LCN2",
                "SERPINE1",
                "TACSTD2",
                "LGALS3",
            ],
        ),
        (
            "stellate_fibroblast",
            [
                "COL1A1",
                "COL1A2",
                "COL3A1",
                "DCN",
                "LUM",
                "PDGFRA",
                "PDGFRB",
                "RGS5",
            ],
        ),
        (
            "smooth_muscle_pericyte",
            [
                "ACTA2",
                "TAGLN",
                "MYH11",
                "MCAM",
                "CSPG4",
                "RGS5",
                "PDGFRB",
            ],
        ),
        (
            "endothelial",
            [
                "PECAM1",
                "VWF",
                "KDR",
                "FLT1",
                "CDH5",
                "EMCN",
                "PLVAP",
                "RAMP2",
                "ESAM",
            ],
        ),
        (
            "immune_pan",
            [
                "PTPRC",
                "LST1",
                "AIF1",
                "CD3D",
                "CD3E",
                "TRAC",
                "NKG7",
                "GNLY",
                "MS4A1",
                "CD79A",
                "MZB1",
            ],
        ),
        (
            "macrophage_monocyte",
            [
                "LST1",
                "AIF1",
                "CD68",
                "CD14",
                "FCGR3A",
                "MS4A7",
                "C1QA",
                "C1QB",
                "C1QC",
                "LYZ",
            ],
        ),
        (
            "t_nk_cell",
            [
                "CD3D",
                "CD3E",
                "TRAC",
                "IL7R",
                "NKG7",
                "GNLY",
                "KLRD1",
                "KLRB1",
            ],
        ),
        (
            "b_plasma_cell",
            [
                "MS4A1",
                "CD79A",
                "CD79B",
                "MZB1",
                "JCHAIN",
                "IGKC",
                "IGHG1",
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
            ],
        ),
        (
            "schwann_neural",
            [
                "SOX10",
                "S100B",
                "PLP1",
                "MPZ",
                "GFRA3",
                "NRXN1",
            ],
        ),
        (
            "erythroid",
            [
                "HBA1",
                "HBA2",
                "HBB",
                "ALAS2",
                "GYPA",
            ],
        ),
        (
            "cycling",
            [
                "MKI67",
                "TOP2A",
                "PCNA",
                "HMGB2",
                "TYMS",
                "CENPF",
                "UBE2C",
            ],
        ),
        (
            "immediate_early_stress",
            [
                "FOS",
                "JUN",
                "JUNB",
                "DUSP1",
                "ATF3",
                "HSPA1A",
                "HSPA1B",
            ],
        ),
        (
            "hypoxia_glycolysis",
            [
                "VEGFA",
                "CA9",
                "SLC2A1",
                "LDHA",
                "BNIP3",
                "ENO1",
            ],
        ),
    ]
)


CONTRASTS = OrderedDict(
    [
        ("beta_vs_other_endocrine", ("beta_cell", ["alpha_cell", "delta_cell", "pp_gamma_cell", "epsilon_cell"])),
        ("alpha_vs_beta", ("alpha_cell", ["beta_cell"])),
        ("delta_vs_alpha_beta", ("delta_cell", ["alpha_cell", "beta_cell"])),
        ("pp_epsilon_minor_endocrine", ("pp_gamma_cell", ["alpha_cell", "beta_cell", "delta_cell"])),
        ("ductal_vs_acinar", ("ductal", ["acinar"])),
        ("endocrine_vs_exocrine", ("endocrine", ["acinar", "ductal"])),
        ("stromal_vs_epithelial", ("stellate_fibroblast", ["endocrine", "acinar", "ductal"])),
        ("vascular_vs_stromal", ("endothelial", ["stellate_fibroblast", "smooth_muscle_pericyte"])),
        ("myeloid_vs_lymphoid", ("macrophage_monocyte", ["t_nk_cell", "b_plasma_cell"])),
        ("ductal_activation_balance", ("activated_ductal", ["ductal"])),
        ("beta_maturation_balance", ("mature_beta", ["endocrine_progenitor_like"])),
        ("injury_vs_secretory_epithelium", ("immediate_early_stress", ["acinar", "ductal", "endocrine"])),
    ]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Embed human pancreas scRNA-seq cells using fixed literature marker counts."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, default=0, help="Run tag; accepted for interface compatibility.")
    return parser.parse_args()


def build_gene_index(var_names: pd.Index) -> dict[str, list[int]]:
    gene_to_indices: dict[str, list[int]] = {}
    for idx, name in enumerate(var_names.astype(str)):
        gene_to_indices.setdefault(name.upper(), []).append(idx)
    return gene_to_indices


def marker_indices(markers: list[str], gene_to_indices: dict[str, list[int]]) -> list[int]:
    indices: list[int] = []
    seen: set[int] = set()
    for marker in markers:
        for idx in gene_to_indices.get(marker.upper(), []):
            if idx not in seen:
                indices.append(idx)
                seen.add(idx)
    return indices


def as_csr_float(matrix) -> sparse.csr_matrix:
    if sparse.issparse(matrix):
        return matrix.tocsr().astype(np.float64)
    return sparse.csr_matrix(np.asarray(matrix, dtype=np.float64))


def log_normalized_marker_mean(
    counts: sparse.csr_matrix,
    indices: list[int],
    per_cell_scale: np.ndarray,
) -> np.ndarray:
    if not indices:
        return np.zeros(counts.shape[0], dtype=np.float64)

    sub = counts[:, indices].astype(np.float64).multiply(per_cell_scale[:, None])
    if sub.nnz:
        sub.data = np.log1p(sub.data)
    return np.asarray(sub.sum(axis=1)).ravel() / float(len(indices))


def robust_scale(features: np.ndarray) -> np.ndarray:
    scaled = np.zeros_like(features, dtype=np.float64)
    for col in range(features.shape[1]):
        values = features[:, col]
        finite = np.isfinite(values)
        if not finite.any():
            continue

        clean = values[finite]
        median = np.median(clean)
        q25, q75 = np.percentile(clean, [25, 75])
        scale = (q75 - q25) / 1.349
        if scale <= 1e-8:
            scale = np.std(clean)
        if scale <= 1e-8:
            scaled[:, col] = 0.0
        else:
            scaled[:, col] = np.clip((values - median) / scale, -8.0, 8.0)
    return scaled


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.input)

    if "counts" not in adata.layers:
        raise KeyError('Input AnnData must contain raw counts in layers["counts"].')

    counts = as_csr_float(adata.layers["counts"])
    if counts.shape[0] != adata.n_obs or counts.shape[1] != adata.n_vars:
        raise ValueError('layers["counts"] shape does not match AnnData obs/var dimensions.')

    library_size = np.asarray(counts.sum(axis=1)).ravel()
    safe_library_size = np.where(library_size > 0.0, library_size, 1.0)
    per_cell_scale = 10000.0 / safe_library_size

    gene_to_indices = build_gene_index(adata.var_names)
    score_columns: OrderedDict[str, np.ndarray] = OrderedDict()

    for signature_name, markers in PANCREAS_SIGNATURES.items():
        indices = marker_indices(markers, gene_to_indices)
        score_columns[signature_name] = log_normalized_marker_mean(counts, indices, per_cell_scale)

    for contrast_name, (positive_name, negative_names) in CONTRASTS.items():
        positive = score_columns[positive_name]
        negative = np.vstack([score_columns[name] for name in negative_names]).max(axis=0)
        score_columns[contrast_name] = positive - negative

    feature_matrix = np.column_stack(list(score_columns.values()))
    embedding = robust_scale(feature_matrix).astype(np.float32)

    output = pd.DataFrame({"cell_id": adata.obs_names.astype(str)})
    for dim in range(embedding.shape[1]):
        output[f"dim_{dim}"] = embedding[:, dim]

    output_dir = os.path.dirname(os.path.abspath(args.output))
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    output.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
