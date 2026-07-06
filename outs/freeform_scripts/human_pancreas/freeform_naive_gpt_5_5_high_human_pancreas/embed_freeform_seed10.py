#!/usr/bin/env python3
"""Literature-grounded freeform embedding for human pancreas scRNA-seq cells.

The embedding is deliberately "naive": each coordinate is a deterministic
count-derived score for a pre-specified marker gene or marker module from the
human pancreas literature.  The script does not train on the input data, infer
labels, or call any pretrained model.
"""

from __future__ import annotations

import argparse
import math
from collections import OrderedDict
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


PANCREAS_MARKER_MODULES: "OrderedDict[str, list[str]]" = OrderedDict(
    [
        # Endocrine islet lineages
        ("beta_cell_insulin", ["INS", "IAPP", "MAFA", "PDX1", "NKX6-1", "SLC30A8", "PCSK1", "G6PC2"]),
        ("alpha_cell_glucagon", ["GCG", "ARX", "MAFB", "TTR", "PCSK2", "IRX2", "LOXL4"]),
        ("delta_cell_somatostatin", ["SST", "HHEX", "RBP4", "LEPR", "BCHE"]),
        ("gamma_pp_cell", ["PPY", "PAX6", "FAP", "MEIS2"]),
        ("epsilon_ghrelin_cell", ["GHRL", "IAPP", "ISL1", "PCSK1N"]),
        ("endocrine_pan_islet", ["CHGA", "CHGB", "SCG2", "SCG3", "PAX6", "NEUROD1", "ISL1", "PCSK1N"]),
        ("endocrine_progenitor_like", ["NEUROG3", "SOX4", "PAX6", "NKX2-2", "INSM1", "RFX6"]),
        ("beta_maturity_function", ["MAFA", "PDX1", "NKX6-1", "SLC30A8", "G6PC2", "UCN3", "PCSK1"]),
        ("beta_stress_dedifferentiation", ["ALDOB", "LDHA", "SOX9", "HES1", "MYC", "DDIT3"]),
        ("islet_secretory_granule", ["PCSK1", "PCSK2", "CPE", "SLC30A8", "SCG5", "SYT7", "VAMP2"]),
        ("glucose_sensing", ["GCK", "SLC2A2", "SLC2A1", "KCNJ11", "ABCC8", "GLIS3"]),
        # Exocrine and epithelial compartments
        ("acinar_digestive_enzyme", ["PRSS1", "PRSS2", "PRSS3", "CTRB1", "CTRB2", "CPA1", "CPA2", "CEL", "PNLIP", "REG1A"]),
        ("acinar_identity", ["PTF1A", "RBPJL", "MIST1", "BHLHA15", "NR5A2", "GATA6"]),
        ("ductal_epithelial", ["KRT19", "KRT7", "KRT8", "KRT18", "SOX9", "SLC4A4", "CFTR", "MUC1", "SPP1"]),
        ("ductal_bicarbonate_fluid", ["CFTR", "SLC4A4", "SLC26A6", "AQP1", "AQP5", "CLDN10"]),
        ("centroacinar_ductal", ["KRT19", "SOX9", "MMP7", "S100P", "ANXA4", "TFF1"]),
        ("epithelial_injury_regeneration", ["REG1A", "REG1B", "REG3A", "SPINK1", "TFF1", "TFF2", "KRT19", "SPP1"]),
        # Stromal and vascular compartments
        ("pancreatic_stellate_fibroblast", ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "COL6A3", "ACTA2"]),
        ("activated_stellate_myofibroblast", ["ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "FN1", "POSTN", "THY1"]),
        ("endothelial", ["PECAM1", "VWF", "KDR", "FLT1", "ESAM", "EMCN", "PLVAP", "CDH5"]),
        ("capillary_endothelium", ["RGCC", "RAMP2", "CA4", "ACKR1", "PLVAP", "CLDN5"]),
        ("pericyte_vascular_smooth_muscle", ["RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3", "ACTA2", "MYH11"]),
        ("schwann_neural", ["S100B", "SOX10", "PLP1", "MPZ", "GFRA3", "NGFR"]),
        # Immune compartments
        ("immune_pan_leukocyte", ["PTPRC", "B2M", "HLA-DRA", "HLA-DRB1", "CXCR4", "LST1"]),
        ("macrophage_monocyte", ["LST1", "AIF1", "C1QA", "C1QB", "C1QC", "CD14", "FCGR3A", "MS4A7"]),
        ("dendritic_antigen_presenting", ["HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CD74", "ITGAX", "CLEC10A", "FCER1A"]),
        ("t_cell", ["CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CCR7", "NKG7"]),
        ("b_cell", ["MS4A1", "CD79A", "CD79B", "BANK1", "CD74", "HLA-DRA"]),
        ("plasma_cell", ["MZB1", "JCHAIN", "XBP1", "SDC1", "IGHG1", "IGKC"]),
        ("mast_cell", ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC"]),
        # Broad biological axes useful within pancreas preparations
        ("cell_cycle_g1s", ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2"]),
        ("cell_cycle_g2m", ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TOP2A", "MKI67", "CENPF"]),
        ("interferon_response", ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "STAT1", "IFI6"]),
        ("hypoxia_stress", ["HIF1A", "VEGFA", "LDHA", "ENO1", "SLC2A1", "BNIP3", "DDIT4"]),
        ("apoptosis_er_stress", ["DDIT3", "ATF4", "HSPA5", "XBP1", "HERPUD1", "DNAJB9", "PPP1R15A"]),
        ("mitochondrial_activity", ["MT-CO1", "MT-CO2", "MT-CO3", "MT-CYB", "MT-ND1", "MT-ND2", "MT-ATP6"]),
        ("ribosomal_translation", ["RPL3", "RPL5", "RPL10", "RPL13A", "RPS3", "RPS6", "RPS18", "RPS27"]),
    ]
)


PANCREAS_SINGLE_MARKERS: list[str] = [
    "INS",
    "IAPP",
    "GCG",
    "SST",
    "PPY",
    "GHRL",
    "CHGA",
    "CHGB",
    "PCSK1",
    "PCSK2",
    "SLC30A8",
    "MAFA",
    "PDX1",
    "NKX6-1",
    "G6PC2",
    "UCN3",
    "ARX",
    "MAFB",
    "HHEX",
    "KRT19",
    "SOX9",
    "CFTR",
    "MUC1",
    "SPP1",
    "PRSS1",
    "PRSS2",
    "CTRB1",
    "CPA1",
    "CPA2",
    "CEL",
    "PNLIP",
    "REG1A",
    "SPINK1",
    "COL1A1",
    "COL3A1",
    "DCN",
    "LUM",
    "ACTA2",
    "PECAM1",
    "VWF",
    "PLVAP",
    "RGS5",
    "PDGFRB",
    "PTPRC",
    "LST1",
    "HLA-DRA",
    "CD3D",
    "MS4A1",
    "MZB1",
    "TPSAB1",
    "S100B",
    "MKI67",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Embed human pancreas cells with literature marker count scores.")
    parser.add_argument("--input", required=True, help="Input .h5ad file with counts in layers['counts'].")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, required=True, help="Run tag; accepted for interface compatibility.")
    return parser.parse_args()


def read_h5ad(path: str):
    try:
        import anndata as ad
    except ImportError as exc:
        raise SystemExit("This script requires anndata to read .h5ad inputs.") from exc
    return ad.read_h5ad(path)


def as_upper_symbol(symbol: object) -> str:
    return str(symbol).strip().upper()


def build_gene_index(var_names: Sequence[object]) -> dict[str, int]:
    gene_to_index: dict[str, int] = {}
    for i, name in enumerate(var_names):
        key = as_upper_symbol(name)
        gene_to_index.setdefault(key, i)
        if "." in key:
            gene_to_index.setdefault(key.split(".", 1)[0], i)
    return gene_to_index


def indices_for_genes(gene_to_index: dict[str, int], genes: Iterable[str]) -> list[int]:
    seen: set[int] = set()
    indices: list[int] = []
    for gene in genes:
        idx = gene_to_index.get(as_upper_symbol(gene))
        if idx is not None and idx not in seen:
            indices.append(idx)
            seen.add(idx)
    return indices


def to_numpy_1d(values) -> np.ndarray:
    array = np.asarray(values)
    if array.ndim == 2 and array.shape[1] == 1:
        array = array[:, 0]
    return array.astype(np.float64, copy=False).ravel()


def sparse_aware_column_sum(matrix, indices: Sequence[int]) -> np.ndarray:
    if not indices:
        return np.zeros(matrix.shape[0], dtype=np.float64)
    return to_numpy_1d(matrix[:, list(indices)].sum(axis=1))


def safe_log1p_rate(count_sum: np.ndarray, library_size: np.ndarray, scale: float = 10000.0) -> np.ndarray:
    rate = np.divide(count_sum * scale, library_size, out=np.zeros_like(count_sum, dtype=np.float64), where=library_size > 0)
    return np.log1p(rate)


def zscore(values: np.ndarray) -> np.ndarray:
    values = values.astype(np.float64, copy=False)
    finite = np.isfinite(values)
    if not finite.any():
        return np.zeros(values.shape[0], dtype=np.float64)
    mean = float(values[finite].mean())
    std = float(values[finite].std())
    if not math.isfinite(std) or std < 1e-8:
        return values - mean
    output = (values - mean) / std
    output[~np.isfinite(output)] = 0.0
    return output


def main() -> None:
    args = parse_args()
    adata = read_h5ad(args.input)

    if "counts" not in adata.layers:
        raise SystemExit("Input AnnData must contain raw counts in layers['counts'].")

    counts = adata.layers["counts"]
    gene_to_index = build_gene_index(list(adata.var_names))
    library_size = to_numpy_1d(counts.sum(axis=1))

    features: "OrderedDict[str, np.ndarray]" = OrderedDict()

    module_raw_scores: dict[str, np.ndarray] = {}
    for module_name, genes in PANCREAS_MARKER_MODULES.items():
        idx = indices_for_genes(gene_to_index, genes)
        score = safe_log1p_rate(sparse_aware_column_sum(counts, idx), library_size)
        if idx:
            score = score / math.sqrt(len(idx))
        module_raw_scores[module_name] = score
        features[f"module_{module_name}"] = zscore(score)

    for gene in PANCREAS_SINGLE_MARKERS:
        idx = indices_for_genes(gene_to_index, [gene])
        score = safe_log1p_rate(sparse_aware_column_sum(counts, idx), library_size)
        features[f"marker_{as_upper_symbol(gene).replace('-', '_')}"] = zscore(score)

    # Fixed literature contrasts emphasize biologically meaningful separations in pancreas.
    contrasts = OrderedDict(
        [
            ("beta_vs_alpha", ("beta_cell_insulin", "alpha_cell_glucagon")),
            ("endocrine_vs_exocrine", ("endocrine_pan_islet", "acinar_digestive_enzyme")),
            ("ductal_vs_acinar", ("ductal_epithelial", "acinar_digestive_enzyme")),
            ("stromal_vs_epithelial", ("pancreatic_stellate_fibroblast", "ductal_epithelial")),
            ("vascular_vs_stromal", ("endothelial", "pancreatic_stellate_fibroblast")),
            ("immune_vs_parenchymal", ("immune_pan_leukocyte", "endocrine_pan_islet")),
            ("activated_stellate_vs_quiescent_matrix", ("activated_stellate_myofibroblast", "pancreatic_stellate_fibroblast")),
            ("beta_maturity_vs_stress", ("beta_maturity_function", "beta_stress_dedifferentiation")),
            ("secretory_granule_vs_er_stress", ("islet_secretory_granule", "apoptosis_er_stress")),
            ("proliferation_g2m_vs_g1s", ("cell_cycle_g2m", "cell_cycle_g1s")),
        ]
    )
    for name, (left, right) in contrasts.items():
        features[f"contrast_{name}"] = zscore(module_raw_scores[left] - module_raw_scores[right])

    data = {"cell_id": np.asarray(adata.obs_names, dtype=object)}
    for dim, (_, values) in enumerate(features.items()):
        data[f"dim_{dim}"] = values.astype(np.float32)

    pd.DataFrame(data).to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
