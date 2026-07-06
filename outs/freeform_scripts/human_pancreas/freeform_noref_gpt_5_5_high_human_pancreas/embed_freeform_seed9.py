#!/usr/bin/env python3
"""Knowledge-built embedding for a healthy human pancreas cell landscape.

The scenario is human pancreas with no age, sex, condition, or cell-type
restriction. The embedding therefore uses explicit marker signatures for
major pancreatic epithelial, endocrine, stromal, vascular, neural, immune,
and shared cell-state axes. Coordinates are not learned; they are robust
module scores from marker genes encoded below.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


# The dimensionality is intentionally biological rather than fixed. Human
# pancreas requires separate axes for exocrine acinar and ductal epithelium,
# islet endocrine identity and principal hormone-producing lineages, rare
# endocrine subtypes/progenitor-like programs, stromal/vascular/mural/neural
# support cells, major immune compartments, and common activation states.
# K = 20 keeps these compartments distinct without splitting every marker
# pair into a separate, noisy coordinate.
SIGNATURES: "OrderedDict[str, dict[str, list[str]]]" = OrderedDict(
    [
        (
            "acinar_digestive_enzyme",
            {
                "pos": [
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
                    "PNLIP",
                    "REG1A",
                    "REG1B",
                    "RBPJL",
                    "SYCN",
                ],
                "neg": ["KRT19", "INS", "GCG", "COL1A1", "PECAM1", "PTPRC"],
            },
        ),
        (
            "ductal_bicarbonate_mucin",
            {
                "pos": [
                    "KRT19",
                    "KRT7",
                    "KRT8",
                    "KRT18",
                    "SOX9",
                    "CFTR",
                    "SLC4A4",
                    "SLC26A6",
                    "MUC1",
                    "MUC6",
                    "SPP1",
                    "ANXA4",
                    "TSPAN8",
                    "EPCAM",
                ],
                "neg": ["PRSS1", "CPA1", "INS", "GCG", "COL1A1", "PTPRC"],
            },
        ),
        (
            "centroacinar_terminal_duct",
            {
                "pos": [
                    "CFTR",
                    "SOX9",
                    "KRT19",
                    "PROM1",
                    "ALDH1A1",
                    "ANXA4",
                    "CLU",
                    "SLC12A2",
                    "AQP1",
                    "TFF1",
                    "TFF2",
                ],
                "neg": ["PRSS1", "CPA1", "INS", "GCG", "PECAM1", "COL1A1"],
            },
        ),
        (
            "pan_endocrine_islet",
            {
                "pos": [
                    "CHGA",
                    "CHGB",
                    "SCG2",
                    "SCG3",
                    "SCG5",
                    "PCSK1N",
                    "PCSK1",
                    "PCSK2",
                    "SYT4",
                    "RAB3A",
                    "SLC18A1",
                    "ISL1",
                    "PAX6",
                    "NEUROD1",
                ],
                "neg": ["PRSS1", "KRT19", "COL1A1", "PECAM1", "PTPRC"],
            },
        ),
        (
            "beta_cell_mature",
            {
                "pos": [
                    "INS",
                    "IAPP",
                    "MAFA",
                    "PDX1",
                    "NKX6-1",
                    "SLC2A2",
                    "G6PC2",
                    "PCSK1",
                    "ABCC8",
                    "KCNJ11",
                    "ADCYAP1",
                    "DLK1",
                    "RBP4",
                    "FFAR1",
                ],
                "neg": ["GCG", "SST", "PPY", "GHRL", "KRT19", "PRSS1"],
            },
        ),
        (
            "alpha_cell_glucagon",
            {
                "pos": [
                    "GCG",
                    "ARX",
                    "IRX2",
                    "MAFB",
                    "TTR",
                    "GC",
                    "TM4SF4",
                    "PCSK2",
                    "SLC7A2",
                    "LOXL4",
                    "FEV",
                ],
                "neg": ["INS", "SST", "PPY", "GHRL", "KRT19", "PRSS1"],
            },
        ),
        (
            "delta_cell_somatostatin",
            {
                "pos": [
                    "SST",
                    "HHEX",
                    "RBP4",
                    "LEPR",
                    "BCHE",
                    "PCSK1",
                    "PCSK2",
                    "UCN3",
                    "PRG4",
                    "RGS2",
                ],
                "neg": ["INS", "GCG", "PPY", "GHRL", "PRSS1", "KRT19"],
            },
        ),
        (
            "pp_epsilon_rare_endocrine",
            {
                "pos": [
                    "PPY",
                    "GHRL",
                    "PYY",
                    "NPW",
                    "IAPP",
                    "PCSK1N",
                    "CHGA",
                    "CHGB",
                    "ISL1",
                    "PAX6",
                ],
                "neg": ["PRSS1", "KRT19", "COL1A1", "PECAM1", "PTPRC"],
            },
        ),
        (
            "endocrine_progenitor_plasticity",
            {
                "pos": [
                    "NEUROG3",
                    "NEUROD1",
                    "PAX4",
                    "PAX6",
                    "ISL1",
                    "INSM1",
                    "RFX6",
                    "NKX2-2",
                    "SOX4",
                    "SOX9",
                    "HES1",
                    "DLL1",
                ],
                "neg": ["PRSS1", "CPA1", "PECAM1", "PTPRC", "COL1A1"],
            },
        ),
        (
            "stellate_fibroblast_ecm",
            {
                "pos": [
                    "COL1A1",
                    "COL1A2",
                    "COL3A1",
                    "COL5A1",
                    "COL6A1",
                    "DCN",
                    "LUM",
                    "FBLN1",
                    "PDGFRA",
                    "FAP",
                    "MMP2",
                    "SPARC",
                    "RGS5",
                    "THY1",
                ],
                "neg": ["EPCAM", "PRSS1", "INS", "PECAM1", "PTPRC"],
            },
        ),
        (
            "activated_stellate_myofibroblast",
            {
                "pos": [
                    "ACTA2",
                    "TAGLN",
                    "MYL9",
                    "TPM2",
                    "COL1A1",
                    "COL1A2",
                    "COL3A1",
                    "POSTN",
                    "FN1",
                    "CTHRC1",
                    "MMP11",
                    "TIMP1",
                    "SERPINE1",
                    "TGFBI",
                ],
                "neg": ["EPCAM", "PECAM1", "PTPRC", "INS", "GCG"],
            },
        ),
        (
            "vascular_endothelial",
            {
                "pos": [
                    "PECAM1",
                    "VWF",
                    "CDH5",
                    "KDR",
                    "FLT1",
                    "EMCN",
                    "RAMP2",
                    "ESAM",
                    "CLDN5",
                    "ENG",
                    "PLVAP",
                    "ADGRL4",
                ],
                "neg": ["EPCAM", "PRSS1", "INS", "COL1A1", "PTPRC"],
            },
        ),
        (
            "lymphatic_endothelial",
            {
                "pos": [
                    "PROX1",
                    "PDPN",
                    "LYVE1",
                    "FLT4",
                    "CCL21",
                    "MMRN1",
                    "TFF3",
                    "CLDN5",
                    "PECAM1",
                ],
                "neg": ["PRSS1", "INS", "COL1A1", "PTPRC", "KRT19"],
            },
        ),
        (
            "pericyte_smooth_muscle",
            {
                "pos": [
                    "RGS5",
                    "PDGFRB",
                    "CSPG4",
                    "MCAM",
                    "NOTCH3",
                    "ACTA2",
                    "TAGLN",
                    "MYH11",
                    "MYL9",
                    "TPM2",
                    "CNN1",
                    "DES",
                ],
                "neg": ["EPCAM", "PECAM1", "PTPRC", "INS", "PRSS1"],
            },
        ),
        (
            "myeloid_macrophage_apc",
            {
                "pos": [
                    "LYZ",
                    "LST1",
                    "AIF1",
                    "TYROBP",
                    "FCER1G",
                    "CD68",
                    "CD14",
                    "CSF1R",
                    "C1QA",
                    "C1QB",
                    "C1QC",
                    "HLA-DRA",
                    "HLA-DRB1",
                    "HLA-DPA1",
                    "HLA-DPB1",
                    "MS4A7",
                ],
                "neg": ["EPCAM", "PRSS1", "INS", "PECAM1", "COL1A1"],
            },
        ),
        (
            "t_nk_lymphocyte",
            {
                "pos": [
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
                    "CTSW",
                    "LCK",
                ],
                "neg": ["EPCAM", "PRSS1", "INS", "PECAM1", "COL1A1"],
            },
        ),
        (
            "b_plasma_cell",
            {
                "pos": [
                    "MS4A1",
                    "CD79A",
                    "CD79B",
                    "BANK1",
                    "CD74",
                    "HLA-DRA",
                    "MZB1",
                    "JCHAIN",
                    "XBP1",
                    "SDC1",
                    "IGHG1",
                    "IGKC",
                    "IGLC2",
                ],
                "neg": ["EPCAM", "PRSS1", "INS", "PECAM1", "COL1A1"],
            },
        ),
        (
            "schwann_neural",
            {
                "pos": [
                    "SOX10",
                    "S100B",
                    "PLP1",
                    "MPZ",
                    "PMP22",
                    "GFRA3",
                    "NGFR",
                    "NRXN1",
                    "NRXN3",
                    "STMN2",
                    "TUBB3",
                    "GAP43",
                ],
                "neg": ["EPCAM", "PRSS1", "INS", "PECAM1", "PTPRC"],
            },
        ),
        (
            "cycling_growth",
            {
                "pos": [
                    "MKI67",
                    "TOP2A",
                    "UBE2C",
                    "CENPF",
                    "NUSAP1",
                    "HMGB2",
                    "PCNA",
                    "MCM2",
                    "MCM5",
                    "TYMS",
                    "STMN1",
                    "CDK1",
                    "CCNB1",
                ],
                "neg": [],
            },
        ),
        (
            "stress_interferon_hypoxia",
            {
                "pos": [
                    "FOS",
                    "JUN",
                    "JUNB",
                    "DUSP1",
                    "ATF3",
                    "HSPA1A",
                    "HSPA1B",
                    "HSP90AA1",
                    "ISG15",
                    "IFIT1",
                    "IFIT3",
                    "MX1",
                    "IFI6",
                    "HIF1A",
                    "VEGFA",
                    "DDIT3",
                ],
                "neg": [],
            },
        ),
    ]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Score a human pancreas AnnData object into an interpretable "
            "knowledge-based cell-type/state embedding."
        )
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file")
    parser.add_argument("--output", required=True, help="Output parquet path")
    parser.add_argument(
        "--seed",
        type=int,
        required=True,
        help="Run tag supplied by orchestrator; accepted for contract compatibility",
    )
    return parser.parse_args()


def _as_csr_or_array(x):
    try:
        from scipy import sparse

        if sparse.issparse(x):
            return x.tocsr().astype(np.float32)
    except Exception:
        pass
    return np.asarray(x, dtype=np.float32)


def _first_gene_index(var_names: Iterable[object]) -> dict[str, int]:
    index: dict[str, int] = {}
    for i, name in enumerate(var_names):
        key = str(name).upper()
        if key not in index:
            index[key] = i
    return index


def _gene_mean_var(x):
    try:
        from scipy import sparse

        if sparse.issparse(x):
            mean = np.asarray(x.mean(axis=0)).ravel().astype(np.float32)
            mean_sq = np.asarray(x.multiply(x).mean(axis=0)).ravel().astype(np.float32)
            var = np.maximum(mean_sq - mean * mean, 0.0)
            return mean, var
    except Exception:
        pass

    mean = np.asarray(x.mean(axis=0), dtype=np.float32).ravel()
    var = np.asarray(x.var(axis=0), dtype=np.float32).ravel()
    return mean, var


def _weighted_sum_by_cell(x, weights: np.ndarray) -> np.ndarray:
    values = x @ weights.astype(np.float32)
    if hasattr(values, "A1"):
        values = values.A1
    return np.asarray(values, dtype=np.float32).ravel()


def _robust_clip_scale(values: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=np.float32)
    finite = np.isfinite(values)
    if not np.any(finite):
        return np.zeros(values.shape, dtype=np.float32)

    clean = np.where(finite, values, np.nan)
    center = np.nanmedian(clean).astype(np.float32)
    q25, q75 = np.nanpercentile(clean, [25, 75]).astype(np.float32)
    iqr = float(q75 - q25)
    if not np.isfinite(iqr) or iqr < 1e-6:
        scale = float(np.nanstd(clean))
    else:
        scale = iqr / 1.349
    if not np.isfinite(scale) or scale < 1e-6:
        return np.zeros(values.shape, dtype=np.float32)

    scaled = (values - center) / scale
    return np.clip(scaled, -8.0, 8.0).astype(np.float32)


def _score_gene_set(
    x,
    gene_to_pos: dict[str, int],
    means: np.ndarray,
    stds: np.ndarray,
    genes: list[str],
) -> np.ndarray | None:
    idx = [gene_to_pos[g.upper()] for g in genes if g.upper() in gene_to_pos]
    if not idx:
        return None

    idx_arr = np.asarray(idx, dtype=np.int64)
    weights = (1.0 / stds[idx_arr]).astype(np.float32) / float(len(idx_arr))
    offset = float(np.sum(means[idx_arr] / stds[idx_arr]) / float(len(idx_arr)))
    return _weighted_sum_by_cell(x[:, idx_arr], weights) - np.float32(offset)


def build_embedding(adata) -> pd.DataFrame:
    x = adata.X
    if x is None:
        raise ValueError("AnnData.X is required and should contain log1p-normalized values")

    x = _as_csr_or_array(x)
    n_cells = int(x.shape[0])
    gene_to_pos = _first_gene_index(adata.var_names)
    means, variances = _gene_mean_var(x)
    stds = np.sqrt(np.maximum(variances, 1e-6)).astype(np.float32)

    columns: dict[str, np.ndarray] = {
        "cell_id": np.asarray(adata.obs_names.astype(str), dtype=object)
    }

    for dim, (_, signature) in enumerate(SIGNATURES.items()):
        pos_score = _score_gene_set(x, gene_to_pos, means, stds, signature["pos"])
        if pos_score is None:
            coord = np.zeros(n_cells, dtype=np.float32)
        else:
            coord = pos_score.astype(np.float32)

        neg_genes = signature.get("neg", [])
        neg_score = _score_gene_set(x, gene_to_pos, means, stds, neg_genes)
        if neg_score is not None:
            coord = coord - 0.35 * neg_score.astype(np.float32)

        columns[f"dim_{dim}"] = _robust_clip_scale(coord)

    out = pd.DataFrame(columns)
    expected = ["cell_id"] + [f"dim_{i}" for i in range(len(SIGNATURES))]
    out = out.loc[:, expected]
    for col in expected[1:]:
        out[col] = out[col].astype(np.float32)
    return out


def main() -> None:
    args = parse_args()
    import anndata as ad

    adata = ad.read_h5ad(args.input)
    embedding = build_embedding(adata)
    output = Path(args.output)
    embedding.to_parquet(output, index=False)


if __name__ == "__main__":
    main()
