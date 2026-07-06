#!/usr/bin/env python3
"""Knowledge-based embedding for a broad healthy human pancreas landscape.

The scenario is human pancreas with unspecified cell type, age, sex, and
condition.  The embedding therefore covers major pancreatic epithelial,
endocrine, stromal, vascular, immune, and common homeostatic/state programs.
No model is trained and no pretrained model is applied; coordinates are module
scores from literal gene signatures.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


@dataclass(frozen=True)
class Signature:
    name: str
    positive: tuple[str, ...]
    negative: tuple[str, ...] = ()


# The order of this list is the order of dim_0 ... dim_22 in the output.
# These dimensions are intentionally biological programs, not learned factors.
SIGNATURES: tuple[Signature, ...] = (
    Signature(
        "pancreatic_epithelial",
        (
            "EPCAM",
            "KRT8",
            "KRT18",
            "KRT19",
            "TACSTD2",
            "CLDN4",
            "CDH1",
            "ONECUT1",
        ),
        ("PTPRC", "PECAM1", "COL1A1", "VWF", "LST1"),
    ),
    Signature(
        "pan_endocrine_islet",
        (
            "CHGA",
            "CHGB",
            "SCG2",
            "SCG3",
            "PCSK1",
            "PCSK2",
            "ISL1",
            "PAX6",
            "SYP",
            "ENO2",
        ),
        ("PRSS1", "KRT19", "COL1A1", "PECAM1", "PTPRC"),
    ),
    Signature(
        "beta_cell",
        (
            "INS",
            "IAPP",
            "SLC30A8",
            "MAFA",
            "NKX6-1",
            "PDX1",
            "PCSK1",
            "ADCYAP1",
            "SIX3",
            "HADH",
        ),
        ("GCG", "SST", "PPY", "GHRL", "PRSS1", "KRT19"),
    ),
    Signature(
        "alpha_cell",
        (
            "GCG",
            "TTR",
            "ARX",
            "IRX2",
            "MAFB",
            "PCSK2",
            "TM4SF4",
            "GC",
            "LOXL4",
        ),
        ("INS", "SST", "PPY", "PRSS1", "KRT19"),
    ),
    Signature(
        "delta_cell",
        (
            "SST",
            "HHEX",
            "RBP4",
            "LEPR",
            "BCHE",
            "GHSR",
            "PCSK1N",
        ),
        ("INS", "GCG", "PPY", "PRSS1", "KRT19"),
    ),
    Signature(
        "pp_epsilon_minor_endocrine",
        (
            "PPY",
            "GHRL",
            "PNOC",
            "TM4SF4",
            "SLC22A18",
            "PCSK1N",
            "CHGA",
        ),
        ("INS", "GCG", "SST", "PRSS1", "KRT19"),
    ),
    Signature(
        "endocrine_progenitor_immature",
        (
            "NEUROG3",
            "NEUROD1",
            "RFX6",
            "PAX4",
            "PAX6",
            "ISL1",
            "INSM1",
            "SOX4",
            "FEV",
            "MYT1",
            "NKX2-2",
        ),
        ("INS", "GCG", "SST", "PRSS1", "COL1A1"),
    ),
    Signature(
        "acinar",
        (
            "PRSS1",
            "PRSS2",
            "CPA1",
            "CPA2",
            "CPB1",
            "CTRB1",
            "CTRB2",
            "PNLIP",
            "CELA3A",
            "CELA3B",
            "REG1A",
            "RBPJL",
            "AMY2A",
        ),
        ("KRT19", "INS", "GCG", "COL1A1", "PTPRC"),
    ),
    Signature(
        "ductal_centroacinar",
        (
            "KRT19",
            "KRT7",
            "SOX9",
            "CFTR",
            "MUC1",
            "ANXA4",
            "SPP1",
            "CLDN10",
            "ONECUT2",
            "HNF1B",
            "LCN2",
            "SLC4A4",
        ),
        ("PRSS1", "INS", "GCG", "COL1A1", "PTPRC"),
    ),
    Signature(
        "rare_ductal_tuft_ionocyte",
        (
            "POU2F3",
            "TRPM5",
            "AVIL",
            "SH2D6",
            "DCLK1",
            "CHAT",
            "ASCL3",
            "FOXI1",
            "BSND",
            "CFTR",
        ),
        ("PRSS1", "INS", "COL1A1", "PTPRC"),
    ),
    Signature(
        "blood_endothelial",
        (
            "PECAM1",
            "VWF",
            "CDH5",
            "KDR",
            "EMCN",
            "ESAM",
            "PLVAP",
            "FLT1",
            "RAMP2",
        ),
        ("EPCAM", "COL1A1", "PTPRC", "RGS5"),
    ),
    Signature(
        "lymphatic_endothelial",
        (
            "PROX1",
            "PDPN",
            "LYVE1",
            "FLT4",
            "CCL21",
            "MMRN1",
            "TFF3",
            "RELN",
        ),
        ("VWF", "EPCAM", "COL1A1", "PTPRC"),
    ),
    Signature(
        "mural_pericyte_smooth_muscle",
        (
            "RGS5",
            "PDGFRB",
            "CSPG4",
            "MCAM",
            "NOTCH3",
            "ACTA2",
            "MYH11",
            "TAGLN",
            "MYLK",
        ),
        ("PECAM1", "EPCAM", "PTPRC", "COL1A1"),
    ),
    Signature(
        "stellate_fibroblast_ecm",
        (
            "COL1A1",
            "COL1A2",
            "COL3A1",
            "DCN",
            "LUM",
            "PDGFRA",
            "PI16",
            "COL6A1",
            "MFAP5",
            "ADIRF",
        ),
        ("EPCAM", "PECAM1", "PTPRC", "RGS5"),
    ),
    Signature(
        "activated_stellate_remodeling",
        (
            "ACTA2",
            "TAGLN",
            "FN1",
            "POSTN",
            "CTHRC1",
            "THY1",
            "MMP2",
            "TIMP1",
            "CXCL12",
            "COL1A1",
        ),
        ("EPCAM", "PECAM1", "PTPRC", "INS", "PRSS1"),
    ),
    Signature(
        "myeloid_macrophage_dendritic",
        (
            "LST1",
            "TYROBP",
            "AIF1",
            "CD68",
            "FCGR3A",
            "C1QA",
            "C1QB",
            "HLA-DRA",
            "CD74",
            "ITGAX",
            "LILRA4",
            "CLEC10A",
        ),
        ("EPCAM", "PECAM1", "COL1A1", "CD3D"),
    ),
    Signature(
        "t_nk_cytotoxic",
        (
            "CD3D",
            "CD3E",
            "TRAC",
            "IL7R",
            "CCR7",
            "NKG7",
            "GNLY",
            "PRF1",
            "KLRD1",
            "TRDC",
        ),
        ("EPCAM", "COL1A1", "MS4A1", "LST1"),
    ),
    Signature(
        "b_plasma",
        (
            "MS4A1",
            "CD79A",
            "CD79B",
            "IGHM",
            "CD27",
            "JCHAIN",
            "MZB1",
            "SDC1",
            "XBP1",
            "BANK1",
        ),
        ("EPCAM", "COL1A1", "CD3D", "LST1"),
    ),
    Signature(
        "mast_granulocyte",
        (
            "TPSAB1",
            "TPSB2",
            "CPA3",
            "MS4A2",
            "KIT",
            "FCER1A",
            "S100A8",
            "S100A9",
            "CSF3R",
            "FCGR3B",
        ),
        ("EPCAM", "COL1A1", "CD3D", "MS4A1"),
    ),
    Signature(
        "cell_cycle",
        (
            "MKI67",
            "TOP2A",
            "UBE2C",
            "HMGB2",
            "PCNA",
            "MCM5",
            "TYMS",
            "STMN1",
            "CENPF",
            "BIRC5",
        ),
    ),
    Signature(
        "interferon_antigen_presentation",
        (
            "IFIT1",
            "IFIT2",
            "IFIT3",
            "ISG15",
            "MX1",
            "OAS1",
            "OAS2",
            "CXCL10",
            "STAT1",
            "IRF7",
            "HLA-DRA",
            "HLA-DRB1",
            "B2M",
        ),
    ),
    Signature(
        "secretory_er_hypoxia_stress",
        (
            "HSPA5",
            "DDIT3",
            "ATF3",
            "XBP1",
            "HERPUD1",
            "DNAJB9",
            "HSP90B1",
            "PDIA4",
            "VEGFA",
            "ENO1",
            "LDHA",
        ),
    ),
    Signature(
        "epithelial_injury_regeneration",
        (
            "KRT17",
            "KRT19",
            "KRT7",
            "TACSTD2",
            "SPP1",
            "LGALS3",
            "VIM",
            "FN1",
            "ITGA5",
            "EPCAM",
            "MMP7",
            "CLU",
        ),
        ("PRSS1", "INS", "PECAM1", "PTPRC", "COL1A1"),
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a knowledge-based freeform embedding for human pancreas "
            "single-cell expression data."
        )
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file")
    parser.add_argument(
        "--output", required=True, help="Output parquet path with cell_id and dim_* columns"
    )
    parser.add_argument(
        "--seed",
        type=int,
        required=True,
        help="Run tag supplied by the orchestrator; not used by this deterministic scorer",
    )
    return parser.parse_args()


def _upper_symbols(symbols: Iterable[str]) -> list[str]:
    return [str(symbol).upper() for symbol in symbols]


def _gene_index(var_names: Iterable[str]) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for i, gene in enumerate(var_names):
        key = str(gene).upper()
        index.setdefault(key, []).append(i)
    return index


def _signature_indices(genes: tuple[str, ...], index: dict[str, list[int]]) -> list[int]:
    seen: set[int] = set()
    indices: list[int] = []
    for gene in _upper_symbols(genes):
        for idx in index.get(gene, ()):
            if idx not in seen:
                seen.add(idx)
                indices.append(idx)
    return indices


def _dense_column_stats(matrix):
    import numpy as np

    arr = np.asarray(matrix, dtype=np.float32)
    means = arr.mean(axis=0, dtype=np.float64)
    sq_means = np.square(arr, dtype=np.float64).mean(axis=0)
    variances = np.maximum(sq_means - means * means, 0.0)
    return means, np.sqrt(variances)


def _sparse_column_stats(matrix):
    import numpy as np

    means = np.asarray(matrix.mean(axis=0)).ravel().astype(np.float64)
    sq_means = np.asarray(matrix.power(2).mean(axis=0)).ravel().astype(np.float64)
    variances = np.maximum(sq_means - means * means, 0.0)
    return means, np.sqrt(variances)


def _mean_gene_z(matrix, indices: list[int], sparse) -> "object":
    import numpy as np

    if not indices:
        return np.zeros(matrix.shape[0], dtype=np.float32)

    sub = matrix[:, indices]
    if sparse.issparse(sub):
        means, stds = _sparse_column_stats(sub)
    else:
        means, stds = _dense_column_stats(sub)

    valid = stds > 1.0e-8
    if not np.any(valid):
        return np.zeros(matrix.shape[0], dtype=np.float32)

    sub = sub[:, valid]
    means = means[valid]
    stds = stds[valid]
    weights = (1.0 / stds).astype(np.float64)

    summed = np.asarray(sub @ weights).ravel().astype(np.float64)
    centered = summed - float(np.sum(means * weights))
    score = centered / float(np.sum(valid))
    return score.astype(np.float32, copy=False)


def _robust_standardize(values) -> "object":
    import numpy as np

    x = np.asarray(values, dtype=np.float32)
    if x.size == 0:
        return x

    center = float(np.nanmedian(x))
    p10, p90 = np.nanpercentile(x, [10.0, 90.0])
    scale = float((p90 - p10) / 2.563)
    if not np.isfinite(scale) or scale < 1.0e-6:
        scale = float(np.nanstd(x))
    if not np.isfinite(scale) or scale < 1.0e-6:
        return np.zeros_like(x, dtype=np.float32)

    z = (x - center) / scale
    return np.clip(z, -8.0, 8.0).astype(np.float32, copy=False)


def _score_signature(matrix, sig: Signature, index: dict[str, list[int]], sparse):
    import numpy as np

    pos = _signature_indices(sig.positive, index)
    neg = _signature_indices(sig.negative, index)

    score = _mean_gene_z(matrix, pos, sparse)
    if neg:
        score = score - _mean_gene_z(matrix, neg, sparse)

    score = _robust_standardize(score)
    return np.asarray(score, dtype=np.float32)


def run_embedding(input_path: Path, output_path: Path) -> None:
    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    adata = ad.read_h5ad(input_path)
    if adata.X is None:
        raise ValueError("Input AnnData has no X matrix; expected log1p-normalized expression.")

    matrix = adata.X
    if sparse.issparse(matrix):
        matrix = matrix.tocsr()
    else:
        matrix = np.asarray(matrix, dtype=np.float32)

    index = _gene_index(adata.var_names)
    columns: dict[str, object] = {
        "cell_id": np.asarray(adata.obs_names.astype(str), dtype=object)
    }
    for dim, sig in enumerate(SIGNATURES):
        columns[f"dim_{dim}"] = _score_signature(matrix, sig, index, sparse)

    out = pd.DataFrame(columns)
    dim_cols = [f"dim_{i}" for i in range(len(SIGNATURES))]
    out = out[["cell_id", *dim_cols]]
    for col in dim_cols:
        out[col] = out[col].astype("float32", copy=False)

    parent = output_path.parent
    if str(parent) not in ("", "."):
        parent.mkdir(parents=True, exist_ok=True)
    out.to_parquet(output_path, index=False)


def main() -> None:
    args = parse_args()
    del args.seed
    run_embedding(Path(args.input), Path(args.output))


if __name__ == "__main__":
    main()
