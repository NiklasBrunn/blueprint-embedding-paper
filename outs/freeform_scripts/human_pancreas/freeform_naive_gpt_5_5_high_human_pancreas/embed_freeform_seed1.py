#!/usr/bin/env python3
"""
Literature-grounded naive cell embedding for human pancreas scRNA-seq.

The embedding uses fixed marker programs for pancreatic endocrine,
exocrine, ductal, stromal, vascular, neural, immune, cycling, stress,
and antigen-presentation states. It does not train on the input data and
does not use pretrained-model inference.
"""

from __future__ import annotations

import argparse
import sys
from collections import OrderedDict
from typing import Iterable, Sequence

import numpy as np
import pandas as pd


MARKER_PROGRAMS: "OrderedDict[str, tuple[str, ...]]" = OrderedDict(
    [
        (
            "endocrine_core",
            (
                "CHGA",
                "CHGB",
                "SCG2",
                "SCG3",
                "SCGN",
                "PCSK1",
                "PCSK2",
                "NEUROD1",
                "PAX6",
                "ISL1",
                "INSM1",
                "SYP",
                "SEMA3D",
            ),
        ),
        (
            "beta_cell",
            (
                "INS",
                "IAPP",
                "PDX1",
                "NKX6-1",
                "MAFA",
                "SLC2A2",
                "G6PC2",
                "PCSK1",
                "ABCC8",
                "KCNJ11",
                "DLK1",
                "RBP4",
            ),
        ),
        (
            "mature_beta",
            (
                "INS",
                "IAPP",
                "MAFA",
                "UCN3",
                "SLC2A2",
                "G6PC2",
                "FFAR1",
                "SYT13",
                "PDX1",
            ),
        ),
        (
            "alpha_cell",
            (
                "GCG",
                "ARX",
                "IRX2",
                "MAFB",
                "TTR",
                "PCSK2",
                "LOX",
                "GC",
                "TM4SF4",
                "CRYBA2",
            ),
        ),
        (
            "delta_cell",
            (
                "SST",
                "HHEX",
                "RBP4",
                "LEPR",
                "PCSK1",
                "PCSK2",
                "SLC2A1",
                "APP",
            ),
        ),
        (
            "pp_gamma_cell",
            (
                "PPY",
                "PNOC",
                "FAP",
                "SLC30A8",
                "PCSK1",
                "PCSK2",
                "CHGB",
                "PAX6",
            ),
        ),
        (
            "epsilon_cell",
            (
                "GHRL",
                "MLN",
                "PAX6",
                "ISL1",
                "CHGA",
                "CHGB",
                "PCSK1",
                "PCSK2",
            ),
        ),
        (
            "enteroendocrine_like",
            (
                "TAC1",
                "TAC3",
                "NPY",
                "PYY",
                "SCT",
                "CCK",
                "GAST",
                "REG4",
                "NEUROD1",
            ),
        ),
        (
            "acinar",
            (
                "PRSS1",
                "PRSS2",
                "PRSS3",
                "CPA1",
                "CPA2",
                "CPB1",
                "CTRB1",
                "CTRB2",
                "CLPS",
                "PNLIP",
                "CELA3A",
                "CELA3B",
                "REG1A",
                "REG1B",
                "AMY2A",
                "AMY2B",
                "SYCN",
            ),
        ),
        (
            "ductal",
            (
                "KRT19",
                "KRT7",
                "KRT8",
                "KRT18",
                "KRT17",
                "EPCAM",
                "TACSTD2",
                "MUC1",
                "CFTR",
                "SLC4A4",
                "SOX9",
                "SPP1",
                "TFF1",
                "TFF2",
                "ANXA4",
                "FXYD2",
            ),
        ),
        (
            "centroacinar_terminal_duct",
            (
                "CLDN10",
                "CLDN18",
                "AQP1",
                "SOX9",
                "SLC4A4",
                "SPP1",
                "KRT19",
                "KRT7",
                "MMP7",
                "MUC6",
            ),
        ),
        (
            "epithelial",
            (
                "EPCAM",
                "TACSTD2",
                "KRT8",
                "KRT18",
                "KRT19",
                "MUC1",
                "CLDN4",
                "CLDN7",
            ),
        ),
        (
            "stellate_fibroblast",
            (
                "COL1A1",
                "COL1A2",
                "COL3A1",
                "COL6A1",
                "DCN",
                "LUM",
                "PDGFRA",
                "SPARC",
                "MMP2",
                "VIM",
                "THY1",
                "C7",
            ),
        ),
        (
            "activated_stellate_myofibroblast",
            (
                "ACTA2",
                "TAGLN",
                "MYL9",
                "TPM2",
                "POSTN",
                "COL1A1",
                "COL5A1",
                "COL12A1",
                "CTHRC1",
                "SPARC",
            ),
        ),
        (
            "endothelial",
            (
                "PECAM1",
                "VWF",
                "KDR",
                "FLT1",
                "CLDN5",
                "ENG",
                "ESAM",
                "RAMP2",
                "PLVAP",
                "EMCN",
                "CDH5",
            ),
        ),
        (
            "pericyte_smooth_muscle",
            (
                "RGS5",
                "PDGFRB",
                "MCAM",
                "CSPG4",
                "NOTCH3",
                "DES",
                "MYH11",
                "ACTA2",
                "TAGLN",
                "MYL9",
            ),
        ),
        (
            "immune_myeloid",
            (
                "PTPRC",
                "LST1",
                "TYROBP",
                "FCER1G",
                "C1QA",
                "C1QB",
                "C1QC",
                "CD68",
                "AIF1",
                "LYZ",
                "CSF1R",
                "FCGR3A",
            ),
        ),
        (
            "t_cell_nk",
            (
                "PTPRC",
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
            ),
        ),
        (
            "b_plasma",
            (
                "PTPRC",
                "MS4A1",
                "CD79A",
                "CD79B",
                "BANK1",
                "MZB1",
                "JCHAIN",
                "XBP1",
                "IGHG1",
                "IGKC",
            ),
        ),
        (
            "mast_cell",
            (
                "TPSAB1",
                "TPSB2",
                "CPA3",
                "KIT",
                "MS4A2",
                "HDC",
                "GATA2",
                "HPGDS",
            ),
        ),
        (
            "neural_schwann",
            (
                "S100B",
                "SOX10",
                "PLP1",
                "MPZ",
                "GFRA3",
                "NGFR",
                "ERBB3",
                "NRXN1",
            ),
        ),
        (
            "mesothelial",
            (
                "MSLN",
                "WT1",
                "UPK3B",
                "KRT19",
                "KRT8",
                "KRT18",
                "CALB2",
                "ITLN1",
            ),
        ),
        (
            "erythroid",
            (
                "HBB",
                "HBA1",
                "HBA2",
                "ALAS2",
                "SLC4A1",
                "GYPA",
                "AHSP",
            ),
        ),
        (
            "cycling",
            (
                "MKI67",
                "TOP2A",
                "UBE2C",
                "PCNA",
                "HMGB2",
                "CENPF",
                "STMN1",
                "TYMS",
            ),
        ),
        (
            "antigen_presentation",
            (
                "HLA-DRA",
                "HLA-DRB1",
                "HLA-DPA1",
                "HLA-DPB1",
                "HLA-DQA1",
                "HLA-DQB1",
                "CD74",
                "B2M",
                "HLA-A",
                "HLA-B",
                "HLA-C",
            ),
        ),
        (
            "interferon_inflammation",
            (
                "ISG15",
                "IFI6",
                "IFIT1",
                "IFIT2",
                "IFIT3",
                "MX1",
                "OAS1",
                "STAT1",
                "CXCL10",
                "IRF7",
            ),
        ),
        (
            "cytokine_chemokine",
            (
                "CXCL8",
                "CXCL2",
                "CXCL3",
                "CCL2",
                "CCL3",
                "CCL4",
                "IL1B",
                "TNF",
                "NFKBIA",
                "JUNB",
            ),
        ),
        (
            "stress_unfolded_protein",
            (
                "FOS",
                "JUN",
                "ATF3",
                "DDIT3",
                "HSPA1A",
                "HSPA1B",
                "HSP90AA1",
                "XBP1",
                "HERPUD1",
                "DNAJB1",
            ),
        ),
    ]
)


ALIASES = {
    "NKX6-1": ("NKX6-1", "NKX6.1", "NKX61"),
    "HLA-DRA": ("HLA-DRA", "HLADRA"),
    "HLA-DRB1": ("HLA-DRB1", "HLADRB1"),
    "HLA-DPA1": ("HLA-DPA1", "HLADPA1"),
    "HLA-DPB1": ("HLA-DPB1", "HLADPB1"),
    "HLA-DQA1": ("HLA-DQA1", "HLADQA1"),
    "HLA-DQB1": ("HLA-DQB1", "HLADQB1"),
    "HLA-A": ("HLA-A", "HLAA"),
    "HLA-B": ("HLA-B", "HLAB"),
    "HLA-C": ("HLA-C", "HLAC"),
}


DERIVED_AXES: tuple[tuple[str, str, str], ...] = (
    ("endocrine_vs_exocrine", "endocrine_core", "exocrine_max"),
    ("acinar_vs_ductal", "acinar", "ductal"),
    ("beta_vs_alpha", "beta_cell", "alpha_cell"),
    ("delta_vs_pp", "delta_cell", "pp_gamma_cell"),
    ("vascular_vs_stromal", "endothelial", "stromal_max"),
    ("immune_vs_epithelial", "immune_max", "epithelial"),
    ("myeloid_vs_lymphoid", "immune_myeloid", "lymphoid_max"),
    ("stellate_activation", "activated_stellate_myofibroblast", "stellate_fibroblast"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Embed human pancreas scRNA-seq cells with fixed literature marker programs."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, default=0, help="Run tag; accepted for interface compatibility.")
    return parser.parse_args()


def import_anndata():
    try:
        import anndata as ad
    except ImportError as exc:
        raise SystemExit("embed_freeform.py requires the 'anndata' package to read .h5ad files.") from exc
    return ad


def is_sparse_matrix(x) -> bool:
    try:
        from scipy import sparse
    except ImportError:
        return False
    return sparse.issparse(x)


def normalize_gene_name(name: object) -> str:
    return str(name).strip().upper()


def build_gene_index(var_names: Sequence[object]) -> dict[str, int]:
    index: dict[str, int] = {}
    for i, name in enumerate(var_names):
        key = normalize_gene_name(name)
        if key and key not in index:
            index[key] = i
    return index


def marker_indices(markers: Iterable[str], gene_index: dict[str, int]) -> list[int]:
    indices: list[int] = []
    seen: set[int] = set()
    for marker in markers:
        candidates = ALIASES.get(marker, (marker,))
        for candidate in candidates:
            idx = gene_index.get(normalize_gene_name(candidate))
            if idx is not None and idx not in seen:
                seen.add(idx)
                indices.append(idx)
                break
    return indices


def dense_vector(values) -> np.ndarray:
    arr = np.asarray(values)
    if arr.ndim == 2:
        arr = arr.reshape(arr.shape[0])
    return arr.astype(np.float32, copy=False)


def row_mean(matrix, indices: Sequence[int], n_obs: int) -> np.ndarray:
    if not indices:
        return np.zeros(n_obs, dtype=np.float32)
    subset = matrix[:, list(indices)]
    return dense_vector(subset.mean(axis=1))


def row_detection_fraction(counts, indices: Sequence[int], n_obs: int) -> np.ndarray:
    if not indices:
        return np.zeros(n_obs, dtype=np.float32)
    subset = counts[:, list(indices)]
    if is_sparse_matrix(subset):
        detected = subset.copy()
        detected.data = (detected.data > 0).astype(np.float32, copy=False)
        return dense_vector(detected.mean(axis=1))
    return (np.asarray(subset) > 0).mean(axis=1).astype(np.float32, copy=False)


def bounded_expression_score(mean_log_expr: np.ndarray) -> np.ndarray:
    mean_log_expr = np.maximum(mean_log_expr.astype(np.float32, copy=False), 0.0)
    return (1.0 - np.exp(-mean_log_expr / np.float32(1.25))).astype(np.float32, copy=False)


def module_score(log_matrix, count_matrix, indices: Sequence[int], n_obs: int) -> np.ndarray:
    mean_expr = row_mean(log_matrix, indices, n_obs)
    detection = row_detection_fraction(count_matrix, indices, n_obs)
    return (0.68 * bounded_expression_score(mean_expr) + 0.32 * detection).astype(np.float32, copy=False)


def row_max(score_map: dict[str, np.ndarray], names: Sequence[str], n_obs: int) -> np.ndarray:
    if not names:
        return np.zeros(n_obs, dtype=np.float32)
    stacked = np.vstack([score_map[name] for name in names])
    return stacked.max(axis=0).astype(np.float32, copy=False)


def safe_signed_difference(left: np.ndarray, right: np.ndarray) -> np.ndarray:
    return np.tanh((left - right).astype(np.float32, copy=False) * np.float32(2.2)).astype(
        np.float32, copy=False
    )


def get_count_matrix(adata):
    if "counts" not in adata.layers:
        raise SystemExit("Input AnnData must contain raw counts in layers['counts'].")
    return adata.layers["counts"]


def main() -> None:
    args = parse_args()
    ad = import_anndata()
    adata = ad.read_h5ad(args.input)

    n_obs = int(adata.n_obs)
    gene_index = build_gene_index(adata.var_names)
    count_matrix = get_count_matrix(adata)
    log_matrix = adata.X

    scores: dict[str, np.ndarray] = {}
    for program_name, markers in MARKER_PROGRAMS.items():
        indices = marker_indices(markers, gene_index)
        scores[program_name] = module_score(log_matrix, count_matrix, indices, n_obs)

    scores["exocrine_max"] = row_max(scores, ("acinar", "ductal", "centroacinar_terminal_duct"), n_obs)
    scores["stromal_max"] = row_max(
        scores, ("stellate_fibroblast", "activated_stellate_myofibroblast", "pericyte_smooth_muscle"), n_obs
    )
    scores["immune_max"] = row_max(scores, ("immune_myeloid", "t_cell_nk", "b_plasma", "mast_cell"), n_obs)
    scores["lymphoid_max"] = row_max(scores, ("t_cell_nk", "b_plasma"), n_obs)

    dims: list[np.ndarray] = []

    for _axis_name, left_name, right_name in DERIVED_AXES:
        dims.append(safe_signed_difference(scores[left_name], scores[right_name]))

    # Direct lineage and state coordinates follow the signed hierarchy axes.
    for program_name in MARKER_PROGRAMS.keys():
        dims.append(scores[program_name].astype(np.float32, copy=False))

    dims.extend(
        [
            row_max(scores, ("beta_cell", "alpha_cell", "delta_cell", "pp_gamma_cell", "epsilon_cell"), n_obs),
            row_max(scores, ("acinar", "ductal", "epithelial"), n_obs),
            row_max(scores, ("endothelial", "pericyte_smooth_muscle"), n_obs),
            row_max(scores, ("antigen_presentation", "interferon_inflammation", "cytokine_chemokine"), n_obs),
        ]
    )

    data: dict[str, object] = {"cell_id": pd.Index(adata.obs_names.astype(str), dtype="object")}
    for i, values in enumerate(dims):
        data[f"dim_{i}"] = values.astype(np.float32, copy=False)

    out = pd.DataFrame(data)
    out.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
