#!/usr/bin/env python3
"""Literature-marker embedding for human pancreas scRNA-seq cells.

The embedding is intentionally naive: each dimension is a fixed marker-program
score derived from published human pancreas cell-type biology. It does not fit
models, cluster cells, run PCA/UMAP, or call pretrained inference.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


def strip_jsonc(text: str) -> str:
    """Remove // and /* */ comments while preserving quoted strings."""
    out = []
    i = 0
    in_string = False
    quote = ""
    while i < len(text):
        ch = text[i]
        nxt = text[i + 1] if i + 1 < len(text) else ""
        if in_string:
            out.append(ch)
            if ch == "\\" and nxt:
                out.append(nxt)
                i += 2
                continue
            if ch == quote:
                in_string = False
            i += 1
            continue
        if ch in ("'", '"'):
            in_string = True
            quote = ch
            out.append(ch)
            i += 1
            continue
        if ch == "/" and nxt == "/":
            i = text.find("\n", i)
            if i == -1:
                break
            out.append("\n")
            i += 1
            continue
        if ch == "/" and nxt == "*":
            end = text.find("*/", i + 2)
            i = len(text) if end == -1 else end + 2
            continue
        out.append(ch)
        i += 1
    return "".join(out)


def read_config() -> dict:
    cfg_path = Path("scenario_config.jsonc")
    if not cfg_path.exists():
        return {}
    return json.loads(strip_jsonc(cfg_path.read_text(encoding="utf-8")))


def norm_symbol(symbol: object) -> str:
    return re.sub(r"\s+", "", str(symbol)).upper()


def build_gene_index(var_names: Iterable[object]) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for pos, symbol in enumerate(var_names):
        key = norm_symbol(symbol)
        index.setdefault(key, []).append(pos)
    return index


def marker_indices(markers: Iterable[str], gene_index: dict[str, list[int]]) -> list[int]:
    seen: set[int] = set()
    idx: list[int] = []
    for marker in markers:
        for gene in gene_index.get(norm_symbol(marker), []):
            if gene not in seen:
                idx.append(gene)
                seen.add(gene)
    return idx


def as_1d(x) -> np.ndarray:
    return np.asarray(x).reshape(-1)


def score_module(counts, indices: list[int], inv_total: np.ndarray) -> np.ndarray:
    """Average log1p(CP10K) over present literature markers."""
    n_cells = counts.shape[0]
    if not indices:
        return np.zeros(n_cells, dtype=np.float32)

    if sparse.issparse(counts):
        sub = counts[:, indices].tocsr(copy=True).astype(np.float64)
        sub = sub.multiply((inv_total * 10_000.0)[:, None])
        sub.data = np.log1p(sub.data)
        return (as_1d(sub.sum(axis=1)) / len(indices)).astype(np.float32)

    sub = np.asarray(counts[:, indices], dtype=np.float64)
    sub = np.log1p(sub * (inv_total * 10_000.0)[:, None])
    return sub.mean(axis=1).astype(np.float32)


# Fixed human pancreas marker modules. The order below is the embedding axis
# order. Programs include endocrine, exocrine, ductal, stromal, vascular,
# immune, and neural compartments expected in adult healthy human pancreas.
PANCREAS_MARKER_MODULES: list[tuple[str, list[str]]] = [
    (
        "pancreatic_beta_cell",
        [
            "INS",
            "IAPP",
            "SLC30A8",
            "G6PC2",
            "MAFA",
            "PDX1",
            "NKX6-1",
            "PCSK1",
            "UCN3",
            "HADH",
            "DLK1",
        ],
    ),
    (
        "pancreatic_alpha_cell",
        [
            "GCG",
            "ARX",
            "IRX2",
            "MAFB",
            "PCSK2",
            "TTR",
            "LOXL4",
            "GC",
            "CHGB",
        ],
    ),
    (
        "pancreatic_delta_cell",
        [
            "SST",
            "HHEX",
            "RBP4",
            "LEPR",
            "GHSR",
            "PCSK1N",
            "NEUROD1",
        ],
    ),
    (
        "pancreatic_pp_cell",
        [
            "PPY",
            "PAX6",
            "ISL1",
            "PCSK1",
            "PCSK2",
            "CHGA",
            "CHGB",
        ],
    ),
    (
        "pancreatic_epsilon_cell",
        [
            "GHRL",
            "MLN",
            "ISL1",
            "PAX6",
            "CHGA",
            "CHGB",
        ],
    ),
    (
        "endocrine_general",
        [
            "CHGA",
            "CHGB",
            "SCG2",
            "SCG3",
            "PCSK1",
            "PCSK2",
            "PAX6",
            "ISL1",
            "NEUROD1",
            "RFX6",
        ],
    ),
    (
        "acinar_digestive_enzyme",
        [
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
            "CLPS",
            "PNLIP",
            "AMY2A",
            "AMY2B",
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
            "MUC1",
            "CFTR",
            "SLC4A4",
            "ANXA4",
            "TSPAN8",
            "SPP1",
        ],
    ),
    (
        "centroacinar_terminal_duct",
        [
            "ALDH1A1",
            "CLDN10",
            "SPP1",
            "ONECUT2",
            "CFTR",
            "KRT19",
            "SOX9",
        ],
    ),
    (
        "epithelial_general",
        [
            "EPCAM",
            "KRT8",
            "KRT18",
            "KRT19",
            "CDH1",
            "TACSTD2",
            "CLDN4",
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
            "MMP2",
            "PDGFRA",
            "FAP",
            "THY1",
        ],
    ),
    (
        "activated_stellate_myofibroblast",
        [
            "ACTA2",
            "TAGLN",
            "MYL9",
            "CNN1",
            "COL1A1",
            "COL1A2",
            "TIMP1",
            "POSTN",
        ],
    ),
    (
        "endothelial_blood_vessel",
        [
            "PECAM1",
            "VWF",
            "KDR",
            "FLT1",
            "ESAM",
            "CDH5",
            "ENG",
            "PLVAP",
            "EMCN",
            "RAMP2",
        ],
    ),
    (
        "lymphatic_endothelial",
        [
            "PROX1",
            "LYVE1",
            "PDPN",
            "FLT4",
            "CCL21",
            "MMRN1",
            "TFF3",
        ],
    ),
    (
        "pericyte_smooth_muscle",
        [
            "RGS5",
            "PDGFRB",
            "MCAM",
            "CSPG4",
            "NOTCH3",
            "ACTA2",
            "TAGLN",
            "MYH11",
            "NDUFA4L2",
        ],
    ),
    (
        "macrophage_monocyte",
        [
            "LYZ",
            "LST1",
            "AIF1",
            "CD14",
            "CD68",
            "MS4A7",
            "FCGR3A",
            "C1QA",
            "C1QB",
            "C1QC",
            "CSF1R",
        ],
    ),
    (
        "dendritic_antigen_presenting",
        [
            "HLA-DRA",
            "HLA-DRB1",
            "HLA-DPA1",
            "HLA-DPB1",
            "CD74",
            "FCER1A",
            "CLEC10A",
            "LILRA4",
            "IRF7",
        ],
    ),
    (
        "t_cell",
        [
            "CD3D",
            "CD3E",
            "CD3G",
            "TRAC",
            "TRBC1",
            "TRBC2",
            "IL7R",
            "CCR7",
            "LCK",
        ],
    ),
    (
        "nk_cell",
        [
            "NKG7",
            "GNLY",
            "PRF1",
            "GZMB",
            "GZMH",
            "KLRD1",
            "KLRF1",
            "FCGR3A",
        ],
    ),
    (
        "b_cell",
        [
            "MS4A1",
            "CD79A",
            "CD79B",
            "BANK1",
            "CD37",
            "PAX5",
            "IGKC",
            "IGHM",
        ],
    ),
    (
        "plasma_cell",
        [
            "MZB1",
            "JCHAIN",
            "XBP1",
            "SDC1",
            "DERL3",
            "IGKC",
            "IGHG1",
            "IGHG3",
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
        "schwann_glial",
        [
            "S100B",
            "SOX10",
            "PLP1",
            "MPZ",
            "MBP",
            "GFRA3",
            "NCAM1",
        ],
    ),
    (
        "neuronal",
        [
            "SNAP25",
            "TUBB3",
            "STMN2",
            "ELAVL4",
            "RBFOX3",
            "NEFM",
            "SCN1A",
        ],
    ),
    (
        "cycling_cells",
        [
            "MKI67",
            "TOP2A",
            "UBE2C",
            "BIRC5",
            "HMGB2",
            "PCNA",
            "TYMS",
        ],
    ),
    (
        "stress_immediate_early",
        [
            "FOS",
            "JUN",
            "JUNB",
            "EGR1",
            "IER2",
            "DUSP1",
            "HSPB1",
            "HSPA1A",
            "HSPA1B",
        ],
    ),
]


def validate_scenario(config: dict) -> None:
    species = {str(x).lower() for x in config.get("species") or []}
    tissues = {str(x).lower() for x in config.get("tissue") or []}
    if species and "human" not in species:
        raise ValueError(
            "This embedding script encodes human pancreas literature markers, "
            f"but scenario_config.jsonc species is {sorted(species)}."
        )
    if tissues and "pancreas" not in tissues:
        raise ValueError(
            "This embedding script encodes human pancreas literature markers, "
            f"but scenario_config.jsonc tissue is {sorted(tissues)}."
        )


def compute_embedding(adata_obj: ad.AnnData) -> pd.DataFrame:
    if "counts" not in adata_obj.layers:
        raise KeyError('Input AnnData must contain raw counts in layers["counts"].')

    counts = adata_obj.layers["counts"]
    if counts.shape[0] != adata_obj.n_obs or counts.shape[1] != adata_obj.n_vars:
        raise ValueError('layers["counts"] must have shape n_obs by n_vars.')

    if sparse.issparse(counts):
        counts = counts.tocsr()
        total_counts = as_1d(counts.sum(axis=1)).astype(np.float64)
    else:
        counts = np.asarray(counts)
        total_counts = counts.sum(axis=1, dtype=np.float64)

    inv_total = np.zeros_like(total_counts, dtype=np.float64)
    np.divide(1.0, total_counts, out=inv_total, where=total_counts > 0)

    gene_index = build_gene_index(adata_obj.var_names)
    dims = {}
    for dim, (_name, markers) in enumerate(PANCREAS_MARKER_MODULES):
        idx = marker_indices(markers, gene_index)
        dims[f"dim_{dim}"] = score_module(counts, idx, inv_total)

    out = pd.DataFrame({"cell_id": adata_obj.obs_names.astype(str)})
    for key, values in dims.items():
        out[key] = values
    return out


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Embed human pancreas scRNA-seq cells using fixed literature marker modules."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file")
    parser.add_argument("--output", required=True, help="Output parquet path")
    parser.add_argument("--seed", required=True, type=int, help="Run tag; accepted but not used")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    config = read_config()
    validate_scenario(config)

    adata_obj = ad.read_h5ad(args.input)
    embedding = compute_embedding(adata_obj)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    embedding.to_parquet(output, index=False)


if __name__ == "__main__":
    main()
