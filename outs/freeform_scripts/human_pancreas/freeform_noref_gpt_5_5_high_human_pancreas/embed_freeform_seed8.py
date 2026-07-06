#!/usr/bin/env python3
"""Marker-signature embedding for baseline human pancreas.

Scenario read from scenario_config.jsonc:
    species: human
    tissue: pancreas
    celltype: unspecified, so cover the full pancreatic cell landscape
    age/sex/condition/notes: unspecified baseline scope

The 16 dimensions below are chosen to cover the expected adult/baseline human
pancreas landscape without inventing additional constraints: endocrine islet
lineages, exocrine epithelium, stromal and vascular support cells, immune
compartments, and two broad state axes commonly visible in dissociated tissue.
"""

from __future__ import annotations

import argparse
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple


SIGNATURES: Mapping[str, Tuple[str, ...]] = {
    # Pan-endocrine and islet-lineage programs.
    "endocrine": (
        "CHGA",
        "CHGB",
        "SCG2",
        "SCG3",
        "SCG5",
        "PCSK1N",
        "PAX6",
        "NEUROD1",
        "ISL1",
        "RFX6",
        "SLC18A2",
        "KCNJ11",
    ),
    "beta": (
        "INS",
        "IAPP",
        "MAFA",
        "PDX1",
        "NKX6-1",
        "PCSK1",
        "SLC30A8",
        "G6PC2",
        "ABCC8",
        "KCNJ11",
        "UCN3",
        "HADH",
    ),
    "alpha": (
        "GCG",
        "ARX",
        "MAFB",
        "IRX2",
        "PCSK2",
        "TTR",
        "LOX",
        "TM4SF4",
        "GC",
        "SLC38A5",
    ),
    "delta": (
        "SST",
        "HHEX",
        "RBP4",
        "LEPR",
        "GHSR",
        "RGS4",
        "PCSK1",
        "NEUROD1",
    ),
    "pp_epsilon": (
        "PPY",
        "PNOC",
        "SERTM1",
        "GHRL",
        "RBP4",
        "TSPAN8",
        "REG4",
        "PAX6",
    ),
    # Exocrine epithelium.
    "acinar": (
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
        "CLPS",
        "REG1A",
        "REG1B",
        "AMY2A",
        "AMY2B",
        "SPINK1",
    ),
    "ductal": (
        "KRT19",
        "KRT7",
        "KRT8",
        "KRT18",
        "SOX9",
        "CFTR",
        "SLC4A4",
        "SLC26A3",
        "MUC1",
        "MMP7",
        "TACSTD2",
        "EPCAM",
        "ANXA4",
        "ONECUT2",
        "AQP1",
    ),
    # Stromal and vascular support compartments.
    "fibro_stellate": (
        "COL1A1",
        "COL1A2",
        "COL3A1",
        "COL6A1",
        "COL6A2",
        "DCN",
        "LUM",
        "PDGFRA",
        "FAP",
        "THY1",
        "VIM",
        "MMP2",
        "TIMP1",
        "RGS5",
    ),
    "endothelial": (
        "PECAM1",
        "VWF",
        "KDR",
        "FLT1",
        "EMCN",
        "CDH5",
        "CLDN5",
        "RAMP2",
        "ESAM",
        "PLVAP",
        "LYVE1",
        "PROX1",
        "PDPN",
        "FLT4",
    ),
    "mural": (
        "ACTA2",
        "TAGLN",
        "MYH11",
        "MCAM",
        "PDGFRB",
        "CSPG4",
        "RGS5",
        "NOTCH3",
        "MYL9",
        "TPM2",
        "CNN1",
    ),
    # Immune cells commonly recovered from human pancreatic samples.
    "t_nk": (
        "CD3D",
        "CD3E",
        "CD3G",
        "TRAC",
        "IL7R",
        "CCR7",
        "NKG7",
        "GNLY",
        "GZMB",
        "PRF1",
        "KLRD1",
        "KLRB1",
    ),
    "b_plasma": (
        "MS4A1",
        "CD79A",
        "CD79B",
        "BANK1",
        "CD37",
        "IGHM",
        "IGKC",
        "JCHAIN",
        "MZB1",
        "XBP1",
        "SDC1",
        "PRDM1",
    ),
    "myeloid_dc": (
        "LYZ",
        "AIF1",
        "C1QA",
        "C1QB",
        "C1QC",
        "FCGR3A",
        "LST1",
        "TYROBP",
        "MS4A7",
        "CD68",
        "CD14",
        "ITGAX",
        "FCER1A",
        "CLEC10A",
        "HLA-DRA",
        "HLA-DPB1",
    ),
    "mast_granulocyte": (
        "TPSAB1",
        "TPSB2",
        "CPA3",
        "KIT",
        "MS4A2",
        "HPGDS",
        "S100A8",
        "S100A9",
        "CSF3R",
        "FCGR3B",
    ),
    # Broad states, included as separate axes rather than extra cell types.
    "cell_cycle": (
        "MKI67",
        "TOP2A",
        "UBE2C",
        "CDK1",
        "CCNB1",
        "CCNB2",
        "HMGB2",
        "NUSAP1",
        "TYMS",
        "PCNA",
        "MCM5",
    ),
    "stress_inflammatory": (
        "HSPA1A",
        "HSPA1B",
        "HSP90AA1",
        "DNAJB1",
        "JUN",
        "FOS",
        "ATF3",
        "DDIT3",
        "XBP1",
        "ISG15",
        "IFIT1",
        "IFIT2",
        "IFIT3",
        "MX1",
        "STAT1",
    ),
}


DIMENSIONS: Sequence[Tuple[str, Sequence[Tuple[str, float]]]] = (
    ("dim_0", (("endocrine", 1.0), ("acinar", -0.15), ("ductal", -0.15))),
    ("dim_1", (("beta", 1.0), ("alpha", -0.20), ("delta", -0.20), ("pp_epsilon", -0.20))),
    ("dim_2", (("alpha", 1.0), ("beta", -0.20), ("delta", -0.20), ("pp_epsilon", -0.20))),
    ("dim_3", (("delta", 1.0), ("beta", -0.15), ("alpha", -0.15))),
    ("dim_4", (("pp_epsilon", 1.0), ("beta", -0.15), ("alpha", -0.15), ("delta", -0.15))),
    ("dim_5", (("acinar", 1.0), ("ductal", -0.25), ("endocrine", -0.10))),
    ("dim_6", (("ductal", 1.0), ("acinar", -0.25), ("endocrine", -0.10))),
    ("dim_7", (("fibro_stellate", 1.0),)),
    ("dim_8", (("endothelial", 1.0),)),
    ("dim_9", (("mural", 1.0), ("endothelial", -0.10), ("fibro_stellate", -0.10))),
    ("dim_10", (("t_nk", 1.0),)),
    ("dim_11", (("b_plasma", 1.0),)),
    ("dim_12", (("myeloid_dc", 1.0),)),
    ("dim_13", (("mast_granulocyte", 1.0),)),
    ("dim_14", (("cell_cycle", 1.0),)),
    ("dim_15", (("stress_inflammatory", 1.0),)),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a biological marker-signature embedding for baseline human "
            "pancreas from a log-normalized AnnData .h5ad file."
        )
    )
    parser.add_argument("--input", required=True, help="Input AnnData .h5ad path.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, required=True, help="Run tag; accepted for contract compatibility.")
    return parser.parse_args()


def _upper_gene_index(var_names: Iterable[str]) -> Dict[str, List[int]]:
    gene_to_indices: Dict[str, List[int]] = {}
    for idx, gene in enumerate(var_names):
        key = str(gene).upper()
        gene_to_indices.setdefault(key, []).append(idx)
    return gene_to_indices


def _mean_expression(matrix, indices: Sequence[int], n_cells: int, np):
    if not indices:
        return np.zeros(n_cells, dtype=np.float32)
    values = matrix[:, list(indices)].mean(axis=1)
    if hasattr(values, "A1"):
        values = values.A1
    else:
        values = np.asarray(values).reshape(-1)
    return values.astype(np.float32, copy=False)


def _signature_scores(adata, np) -> Dict[str, object]:
    gene_to_indices = _upper_gene_index(adata.var_names)
    n_cells = int(adata.n_obs)
    matrix = adata.X
    scores: Dict[str, object] = {}

    for name, genes in SIGNATURES.items():
        indices: List[int] = []
        for gene in genes:
            indices.extend(gene_to_indices.get(gene.upper(), ()))
        scores[name] = _mean_expression(matrix, indices, n_cells, np)

    return scores


def _robust_scale(values, np):
    values = np.asarray(values, dtype=np.float32)
    finite = np.isfinite(values)
    if not finite.any():
        return np.zeros(values.shape[0], dtype=np.float32)

    valid = values[finite]
    median = np.median(valid)
    q25, q75 = np.percentile(valid, [25.0, 75.0])
    scale = (q75 - q25) / 1.349
    if not np.isfinite(scale) or scale < 1.0e-6:
        mad = np.median(np.abs(valid - median)) * 1.4826
        scale = mad if np.isfinite(mad) and mad >= 1.0e-6 else 1.0

    scaled = (values - median) / scale
    scaled[~finite] = 0.0
    # Bound very high hormone or immunoglobulin expression while preserving rank.
    return np.tanh(scaled / 3.0).astype(np.float32, copy=False)


def _combine_dimensions(scores: Mapping[str, object], np):
    dims = []
    for _, terms in DIMENSIONS:
        combined = None
        for name, weight in terms:
            term = np.asarray(scores[name], dtype=np.float32) * np.float32(weight)
            combined = term.copy() if combined is None else combined + term
        dims.append(_robust_scale(combined, np))
    return np.column_stack(dims).astype(np.float32, copy=False)


def main() -> None:
    args = parse_args()

    import anndata as ad
    import numpy as np
    import pandas as pd

    adata = ad.read_h5ad(args.input)
    scores = _signature_scores(adata, np)
    embedding = _combine_dimensions(scores, np)

    dim_columns = [f"dim_{idx}" for idx in range(embedding.shape[1])]
    output = pd.DataFrame(embedding, columns=dim_columns)
    output.insert(0, "cell_id", adata.obs_names.astype(str).to_numpy())

    # Enforce the exact output schema and float32 coordinate dtype.
    output = output[["cell_id", *dim_columns]]
    for column in dim_columns:
        output[column] = output[column].astype("float32")
    output["cell_id"] = output["cell_id"].astype(str)

    output.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
