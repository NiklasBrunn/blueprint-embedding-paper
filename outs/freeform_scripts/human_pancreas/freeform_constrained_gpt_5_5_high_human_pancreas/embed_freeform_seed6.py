#!/usr/bin/env python3
"""Literature-derived embedding for healthy human pancreas single-cell data.

The scenario is unrestricted healthy human pancreas. A compact embedding needs
to resolve the major epithelial endocrine/exocrine lineages, vascular and
stromal compartments, rare neural/mast populations, and common non-disease
states such as cell cycling and unfolded-protein/secretory stress. Therefore
K=17 dimensions are used: one robust-scaled score for each compartment/state
module below. No model is trained and no pretrained model is applied.
"""

from __future__ import annotations

import argparse
import sys
from typing import Dict, Iterable, List, Mapping, Sequence


# Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020; Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
BETA_CELL_GENES = [
    "INS",
    "IAPP",
    "MAFA",
    "PDX1",
    "NKX6-1",
    "SLC2A2",
    "G6PC2",
    "PCSK1",
    "HADH",
    "UCN3",
    "RBP4",
    "DLK1",
]

# Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020; Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
ALPHA_CELL_GENES = [
    "GCG",
    "ARX",
    "MAFB",
    "IRX2",
    "PCSK2",
    "TTR",
    "GC",
    "LOXL4",
    "TM4SF4",
    "SLC38A5",
]

# Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020; Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
DELTA_CELL_GENES = [
    "SST",
    "HHEX",
    "LEPR",
    "RBP4",
    "GHSR",
    "RGS4",
    "SEC11C",
    "BCHE",
]

# Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020; Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
PP_GAMMA_CELL_GENES = [
    "PPY",
    "PYY",
    "SLC18A2",
    "MEIS2",
    "TSPAN8",
    "SLC4A4",
    "SLC30A8",
]

# Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020; Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
EPSILON_CELL_GENES = [
    "GHRL",
    "MLN",
    "PCSK1N",
    "ISL1",
    "NEUROD1",
]

# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
ACINAR_CELL_GENES = [
    "PRSS1",
    "PRSS2",
    "CTRB1",
    "CTRB2",
    "CPA1",
    "CPA2",
    "CPB1",
    "PNLIP",
    "CELA3A",
    "CELA3B",
    "CEL",
    "CLPS",
    "AMY2A",
    "AMY2B",
    "REG1A",
    "REG1B",
]

# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
DUCTAL_CELL_GENES = [
    "KRT19",
    "KRT7",
    "KRT8",
    "KRT18",
    "SOX9",
    "EPCAM",
    "TACSTD2",
    "MUC1",
    "CFTR",
    "SLC4A4",
    "ANXA4",
    "CLDN4",
    "SPP1",
    "MSLN",
]

# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896; GO angiogenesis gene set GO:0001525.
ENDOTHELIAL_CELL_GENES = [
    "PECAM1",
    "VWF",
    "KDR",
    "FLT1",
    "ESAM",
    "EMCN",
    "CDH5",
    "CLDN5",
    "RAMP2",
    "ENG",
    "PLVAP",
    "CLEC14A",
]

# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896; Reactome extracellular matrix organization R-HSA-1474244.
FIBROBLAST_STELLATE_GENES = [
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "COL6A1",
    "DCN",
    "LUM",
    "BGN",
    "PDGFRA",
    "FAP",
    "MMP2",
    "SPARC",
    "TIMP1",
]

# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896; Reactome smooth muscle contraction R-HSA-445355.
PERICYTE_SMOOTH_MUSCLE_GENES = [
    "RGS5",
    "PDGFRB",
    "CSPG4",
    "MCAM",
    "NOTCH3",
    "ACTA2",
    "TAGLN",
    "MYH11",
    "MYL9",
    "COX4I2",
]

# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896; GO immune response gene set GO:0006955; GO complement activation GO:0006956.
MYELOID_MACROPHAGE_GENES = [
    "LST1",
    "AIF1",
    "C1QA",
    "C1QB",
    "C1QC",
    "CD68",
    "FCGR3A",
    "MS4A7",
    "HLA-DRA",
    "CD14",
    "MRC1",
    "CSF1R",
]

# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896; GO adaptive immune response gene set GO:0002250.
T_NK_LYMPHOID_GENES = [
    "CD3D",
    "CD3E",
    "TRAC",
    "IL7R",
    "CCL5",
    "NKG7",
    "GNLY",
    "KLRD1",
    "GZMB",
    "PRF1",
]

# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896; GO B cell activation gene set GO:0042113.
B_PLASMA_CELL_GENES = [
    "MS4A1",
    "CD79A",
    "CD79B",
    "BANK1",
    "CD74",
    "IGHM",
    "IGKC",
    "JCHAIN",
    "MZB1",
    "XBP1",
]

# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896; GO mast cell activation gene set GO:0097533.
MAST_CELL_GENES = [
    "TPSAB1",
    "TPSB2",
    "CPA3",
    "KIT",
    "MS4A2",
    "HDC",
    "HPGDS",
]

# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896; GO ensheathment of neurons gene set GO:0007272.
SCHWANN_NEURAL_GENES = [
    "MPZ",
    "PLP1",
    "SOX10",
    "S100B",
    "MBP",
    "PMP22",
    "ERBB3",
    "NGFR",
]

# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896; MSigDB Hallmark gene sets HALLMARK_E2F_TARGETS and HALLMARK_G2M_CHECKPOINT.
CYCLING_GENES = [
    "MKI67",
    "TOP2A",
    "CENPF",
    "UBE2C",
    "HMGB2",
    "TYMS",
    "PCNA",
    "MCM5",
    "STMN1",
]

# Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020; GO response to endoplasmic reticulum stress gene set GO:0034976; MSigDB Hallmark gene set HALLMARK_UNFOLDED_PROTEIN_RESPONSE.
UPR_SECRETORY_STRESS_GENES = [
    "HSPA5",
    "HSP90B1",
    "XBP1",
    "ATF4",
    "DDIT3",
    "HERPUD1",
    "DNAJB9",
    "PPP1R15A",
    "ERN1",
]


MODULES: Mapping[str, Sequence[str]] = {
    "beta_cell": BETA_CELL_GENES,
    "alpha_cell": ALPHA_CELL_GENES,
    "delta_cell": DELTA_CELL_GENES,
    "pp_gamma_cell": PP_GAMMA_CELL_GENES,
    "epsilon_cell": EPSILON_CELL_GENES,
    "acinar_cell": ACINAR_CELL_GENES,
    "ductal_cell": DUCTAL_CELL_GENES,
    "endothelial_cell": ENDOTHELIAL_CELL_GENES,
    "fibroblast_stellate": FIBROBLAST_STELLATE_GENES,
    "pericyte_smooth_muscle": PERICYTE_SMOOTH_MUSCLE_GENES,
    "myeloid_macrophage": MYELOID_MACROPHAGE_GENES,
    "t_nk_lymphoid": T_NK_LYMPHOID_GENES,
    "b_plasma_cell": B_PLASMA_CELL_GENES,
    "mast_cell": MAST_CELL_GENES,
    "schwann_neural": SCHWANN_NEURAL_GENES,
    "cycling": CYCLING_GENES,
    "upr_secretory_stress": UPR_SECRETORY_STRESS_GENES,
}

DIMENSION_MODULES: Sequence[str] = (
    "beta_cell",
    "alpha_cell",
    "delta_cell",
    "pp_gamma_cell",
    "epsilon_cell",
    "acinar_cell",
    "ductal_cell",
    "endothelial_cell",
    "fibroblast_stellate",
    "pericyte_smooth_muscle",
    "myeloid_macrophage",
    "t_nk_lymphoid",
    "b_plasma_cell",
    "mast_cell",
    "schwann_neural",
    "cycling",
    "upr_secretory_stress",
)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a literature-derived freeform embedding for healthy human pancreas cells."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, required=True, help="Run tag; accepted for orchestrator compatibility.")
    return parser.parse_args(argv)


def as_numpy_vector(matrix) -> "object":
    import numpy as np

    if hasattr(matrix, "A1"):
        return matrix.A1
    return np.asarray(matrix).reshape(-1)


def build_gene_index(var_names: Iterable[object]) -> Dict[str, List[int]]:
    index: Dict[str, List[int]] = {}
    for i, name in enumerate(var_names):
        symbol = str(name)
        index.setdefault(symbol, []).append(i)
        upper = symbol.upper()
        if upper != symbol:
            index.setdefault(upper, []).append(i)
    return index


def present_indices(genes: Sequence[str], gene_index: Mapping[str, Sequence[int]]) -> List[int]:
    indices: List[int] = []
    seen = set()
    for gene in genes:
        for idx in gene_index.get(gene, gene_index.get(gene.upper(), ())):
            if idx not in seen:
                indices.append(idx)
                seen.add(idx)
    return indices


def score_module(x, indices: Sequence[int], n_cells: int):
    import numpy as np

    if not indices:
        return np.zeros(n_cells, dtype=np.float64)
    values = x[:, list(indices)].mean(axis=1)
    return as_numpy_vector(values).astype(np.float64, copy=False)


def robust_zscore(values):
    import numpy as np

    values = np.asarray(values, dtype=np.float64)
    median = float(np.median(values))
    mad = float(np.median(np.abs(values - median)) * 1.4826)
    if not np.isfinite(mad) or mad < 1e-6:
        mad = float(np.std(values))
    if not np.isfinite(mad) or mad < 1e-6:
        return np.zeros(values.shape[0], dtype=np.float32)
    z = (values - median) / mad
    return np.clip(z, -8.0, 8.0).astype(np.float32)


def build_embedding(input_path: str):
    import anndata as ad
    import numpy as np
    import pandas as pd

    adata = ad.read_h5ad(input_path)
    x = adata.X
    n_cells = int(adata.n_obs)
    gene_index = build_gene_index(adata.var_names)

    raw_scores = {}
    for module_name, genes in MODULES.items():
        raw_scores[module_name] = score_module(x, present_indices(genes, gene_index), n_cells)

    data = {"cell_id": np.asarray(adata.obs_names, dtype=str)}
    for dim, module_name in enumerate(DIMENSION_MODULES):
        data[f"dim_{dim}"] = robust_zscore(raw_scores[module_name])

    return pd.DataFrame(data, columns=["cell_id"] + [f"dim_{i}" for i in range(len(DIMENSION_MODULES))])


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    embedding = build_embedding(args.input)
    embedding.to_parquet(args.output, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
