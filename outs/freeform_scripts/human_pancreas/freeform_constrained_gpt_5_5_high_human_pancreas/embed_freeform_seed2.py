#!/usr/bin/env python3
"""Literature-derived embedding for the baseline human pancreas."""

from __future__ import annotations

import argparse
from collections.abc import Iterable


# Baron et al., Cell Systems 2016, human/mouse pancreas single-cell atlas, 10.1016/j.cels.2016.08.011;
# Muraro et al., Cell Systems 2016, human pancreas single-cell atlas, 10.1016/j.cels.2016.09.002;
# Tosti et al., Gastroenterology 2021, human pancreas single-nucleus/spatial atlas, 10.1053/j.gastro.2020.11.010;
# KEGG pancreatic secretion pathway hsa04972.
ACINAR_ZYMOGEN = [
    "PRSS1",
    "PRSS2",
    "PRSS3",
    "CTRB1",
    "CTRB2",
    "CPA1",
    "CPA2",
    "CPB1",
    "CELA2A",
    "CELA3A",
    "CELA3B",
    "PNLIP",
    "CLPS",
    "CEL",
    "REG1A",
    "REG1B",
    "SYCN",
    "CUZD1",
]

# Tosti et al., Gastroenterology 2021, acinar heterogeneity in human pancreas, 10.1053/j.gastro.2020.11.010;
# Enge et al., Cell 2017, human pancreas single-cell aging atlas, 10.1016/j.cell.2017.09.004.
ACINAR_RIBO_STRESS = [
    "RPLP1",
    "RPLP2",
    "RPL3",
    "RPS3",
    "RPS6",
    "RPS8",
    "EEF1A1",
    "HSPA1A",
    "HSPA1B",
    "HSP90AA1",
    "FOS",
    "JUN",
]

# Baron et al., Cell Systems 2016, human/mouse pancreas single-cell atlas, 10.1016/j.cels.2016.08.011;
# Muraro et al., Cell Systems 2016, human pancreas single-cell atlas, 10.1016/j.cels.2016.09.002;
# Tosti et al., Gastroenterology 2021, human pancreas single-nucleus/spatial atlas, 10.1053/j.gastro.2020.11.010;
# KEGG pancreatic secretion pathway hsa04972.
DUCTAL_EPITHELIUM = [
    "KRT7",
    "KRT8",
    "KRT18",
    "KRT19",
    "SOX9",
    "MUC1",
    "MUC6",
    "CFTR",
    "SLC4A4",
    "SLC26A6",
    "ANXA4",
    "SPP1",
    "TFF1",
    "TFF2",
]

# Segerstolpe et al., Cell Metabolism 2016, human pancreatic islet scRNA-seq, 10.1016/j.cmet.2016.08.020;
# Baron et al., Cell Systems 2016, human/mouse pancreas single-cell atlas, 10.1016/j.cels.2016.08.011;
# KEGG insulin secretion pathway hsa04911.
BETA_CELL = [
    "INS",
    "IAPP",
    "PDX1",
    "NKX6-1",
    "MAFA",
    "PCSK1",
    "SLC30A8",
    "G6PC2",
    "ADCYAP1",
    "HADH",
    "RBP4",
]

# Segerstolpe et al., Cell Metabolism 2016, human pancreatic islet scRNA-seq, 10.1016/j.cmet.2016.08.020;
# Baron et al., Cell Systems 2016, human/mouse pancreas single-cell atlas, 10.1016/j.cels.2016.08.011.
ALPHA_CELL = [
    "GCG",
    "ARX",
    "IRX2",
    "PCSK2",
    "TTR",
    "TM4SF4",
    "LOXL4",
    "MAFB",
    "CHGB",
]

# Segerstolpe et al., Cell Metabolism 2016, human pancreatic islet scRNA-seq, 10.1016/j.cmet.2016.08.020;
# Muraro et al., Cell Systems 2016, human pancreas single-cell atlas, 10.1016/j.cels.2016.09.002.
DELTA_CELL = [
    "SST",
    "HHEX",
    "RBP4",
    "LEPR",
    "BCHE",
    "ERBB4",
    "RGS2",
    "PCSK1N",
]

# Segerstolpe et al., Cell Metabolism 2016, human pancreatic islet scRNA-seq, 10.1016/j.cmet.2016.08.020;
# Baron et al., Cell Systems 2016, human/mouse pancreas single-cell atlas, 10.1016/j.cels.2016.08.011.
PP_GAMMA_CELL = [
    "PPY",
    "PYY",
    "SERTM1",
    "SLC4A10",
    "DPP4",
    "ID4",
    "CHGA",
    "CHGB",
]

# Segerstolpe et al., Cell Metabolism 2016, human pancreatic islet scRNA-seq, 10.1016/j.cmet.2016.08.020;
# Muraro et al., Cell Systems 2016, human pancreas single-cell atlas, 10.1016/j.cels.2016.09.002.
EPSILON_CELL = [
    "GHRL",
    "MLN",
    "NEUROD1",
    "ISL1",
    "CHGA",
    "CHGB",
    "PCSK1",
]

# Tosti et al., Gastroenterology 2021, human pancreas single-nucleus/spatial atlas, 10.1053/j.gastro.2020.11.010;
# The Tabula Sapiens Consortium, Science 2022, multi-organ human cell atlas, 10.1126/science.abl4896;
# MSigDB Hallmark gene set HALLMARK_ANGIOGENESIS.
ENDOTHELIAL = [
    "PECAM1",
    "VWF",
    "CDH5",
    "KDR",
    "FLT1",
    "RAMP2",
    "EMCN",
    "ESAM",
    "PLVAP",
    "CLEC14A",
    "ENG",
    "ACKR1",
]

# Tosti et al., Gastroenterology 2021, human pancreas single-nucleus/spatial atlas, 10.1053/j.gastro.2020.11.010;
# The Tabula Sapiens Consortium, Science 2022, multi-organ human cell atlas, 10.1126/science.abl4896;
# Reactome extracellular matrix organization R-HSA-1474244.
STELLATE_FIBROBLAST = [
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "FBLN1",
    "PDGFRA",
    "MMP2",
    "SPARC",
    "VIM",
    "THY1",
    "CXCL12",
]

# Tosti et al., Gastroenterology 2021, human pancreas single-nucleus/spatial atlas, 10.1053/j.gastro.2020.11.010;
# The Tabula Sapiens Consortium, Science 2022, multi-organ human cell atlas, 10.1126/science.abl4896;
# MSigDB Hallmark gene set HALLMARK_SMOOTH_MUSCLE_CONTRACTION is absent, so genes are anchored by atlas markers and Reactome muscle contraction R-HSA-397014.
MURAL_SMOOTH_MUSCLE = [
    "RGS5",
    "PDGFRB",
    "CSPG4",
    "MCAM",
    "ACTA2",
    "TAGLN",
    "MYH11",
    "MYL9",
    "NOTCH3",
    "COX4I2",
]

# Tosti et al., Gastroenterology 2021, human pancreas single-nucleus/spatial atlas, 10.1053/j.gastro.2020.11.010;
# The Tabula Sapiens Consortium, Science 2022, multi-organ human cell atlas, 10.1126/science.abl4896.
SCHWANN_NEURAL = [
    "SOX10",
    "S100B",
    "PLP1",
    "MPZ",
    "MBP",
    "ERBB3",
    "GFRA3",
    "NGFR",
    "NCAM1",
    "NRXN1",
]

# Tosti et al., Gastroenterology 2021, human pancreas single-nucleus/spatial atlas, 10.1053/j.gastro.2020.11.010;
# The Tabula Sapiens Consortium, Science 2022, multi-organ human cell atlas, 10.1126/science.abl4896;
# MSigDB Hallmark gene set HALLMARK_INFLAMMATORY_RESPONSE.
MYELOID_ANTIGEN_PRESENTING = [
    "PTPRC",
    "LST1",
    "TYROBP",
    "FCER1G",
    "C1QA",
    "C1QB",
    "C1QC",
    "CD68",
    "CD14",
    "FCGR3A",
    "HLA-DRA",
    "HLA-DRB1",
    "AIF1",
]

# Tosti et al., Gastroenterology 2021, human pancreas single-nucleus/spatial atlas, 10.1053/j.gastro.2020.11.010;
# The Tabula Sapiens Consortium, Science 2022, multi-organ human cell atlas, 10.1126/science.abl4896;
# GO adaptive immune response GO:0002250.
LYMPHOID = [
    "PTPRC",
    "CD3D",
    "CD3E",
    "TRAC",
    "TRBC1",
    "IL7R",
    "NKG7",
    "GNLY",
    "KLRD1",
    "MS4A1",
    "CD79A",
    "JCHAIN",
]


SIGNATURES = {
    "acinar_zymogen": ACINAR_ZYMOGEN,
    "acinar_ribo_stress": ACINAR_RIBO_STRESS,
    "ductal_epithelium": DUCTAL_EPITHELIUM,
    "beta_cell": BETA_CELL,
    "alpha_cell": ALPHA_CELL,
    "delta_cell": DELTA_CELL,
    "pp_gamma_cell": PP_GAMMA_CELL,
    "epsilon_cell": EPSILON_CELL,
    "endothelial": ENDOTHELIAL,
    "stellate_fibroblast": STELLATE_FIBROBLAST,
    "mural_smooth_muscle": MURAL_SMOOTH_MUSCLE,
    "schwann_neural": SCHWANN_NEURAL,
    "myeloid_antigen_presenting": MYELOID_ANTIGEN_PRESENTING,
    "lymphoid": LYMPHOID,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a fixed, literature-derived embedding for baseline human "
            "pancreas single-cell RNA-seq data."
        )
    )
    parser.add_argument("--input", required=True, help="Input .h5ad file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, required=True, help="Run tag; ignored.")
    return parser.parse_args()


def unique_preserving_order(items: Iterable[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out


def build_gene_lookup(var_names: Iterable[str]) -> dict[str, str]:
    lookup: dict[str, str] = {}
    for gene in var_names:
        gene_str = str(gene)
        lookup.setdefault(gene_str.upper(), gene_str)
    return lookup


def as_dense_float32(matrix):
    import numpy as np
    from scipy import sparse

    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=np.float32)


def signature_score(adata, genes: list[str], gene_lookup: dict[str, str], min_genes: int = 2):
    import numpy as np

    present = unique_preserving_order(
        gene_lookup[g.upper()] for g in genes if g.upper() in gene_lookup
    )
    if len(present) < min_genes:
        return np.zeros(adata.n_obs, dtype=np.float32)

    expr = as_dense_float32(adata[:, present].X)
    means = expr.mean(axis=0, dtype=np.float64)
    centered = expr - means.astype(np.float32, copy=False)
    variances = np.mean(centered * centered, axis=0, dtype=np.float64)
    stds = np.sqrt(np.maximum(variances, 1e-6)).astype(np.float32, copy=False)
    z = centered / stds
    z = np.clip(z, -5.0, 5.0)
    return z.mean(axis=1, dtype=np.float32)


def robust_scale(vector):
    import numpy as np

    vector = np.asarray(vector, dtype=np.float32)
    median = np.median(vector)
    q25, q75 = np.percentile(vector, [25, 75])
    scale = (q75 - q25) / 1.349
    if not np.isfinite(scale) or scale < 1e-6:
        scale = float(np.std(vector))
    if not np.isfinite(scale) or scale < 1e-6:
        return np.zeros_like(vector, dtype=np.float32)
    out = (vector - median) / scale
    return np.clip(out, -6.0, 6.0).astype(np.float32)


def mean_scores(scores: dict[str, object], names: list[str]):
    import numpy as np

    return np.mean([scores[name] for name in names], axis=0, dtype=np.float32)


def build_embedding(scores: dict[str, object]):
    import numpy as np

    epithelial = mean_scores(
        scores,
        [
            "acinar_zymogen",
            "ductal_epithelium",
            "beta_cell",
            "alpha_cell",
            "delta_cell",
            "pp_gamma_cell",
            "epsilon_cell",
        ],
    )
    endocrine = mean_scores(
        scores,
        ["beta_cell", "alpha_cell", "delta_cell", "pp_gamma_cell", "epsilon_cell"],
    )
    exocrine = mean_scores(scores, ["acinar_zymogen", "ductal_epithelium"])
    mesenchymal = mean_scores(scores, ["stellate_fibroblast", "mural_smooth_muscle"])
    immune = mean_scores(scores, ["myeloid_antigen_presenting", "lymphoid"])
    vascular_neural = mean_scores(scores, ["endothelial", "schwann_neural"])

    raw_dims = [
        scores["acinar_zymogen"] - scores["ductal_epithelium"],
        endocrine - exocrine,
        scores["beta_cell"] - scores["alpha_cell"],
        scores["delta_cell"] - mean_scores(scores, ["pp_gamma_cell", "epsilon_cell"]),
        scores["pp_gamma_cell"] - scores["epsilon_cell"],
        scores["endothelial"] - epithelial,
        scores["stellate_fibroblast"] - scores["mural_smooth_muscle"],
        immune - mean_scores(scores, ["epithelial_proxy"]),
        scores["myeloid_antigen_presenting"] - scores["lymphoid"],
        scores["schwann_neural"] - mean_scores(scores, ["endothelial", "mesenchymal_proxy"]),
        scores["acinar_ribo_stress"] - scores["acinar_zymogen"],
        epithelial - mean_scores(scores, ["immune_proxy", "vascular_neural_proxy", "mesenchymal_proxy"]),
    ]

    return np.column_stack([robust_scale(dim) for dim in raw_dims]).astype(np.float32)


def main() -> None:
    args = parse_args()

    import anndata as ad
    import pandas as pd

    adata = ad.read_h5ad(args.input)
    gene_lookup = build_gene_lookup(adata.var_names)

    scores = {
        name: signature_score(adata, genes, gene_lookup)
        for name, genes in SIGNATURES.items()
    }
    scores["epithelial_proxy"] = mean_scores(
        scores,
        [
            "acinar_zymogen",
            "ductal_epithelium",
            "beta_cell",
            "alpha_cell",
            "delta_cell",
            "pp_gamma_cell",
            "epsilon_cell",
        ],
    )
    scores["mesenchymal_proxy"] = mean_scores(
        scores, ["stellate_fibroblast", "mural_smooth_muscle"]
    )
    scores["immune_proxy"] = mean_scores(
        scores, ["myeloid_antigen_presenting", "lymphoid"]
    )
    scores["vascular_neural_proxy"] = mean_scores(scores, ["endothelial", "schwann_neural"])

    embedding = build_embedding(scores)
    columns = {f"dim_{idx}": embedding[:, idx].astype("float32") for idx in range(embedding.shape[1])}
    output = pd.DataFrame({"cell_id": adata.obs_names.astype(str), **columns})
    output.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
