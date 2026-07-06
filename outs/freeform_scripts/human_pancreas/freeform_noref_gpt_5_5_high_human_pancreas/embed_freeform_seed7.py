#!/usr/bin/env python3
"""Biology-defined embedding for a baseline human pancreas cell landscape.

The scenario is human pancreas, with no requested age, sex, condition, or
cell-type restriction. The embedding therefore covers the major pancreatic
epithelial, endocrine, stromal, vascular, neural, and immune populations, plus
common healthy-tissue state axes that help separate maturation, cycling, stress,
and inflammatory activation without invoking any learned model.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict


# Dimensionality is set by the biological surface to resolve: exocrine and
# islet epithelial compartments, age-unspecified developmental programs,
# endothelial/lymphatic/perivascular/stromal neighborhoods, resident immune
# lineages, neural/glial cells, and shared cell states. Each dimension is a
# direct marker-module coordinate, not a trained latent factor.
SIGNATURES = OrderedDict(
    [
        (
            "acinar_exocrine",
            [
                "PRSS1",
                "PRSS2",
                "PRSS3",
                "CPA1",
                "CPA2",
                "CPB1",
                "CTRB1",
                "CTRB2",
                "CELA3A",
                "CELA3B",
                "CLPS",
                "PNLIP",
                "PNLIPRP1",
                "AMY2A",
                "AMY2B",
                "GP2",
                "REG1A",
                "REG1B",
                "CTRL",
            ],
        ),
        (
            "ductal_centroacinar",
            [
                "KRT19",
                "KRT7",
                "KRT8",
                "KRT18",
                "SOX9",
                "HNF1B",
                "CFTR",
                "SLC4A4",
                "SLC9A3R1",
                "ANXA4",
                "MMP7",
                "SPP1",
                "MUC1",
                "CLDN18",
                "TSPAN8",
                "EPCAM",
            ],
        ),
        (
            "beta_cell",
            [
                "INS",
                "IAPP",
                "MAFA",
                "PDX1",
                "NKX6-1",
                "SLC2A2",
                "PCSK1",
                "PCSK2",
                "ABCC8",
                "KCNJ11",
                "G6PC2",
                "UCN3",
                "NPTX2",
                "RBP4",
            ],
        ),
        (
            "alpha_cell",
            [
                "GCG",
                "ARX",
                "MAFB",
                "IRX2",
                "PCSK2",
                "TTR",
                "GC",
                "ISL1",
                "LOXL4",
                "TM4SF4",
                "SLC38A5",
            ],
        ),
        (
            "delta_cell",
            [
                "SST",
                "HHEX",
                "RBP4",
                "LEPR",
                "PCSK1",
                "PCSK2",
                "GHSR",
                "CRH",
                "LY6H",
            ],
        ),
        (
            "pp_gamma_cell",
            [
                "PPY",
                "PNMAL1",
                "SLC18A2",
                "PCSK1N",
                "CHGA",
                "CHGB",
                "TM4SF4",
            ],
        ),
        (
            "epsilon_ghrelin_cell",
            [
                "GHRL",
                "MLN",
                "NEUROD1",
                "CHGA",
                "CHGB",
                "PCSK1",
                "RFX6",
            ],
        ),
        (
            "pan_endocrine_neurosecretory",
            [
                "CHGA",
                "CHGB",
                "SCG2",
                "SCG3",
                "SYP",
                "PCSK1N",
                "SLC18A1",
                "SLC18A2",
                "ENO2",
                "INSM1",
                "PAX6",
            ],
        ),
        (
            "pancreatic_epithelial_progenitor",
            [
                "PDX1",
                "SOX9",
                "ONECUT1",
                "ONECUT2",
                "HNF1B",
                "FOXA2",
                "GATA4",
                "GATA6",
                "MNX1",
                "NKX6-1",
                "NKX6-2",
            ],
        ),
        (
            "endocrine_progenitor_immature",
            [
                "NEUROG3",
                "NEUROD1",
                "INSM1",
                "NKX2-2",
                "PAX6",
                "RFX6",
                "ISL1",
                "PAX4",
                "FEV",
                "SOX4",
                "MYT1",
            ],
        ),
        (
            "blood_endothelial",
            [
                "PECAM1",
                "CDH5",
                "VWF",
                "KDR",
                "FLT1",
                "ESAM",
                "EMCN",
                "RAMP2",
                "ENG",
                "CLDN5",
                "PLVAP",
                "RGCC",
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
                "TBX1",
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
                "DES",
                "MYL9",
                "TPM2",
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
                "FBLN1",
                "SPARC",
            ],
        ),
        (
            "mesothelial",
            [
                "MSLN",
                "WT1",
                "UPK3B",
                "KRT19",
                "KRT7",
                "CALB1",
                "LRRN4",
                "ITLN1",
                "C3",
            ],
        ),
        (
            "macrophage_monocyte_dendritic",
            [
                "LYZ",
                "AIF1",
                "C1QA",
                "C1QB",
                "C1QC",
                "TYROBP",
                "FCER1G",
                "CD68",
                "MSR1",
                "CD14",
                "FCGR3A",
                "HLA-DRA",
                "CD74",
                "FCN1",
                "ITGAX",
                "LST1",
            ],
        ),
        (
            "t_nk_cell",
            [
                "CD3D",
                "CD3E",
                "CD3G",
                "TRAC",
                "CD2",
                "IL7R",
                "CCR7",
                "LCK",
                "NKG7",
                "GNLY",
                "PRF1",
                "GZMB",
                "KLRD1",
                "CTLA4",
            ],
        ),
        (
            "b_plasma_cell",
            [
                "MS4A1",
                "CD79A",
                "CD79B",
                "CD19",
                "CD22",
                "IGHM",
                "IGKC",
                "JCHAIN",
                "MZB1",
                "XBP1",
                "SDC1",
                "PRDM1",
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
                "HPGDS",
                "FCER1A",
            ],
        ),
        (
            "neural_schwann_glial",
            [
                "SOX10",
                "PLP1",
                "MPZ",
                "MBP",
                "S100B",
                "GFRA3",
                "NRXN1",
                "PHOX2B",
                "FOXD3",
                "GAP43",
            ],
        ),
        (
            "cycling",
            [
                "MKI67",
                "TOP2A",
                "HMGB2",
                "UBE2C",
                "CCNB1",
                "CCNB2",
                "CDK1",
                "PCNA",
                "STMN1",
                "TUBA1B",
                "TYMS",
                "MCM5",
            ],
        ),
        (
            "stress_interferon_inflammatory",
            [
                "HSPA1A",
                "HSPA1B",
                "DNAJB1",
                "HSPH1",
                "FOS",
                "JUN",
                "JUNB",
                "EGR1",
                "ATF3",
                "DUSP1",
                "DDIT3",
                "ATF4",
                "XBP1",
                "HERPUD1",
                "ISG15",
                "IFIT1",
                "IFIT2",
                "IFI6",
                "MX1",
                "OAS1",
                "CXCL8",
                "NFKBIA",
            ],
        ),
    ]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a marker-module embedding for human pancreas AnnData."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData path.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument(
        "--seed",
        type=int,
        required=True,
        help="Run tag accepted for orchestrator compatibility; not used.",
    )
    return parser.parse_args()


def build_gene_lookup(var_names) -> dict[str, int]:
    lookup: dict[str, int] = {}
    for idx, gene in enumerate(map(str, var_names)):
        if gene not in lookup:
            lookup[gene] = idx
        upper = gene.upper()
        if upper not in lookup:
            lookup[upper] = idx
    return lookup


def present_gene_indices(genes: list[str], lookup: dict[str, int]) -> list[int]:
    seen: set[int] = set()
    indices: list[int] = []
    for gene in genes:
        idx = lookup.get(gene)
        if idx is None:
            idx = lookup.get(gene.upper())
        if idx is not None and idx not in seen:
            seen.add(idx)
            indices.append(idx)
    return indices


def dense_float32(matrix):
    import numpy as np
    from scipy import sparse

    if sparse.issparse(matrix):
        return matrix.toarray().astype(np.float32, copy=False)
    return np.asarray(matrix, dtype=np.float32)


def gene_zscore(matrix):
    import numpy as np

    means = matrix.mean(axis=0, dtype=np.float64)
    centered = matrix - means.astype(np.float32)
    variances = np.mean(centered.astype(np.float64) ** 2, axis=0)
    scales = np.sqrt(variances).astype(np.float32)
    valid = scales > 1e-6
    if not np.any(valid):
        return np.zeros_like(matrix, dtype=np.float32)
    centered[:, valid] /= scales[valid]
    centered[:, ~valid] = 0.0
    return np.clip(centered, -6.0, 6.0).astype(np.float32, copy=False)


def robust_standardize(values):
    import numpy as np

    values = np.asarray(values, dtype=np.float32)
    median = np.median(values).astype(np.float32)
    centered = values - median
    mad = np.median(np.abs(centered)).astype(np.float32)
    scale = float(1.4826 * mad)
    if scale <= 1e-6:
        scale = float(np.std(values))
    if scale <= 1e-6:
        return np.zeros_like(values, dtype=np.float32)
    return np.clip(centered / scale, -8.0, 8.0).astype(np.float32, copy=False)


def score_signature(expression, indices: list[int], n_cells: int):
    import numpy as np

    if not indices:
        return np.zeros(n_cells, dtype=np.float32)
    subset = dense_float32(expression[:, indices])
    z = gene_zscore(subset)
    raw_score = z.mean(axis=1, dtype=np.float32)
    return robust_standardize(raw_score)


def main() -> None:
    args = parse_args()

    import anndata as ad
    import numpy as np
    import pandas as pd

    adata = ad.read_h5ad(args.input)
    n_cells = int(adata.n_obs)
    lookup = build_gene_lookup(adata.var_names)
    coords = np.zeros((n_cells, len(SIGNATURES)), dtype=np.float32)

    for dim_idx, genes in enumerate(SIGNATURES.values()):
        indices = present_gene_indices(genes, lookup)
        coords[:, dim_idx] = score_signature(adata.X, indices, n_cells)

    data = {"cell_id": pd.Series(adata.obs_names.astype(str), dtype="string")}
    for dim_idx in range(coords.shape[1]):
        data[f"dim_{dim_idx}"] = coords[:, dim_idx].astype(np.float32, copy=False)

    pd.DataFrame(data, columns=["cell_id"] + [f"dim_{i}" for i in range(coords.shape[1])]).to_parquet(
        args.output, index=False
    )


if __name__ == "__main__":
    main()
