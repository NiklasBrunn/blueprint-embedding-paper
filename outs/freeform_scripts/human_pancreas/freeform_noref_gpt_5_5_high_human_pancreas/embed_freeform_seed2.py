#!/usr/bin/env python3
"""
Biology-defined embedding for broad human pancreas single-cell data.

The scenario_config.jsonc for this run is unrestricted human pancreas:
all expected pancreatic epithelial, endocrine, stromal, vascular, neural,
and resident/infiltrating immune compartments are in scope. A compact
whole-tissue representation needs more than the endocrine/exocrine axes,
because pancreas dissociation commonly contains endothelial, stellate,
pericyte/smooth muscle, neural/glial, myeloid, lymphoid, mast, and plasma
cell populations, plus state programs such as inflammation, stress, and
cycling. The dimensions below are therefore interpretable signature-score
coordinates rather than fitted latent factors.
"""

import argparse
from pathlib import Path


MODULES = [
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
            "CEL",
            "PNLIP",
            "PNLIPRP1",
            "PNLIPRP2",
            "CLPS",
            "REG1A",
            "REG1B",
            "REG3A",
            "AMY2A",
            "AMY2B",
            "CELA2A",
            "CELA3A",
            "CELA3B",
            "SYCN",
            "GP2",
            "CTRL",
            "CUZD1",
        ],
    ),
    (
        "ductal_centroacinar_epithelium",
        [
            "KRT19",
            "KRT7",
            "KRT8",
            "KRT18",
            "EPCAM",
            "MUC1",
            "MUC5B",
            "CFTR",
            "SLC4A4",
            "SLC26A6",
            "AQP1",
            "AQP5",
            "SOX9",
            "ONECUT1",
            "HNF1B",
            "TSPAN8",
            "ANXA4",
            "SPP1",
            "MMP7",
            "TFF1",
            "TFF2",
            "TFF3",
        ],
    ),
    (
        "pan_endocrine_islet",
        [
            "CHGA",
            "CHGB",
            "SCG2",
            "SCG3",
            "SCG5",
            "PCSK1",
            "PCSK2",
            "PAX6",
            "ISL1",
            "NEUROD1",
            "INSM1",
            "RFX6",
            "SLC18A2",
            "SYP",
            "SNAP25",
            "VGF",
            "IAPP",
        ],
    ),
    (
        "beta_cell_mature_insulin_secretory",
        [
            "INS",
            "IAPP",
            "MAFA",
            "PDX1",
            "NKX6-1",
            "SLC30A8",
            "PCSK1",
            "PCSK1N",
            "G6PC2",
            "SLC2A2",
            "ABCC8",
            "KCNJ11",
            "RBP4",
            "NPTX2",
            "DLK1",
            "FFAR4",
            "ADCYAP1",
            "ERO1B",
            "HADH",
        ],
    ),
    (
        "alpha_cell_glucagon",
        [
            "GCG",
            "ARX",
            "MAFB",
            "IRX1",
            "IRX2",
            "TTR",
            "LOXL4",
            "TM4SF4",
            "PCSK2",
            "DPP4",
            "GC",
            "CHGB",
            "CRYBA2",
        ],
    ),
    (
        "delta_cell_somatostatin",
        [
            "SST",
            "HES1",
            "RBP4",
            "LEPR",
            "BCHE",
            "GHSR",
            "PCSK1",
            "PCSK2",
            "CHGA",
            "CHGB",
            "LY6H",
        ],
    ),
    (
        "gamma_pp_epsilon_endocrine",
        [
            "PPY",
            "PNOC",
            "GHRL",
            "TMEM176A",
            "TMEM176B",
            "SLC6A4",
            "ALDH1A1",
            "PAX6",
            "ISL1",
            "PCSK1",
            "CHGA",
            "CHGB",
        ],
    ),
    (
        "developmental_endocrine_progenitor",
        [
            "PDX1",
            "NKX6-1",
            "NKX2-2",
            "PAX6",
            "NEUROG3",
            "NEUROD1",
            "INSM1",
            "RFX6",
            "SOX9",
            "ONECUT1",
            "HES1",
            "DLL1",
            "MNX1",
            "PAX4",
            "ISL1",
            "FEV",
        ],
    ),
    (
        "blood_endothelium",
        [
            "PECAM1",
            "VWF",
            "KDR",
            "FLT1",
            "ESAM",
            "CDH5",
            "EMCN",
            "RAMP2",
            "SOX18",
            "ERG",
            "CLEC14A",
            "PLVAP",
            "IFI27",
            "AQP1",
        ],
    ),
    (
        "lymphatic_endothelium",
        [
            "PROX1",
            "PDPN",
            "LYVE1",
            "FLT4",
            "CCL21",
            "MMRN1",
            "RELN",
            "TFF3",
            "ACKR4",
            "TBX1",
        ],
    ),
    (
        "stellate_fibroblast_matrix",
        [
            "COL1A1",
            "COL1A2",
            "COL3A1",
            "COL6A1",
            "COL6A2",
            "COL6A3",
            "DCN",
            "LUM",
            "BGN",
            "MMP2",
            "PDGFRA",
            "PDGFRB",
            "THY1",
            "FAP",
            "CFD",
            "APOD",
            "RGS5",
            "SPARC",
        ],
    ),
    (
        "activated_stellate_myofibroblast",
        [
            "ACTA2",
            "TAGLN",
            "MYL9",
            "TPM2",
            "CNN1",
            "COL1A1",
            "COL1A2",
            "COL5A1",
            "POSTN",
            "CTHRC1",
            "FN1",
            "TNC",
            "LOX",
            "SERPINE1",
            "INHBA",
            "TIMP1",
        ],
    ),
    (
        "pericyte_smooth_muscle",
        [
            "RGS5",
            "MCAM",
            "CSPG4",
            "NOTCH3",
            "PDGFRB",
            "ACTA2",
            "MYH11",
            "TAGLN",
            "CNN1",
            "MYLK",
            "DES",
            "NDUFA4L2",
            "COX4I2",
            "GUCY1A2",
        ],
    ),
    (
        "schwann_neural_glial",
        [
            "SOX10",
            "PLP1",
            "MPZ",
            "MBP",
            "PMP22",
            "S100B",
            "GPM6B",
            "NRXN1",
            "NRXN3",
            "NCAM1",
            "GFRA3",
            "NTRK2",
            "CHL1",
            "PRPH",
            "STMN2",
        ],
    ),
    (
        "myeloid_macrophage_monocyte",
        [
            "LYZ",
            "C1QA",
            "C1QB",
            "C1QC",
            "CSF1R",
            "CD68",
            "CD14",
            "FCGR3A",
            "MS4A7",
            "LST1",
            "AIF1",
            "TYROBP",
            "CTSS",
            "MRC1",
            "MARCO",
            "VSIG4",
        ],
    ),
    (
        "dendritic_antigen_presenting",
        [
            "HLA-DRA",
            "HLA-DRB1",
            "HLA-DPA1",
            "HLA-DPB1",
            "HLA-DQA1",
            "HLA-DQB1",
            "CD74",
            "ITGAX",
            "CLEC10A",
            "CLEC9A",
            "FCER1A",
            "LAMP3",
            "CCR7",
            "IRF8",
            "BATF3",
        ],
    ),
    (
        "t_nk_lymphocyte",
        [
            "CD3D",
            "CD3E",
            "CD3G",
            "TRAC",
            "IL7R",
            "CCR7",
            "LTB",
            "NKG7",
            "GNLY",
            "GZMB",
            "GZMK",
            "PRF1",
            "KLRD1",
            "KLRB1",
            "NCR1",
            "CD8A",
            "CD4",
        ],
    ),
    (
        "b_plasma_cell",
        [
            "MS4A1",
            "CD79A",
            "CD79B",
            "BANK1",
            "CD37",
            "BLK",
            "PAX5",
            "MZB1",
            "JCHAIN",
            "XBP1",
            "SDC1",
            "PRDM1",
            "IGHM",
            "IGKC",
            "IGLC2",
            "TNFRSF17",
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
            "GATA2",
            "VWA5A",
            "SLC18A2",
            "RGS13",
        ],
    ),
    (
        "inflammatory_interferon_chemokine",
        [
            "ISG15",
            "IFIT1",
            "IFIT2",
            "IFIT3",
            "IFI6",
            "IFI27",
            "IFI44L",
            "MX1",
            "MX2",
            "OAS1",
            "OAS2",
            "STAT1",
            "IRF7",
            "CXCL9",
            "CXCL10",
            "CXCL11",
            "CCL2",
            "CCL3",
            "CCL4",
            "TNFAIP3",
            "NFKBIA",
        ],
    ),
    (
        "stress_unfolded_protein_hypoxia",
        [
            "HSPA1A",
            "HSPA1B",
            "HSP90AA1",
            "DNAJB1",
            "JUN",
            "FOS",
            "ATF3",
            "DDIT3",
            "XBP1",
            "HSPA5",
            "HERPUD1",
            "ERN1",
            "ATF4",
            "EDEM1",
            "PDIA4",
            "HIF1A",
            "VEGFA",
            "LDHA",
            "SLC2A1",
        ],
    ),
    (
        "cell_cycle",
        [
            "MKI67",
            "TOP2A",
            "UBE2C",
            "BIRC5",
            "HMGB2",
            "CENPF",
            "NUSAP1",
            "TYMS",
            "PCNA",
            "MCM2",
            "MCM3",
            "MCM4",
            "MCM5",
            "MCM6",
            "MCM7",
            "CDK1",
            "CCNB1",
            "CCNB2",
        ],
    ),
]


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Create a biology-defined module-score embedding for human "
            "pancreas AnnData."
        )
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData path")
    parser.add_argument("--output", required=True, help="Output parquet path")
    parser.add_argument("--seed", type=int, default=0, help="Accepted run tag")
    return parser.parse_args()


def build_gene_index(var_names):
    exact = {}
    upper = {}
    for idx, symbol in enumerate(var_names):
        text = str(symbol)
        exact.setdefault(text, idx)
        upper.setdefault(text.upper(), idx)
    return exact, upper


def module_indices(genes, exact_index, upper_index):
    indices = []
    seen = set()
    for gene in genes:
        idx = exact_index.get(gene)
        if idx is None:
            idx = upper_index.get(gene.upper())
        if idx is not None and idx not in seen:
            seen.add(idx)
            indices.append(idx)
    return indices


def dense_vector(values):
    import numpy as np

    if hasattr(values, "A1"):
        return values.A1
    return np.asarray(values).reshape(-1)


def module_mean(matrix, indices):
    import numpy as np

    if not indices:
        return np.zeros(matrix.shape[0], dtype=np.float32)
    values = matrix[:, indices].mean(axis=1)
    return dense_vector(values).astype(np.float32, copy=False)


def robust_center_scale(values):
    import numpy as np

    vector = np.asarray(values, dtype=np.float32)
    finite = np.isfinite(vector)
    if not finite.any():
        return np.zeros_like(vector, dtype=np.float32)

    valid = vector[finite]
    median = np.median(valid)
    q25, q75 = np.percentile(valid, [25, 75])
    scale = (q75 - q25) / 1.349
    if scale < 1e-6:
        scale = np.std(valid)
    if scale < 1e-6:
        out = np.zeros_like(vector, dtype=np.float32)
    else:
        out = (vector - median) / scale
    out[~finite] = 0.0
    return np.clip(out, -8.0, 8.0).astype(np.float32)


def compute_embedding(adata):
    import numpy as np

    matrix = adata.X
    var_names = list(adata.var_names)
    exact_index, upper_index = build_gene_index(var_names)

    coords = []
    for _, genes in MODULES:
        indices = module_indices(genes, exact_index, upper_index)
        score = module_mean(matrix, indices)
        coords.append(robust_center_scale(score))

    return np.column_stack(coords).astype(np.float32)


def main():
    args = parse_args()

    import anndata as ad
    import pandas as pd

    adata = ad.read_h5ad(args.input)
    embedding = compute_embedding(adata)

    columns = {"cell_id": pd.Series(adata.obs_names.astype(str), dtype="string")}
    for dim in range(embedding.shape[1]):
        columns[f"dim_{dim}"] = embedding[:, dim].astype("float32", copy=False)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(columns).to_parquet(output, index=False)


if __name__ == "__main__":
    main()
