#!/usr/bin/env python3
"""Knowledge-defined embedding for the human pancreas cell landscape.

The scenario is unrestricted human pancreas under baseline/healthy
conditions. The dimensionality below is intentionally tied to the biology:
major epithelial compartments, endocrine hormone subtypes, stromal/vascular/
immune/neural lineages, and normal tissue state programs that commonly
separate cells in pancreas single-cell data.
"""

from __future__ import annotations

import argparse
from pathlib import Path


DIMENSIONS = [
    {
        "name": "acinar_zymogen_secretory",
        "genes": [
            "PRSS1",
            "PRSS2",
            "CPA1",
            "CPA2",
            "CPB1",
            "CTRB1",
            "CTRB2",
            "CELA3A",
            "CELA3B",
            "CLPS",
            "PNLIP",
            "GP2",
            "REG1A",
            "REG1B",
            "RBPJL",
            "AMY2A",
            "AMY2B",
        ],
    },
    {
        "name": "acinar_protective_regeneration",
        "genes": [
            "SPINK1",
            "REG1A",
            "REG1B",
            "REG3A",
            "LCN2",
            "AGR2",
            "TFF2",
            "MMP7",
            "KRT19",
            "SOX9",
            "SERPINA1",
        ],
    },
    {
        "name": "ductal_bicarbonate_mucin",
        "genes": [
            "KRT19",
            "KRT7",
            "KRT8",
            "KRT18",
            "EPCAM",
            "SOX9",
            "CFTR",
            "SLC4A4",
            "SLC26A3",
            "SLC26A6",
            "MUC1",
            "MUC6",
            "TFF1",
            "TFF2",
            "CLDN10",
        ],
    },
    {
        "name": "centroacinar_terminal_duct",
        "genes": [
            "CLDN10",
            "PROM1",
            "ALDH1A1",
            "SLC4A4",
            "CFTR",
            "SPP1",
            "CA2",
            "SOX9",
            "KRT19",
            "MMP7",
        ],
    },
    {
        "name": "pan_endocrine_neuroendocrine",
        "genes": [
            "CHGA",
            "CHGB",
            "SCG2",
            "SCG3",
            "PCSK1",
            "PCSK2",
            "PAX6",
            "ISL1",
            "NEUROD1",
            "INSM1",
            "NKX2-2",
            "SYP",
        ],
    },
    {
        "name": "beta_cell_insulin_mature",
        "genes": [
            "INS",
            "IAPP",
            "MAFA",
            "PDX1",
            "NKX6-1",
            "SLC2A2",
            "G6PC2",
            "SLC30A8",
            "PCSK1",
            "UCN3",
            "HADH",
            "RBP4",
        ],
    },
    {
        "name": "alpha_cell_glucagon",
        "genes": [
            "GCG",
            "ARX",
            "IRX2",
            "IRX1",
            "MAFB",
            "PCSK2",
            "TTR",
            "GC",
            "LOXL4",
            "CRYBA2",
        ],
    },
    {
        "name": "delta_cell_somatostatin",
        "genes": [
            "SST",
            "HHEX",
            "RBP4",
            "LEPR",
            "GHSR",
            "PCSK1",
            "PCSK2",
            "PAX6",
            "SEC11C",
        ],
    },
    {
        "name": "pp_gamma_cell",
        "genes": [
            "PPY",
            "PRLHR",
            "TSPAN8",
            "SLC30A8",
            "PCSK1",
            "PCSK2",
            "PAX6",
            "CHGA",
            "CHGB",
        ],
    },
    {
        "name": "epsilon_ghrelin_cell",
        "genes": [
            "GHRL",
            "MLN",
            "PCSK1N",
            "CHGA",
            "CHGB",
            "PAX6",
            "ISL1",
            "INSM1",
        ],
    },
    {
        "name": "endocrine_commitment_immature",
        "genes": [
            "NEUROG3",
            "SOX4",
            "SOX11",
            "MYT1",
            "PAX4",
            "HES6",
            "DLL1",
            "NKX6-1",
            "PDX1",
            "ISL1",
            "INSM1",
        ],
    },
    {
        "name": "beta_er_secretory_stress",
        "genes": [
            "HSPA5",
            "HSP90B1",
            "XBP1",
            "ERN1",
            "HERPUD1",
            "DNAJB9",
            "DDIT3",
            "ATF3",
            "ATF4",
            "PPP1R15A",
            "PDIA4",
        ],
    },
    {
        "name": "fibroblast_stellate_matrix",
        "genes": [
            "COL1A1",
            "COL1A2",
            "COL3A1",
            "COL6A1",
            "COL6A2",
            "DCN",
            "LUM",
            "FBLN1",
            "PDGFRA",
            "CXCL12",
            "MMP2",
            "THY1",
        ],
    },
    {
        "name": "myofibroblast_activated_stellate",
        "genes": [
            "ACTA2",
            "TAGLN",
            "MYL9",
            "TPM2",
            "POSTN",
            "FN1",
            "COL1A1",
            "COL1A2",
            "CTHRC1",
            "SPARC",
            "THY1",
        ],
    },
    {
        "name": "pericyte_smooth_muscle",
        "genes": [
            "RGS5",
            "PDGFRB",
            "CSPG4",
            "MCAM",
            "NOTCH3",
            "ACTA2",
            "MYH11",
            "TAGLN",
            "PLN",
            "DES",
        ],
    },
    {
        "name": "blood_endothelial",
        "genes": [
            "PECAM1",
            "VWF",
            "CDH5",
            "KDR",
            "FLT1",
            "EMCN",
            "ESAM",
            "RAMP2",
            "CLDN5",
            "ENG",
        ],
    },
    {
        "name": "arterial_endothelial",
        "genes": [
            "GJA5",
            "EFNB2",
            "SOX17",
            "HEY1",
            "DLL4",
            "NOTCH4",
            "CXCR4",
            "SEMA3G",
        ],
    },
    {
        "name": "venous_capillary_endothelial",
        "genes": [
            "ACKR1",
            "SELE",
            "VCAM1",
            "PLVAP",
            "RGCC",
            "CA4",
            "RAMP3",
            "VWF",
            "EMCN",
        ],
    },
    {
        "name": "lymphatic_endothelial",
        "genes": [
            "PROX1",
            "PDPN",
            "LYVE1",
            "FLT4",
            "CCL21",
            "MMRN1",
            "RELN",
            "TFF3",
        ],
    },
    {
        "name": "mesothelial_capsular",
        "genes": [
            "MSLN",
            "UPK3B",
            "WT1",
            "KRT8",
            "KRT18",
            "PDPN",
            "LUM",
            "C3",
            "ITLN1",
        ],
    },
    {
        "name": "schwann_neural",
        "genes": [
            "SOX10",
            "S100B",
            "PLP1",
            "MPZ",
            "MBP",
            "GFRA3",
            "CHL1",
            "NRXN1",
            "SNAP25",
            "GAP43",
        ],
    },
    {
        "name": "macrophage_monocyte_resident",
        "genes": [
            "LST1",
            "TYROBP",
            "AIF1",
            "FCGR3A",
            "CSF1R",
            "CD68",
            "CD163",
            "MRC1",
            "C1QA",
            "C1QB",
            "C1QC",
            "MS4A7",
        ],
    },
    {
        "name": "dendritic_antigen_presenting",
        "genes": [
            "FCER1A",
            "CD1C",
            "CLEC10A",
            "CLEC9A",
            "XCR1",
            "BATF3",
            "LILRA4",
            "IL3RA",
            "HLA-DRA",
            "HLA-DPB1",
        ],
    },
    {
        "name": "t_nk_cytotoxic",
        "genes": [
            "CD3D",
            "CD3E",
            "TRAC",
            "IL7R",
            "CCR7",
            "NKG7",
            "GNLY",
            "PRF1",
            "KLRD1",
            "GZMB",
        ],
    },
    {
        "name": "b_plasma_humoral",
        "genes": [
            "MS4A1",
            "CD79A",
            "CD79B",
            "BANK1",
            "CD74",
            "MZB1",
            "JCHAIN",
            "IGKC",
            "IGHG1",
            "XBP1",
        ],
    },
    {
        "name": "mast_cell",
        "genes": [
            "TPSAB1",
            "TPSB2",
            "CPA3",
            "KIT",
            "MS4A2",
            "HDC",
            "HPGDS",
            "GATA2",
        ],
    },
    {
        "name": "cell_cycle_proliferation",
        "genes": [
            "MKI67",
            "TOP2A",
            "PCNA",
            "TYMS",
            "CENPF",
            "UBE2C",
            "BIRC5",
            "HMGB2",
            "STMN1",
            "MCM5",
        ],
    },
    {
        "name": "interferon_antiviral",
        "genes": [
            "ISG15",
            "IFIT1",
            "IFIT2",
            "IFIT3",
            "MX1",
            "OAS1",
            "OAS2",
            "STAT1",
            "IRF7",
            "IFI6",
            "IFI27",
        ],
    },
    {
        "name": "hypoxia_glycolytic_stress",
        "genes": [
            "VEGFA",
            "CA9",
            "SLC2A1",
            "LDHA",
            "BNIP3",
            "PGK1",
            "ENO1",
            "ALDOA",
            "NDRG1",
        ],
    },
    {
        "name": "epithelial_antigen_interferon_response",
        "genes": [
            "HLA-A",
            "HLA-B",
            "HLA-C",
            "B2M",
            "TAP1",
            "TAP2",
            "NLRC5",
            "CXCL10",
            "CXCL11",
            "STAT1",
        ],
    },
]


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Create a fixed, biologically defined embedding for human pancreas "
            "AnnData input."
        )
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file")
    parser.add_argument("--output", required=True, help="Output parquet path")
    parser.add_argument(
        "--seed",
        type=int,
        required=True,
        help="Run tag supplied by the orchestrator; not used by this deterministic scorer",
    )
    return parser


def make_gene_index(var_names) -> dict[str, int]:
    """Map case-normalized gene symbols to their first matrix column index."""
    gene_to_idx: dict[str, int] = {}
    for idx, gene in enumerate(var_names):
        key = str(gene).strip().upper()
        if key and key not in gene_to_idx:
            gene_to_idx[key] = idx
    return gene_to_idx


def dense_subset(matrix, indices):
    import numpy as np
    from scipy import sparse

    subset = matrix[:, indices]
    if sparse.issparse(subset):
        subset = subset.toarray()
    else:
        subset = np.asarray(subset)
    if subset.ndim == 1:
        subset = subset.reshape(-1, 1)
    return subset.astype(np.float32, copy=False)


def robust_scale(values):
    import numpy as np

    values = np.asarray(values, dtype=np.float32)
    if values.size == 0:
        return values

    median = np.nanmedian(values)
    q25, q75 = np.nanpercentile(values, [25.0, 75.0])
    scale = q75 - q25
    if not np.isfinite(scale) or scale < 1e-6:
        scale = np.nanstd(values)
    if not np.isfinite(scale) or scale < 1e-6:
        return np.zeros(values.shape, dtype=np.float32)

    scaled = (values - median) / scale
    scaled = np.clip(scaled, -8.0, 8.0)
    return np.nan_to_num(scaled, nan=0.0, posinf=8.0, neginf=-8.0).astype(
        np.float32
    )


def score_signature(matrix, gene_to_idx, genes):
    import numpy as np

    present_indices = []
    seen = set()
    for gene in genes:
        key = gene.upper()
        idx = gene_to_idx.get(key)
        if idx is not None and idx not in seen:
            present_indices.append(idx)
            seen.add(idx)

    n_cells = matrix.shape[0]
    if not present_indices:
        return np.zeros(n_cells, dtype=np.float32)

    expression = dense_subset(matrix, present_indices)
    raw_score = np.asarray(expression.mean(axis=1)).reshape(-1)
    return robust_scale(raw_score)


def compute_embedding(adata):
    import numpy as np

    matrix = adata.X
    gene_to_idx = make_gene_index(adata.var_names)
    columns = []
    for dimension in DIMENSIONS:
        columns.append(score_signature(matrix, gene_to_idx, dimension["genes"]))
    return np.column_stack(columns).astype(np.float32, copy=False)


def run(input_path: str, output_path: str) -> None:
    import anndata as ad
    import pandas as pd

    adata = ad.read_h5ad(input_path)
    embedding = compute_embedding(adata)

    data = {"cell_id": [str(cell_id) for cell_id in adata.obs_names]}
    for dim_idx in range(embedding.shape[1]):
        data[f"dim_{dim_idx}"] = embedding[:, dim_idx].astype("float32", copy=False)

    out = Path(output_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(data).to_parquet(out, index=False)


def main() -> None:
    args = build_arg_parser().parse_args()
    _ = args.seed
    run(args.input, args.output)


if __name__ == "__main__":
    main()
