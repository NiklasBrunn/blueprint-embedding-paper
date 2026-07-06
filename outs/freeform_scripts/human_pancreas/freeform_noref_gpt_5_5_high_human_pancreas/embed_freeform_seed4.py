#!/usr/bin/env python3
"""Knowledge-based embedding for a human pancreas single-cell landscape.

The scenario is unconstrained human pancreas, so the coordinates below cover
major pancreatic epithelial, endocrine, stromal, vascular, neural, immune, and
common physiologic state programs using fixed gene signatures only.
"""

from __future__ import annotations

import argparse
from collections import OrderedDict
from pathlib import Path


SIGNATURES = OrderedDict(
    [
        (
            "epithelial_identity",
            [
                "EPCAM",
                "KRT8",
                "KRT18",
                "KRT19",
                "CDH1",
                "CLDN4",
                "MUC1",
                "TACSTD2",
            ],
        ),
        (
            "acinar_digestive_enzyme",
            [
                "PRSS1",
                "PRSS2",
                "PRSS3",
                "CPA1",
                "CPA2",
                "CPB1",
                "CTRB1",
                "CTRB2",
                "CELA2A",
                "CELA3A",
                "CELA3B",
                "PNLIP",
                "PNLIPRP1",
                "AMY2A",
                "AMY2B",
                "REG1A",
                "REG1B",
                "PLA2G1B",
                "CLPS",
                "CEL",
            ],
        ),
        (
            "acinar_secretory_granule",
            [
                "GP2",
                "SYCN",
                "RBPJL",
                "CUZD1",
                "SPINK1",
                "PDIA2",
                "CELA2B",
                "CTRL",
                "ZG16B",
                "KRTCAP3",
            ],
        ),
        (
            "ductal_ion_mucin",
            [
                "KRT19",
                "KRT7",
                "SOX9",
                "CFTR",
                "SLC4A4",
                "SLC26A3",
                "SLC26A6",
                "MUC1",
                "MUC6",
                "SPP1",
                "ANXA4",
                "TFF1",
                "TFF2",
                "AQP1",
                "LCN2",
                "KRT17",
            ],
        ),
        (
            "centroacinar_ductal",
            [
                "SOX9",
                "SPINK1",
                "KRT19",
                "MMP7",
                "SLC4A4",
                "TSPAN8",
                "CLDN10",
                "CEACAM6",
                "FXYD2",
            ],
        ),
        (
            "endocrine_pan_islet",
            [
                "CHGA",
                "CHGB",
                "SCG2",
                "SCG3",
                "PCSK1N",
                "PCSK1",
                "PCSK2",
                "NEUROD1",
                "PAX6",
                "ISL1",
                "SYP",
                "KCNJ11",
            ],
        ),
        (
            "beta_cell_insulin",
            [
                "INS",
                "IAPP",
                "MAFA",
                "PDX1",
                "NKX6-1",
                "SLC30A8",
                "PCSK1",
                "G6PC2",
                "UCN3",
                "ADCYAP1",
                "HADH",
                "DLK1",
                "RBP4",
                "SLC2A2",
            ],
        ),
        (
            "alpha_cell_glucagon",
            [
                "GCG",
                "ARX",
                "IRX2",
                "IRX1",
                "MAFB",
                "TTR",
                "LOXL4",
                "PCSK2",
                "GC",
                "TM4SF4",
                "SLC38A5",
            ],
        ),
        (
            "delta_cell_somatostatin",
            [
                "SST",
                "HHEX",
                "RBP4",
                "LEPR",
                "GHSR",
                "SLC17A6",
                "BCHE",
                "PCSK1",
                "PCSK2",
            ],
        ),
        (
            "pp_cell_pancreatic_polypeptide",
            [
                "PPY",
                "FAP",
                "SERTM1",
                "ID4",
                "PAX6",
                "ISL1",
                "PCSK1N",
                "CHGA",
                "CHGB",
            ],
        ),
        (
            "epsilon_cell_ghrelin",
            [
                "GHRL",
                "MLN",
                "ISL1",
                "PAX6",
                "SCG2",
                "CHGA",
                "CHGB",
                "PCSK1N",
            ],
        ),
        (
            "islet_maturation_function",
            [
                "MAFA",
                "UCN3",
                "SLC30A8",
                "G6PC2",
                "ABCC8",
                "KCNJ11",
                "PDX1",
                "NKX6-1",
                "GLIS3",
                "PAX6",
                "RFX6",
            ],
        ),
        (
            "fibroblast_stellate_quiescent",
            [
                "COL1A1",
                "COL1A2",
                "COL3A1",
                "DCN",
                "LUM",
                "PDGFRA",
                "COL6A1",
                "COL6A2",
                "COL6A3",
                "MGP",
                "RGS5",
                "RARRES2",
                "CFD",
            ],
        ),
        (
            "activated_stellate_myofibroblast",
            [
                "ACTA2",
                "TAGLN",
                "MYL9",
                "TPM2",
                "COL1A1",
                "COL1A2",
                "POSTN",
                "SPARC",
                "FN1",
                "THY1",
                "CTHRC1",
                "SERPINE1",
                "TIMP1",
            ],
        ),
        (
            "endothelial_blood_vessel",
            [
                "PECAM1",
                "VWF",
                "CDH5",
                "KDR",
                "FLT1",
                "ESAM",
                "EMCN",
                "PLVAP",
                "RAMP2",
                "CLEC14A",
                "ENG",
                "ERG",
            ],
        ),
        (
            "endothelial_arterial_capillary",
            [
                "GJA5",
                "HEY1",
                "EFNB2",
                "SOX17",
                "DLL4",
                "ACKR1",
                "RGCC",
                "CA4",
                "IFI27",
                "EDN1",
            ],
        ),
        (
            "lymphatic_endothelial",
            [
                "PROX1",
                "PDPN",
                "LYVE1",
                "FLT4",
                "CCL21",
                "MMRN1",
                "TFF3",
                "RELN",
                "TBX1",
            ],
        ),
        (
            "mural_pericyte_smooth_muscle",
            [
                "RGS5",
                "PDGFRB",
                "CSPG4",
                "MCAM",
                "NOTCH3",
                "ACTA2",
                "MYH11",
                "TAGLN",
                "MYLK",
                "CNN1",
                "DES",
            ],
        ),
        (
            "schwann_neural",
            [
                "SOX10",
                "S100B",
                "PLP1",
                "MPZ",
                "GFRA3",
                "NGFR",
                "ERBB3",
                "NRXN1",
                "SCN7A",
                "GAP43",
                "NRG1",
            ],
        ),
        (
            "immune_leukocyte",
            [
                "PTPRC",
                "LST1",
                "TYROBP",
                "HLA-DRA",
                "HLA-DRB1",
                "CXCR4",
                "CD52",
                "CORO1A",
                "ARHGDIB",
            ],
        ),
        (
            "t_cell_nk",
            [
                "CD3D",
                "CD3E",
                "CD3G",
                "TRAC",
                "TRBC1",
                "TRBC2",
                "IL7R",
                "CCR7",
                "NKG7",
                "GNLY",
                "PRF1",
                "GZMB",
                "KLRD1",
                "KLRB1",
            ],
        ),
        (
            "b_cell",
            [
                "MS4A1",
                "CD79A",
                "CD79B",
                "BANK1",
                "CD74",
                "HLA-DRA",
                "CD37",
                "VPREB3",
                "TCL1A",
                "PAX5",
            ],
        ),
        (
            "plasma_cell",
            [
                "MZB1",
                "XBP1",
                "JCHAIN",
                "SDC1",
                "DERL3",
                "PRDM1",
                "TNFRSF17",
                "IGHG1",
                "IGHG2",
                "IGHA1",
                "IGKC",
            ],
        ),
        (
            "monocyte_macrophage",
            [
                "LYZ",
                "LST1",
                "S100A8",
                "S100A9",
                "FCGR3A",
                "FCER1G",
                "C1QA",
                "C1QB",
                "C1QC",
                "MS4A7",
                "CD68",
                "MRC1",
                "MARCO",
            ],
        ),
        (
            "dendritic_antigen_presenting",
            [
                "FCER1A",
                "CLEC10A",
                "CD1C",
                "CLEC9A",
                "XCR1",
                "LILRA4",
                "IRF7",
                "HLA-DRA",
                "HLA-DPA1",
                "HLA-DPB1",
                "CCR7",
                "LAMP3",
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
                "HPGDS",
                "RGS13",
            ],
        ),
        (
            "erythroid_blood",
            [
                "HBB",
                "HBA1",
                "HBA2",
                "ALAS2",
                "SLC4A1",
                "GYPA",
                "AHSP",
                "KLF1",
            ],
        ),
        (
            "cell_cycle_g1s",
            [
                "MKI67",
                "TOP2A",
                "PCNA",
                "MCM2",
                "MCM3",
                "MCM4",
                "MCM5",
                "MCM6",
                "MCM7",
                "TYMS",
                "RRM2",
            ],
        ),
        (
            "cell_cycle_g2m",
            [
                "MKI67",
                "TOP2A",
                "CDK1",
                "CCNB1",
                "CCNB2",
                "BIRC5",
                "UBE2C",
                "AURKB",
                "NUSAP1",
                "CENPF",
                "HMGB2",
            ],
        ),
        (
            "stress_immediate_early",
            [
                "FOS",
                "FOSB",
                "JUN",
                "JUNB",
                "JUND",
                "ATF3",
                "IER2",
                "EGR1",
                "EGR2",
                "DUSP1",
                "BTG2",
            ],
        ),
        (
            "unfolded_protein_secretory_stress",
            [
                "HSPA5",
                "HSP90B1",
                "XBP1",
                "ATF4",
                "ATF6",
                "DDIT3",
                "HERPUD1",
                "DNAJB9",
                "PDIA4",
                "PDIA6",
                "ERP27",
            ],
        ),
        (
            "interferon_response",
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
                "CXCL10",
            ],
        ),
        (
            "hypoxia_angiogenic",
            [
                "VEGFA",
                "HILPDA",
                "NDRG1",
                "BNIP3",
                "SLC2A1",
                "ENO1",
                "LDHA",
                "PGK1",
                "CA9",
                "ANGPT2",
                "ADM",
            ],
        ),
        (
            "inflammatory_cytokine_nfkb",
            [
                "NFKBIA",
                "TNFAIP3",
                "CXCL1",
                "CXCL2",
                "CXCL3",
                "CXCL8",
                "IL6",
                "CCL2",
                "ICAM1",
                "PTGS2",
                "BIRC3",
            ],
        ),
    ]
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a fixed, knowledge-based embedding for human pancreas "
            "single-cell data from predefined gene signatures."
        )
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument(
        "--seed",
        required=True,
        type=int,
        help="Run tag supplied by the orchestrator; accepted for provenance.",
    )
    return parser.parse_args()


def _as_1d(array_like):
    import numpy as np

    return np.asarray(array_like).reshape(-1)


def _mean_expression(matrix, indices):
    import numpy as np
    from scipy import sparse

    if not indices:
        return np.zeros(matrix.shape[0], dtype=np.float32)
    sub = matrix[:, indices]
    if sparse.issparse(sub):
        return _as_1d(sub.mean(axis=1)).astype(np.float32)
    return np.asarray(sub, dtype=np.float32).mean(axis=1)


def _detection_fraction(matrix, indices):
    import numpy as np
    from scipy import sparse

    if not indices:
        return np.zeros(matrix.shape[0], dtype=np.float32)
    sub = matrix[:, indices]
    if sparse.issparse(sub):
        return (sub.getnnz(axis=1) / float(len(indices))).astype(np.float32)
    return (np.asarray(sub) > 0).mean(axis=1).astype(np.float32)


def _robust_scale(values):
    import numpy as np

    values = np.asarray(values, dtype=np.float32)
    finite = np.isfinite(values)
    if not finite.any():
        return np.zeros_like(values, dtype=np.float32)

    clean = values[finite]
    center = np.median(clean)
    mad = np.median(np.abs(clean - center)) * 1.4826
    if mad < 1e-6:
        mad = clean.std()
    if mad < 1e-6:
        out = np.zeros_like(values, dtype=np.float32)
        out[finite] = values[finite] - center
        return out.astype(np.float32)

    scaled = (values - center) / mad
    scaled[~finite] = 0.0
    return np.clip(scaled, -6.0, 6.0).astype(np.float32)


def _gene_index(var_names):
    index = {}
    for i, gene in enumerate(var_names):
        key = str(gene).upper()
        index.setdefault(key, []).append(i)
    return index


def _signature_indices(gene_to_indices, genes):
    indices = []
    seen = set()
    for gene in genes:
        for idx in gene_to_indices.get(gene.upper(), []):
            if idx not in seen:
                seen.add(idx)
                indices.append(idx)
    return indices


def build_embedding(input_path):
    import anndata as ad
    import numpy as np

    adata = ad.read_h5ad(input_path)
    matrix = adata.X
    gene_to_indices = _gene_index(adata.var_names)
    embedding = np.zeros((adata.n_obs, len(SIGNATURES)), dtype=np.float32)

    for dim, genes in enumerate(SIGNATURES.values()):
        indices = _signature_indices(gene_to_indices, genes)
        mean_expr = _mean_expression(matrix, indices)
        detection = _detection_fraction(matrix, indices)
        raw_score = mean_expr + 0.35 * detection
        embedding[:, dim] = _robust_scale(raw_score)

    cell_ids = np.asarray(adata.obs_names.astype(str))
    return cell_ids, embedding


def write_embedding(cell_ids, embedding, output_path):
    import pandas as pd

    data = {"cell_id": cell_ids.astype(str)}
    for i in range(embedding.shape[1]):
        data[f"dim_{i}"] = embedding[:, i].astype("float32", copy=False)
    frame = pd.DataFrame(data)
    frame.to_parquet(Path(output_path), index=False)


def main() -> None:
    args = parse_args()
    cell_ids, embedding = build_embedding(args.input)
    write_embedding(cell_ids, embedding, args.output)


if __name__ == "__main__":
    main()
