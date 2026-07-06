#!/usr/bin/env python3
"""Literature-marker embedding for human pancreas scRNA-seq cells.

This intentionally avoids model fitting, clustering, PCA, pretrained models, or
data-derived gene selection. Coordinates are fixed marker-program summaries and
fixed biological contrasts for human pancreas cell identities and states.
"""

from __future__ import annotations

import argparse
import re
import sys
from collections import OrderedDict
from typing import Iterable

import numpy as np
import pandas as pd


MARKER_PANELS: "OrderedDict[str, tuple[str, ...]]" = OrderedDict(
    [
        (
            "pancreas_epithelial",
            (
                "KRT8",
                "KRT18",
                "KRT19",
                "EPCAM",
                "CDH1",
                "CLDN4",
                "TACSTD2",
            ),
        ),
        (
            "endocrine_core",
            (
                "CHGA",
                "CHGB",
                "SCG2",
                "SCG3",
                "PCSK1N",
                "PCSK1",
                "PCSK2",
                "INSM1",
                "ISL1",
                "NEUROD1",
                "PAX6",
            ),
        ),
        (
            "beta_cell",
            (
                "INS",
                "IAPP",
                "SLC30A8",
                "G6PC2",
                "PDX1",
                "NKX6-1",
                "MAFA",
                "PCSK1",
                "DLK1",
            ),
        ),
        (
            "beta_mature",
            (
                "UCN3",
                "MAFA",
                "SLC30A8",
                "G6PC2",
                "ABCC8",
                "KCNJ11",
                "PCSK1",
            ),
        ),
        (
            "beta_functional_glucose_sensing",
            (
                "GCK",
                "SLC2A2",
                "SLC30A8",
                "ABCC8",
                "KCNJ11",
                "PCSK1",
                "PCSK2",
            ),
        ),
        (
            "alpha_cell",
            (
                "GCG",
                "ARX",
                "IRX1",
                "IRX2",
                "TTR",
                "LOXL4",
                "PCSK2",
                "MAFB",
            ),
        ),
        (
            "delta_cell",
            (
                "SST",
                "HHEX",
                "RBP4",
                "LEPR",
                "GHSR",
            ),
        ),
        (
            "gamma_pp_cell",
            (
                "PPY",
                "SLC30A8",
                "PAX6",
                "ISL1",
            ),
        ),
        (
            "epsilon_cell",
            (
                "GHRL",
                "MLN",
                "PAX6",
                "ISL1",
            ),
        ),
        (
            "endocrine_progenitor_like",
            (
                "NEUROG3",
                "SOX4",
                "PAX4",
                "NKX2-2",
                "NKX6-1",
                "RFX6",
                "FEV",
                "ONECUT1",
            ),
        ),
        (
            "acinar_cell",
            (
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
                "REG1A",
                "REG1B",
                "PNLIP",
                "CLPS",
                "AMY2A",
                "AMY2B",
            ),
        ),
        (
            "acinar_zymogen_secretory",
            (
                "PRSS1",
                "PRSS2",
                "CPA1",
                "CTRB1",
                "CTRB2",
                "CELA3A",
                "GP2",
                "SYCN",
                "CEL",
                "PNLIP",
            ),
        ),
        (
            "ductal_cell",
            (
                "KRT19",
                "KRT7",
                "KRT8",
                "KRT18",
                "SOX9",
                "CFTR",
                "MUC1",
                "MUC6",
                "SLC4A4",
                "ANXA4",
                "SPP1",
                "TACSTD2",
            ),
        ),
        (
            "centroacinar_terminal_duct",
            (
                "SOX9",
                "KRT19",
                "CFTR",
                "SLC4A4",
                "CLDN10",
                "PROM1",
                "ALDH1A1",
            ),
        ),
        (
            "tuft_like_ductal",
            (
                "POU2F3",
                "TRPM5",
                "DCLK1",
                "AVIL",
                "SH2D6",
                "GFI1B",
            ),
        ),
        (
            "mesenchymal_stellate",
            (
                "COL1A1",
                "COL1A2",
                "COL3A1",
                "DCN",
                "LUM",
                "PDGFRB",
                "RGS5",
                "ACTA2",
                "TAGLN",
                "DES",
            ),
        ),
        (
            "fibroblast_matrix",
            (
                "COL1A1",
                "COL1A2",
                "COL6A1",
                "COL6A2",
                "DCN",
                "LUM",
                "FBLN1",
                "MMP2",
            ),
        ),
        (
            "activated_stellate_myofibroblast",
            (
                "ACTA2",
                "TAGLN",
                "MYL9",
                "COL1A1",
                "COL1A2",
                "POSTN",
                "TIMP1",
                "SPARC",
            ),
        ),
        (
            "pericyte_smooth_muscle",
            (
                "RGS5",
                "PDGFRB",
                "CSPG4",
                "MCAM",
                "NOTCH3",
                "ACTA2",
                "MYH11",
                "TAGLN",
            ),
        ),
        (
            "blood_endothelial",
            (
                "PECAM1",
                "VWF",
                "KDR",
                "FLT1",
                "CDH5",
                "ESAM",
                "ENG",
                "PLVAP",
            ),
        ),
        (
            "lymphatic_endothelial",
            (
                "PROX1",
                "LYVE1",
                "PDPN",
                "FLT4",
                "CCL21",
                "MMRN1",
            ),
        ),
        (
            "immune_pan_leukocyte",
            (
                "PTPRC",
                "LST1",
                "TYROBP",
                "FCER1G",
                "CD74",
                "HLA-DRA",
            ),
        ),
        (
            "macrophage_monocyte",
            (
                "LST1",
                "TYROBP",
                "FCER1G",
                "AIF1",
                "C1QA",
                "C1QB",
                "C1QC",
                "CD68",
                "CSF1R",
                "LYZ",
            ),
        ),
        (
            "dendritic_antigen_presenting",
            (
                "FCER1A",
                "CLEC10A",
                "CD1C",
                "ITGAX",
                "LILRA4",
                "CLEC4C",
                "HLA-DRA",
                "HLA-DPA1",
                "HLA-DPB1",
            ),
        ),
        (
            "t_nk_cell",
            (
                "CD3D",
                "CD3E",
                "TRAC",
                "IL7R",
                "NKG7",
                "GNLY",
                "KLRD1",
                "PRF1",
            ),
        ),
        (
            "b_plasma_cell",
            (
                "MS4A1",
                "CD79A",
                "CD79B",
                "BANK1",
                "MZB1",
                "JCHAIN",
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
            ),
        ),
        (
            "schwann_neural",
            (
                "S100B",
                "SOX10",
                "PLP1",
                "MPZ",
                "GFRA3",
                "NRXN1",
            ),
        ),
        (
            "erythroid",
            (
                "HBB",
                "HBA1",
                "HBA2",
                "ALAS2",
                "AHSP",
            ),
        ),
        (
            "cycling",
            (
                "MKI67",
                "TOP2A",
                "UBE2C",
                "BIRC5",
                "PCNA",
                "MCM2",
                "MCM5",
                "TYMS",
            ),
        ),
        (
            "stress_unfolded_protein_response",
            (
                "HSPA5",
                "HSP90B1",
                "XBP1",
                "DDIT3",
                "ATF4",
                "HERPUD1",
                "DNAJB9",
            ),
        ),
        (
            "interferon_response",
            (
                "ISG15",
                "IFIT1",
                "IFIT2",
                "IFIT3",
                "MX1",
                "OAS1",
                "STAT1",
                "IRF7",
            ),
        ),
        (
            "hypoxia_glycolysis",
            (
                "HIF1A",
                "VEGFA",
                "SLC2A1",
                "LDHA",
                "ENO1",
                "PGK1",
                "CA9",
            ),
        ),
        (
            "apoptosis_injury",
            (
                "BAX",
                "BBC3",
                "PMAIP1",
                "CASP3",
                "GADD45A",
                "FOS",
                "JUN",
            ),
        ),
        (
            "epithelial_mesenchymal_transition",
            (
                "VIM",
                "FN1",
                "ITGA5",
                "ITGB1",
                "SNAI2",
                "ZEB1",
                "SERPINE1",
            ),
        ),
        (
            "mhc_class_ii",
            (
                "HLA-DRA",
                "HLA-DRB1",
                "HLA-DPA1",
                "HLA-DPB1",
                "HLA-DQA1",
                "HLA-DQB1",
                "CD74",
            ),
        ),
    ]
)


CONTRASTS: tuple[tuple[str, str, str], ...] = (
    ("endocrine_vs_exocrine", "endocrine_core", "acinar_cell"),
    ("endocrine_vs_ductal", "endocrine_core", "ductal_cell"),
    ("acinar_vs_ductal", "acinar_cell", "ductal_cell"),
    ("epithelial_vs_stromal", "pancreas_epithelial", "mesenchymal_stellate"),
    ("epithelial_vs_immune", "pancreas_epithelial", "immune_pan_leukocyte"),
    ("beta_vs_alpha", "beta_cell", "alpha_cell"),
    ("beta_vs_delta", "beta_cell", "delta_cell"),
    ("alpha_vs_delta", "alpha_cell", "delta_cell"),
    ("pp_vs_other_endocrine", "gamma_pp_cell", "endocrine_core"),
    ("epsilon_vs_other_endocrine", "epsilon_cell", "endocrine_core"),
    ("beta_maturity_vs_progenitor", "beta_mature", "endocrine_progenitor_like"),
    ("beta_glucose_sensing_vs_stress", "beta_functional_glucose_sensing", "stress_unfolded_protein_response"),
    ("zymogen_vs_ductal_secretory", "acinar_zymogen_secretory", "ductal_cell"),
    ("centroacinar_vs_acinar", "centroacinar_terminal_duct", "acinar_cell"),
    ("tuft_ductal_vs_canonical_ductal", "tuft_like_ductal", "ductal_cell"),
    ("stellate_vs_endothelial", "mesenchymal_stellate", "blood_endothelial"),
    ("myofibroblast_vs_quiescent_matrix", "activated_stellate_myofibroblast", "fibroblast_matrix"),
    ("pericyte_vs_fibroblast", "pericyte_smooth_muscle", "fibroblast_matrix"),
    ("blood_vs_lymphatic_endothelium", "blood_endothelial", "lymphatic_endothelial"),
    ("myeloid_vs_lymphoid", "macrophage_monocyte", "t_nk_cell"),
    ("antigen_presentation_vs_pan_immune", "mhc_class_ii", "immune_pan_leukocyte"),
    ("b_plasma_vs_t_nk", "b_plasma_cell", "t_nk_cell"),
    ("mast_vs_myeloid", "mast_cell", "macrophage_monocyte"),
    ("injury_stress_vs_endocrine", "stress_unfolded_protein_response", "endocrine_core"),
    ("interferon_vs_mhc_ii", "interferon_response", "mhc_class_ii"),
    ("emt_vs_epithelial", "epithelial_mesenchymal_transition", "pancreas_epithelial"),
)


ALIASES: dict[str, tuple[str, ...]] = {
    "NKX6-1": ("NKX6-1", "NKX6_1", "NKX61"),
    "NKX2-2": ("NKX2-2", "NKX2_2", "NKX22"),
    "HLA-DRA": ("HLA-DRA", "HLADRA"),
    "HLA-DRB1": ("HLA-DRB1", "HLADRB1"),
    "HLA-DPA1": ("HLA-DPA1", "HLADPA1"),
    "HLA-DPB1": ("HLA-DPB1", "HLADPB1"),
    "HLA-DQA1": ("HLA-DQA1", "HLADQA1"),
    "HLA-DQB1": ("HLA-DQB1", "HLADQB1"),
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a fixed literature-marker embedding for human pancreas scRNA-seq cells."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad file")
    parser.add_argument("--output", required=True, help="Output parquet path")
    parser.add_argument("--seed", required=True, type=int, help="Run tag; accepted for interface compatibility")
    return parser.parse_args()


def compact_symbol(symbol: str) -> str:
    return re.sub(r"[^A-Z0-9]", "", symbol.upper())


def build_gene_lookup(var_names: Iterable[object]) -> tuple[dict[str, int], dict[str, int]]:
    exact: dict[str, int] = {}
    compact: dict[str, int] = {}
    for idx, raw in enumerate(var_names):
        name = str(raw).strip()
        upper = name.upper()
        exact.setdefault(upper, idx)
        compact.setdefault(compact_symbol(upper), idx)
    return exact, compact


def resolve_genes(
    genes: Iterable[str], exact_lookup: dict[str, int], compact_lookup: dict[str, int]
) -> list[int]:
    indices: list[int] = []
    seen: set[int] = set()
    for gene in genes:
        candidates = ALIASES.get(gene, (gene,))
        for candidate in candidates:
            idx = exact_lookup.get(candidate.upper())
            if idx is None:
                idx = compact_lookup.get(compact_symbol(candidate))
            if idx is not None and idx not in seen:
                indices.append(idx)
                seen.add(idx)
                break
    return indices


def row_sums(matrix) -> np.ndarray:
    values = np.asarray(matrix.sum(axis=1)).reshape(-1)
    return values.astype(np.float64, copy=False)


def marker_score(counts, gene_indices: list[int], library_size: np.ndarray) -> np.ndarray:
    if not gene_indices:
        return np.zeros(library_size.shape[0], dtype=np.float32)
    marker_counts = row_sums(counts[:, gene_indices])
    safe_library = np.maximum(library_size, 1.0)
    mean_cptt = 10000.0 * marker_counts / safe_library / float(len(gene_indices))
    return np.log1p(mean_cptt).astype(np.float32, copy=False)


def main() -> None:
    args = parse_args()

    try:
        import anndata as ad
    except ImportError as exc:
        raise SystemExit("embed_freeform.py requires anndata to read .h5ad input") from exc

    adata = ad.read_h5ad(args.input)
    if "counts" not in adata.layers:
        raise SystemExit('Input AnnData must contain raw counts in layers["counts"]')

    counts = adata.layers["counts"]
    n_cells = adata.n_obs
    library_size = row_sums(counts)

    exact_lookup, compact_lookup = build_gene_lookup(adata.var_names)
    panel_indices = {
        panel_name: resolve_genes(genes, exact_lookup, compact_lookup)
        for panel_name, genes in MARKER_PANELS.items()
    }

    matched_genes = sum(len(indices) for indices in panel_indices.values())
    if matched_genes == 0:
        print(
            "Warning: no marker genes matched var_names; writing all-zero coordinates.",
            file=sys.stderr,
        )

    panel_scores: dict[str, np.ndarray] = {}
    for panel_name, indices in panel_indices.items():
        panel_scores[panel_name] = marker_score(counts, indices, library_size)

    feature_arrays: list[np.ndarray] = []
    for panel_name in MARKER_PANELS:
        feature_arrays.append(panel_scores[panel_name])

    for _contrast_name, positive_panel, negative_panel in CONTRASTS:
        feature_arrays.append(panel_scores[positive_panel] - panel_scores[negative_panel])

    embedding = np.column_stack(feature_arrays).astype(np.float32, copy=False)
    output = pd.DataFrame({"cell_id": adata.obs_names.astype(str).to_numpy()})
    for dim in range(embedding.shape[1]):
        output[f"dim_{dim}"] = embedding[:, dim]

    output.to_parquet(args.output, index=False)
    print(
        f"Wrote {n_cells} cells x {embedding.shape[1]} literature-marker dimensions to {args.output}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
