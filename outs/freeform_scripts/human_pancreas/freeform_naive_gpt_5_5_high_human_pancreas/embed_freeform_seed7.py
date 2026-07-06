#!/usr/bin/env python3
"""Literature-marker embedding for human pancreas scRNA-seq cells."""

from __future__ import annotations

import argparse
import sys
from collections import OrderedDict

import numpy as np
import pandas as pd


def _require_runtime_packages():
    try:
        import anndata as ad  # noqa: F401
        from scipy import sparse  # noqa: F401
    except Exception as exc:  # pragma: no cover - depends on runtime env
        raise SystemExit(
            "embed_freeform.py requires anndata and scipy to read .h5ad inputs."
        ) from exc


MARKER_MODULES = OrderedDict(
    [
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
                "AMY2A",
                "AMY2B",
                "PNLIP",
                "PNLIPRP1",
                "CLPS",
                "REG1A",
                "REG1B",
                "REG3A",
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
                "MUC1",
                "CFTR",
                "SOX9",
                "HNF1B",
                "ONECUT1",
                "CLDN4",
                "CLDN10",
                "ANXA4",
                "SLC4A4",
                "AQP1",
                "SLC26A6",
            ],
        ),
        (
            "ductal_activated_spp1",
            ["SPP1", "TFF1", "TFF2", "MMP7", "MUC5B", "LCN2", "KRT17", "CEACAM6"],
        ),
        (
            "pan_endocrine",
            [
                "CHGA",
                "CHGB",
                "SCG2",
                "SCG3",
                "SCG5",
                "SCGN",
                "PCSK1N",
                "PAX6",
                "ISL1",
                "NEUROD1",
                "RFX6",
                "SLC18A2",
            ],
        ),
        (
            "beta_cell",
            [
                "INS",
                "IAPP",
                "SLC30A8",
                "PCSK1",
                "PCSK2",
                "PDX1",
                "NKX6-1",
                "G6PC2",
                "RBP4",
                "DLK1",
                "HADH",
            ],
        ),
        (
            "beta_mature",
            ["MAFA", "UCN3", "SLC30A8", "G6PC2", "PDX1", "NKX6-1", "ABCC8", "KCNJ11"],
        ),
        (
            "alpha_cell",
            ["GCG", "ARX", "IRX1", "IRX2", "MAFB", "TTR", "GC", "LOXL4"],
        ),
        (
            "delta_cell",
            ["SST", "HHEX", "RBP4", "LEPR", "BCHE", "LY6H"],
        ),
        ("gamma_pp_cell", ["PPY", "PNMAL1", "SERTM1", "SLITRK6"]),
        ("epsilon_cell", ["GHRL", "MLN", "IAPP", "RFX6"]),
        (
            "enterochromaffin_like",
            ["TPH1", "DDC", "SLC18A1", "TAC1", "LMX1A", "CHGA", "CHGB"],
        ),
        (
            "endocrine_progenitor_like",
            ["NEUROG3", "SOX9", "PAX6", "RFX6", "INSM1", "NKX2-2", "PAX4", "MNX1"],
        ),
        (
            "fibroblast_matrix",
            [
                "COL1A1",
                "COL1A2",
                "COL3A1",
                "COL5A1",
                "COL6A1",
                "COL6A2",
                "COL6A3",
                "DCN",
                "LUM",
                "FBLN1",
                "MMP2",
                "FN1",
                "THY1",
                "PDGFRA",
            ],
        ),
        (
            "stellate_myofibroblast",
            ["ACTA2", "TAGLN", "MYL9", "COL1A1", "COL3A1", "RGS5", "PDGFRB", "MMP11"],
        ),
        (
            "endothelial",
            ["PECAM1", "VWF", "KDR", "FLT1", "CDH5", "ESAM", "EMCN", "CLDN5", "RAMP2"],
        ),
        (
            "capillary_venous_endothelium",
            ["PLVAP", "ACKR1", "SELE", "ICAM1", "VCAM1", "RGCC", "RAMP3"],
        ),
        (
            "pericyte_smooth_muscle",
            ["RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3", "ACTA2", "MYH11", "CNN1"],
        ),
        (
            "myeloid_macrophage",
            [
                "LYZ",
                "LST1",
                "AIF1",
                "C1QA",
                "C1QB",
                "C1QC",
                "CD14",
                "FCER1G",
                "MS4A7",
                "TYROBP",
                "HLA-DRA",
                "HLA-DRB1",
            ],
        ),
        ("dendritic_cell", ["FCER1A", "CLEC10A", "CD1C", "LILRA4", "CLEC4C", "IRF7"]),
        ("t_cell", ["CD3D", "CD3E", "CD3G", "TRAC", "TRBC1", "TRBC2", "IL7R", "CCR7"]),
        ("nk_cytotoxic", ["NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "KLRD1", "KLRF1"]),
        ("b_cell", ["MS4A1", "CD79A", "CD79B", "BANK1", "CD74", "HLA-DRA"]),
        ("plasma_cell", ["MZB1", "JCHAIN", "XBP1", "SDC1", "DERL3", "IGKC", "IGHG1"]),
        ("mast_cell", ["TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC"]),
        ("schwann_neural", ["SOX10", "PLP1", "MPZ", "S100B", "GFRA3", "NRXN1", "NCAM1"]),
        ("erythroid", ["HBB", "HBA1", "HBA2", "ALAS2", "SLC4A1", "GYPA"]),
        (
            "interferon_inflammatory",
            ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OAS2", "STAT1", "CXCL10"],
        ),
        (
            "stress_unfolded_protein",
            ["HSPA5", "HSP90B1", "XBP1", "ATF4", "DDIT3", "HERPUD1", "DNAJB9", "ERN1"],
        ),
        ("immediate_early", ["FOS", "JUN", "JUNB", "EGR1", "IER2", "IER3", "BTG2"]),
        (
            "cell_cycle",
            ["MKI67", "TOP2A", "CENPF", "PCNA", "MCM2", "MCM5", "UBE2C", "CCNB1"],
        ),
    ]
)


CONTRASTS = [
    ("endocrine_vs_exocrine", ["pan_endocrine"], ["acinar_digestive_enzyme", "ductal_epithelial"]),
    ("beta_vs_other_islet", ["beta_cell"], ["alpha_cell", "delta_cell", "gamma_pp_cell", "epsilon_cell"]),
    ("alpha_vs_beta", ["alpha_cell"], ["beta_cell"]),
    ("delta_vs_beta_alpha", ["delta_cell"], ["beta_cell", "alpha_cell"]),
    ("acinar_vs_ductal", ["acinar_digestive_enzyme"], ["ductal_epithelial", "ductal_activated_spp1"]),
    ("ductal_activation_vs_homeostatic", ["ductal_activated_spp1"], ["ductal_epithelial"]),
    ("epithelial_vs_stromal", ["acinar_digestive_enzyme", "ductal_epithelial", "pan_endocrine"], ["fibroblast_matrix", "endothelial", "pericyte_smooth_muscle"]),
    ("vascular_vs_matrix", ["endothelial", "capillary_venous_endothelium"], ["fibroblast_matrix", "stellate_myofibroblast"]),
    ("mural_vs_endothelial", ["pericyte_smooth_muscle"], ["endothelial"]),
    ("immune_vs_parenchymal", ["myeloid_macrophage", "t_cell", "b_cell", "nk_cytotoxic"], ["acinar_digestive_enzyme", "ductal_epithelial", "pan_endocrine"]),
    ("lymphoid_vs_myeloid", ["t_cell", "b_cell", "nk_cytotoxic"], ["myeloid_macrophage", "dendritic_cell"]),
    ("plasma_vs_b_cell", ["plasma_cell"], ["b_cell"]),
    ("cytotoxic_vs_t_cell", ["nk_cytotoxic"], ["t_cell"]),
    ("beta_maturity_vs_progenitor", ["beta_mature"], ["endocrine_progenitor_like"]),
    ("secretory_granule_vs_ductal", ["pan_endocrine"], ["ductal_epithelial"]),
    ("stress_vs_identity", ["stress_unfolded_protein", "immediate_early"], ["acinar_digestive_enzyme", "ductal_epithelial", "pan_endocrine"]),
]


def parse_args(argv: list[str]) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Embed human pancreas scRNA-seq cells with fixed literature marker modules."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file")
    parser.add_argument("--output", required=True, help="Output parquet path")
    parser.add_argument("--seed", required=True, type=int, help="Run tag; accepted for interface compatibility")
    return parser.parse_args(argv)


def _gene_lookup(var_names) -> dict[str, int]:
    lookup: dict[str, int] = {}
    for idx, name in enumerate(var_names):
        key = str(name).strip().upper()
        if key and key not in lookup:
            lookup[key] = idx
    return lookup


def _selected_log_cpm(adata, selected_indices: list[int]):
    from scipy import sparse

    counts = adata.layers["counts"]
    if counts.shape != adata.shape:
        raise ValueError(
            f'layers["counts"] shape {counts.shape} does not match AnnData shape {adata.shape}'
        )

    if sparse.issparse(counts):
        matrix = counts.tocsr().astype(np.float32, copy=False)
        library = np.asarray(matrix.sum(axis=1)).ravel().astype(np.float32)
        scale = np.divide(
            10000.0,
            library,
            out=np.zeros_like(library, dtype=np.float32),
            where=library > 0,
        )
        selected = matrix[:, selected_indices].astype(np.float32, copy=True)
        selected = selected.multiply(scale[:, None]).tocsr()
        np.log1p(selected.data, out=selected.data)
        return selected

    matrix = np.asarray(counts, dtype=np.float32)
    library = matrix.sum(axis=1, dtype=np.float32)
    scale = np.divide(
        10000.0,
        library,
        out=np.zeros_like(library, dtype=np.float32),
        where=library > 0,
    )
    return np.log1p(matrix[:, selected_indices] * scale[:, None]).astype(np.float32, copy=False)


def _column_mean(expr, columns: list[int], n_cells: int) -> np.ndarray:
    if not columns:
        return np.zeros(n_cells, dtype=np.float32)
    values = expr[:, columns]
    mean = np.asarray(values.mean(axis=1)).ravel()
    return mean.astype(np.float32, copy=False)


def _module_average(module_scores: dict[str, np.ndarray], names: list[str], n_cells: int) -> np.ndarray:
    present = [module_scores[name] for name in names if name in module_scores]
    if not present:
        return np.zeros(n_cells, dtype=np.float32)
    return np.mean(np.vstack(present), axis=0).astype(np.float32, copy=False)


def build_embedding(adata) -> np.ndarray:
    n_cells = adata.n_obs
    lookup = _gene_lookup(adata.var_names)

    ordered_genes: list[str] = []
    seen: set[str] = set()
    for genes in MARKER_MODULES.values():
        for gene in genes:
            key = gene.upper()
            if key in lookup and key not in seen:
                ordered_genes.append(key)
                seen.add(key)

    if not ordered_genes:
        raise ValueError("None of the configured human pancreas marker genes were found in var_names.")

    selected_indices = [lookup[gene] for gene in ordered_genes]
    selected_position = {gene: idx for idx, gene in enumerate(ordered_genes)}
    log_cpm = _selected_log_cpm(adata, selected_indices)

    module_scores: OrderedDict[str, np.ndarray] = OrderedDict()
    for module_name, genes in MARKER_MODULES.items():
        columns = [selected_position[gene.upper()] for gene in genes if gene.upper() in selected_position]
        module_scores[module_name] = _column_mean(log_cpm, columns, n_cells)

    dims: list[np.ndarray] = list(module_scores.values())
    for _contrast_name, positive, negative in CONTRASTS:
        pos = _module_average(module_scores, positive, n_cells)
        neg = _module_average(module_scores, negative, n_cells)
        dims.append((pos - neg).astype(np.float32, copy=False))

    return np.column_stack(dims).astype(np.float32, copy=False)


def main(argv: list[str]) -> int:
    args = parse_args(argv)
    _require_runtime_packages()

    import anndata as ad

    adata = ad.read_h5ad(args.input)
    if "counts" not in adata.layers:
        raise ValueError('Input AnnData is missing layers["counts"].')

    embedding = build_embedding(adata)
    output = pd.DataFrame({"cell_id": adata.obs_names.astype(str).to_numpy()})
    for idx in range(embedding.shape[1]):
        output[f"dim_{idx}"] = embedding[:, idx]
    output.to_parquet(args.output, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
