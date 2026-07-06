#!/usr/bin/env python3
"""
Knowledge-defined embedding for healthy human pancreas single-cell data.

The scenario is human pancreas with no narrower age, sex, condition, or
cell-type restriction. The dimensions therefore cover the major pancreatic
compartments and common baseline/sample-handling states:

0 endocrine islet, 1 acinar enzyme, 2 ductal/centroacinar, 3 beta,
4 alpha, 5 delta/PP/epsilon, 6 endocrine maturity/processing,
7 stellate/fibroblast, 8 activated ECM/myofibroblast, 9 endothelial,
10 mural vascular, 11 T/NK lymphoid, 12 B/plasma, 13 myeloid/APC,
14 mast, 15 Schwann/glial, 16 cycling, 17 interferon/inflammatory,
18 hypoxia/acute stress.
"""

from __future__ import annotations

import argparse


SIGNATURES = {
    "epithelial": [
        "EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "TACSTD2", "CLDN4",
        "MUC1",
    ],
    "acinar": [
        "PRSS1", "PRSS2", "PRSS3", "CPA1", "CPA2", "CPB1", "CTRB1",
        "CTRB2", "CTRL", "CELA2A", "CELA3A", "CELA3B", "PNLIP", "PNLIPRP1",
        "REG1A", "REG1B", "SYCN", "RBPJL", "AMY2A", "AMY2B",
    ],
    "ductal": [
        "KRT7", "KRT19", "SOX9", "CFTR", "SLC4A4", "SLC26A6", "ANXA4",
        "MMP7", "MUC1", "MUC6", "TSPAN8", "KLF5", "SPP1", "CLDN10",
    ],
    "centroacinar": [
        "SOX9", "CFTR", "SLC4A4", "SLC26A6", "CLDN10", "AQP1", "ANXA4",
        "MMP7",
    ],
    "endocrine": [
        "CHGA", "CHGB", "SCG2", "SCG3", "SCG5", "PCSK1N", "SYP",
        "ENO2", "ISL1", "NEUROD1", "PAX6", "NKX2-2", "RFX6",
    ],
    "beta": [
        "INS", "IAPP", "MAFA", "PDX1", "NKX6-1", "SLC30A8", "PCSK1",
        "PCSK2", "G6PC2", "ADCYAP1", "DLK1", "HADH", "FFAR1", "ABCC8",
        "KCNJ11",
    ],
    "alpha": [
        "GCG", "ARX", "IRX2", "IRX1", "PCSK2", "TTR", "TM4SF4", "LOXL4",
        "GC", "MAFB",
    ],
    "delta": [
        "SST", "HHEX", "RBP4", "LEPR", "BCHE", "PCSK1", "PCSK2",
    ],
    "pp_epsilon": [
        "PPY", "GHRL", "PAX6", "ISL1", "IAPP", "SLC18A2",
    ],
    "endocrine_maturity": [
        "PCSK1", "PCSK2", "SCG5", "SLC30A8", "MAFA", "PAX6", "NEUROD1",
        "RFX6", "ABCC8", "KCNJ11", "G6PC2", "CHGA", "CHGB",
    ],
    "stellate_fibroblast": [
        "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "C7", "PDGFRA",
        "CXCL12", "RGS5", "ACTA2", "VIM", "MMP2", "SPARC", "COL6A1",
        "COL6A2", "COL6A3",
    ],
    "activated_ecm": [
        "ACTA2", "TAGLN", "MYL9", "COL1A1", "COL1A2", "COL3A1",
        "FN1", "POSTN", "THBS2", "CTGF", "SERPINE1", "TIMP1", "MMP2",
        "MMP14", "LOX",
    ],
    "endothelial": [
        "PECAM1", "VWF", "KDR", "FLT1", "CDH5", "ESAM", "CLDN5",
        "RAMP2", "ENG", "EMCN", "PLVAP", "ACKR1", "SELE", "ICAM2",
    ],
    "mural": [
        "RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3", "ACTA2", "TAGLN",
        "MYH11", "MYL9", "CNN1", "DES", "TPM2",
    ],
    "t_nk": [
        "PTPRC", "CD3D", "CD3E", "CD2", "TRAC", "IL7R", "CCR7", "LCK",
        "NKG7", "GNLY", "PRF1", "GZMB", "KLRD1", "KLRB1",
    ],
    "b_plasma": [
        "PTPRC", "MS4A1", "CD79A", "CD79B", "BANK1", "CD74", "HLA-DRA",
        "MZB1", "JCHAIN", "XBP1", "SDC1", "IGKC", "IGHG1",
    ],
    "myeloid": [
        "PTPRC", "LYZ", "LST1", "FCGR3A", "FCER1G", "TYROBP", "AIF1",
        "CD14", "C1QA", "C1QB", "C1QC", "MS4A7", "HLA-DRA", "HLA-DPA1",
        "HLA-DPB1", "S100A8", "S100A9",
    ],
    "mast": [
        "TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC", "GATA2",
        "HPGDS", "RGS13",
    ],
    "schwann_glial": [
        "SOX10", "S100B", "PLP1", "MPZ", "MBP", "GFRA3", "NCAM1",
        "NRXN1", "NGFR", "ERBB3", "PMP22",
    ],
    "cycling": [
        "MKI67", "TOP2A", "UBE2C", "BIRC5", "CENPF", "NUSAP1", "MCM2",
        "MCM5", "PCNA", "TYMS", "STMN1", "HMGB2",
    ],
    "interferon_inflammatory": [
        "ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "IFI27", "MX1",
        "MX2", "OAS1", "OAS2", "STAT1", "IRF7", "CXCL10", "CXCL8",
        "CCL2", "NFKBIA", "TNFAIP3",
    ],
    "hypoxia_stress": [
        "HIF1A", "VEGFA", "SLC2A1", "LDHA", "BNIP3", "DDIT4", "ADM",
        "JUN", "FOS", "ATF3", "EGR1", "HSPA1A", "HSPA1B", "HSP90AA1",
        "DNAJB1",
    ],
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a marker-signature embedding for human pancreas cells."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument(
        "--seed",
        type=int,
        default=0,
        help="Run tag supplied by the orchestrator; not used by this deterministic scoring script.",
    )
    return parser.parse_args()


def _runtime_imports():
    global np, pd, ad, issparse
    import numpy as np  # noqa: PLC0415
    import pandas as pd  # noqa: PLC0415
    import anndata as ad  # noqa: PLC0415
    from scipy.sparse import issparse  # noqa: PLC0415


def _gene_lookup(var_names):
    exact = {}
    upper = {}
    for idx, name in enumerate(var_names):
        symbol = str(name)
        exact.setdefault(symbol, idx)
        upper.setdefault(symbol.upper(), idx)
    return exact, upper


def _indices_for_genes(genes, exact, upper):
    seen = set()
    indices = []
    for gene in genes:
        idx = exact.get(gene)
        if idx is None:
            idx = upper.get(gene.upper())
        if idx is not None and idx not in seen:
            seen.add(idx)
            indices.append(idx)
    return indices


def _as_dense_float32(matrix):
    if issparse(matrix):
        return matrix.toarray().astype("float32", copy=False)
    return np.asarray(matrix, dtype="float32")


def _robust_scale(values):
    values = np.asarray(values, dtype="float32")
    if values.size == 0:
        return values
    median = np.nanmedian(values)
    q25, q75 = np.nanpercentile(values, [25.0, 75.0])
    scale = (q75 - q25) / 1.349
    if not np.isfinite(scale) or scale < 1e-6:
        scale = np.nanstd(values)
    if not np.isfinite(scale) or scale < 1e-6:
        return np.zeros(values.shape, dtype="float32")
    scaled = (values - median) / scale
    return np.clip(scaled, -8.0, 8.0).astype("float32", copy=False)


def _module_score(x_matrix, indices, n_cells):
    if not indices:
        return np.zeros(n_cells, dtype="float32")

    sub = _as_dense_float32(x_matrix[:, indices])
    if sub.ndim == 1:
        sub = sub.reshape(n_cells, 1)

    gene_sd = sub.std(axis=0)
    usable = gene_sd > 1e-6
    if not np.any(usable):
        return np.zeros(n_cells, dtype="float32")

    sub = sub[:, usable]
    gene_mean = sub.mean(axis=0)
    gene_sd = gene_sd[usable]
    z = (sub - gene_mean) / gene_sd
    return _robust_scale(z.mean(axis=1))


def _build_scores(adata_obj):
    n_cells = adata_obj.n_obs
    exact, upper = _gene_lookup(adata_obj.var_names)
    scores = {}
    for name, genes in SIGNATURES.items():
        scores[name] = _module_score(
            adata_obj.X,
            _indices_for_genes(genes, exact, upper),
            n_cells,
        )
    return scores


def _mean_scores(*arrays):
    return np.mean(np.vstack(arrays), axis=0).astype("float32", copy=False)


def _dims_from_scores(s):
    epithelial_exocrine = _mean_scores(s["acinar"], s["ductal"])
    non_beta_endocrine = _mean_scores(s["alpha"], s["delta"], s["pp_epsilon"])

    raw_dims = [
        s["endocrine"] - 0.35 * epithelial_exocrine,
        s["acinar"] - 0.30 * s["ductal"] - 0.20 * s["endocrine"],
        _mean_scores(s["ductal"], s["centroacinar"]) - 0.25 * s["acinar"],
        s["beta"] + 0.35 * s["endocrine_maturity"] - 0.25 * non_beta_endocrine,
        s["alpha"] - 0.35 * s["beta"],
        _mean_scores(s["delta"], s["pp_epsilon"]) - 0.20 * _mean_scores(s["beta"], s["alpha"]),
        s["endocrine_maturity"] + 0.20 * s["endocrine"] - 0.20 * epithelial_exocrine,
        s["stellate_fibroblast"] - 0.25 * s["epithelial"],
        s["activated_ecm"] - 0.20 * s["stellate_fibroblast"],
        s["endothelial"] - 0.25 * s["epithelial"],
        s["mural"] - 0.30 * s["endothelial"],
        s["t_nk"] - 0.20 * s["myeloid"],
        s["b_plasma"] - 0.20 * s["t_nk"],
        s["myeloid"] - 0.20 * _mean_scores(s["t_nk"], s["b_plasma"]),
        s["mast"] - 0.15 * s["myeloid"],
        s["schwann_glial"] - 0.20 * s["mural"],
        s["cycling"],
        s["interferon_inflammatory"] - 0.15 * s["myeloid"],
        s["hypoxia_stress"] - 0.10 * s["cycling"],
    ]

    return np.column_stack([_robust_scale(dim) for dim in raw_dims]).astype("float32")


def run(input_path, output_path):
    _runtime_imports()
    adata_obj = ad.read_h5ad(input_path)
    if adata_obj.n_obs == 0:
        coords = np.zeros((0, 19), dtype="float32")
    else:
        coords = _dims_from_scores(_build_scores(adata_obj))

    out = pd.DataFrame({"cell_id": adata_obj.obs_names.astype(str).to_numpy()})
    for dim_idx in range(coords.shape[1]):
        out[f"dim_{dim_idx}"] = coords[:, dim_idx].astype("float32", copy=False)
    out.to_parquet(output_path, index=False)


def main():
    args = parse_args()
    run(args.input, args.output)


if __name__ == "__main__":
    main()
