#!/usr/bin/env python3
"""Literature-grounded freeform embedding for healthy human pancreas."""

from __future__ import annotations

import argparse


# Human pancreas acinar/exocrine enzyme signature.
# Atlas/paper support: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004.
ACINAR_GENES = [
    "PRSS1",
    "PRSS2",
    "CPA1",
    "CPA2",
    "CPB1",
    "CTRB1",
    "CTRB2",
    "CEL",
    "PNLIP",
    "CLPS",
    "REG1A",
    "REG1B",
]

# Human pancreas ductal/centroacinar epithelial signature.
# Atlas/paper support: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004.
DUCTAL_GENES = [
    "KRT19",
    "KRT7",
    "KRT8",
    "KRT18",
    "EPCAM",
    "MUC1",
    "CFTR",
    "SOX9",
    "SPP1",
    "ANXA4",
    "TSPAN8",
]

# Pancreatic endocrine/islet core signature.
# Atlas/paper support: Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004.
ENDOCRINE_CORE_GENES = [
    "CHGA",
    "CHGB",
    "SCGN",
    "PCSK1N",
    "NEUROD1",
    "PAX6",
    "ISL1",
    "RFX6",
]

# Human pancreatic beta-cell identity and mature insulin-secretory signature.
# Paper support: Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004.
BETA_GENES = [
    "INS",
    "IAPP",
    "MAFA",
    "PDX1",
    "NKX6-1",
    "PCSK1",
    "SLC30A8",
    "UCN3",
    "G6PC2",
    "DLK1",
]

# Human pancreatic alpha-cell identity signature.
# Paper support: Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004.
ALPHA_GENES = [
    "GCG",
    "ARX",
    "IRX2",
    "TTR",
    "GC",
    "TM4SF4",
    "LOXL4",
    "PCSK2",
]

# Human pancreatic delta-cell identity signature.
# Paper support: Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004.
DELTA_GENES = [
    "SST",
    "HHEX",
    "RBP4",
    "LEPR",
    "BCHE",
    "PCSK1",
    "RGS4",
]

# Human pancreatic gamma/PP-cell identity signature.
# Paper support: Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004.
PP_GAMMA_GENES = [
    "PPY",
    "PAX6",
    "SLC30A8",
    "PCSK1N",
    "TM4SF4",
    "TTR",
]

# Human pancreatic epsilon-cell/ghrelin endocrine signature.
# Paper support: Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004.
EPSILON_GENES = [
    "GHRL",
    "MLN",
    "PAX6",
    "PCSK1N",
    "RBP4",
]

# Pancreatic fibroblast/stellate extracellular-matrix signature.
# Atlas/paper support: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004.
# Database support: MSigDB HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.
FIBROBLAST_STELLATE_GENES = [
    "COL1A1",
    "COL1A2",
    "COL3A1",
    "DCN",
    "LUM",
    "MMP2",
    "PDGFRA",
    "FN1",
    "THY1",
    "FAP",
]

# Pancreatic endothelial vascular signature.
# Atlas/paper support: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004.
# Database support: GO:0001938 positive regulation of endothelial cell proliferation.
ENDOTHELIAL_GENES = [
    "PECAM1",
    "VWF",
    "CDH5",
    "KDR",
    "EMCN",
    "CLDN5",
    "RAMP2",
    "ESAM",
    "ENG",
]

# Pancreatic mural/pericyte/smooth-muscle vascular support signature.
# Atlas/paper support: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004.
MURAL_GENES = [
    "RGS5",
    "NOTCH3",
    "MCAM",
    "CSPG4",
    "PDGFRB",
    "ACTA2",
    "MYH11",
    "TAGLN",
]

# Resident pancreatic myeloid/APC signature.
# Atlas support: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
# Database support: GO:0002376 immune system process; Reactome R-HSA-2132295 MHC class II antigen presentation.
MYELOID_GENES = [
    "LYZ",
    "LST1",
    "FCGR3A",
    "CD14",
    "CD68",
    "MS4A7",
    "HLA-DRA",
    "HLA-DPA1",
    "C1QA",
    "C1QB",
]

# Pancreatic T/NK lymphoid signature.
# Atlas support: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
# Database support: GO:0002376 immune system process.
T_NK_GENES = [
    "CD3D",
    "CD3E",
    "TRAC",
    "IL7R",
    "NKG7",
    "GNLY",
    "PRF1",
    "KLRD1",
]

# Pancreatic B/plasma-cell signature.
# Atlas support: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
# Database support: GO:0002376 immune system process; Reactome R-HSA-983705 antigen activates B-cell receptor.
B_PLASMA_GENES = [
    "MS4A1",
    "CD79A",
    "CD79B",
    "BANK1",
    "MZB1",
    "JCHAIN",
    "IGKC",
    "IGHG1",
]

# Pancreatic mast-cell signature.
# Atlas support: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
# Database support: GO:0002448 mast cell mediated immunity.
MAST_GENES = [
    "TPSAB1",
    "TPSB2",
    "CPA3",
    "KIT",
    "HDC",
    "MS4A2",
]


SIGNATURES = {
    "acinar": ACINAR_GENES,
    "ductal": DUCTAL_GENES,
    "endocrine_core": ENDOCRINE_CORE_GENES,
    "beta": BETA_GENES,
    "alpha": ALPHA_GENES,
    "delta": DELTA_GENES,
    "pp_gamma": PP_GAMMA_GENES,
    "epsilon": EPSILON_GENES,
    "fibroblast_stellate": FIBROBLAST_STELLATE_GENES,
    "endothelial": ENDOTHELIAL_GENES,
    "mural": MURAL_GENES,
    "myeloid": MYELOID_GENES,
    "t_nk": T_NK_GENES,
    "b_plasma": B_PLASMA_GENES,
    "mast": MAST_GENES,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a literature-signature embedding for human pancreas h5ad data."
    )
    parser.add_argument("--input", required=True, help="Input AnnData .h5ad path.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, default=0, help="Run tag; not used by scoring.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    import anndata as ad
    import numpy as np
    import pandas as pd
    from scipy import sparse

    adata = ad.read_h5ad(args.input)
    x = adata.X
    n_cells = adata.n_obs

    var_names = [str(gene) for gene in adata.var_names]
    exact_index = {}
    upper_index = {}
    for idx, gene in enumerate(var_names):
        exact_index.setdefault(gene, idx)
        upper_index.setdefault(gene.upper(), idx)

    def present_indices(genes: list[str]) -> list[int]:
        indices = []
        seen = set()
        for gene in genes:
            idx = exact_index.get(gene)
            if idx is None:
                idx = upper_index.get(gene.upper())
            if idx is not None and idx not in seen:
                indices.append(idx)
                seen.add(idx)
        return indices

    def module_mean(genes: list[str]) -> np.ndarray:
        idx = present_indices(genes)
        if not idx:
            return np.zeros(n_cells, dtype=np.float32)
        sub = x[:, idx]
        if sparse.issparse(sub):
            values = np.asarray(sub.mean(axis=1)).ravel()
        else:
            values = np.asarray(sub, dtype=np.float32).mean(axis=1)
        return values.astype(np.float32, copy=False)

    def robust_scale(values: np.ndarray) -> np.ndarray:
        values = np.asarray(values, dtype=np.float32)
        finite = np.isfinite(values)
        if not finite.any():
            return np.zeros_like(values, dtype=np.float32)
        clean = values[finite]
        median = float(np.median(clean))
        q25, q75 = np.percentile(clean, [25, 75])
        scale = float((q75 - q25) / 1.349)
        if scale < 1e-6:
            scale = float(np.std(clean))
        if scale < 1e-6:
            out = np.zeros_like(values, dtype=np.float32)
        else:
            out = (values - median) / scale
        out = np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)
        return np.clip(out, -8.0, 8.0).astype(np.float32)

    scores = {name: module_mean(genes) for name, genes in SIGNATURES.items()}

    exocrine = 0.5 * (scores["acinar"] + scores["ductal"])
    epithelial = (scores["acinar"] + scores["ductal"] + scores["endocrine_core"]) / 3.0
    mesenchymal = 0.5 * (scores["fibroblast_stellate"] + scores["mural"])
    vascular = 0.5 * (scores["endothelial"] + scores["mural"])
    immune = (
        scores["myeloid"] + scores["t_nk"] + scores["b_plasma"] + scores["mast"]
    ) / 4.0
    nonimmune_tissue = (epithelial + mesenchymal + vascular) / 3.0
    rare_islet = (scores["delta"] + scores["pp_gamma"] + scores["epsilon"]) / 3.0

    dims = [
        robust_scale(scores["endocrine_core"] - exocrine),
        robust_scale(scores["acinar"] - scores["ductal"]),
        robust_scale(scores["beta"] - scores["alpha"]),
        robust_scale(scores["delta"] - (scores["beta"] + scores["alpha"] + rare_islet) / 3.0),
        robust_scale(scores["pp_gamma"] - scores["epsilon"]),
        robust_scale(epithelial - mesenchymal),
        robust_scale(scores["endothelial"] - scores["mural"]),
        robust_scale(scores["fibroblast_stellate"] - scores["mural"]),
        robust_scale(immune - nonimmune_tissue),
        robust_scale((scores["t_nk"] + scores["b_plasma"]) / 2.0 - (scores["myeloid"] + scores["mast"]) / 2.0),
        robust_scale(scores["b_plasma"] - scores["t_nk"]),
        robust_scale(scores["mast"] - scores["myeloid"]),
    ]

    coords = np.column_stack(dims).astype(np.float32, copy=False)
    data = {"cell_id": pd.Index(adata.obs_names.astype(str), dtype="string")}
    for dim_idx in range(coords.shape[1]):
        data[f"dim_{dim_idx}"] = coords[:, dim_idx]
    pd.DataFrame(data).to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
