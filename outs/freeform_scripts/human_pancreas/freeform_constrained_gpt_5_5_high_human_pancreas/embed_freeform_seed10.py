#!/usr/bin/env python3
"""Literature-derived embedding for healthy human pancreas scRNA-seq."""

from __future__ import annotations

import argparse
from collections import OrderedDict
from typing import Iterable


SIGNATURES = OrderedDict(
    [
        (
            "endocrine_common",
            # Baron et al., A Single-Cell Transcriptomic Map of the Human and Mouse Pancreas, 10.1016/j.cels.2016.08.011; Muraro et al., A Single-Cell Transcriptome Atlas of the Human Pancreas, 10.1016/j.cels.2016.09.002
            [
                "CHGA",
                "CHGB",
                "SCG2",
                "SCG3",
                "PCSK1N",
                "PCSK2",
                "SLC18A1",
                "NEUROD1",
                "ISL1",
                "PAX6",
                "RFX6",
            ],
        ),
        (
            "beta_cell",
            # Segerstolpe et al., Single-Cell Transcriptome Profiling of Human Pancreatic Islets in Health and Type 2 Diabetes, 10.1016/j.cmet.2016.08.020; Baron et al., 10.1016/j.cels.2016.08.011
            [
                "INS",
                "IAPP",
                "SLC30A8",
                "PCSK1",
                "NKX6-1",
                "PDX1",
                "MAFA",
                "G6PC2",
                "HADH",
                "DLK1",
            ],
        ),
        (
            "beta_mature",
            # Segerstolpe et al., 10.1016/j.cmet.2016.08.020; Enge et al., Single-Cell Analysis of Human Pancreas Reveals Transcriptional Signatures of Aging and Somatic Mutation Patterns, 10.1016/j.cell.2017.09.004
            [
                "INS",
                "IAPP",
                "SLC30A8",
                "UCN3",
                "MAFA",
                "PDX1",
                "NKX6-1",
                "G6PC2",
                "FFAR4",
                "HADH",
                "PCSK1",
            ],
        ),
        (
            "alpha_cell",
            # Segerstolpe et al., 10.1016/j.cmet.2016.08.020; Muraro et al., 10.1016/j.cels.2016.09.002
            [
                "GCG",
                "ARX",
                "IRX2",
                "TTR",
                "LOXL4",
                "TM4SF4",
                "CRYBA2",
                "GC",
                "PCSK2",
                "RGS4",
            ],
        ),
        (
            "delta_cell",
            # Segerstolpe et al., 10.1016/j.cmet.2016.08.020; Baron et al., 10.1016/j.cels.2016.08.011
            [
                "SST",
                "HHEX",
                "RBP4",
                "LEPR",
                "BCHE",
                "GHSR",
                "PCSK1",
                "PCSK2",
                "SEC11C",
            ],
        ),
        (
            "pp_gamma_cell",
            # Muraro et al., 10.1016/j.cels.2016.09.002; Segerstolpe et al., 10.1016/j.cmet.2016.08.020
            [
                "PPY",
                "PYY",
                "SLC18A2",
                "MEIS2",
                "ID4",
                "ABCC8",
                "PCSK1N",
                "SCGN",
            ],
        ),
        (
            "epsilon_cell",
            # Baron et al., 10.1016/j.cels.2016.08.011; Muraro et al., 10.1016/j.cels.2016.09.002
            [
                "GHRL",
                "MLN",
                "IAPP",
                "PCSK1",
                "PCSK2",
                "CHGA",
                "CHGB",
                "RBP4",
            ],
        ),
        (
            "acinar",
            # Muraro et al., 10.1016/j.cels.2016.09.002; Baron et al., 10.1016/j.cels.2016.08.011
            [
                "PRSS1",
                "PRSS2",
                "CPA1",
                "CPA2",
                "CTRB1",
                "CTRB2",
                "CELA3A",
                "CELA3B",
                "PNLIP",
                "CLPS",
                "REG1A",
                "AMY2A",
                "GP2",
            ],
        ),
        (
            "ductal",
            # Muraro et al., 10.1016/j.cels.2016.09.002; Tabula Sapiens Consortium, The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans, 10.1126/science.abl4896
            [
                "KRT19",
                "KRT7",
                "SOX9",
                "MUC1",
                "CFTR",
                "SPP1",
                "ANXA4",
                "EPCAM",
                "TSPAN8",
                "MMP7",
                "CLDN10",
                "LCN2",
            ],
        ),
        (
            "stellate_quiescent_ecm",
            # Baron et al., 10.1016/j.cels.2016.08.011; Tabula Sapiens Consortium, 10.1126/science.abl4896
            [
                "DCN",
                "LUM",
                "COL1A1",
                "COL1A2",
                "COL3A1",
                "COL6A1",
                "COL6A2",
                "MGP",
                "VIM",
                "PDGFRA",
            ],
        ),
        (
            "stellate_activated_myofibroblast",
            # Baron et al., 10.1016/j.cels.2016.08.011; Enge et al., 10.1016/j.cell.2017.09.004
            [
                "ACTA2",
                "TAGLN",
                "MYL9",
                "CNN1",
                "COL1A1",
                "COL1A2",
                "SPARC",
                "FN1",
                "CTGF",
                "TIMP1",
                "MMP2",
            ],
        ),
        (
            "endothelial",
            # Baron et al., 10.1016/j.cels.2016.08.011; Tabula Sapiens Consortium, 10.1126/science.abl4896
            [
                "PECAM1",
                "VWF",
                "KDR",
                "ESAM",
                "RAMP2",
                "CLDN5",
                "EMCN",
                "PLVAP",
                "CDH5",
                "FLT1",
                "ENG",
            ],
        ),
        (
            "pericyte_smooth_muscle",
            # Tabula Sapiens Consortium, 10.1126/science.abl4896; Baron et al., 10.1016/j.cels.2016.08.011
            [
                "RGS5",
                "PDGFRB",
                "CSPG4",
                "MCAM",
                "NOTCH3",
                "NDUFA4L2",
                "MYH11",
                "ACTA2",
                "TAGLN",
                "MYL9",
            ],
        ),
        (
            "t_nk_cell",
            # Baron et al., 10.1016/j.cels.2016.08.011; Tabula Sapiens Consortium, 10.1126/science.abl4896
            [
                "CD3D",
                "CD3E",
                "CD2",
                "TRAC",
                "IL7R",
                "NKG7",
                "GNLY",
                "KLRD1",
                "GZMB",
                "PRF1",
            ],
        ),
        (
            "b_plasma_cell",
            # Baron et al., 10.1016/j.cels.2016.08.011; Tabula Sapiens Consortium, 10.1126/science.abl4896
            [
                "MS4A1",
                "CD79A",
                "CD79B",
                "BANK1",
                "IGHM",
                "IGKC",
                "JCHAIN",
                "MZB1",
                "XBP1",
                "SDC1",
            ],
        ),
        (
            "myeloid_macrophage",
            # Baron et al., 10.1016/j.cels.2016.08.011; Tabula Sapiens Consortium, 10.1126/science.abl4896
            [
                "LST1",
                "TYROBP",
                "C1QA",
                "C1QB",
                "C1QC",
                "AIF1",
                "FCGR3A",
                "CD14",
                "MSR1",
                "MRC1",
                "LYZ",
            ],
        ),
        (
            "mast_cell",
            # Baron et al., 10.1016/j.cels.2016.08.011; Tabula Sapiens Consortium, 10.1126/science.abl4896
            [
                "TPSAB1",
                "TPSB2",
                "CPA3",
                "KIT",
                "MS4A2",
                "HPGDS",
                "GATA2",
                "HDC",
            ],
        ),
        (
            "schwann_neural",
            # Baron et al., 10.1016/j.cels.2016.08.011; Muraro et al., 10.1016/j.cels.2016.09.002
            [
                "SOX10",
                "PLP1",
                "MPZ",
                "MBP",
                "S100B",
                "NRXN1",
                "NCAM1",
                "GPM6B",
                "ERBB3",
            ],
        ),
        (
            "cell_cycle",
            # MSigDB HALLMARK_G2M_CHECKPOINT and HALLMARK_E2F_TARGETS; Segerstolpe et al. reports proliferating islet cells, 10.1016/j.cmet.2016.08.020
            [
                "MKI67",
                "TOP2A",
                "CENPF",
                "BIRC5",
                "CDK1",
                "CCNB1",
                "CCNB2",
                "UBE2C",
                "AURKB",
                "HMGB2",
                "TYMS",
                "MCM5",
            ],
        ),
    ]
)


DIMENSION_NAMES = [
    "endocrine_vs_exocrine",
    "acinar_vs_ductal",
    "beta_vs_alpha",
    "delta_vs_pp_epsilon",
    "mesenchymal_ecm",
    "vascular_wall",
    "lymphoid",
    "myeloid_mast",
    "neural",
    "beta_maturation",
    "stellate_activation",
    "cell_cycle",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Score fixed, literature-curated human pancreas signatures and write "
            "a 12-dimensional parquet embedding."
        )
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file")
    parser.add_argument("--output", required=True, help="Output parquet path")
    parser.add_argument("--seed", type=int, default=0, help="Run tag; not used")
    return parser.parse_args()


def _as_1d_float(values):
    import numpy as np

    if hasattr(values, "A1"):
        values = values.A1
    else:
        values = np.asarray(values).reshape(-1)
    return values.astype(np.float32, copy=False)


def _mean_expression(matrix, indices):
    import numpy as np

    if not indices:
        return np.zeros(matrix.shape[0], dtype=np.float32)
    block = matrix[:, indices]
    return _as_1d_float(block.mean(axis=1))


def _robust_scale(scores):
    import numpy as np

    scores = scores.astype(np.float32, copy=False)
    median = np.nanmedian(scores, axis=0)
    q25 = np.nanpercentile(scores, 25, axis=0)
    q75 = np.nanpercentile(scores, 75, axis=0)
    scale = q75 - q25
    std = np.nanstd(scores, axis=0)
    scale = np.where(scale > 1e-6, scale, np.where(std > 1e-6, std, 1.0))
    z = (scores - median) / scale
    return np.clip(z, -8.0, 8.0).astype(np.float32)


def _gene_index(var_names: Iterable[str]) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for pos, name in enumerate(var_names):
        key = str(name).upper()
        index.setdefault(key, []).append(pos)
    return index


def _signature_indices(var_names: Iterable[str]) -> dict[str, list[int]]:
    gene_to_pos = _gene_index(var_names)
    out: dict[str, list[int]] = {}
    for sig_name, genes in SIGNATURES.items():
        positions: list[int] = []
        seen: set[int] = set()
        for gene in genes:
            for pos in gene_to_pos.get(gene.upper(), []):
                if pos not in seen:
                    positions.append(pos)
                    seen.add(pos)
        out[sig_name] = positions
    return out


def _score_signatures(adata):
    import numpy as np

    matrix = adata.X
    sig_to_indices = _signature_indices(adata.var_names)
    raw = np.column_stack(
        [_mean_expression(matrix, sig_to_indices[name]) for name in SIGNATURES]
    ).astype(np.float32)
    scaled = _robust_scale(raw)
    return {name: scaled[:, i] for i, name in enumerate(SIGNATURES)}


def _embedding_from_scores(score):
    import numpy as np

    endocrine = score["endocrine_common"]
    exocrine = 0.5 * (score["acinar"] + score["ductal"])
    epithelial = (endocrine + exocrine) / 2.0
    immune = (
        score["t_nk_cell"]
        + score["b_plasma_cell"]
        + score["myeloid_macrophage"]
        + score["mast_cell"]
    ) / 4.0
    mesenchymal = (
        score["stellate_quiescent_ecm"]
        + score["stellate_activated_myofibroblast"]
        + score["pericyte_smooth_muscle"]
    ) / 3.0

    dims = [
        endocrine - exocrine,
        score["acinar"] - score["ductal"],
        score["beta_cell"] - score["alpha_cell"],
        score["delta_cell"] - 0.5 * (score["pp_gamma_cell"] + score["epsilon_cell"]),
        mesenchymal - epithelial,
        score["endothelial"] + 0.5 * score["pericyte_smooth_muscle"] - epithelial,
        0.5 * (score["t_nk_cell"] + score["b_plasma_cell"]) - epithelial,
        0.5 * (score["myeloid_macrophage"] + score["mast_cell"]) - epithelial,
        score["schwann_neural"] - epithelial,
        score["beta_mature"] - 0.5 * (score["alpha_cell"] + score["delta_cell"]),
        score["stellate_activated_myofibroblast"]
        - score["stellate_quiescent_ecm"],
        score["cell_cycle"] - 0.25 * (epithelial + mesenchymal + immune),
    ]
    emb = np.column_stack(dims).astype(np.float32)
    return np.clip(emb, -12.0, 12.0).astype(np.float32)


def main() -> None:
    args = parse_args()

    import anndata as ad
    import numpy as np
    import pandas as pd

    adata = ad.read_h5ad(args.input)
    score = _score_signatures(adata)
    emb = _embedding_from_scores(score)

    data = {"cell_id": np.asarray(adata.obs_names, dtype=str)}
    for i in range(emb.shape[1]):
        data[f"dim_{i}"] = emb[:, i].astype(np.float32, copy=False)
    out = pd.DataFrame(data, columns=["cell_id"] + [f"dim_{i}" for i in range(emb.shape[1])])
    out.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
