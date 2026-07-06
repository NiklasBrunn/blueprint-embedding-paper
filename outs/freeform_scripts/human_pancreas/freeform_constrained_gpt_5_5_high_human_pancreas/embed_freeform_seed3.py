#!/usr/bin/env python3
"""Literature-derived embedding for healthy human pancreas single-cell data.

The scenario is unrestricted adult/healthy human pancreas, so the embedding
covers endocrine, exocrine, ductal, vascular, stromal, neural, immune, and
common non-disease cell-state variation expected in normal pancreas atlases.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Dict, Iterable, List, Mapping, Sequence


@dataclass(frozen=True)
class Signature:
    name: str
    genes: Sequence[str]
    min_genes: int = 2


SIGNATURES: Sequence[Signature] = (
    Signature(
        "pan_endocrine",
        # Segerstolpe et al., Cell Metabolism 2016, 10.1016/j.cmet.2016.08.020;
        # Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004;
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896.
        genes=[
            "CHGA",
            "CHGB",
            "SCG2",
            "SCG3",
            "SCG5",
            "PCSK1",
            "PCSK2",
            "NEUROD1",
            "ISL1",
            "PAX6",
        ],
    ),
    Signature(
        "beta_cell",
        # Segerstolpe et al., Cell Metabolism 2016, 10.1016/j.cmet.2016.08.020;
        # Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
        genes=[
            "INS",
            "IAPP",
            "MAFA",
            "PDX1",
            "NKX6-1",
            "UCN3",
            "HADH",
            "SLC2A2",
            "PCSK1",
            "RBP4",
        ],
    ),
    Signature(
        "alpha_cell",
        # Segerstolpe et al., Cell Metabolism 2016, 10.1016/j.cmet.2016.08.020;
        # Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
        genes=[
            "GCG",
            "TTR",
            "ARX",
            "IRX2",
            "MAFB",
            "PCSK2",
            "GC",
            "LOXL4",
            "CRYBA2",
            "PLCE1",
        ],
    ),
    Signature(
        "delta_cell",
        # Segerstolpe et al., Cell Metabolism 2016, 10.1016/j.cmet.2016.08.020;
        # Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
        genes=[
            "SST",
            "HHEX",
            "LEPR",
            "GHSR",
            "RBP4",
            "PCSK1",
            "PAX6",
            "ISL1",
        ],
    ),
    Signature(
        "gamma_pp_cell",
        # Segerstolpe et al., Cell Metabolism 2016, 10.1016/j.cmet.2016.08.020;
        # Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004.
        genes=[
            "PPY",
            "PYY",
            "TSPAN8",
            "SERTM1",
            "SLITRK6",
            "PCSK1",
            "CHGA",
        ],
    ),
    Signature(
        "epsilon_cell",
        # Segerstolpe et al., Cell Metabolism 2016, 10.1016/j.cmet.2016.08.020;
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896.
        genes=[
            "GHRL",
            "MLN",
            "GAST",
            "CHGA",
            "CHGB",
            "ISL1",
        ],
        min_genes=1,
    ),
    Signature(
        "acinar",
        # Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004;
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896.
        genes=[
            "PRSS1",
            "PRSS2",
            "CTRB1",
            "CTRB2",
            "CPA1",
            "CPA2",
            "CPB1",
            "PNLIP",
            "CEL",
            "CLPS",
            "REG1A",
            "REG3A",
            "SYCN",
            "RBPJL",
        ],
    ),
    Signature(
        "ductal",
        # Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004;
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896.
        genes=[
            "KRT19",
            "KRT7",
            "KRT8",
            "KRT18",
            "SOX9",
            "SLC4A4",
            "CFTR",
            "MSLN",
            "MUC1",
            "EPCAM",
            "SPP1",
            "ANXA4",
        ],
    ),
    Signature(
        "endothelial",
        # Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004;
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896.
        genes=[
            "PECAM1",
            "VWF",
            "KDR",
            "CDH5",
            "EMCN",
            "RAMP2",
            "ESAM",
            "CLDN5",
            "FLT1",
        ],
    ),
    Signature(
        "mural_smooth_muscle",
        # Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004;
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896.
        genes=[
            "RGS5",
            "PDGFRB",
            "CSPG4",
            "MCAM",
            "NOTCH3",
            "ACTA2",
            "TAGLN",
            "MYH11",
            "MYL9",
        ],
    ),
    Signature(
        "fibroblast_stellate",
        # Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004;
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896.
        genes=[
            "COL1A1",
            "COL1A2",
            "COL3A1",
            "DCN",
            "LUM",
            "COL6A1",
            "COL6A2",
            "PDGFRA",
            "CXCL14",
            "MMP2",
        ],
    ),
    Signature(
        "activated_ecm_stellate",
        # Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004;
        # Reactome extracellular matrix organization, R-HSA-1474244.
        genes=[
            "ACTA2",
            "TAGLN",
            "COL1A1",
            "COL3A1",
            "FN1",
            "SPARC",
            "THY1",
            "POSTN",
            "TIMP1",
            "MMP2",
        ],
    ),
    Signature(
        "myeloid_antigen_presenting",
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896;
        # GO antigen processing and presentation of peptide antigen, GO:0048002.
        genes=[
            "LST1",
            "C1QA",
            "C1QB",
            "AIF1",
            "CD68",
            "MS4A7",
            "FCGR3A",
            "HLA-DRA",
            "LYZ",
            "ITGAX",
            "FCER1A",
            "CLEC10A",
        ],
    ),
    Signature(
        "t_nk_lymphoid",
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896;
        # GO T cell activation, GO:0042110; GO natural killer cell activation, GO:0030101.
        genes=[
            "PTPRC",
            "CD3D",
            "CD3E",
            "TRAC",
            "IL7R",
            "NKG7",
            "GNLY",
            "KLRD1",
            "GZMB",
        ],
    ),
    Signature(
        "b_plasma",
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896;
        # GO B cell activation, GO:0042113.
        genes=[
            "MS4A1",
            "CD79A",
            "CD74",
            "BANK1",
            "MZB1",
            "JCHAIN",
            "XBP1",
            "IGKC",
        ],
    ),
    Signature(
        "mast_cell",
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896;
        # GO mast cell activation, GO:0045576.
        genes=[
            "TPSAB1",
            "TPSB2",
            "CPA3",
            "KIT",
            "HDC",
            "MS4A2",
        ],
    ),
    Signature(
        "schwann_neural",
        # Enge et al., Cell 2017, 10.1016/j.cell.2017.09.004;
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896.
        genes=[
            "S100B",
            "SOX10",
            "PLP1",
            "MPZ",
            "GFRA3",
            "NGFR",
            "SCN7A",
            "NRXN1",
        ],
    ),
    Signature(
        "cell_cycle",
        # Segerstolpe et al., Cell Metabolism 2016, 10.1016/j.cmet.2016.08.020;
        # MSigDB HALLMARK_G2M_CHECKPOINT; MSigDB HALLMARK_E2F_TARGETS.
        genes=[
            "MKI67",
            "TOP2A",
            "CENPF",
            "CDK1",
            "BIRC5",
            "UBE2C",
            "PCNA",
            "MCM5",
            "HMGB2",
        ],
    ),
    Signature(
        "interferon_inflammatory",
        # Tabula Sapiens Consortium, Science 2022, 10.1126/science.abl4896;
        # MSigDB HALLMARK_INTERFERON_ALPHA_RESPONSE; MSigDB HALLMARK_INTERFERON_GAMMA_RESPONSE.
        genes=[
            "IFIT1",
            "IFIT2",
            "IFIT3",
            "ISG15",
            "MX1",
            "OAS1",
            "STAT1",
            "CXCL10",
            "IRF7",
        ],
    ),
)


DIMENSION_DESIGN = (
    "dim_0 pan-endocrine separation",
    "dim_1 beta-to-alpha endocrine contrast",
    "dim_2 minor endocrine states",
    "dim_3 acinar exocrine identity",
    "dim_4 ductal epithelial identity",
    "dim_5 acinar-to-ductal epithelial contrast",
    "dim_6 endothelial vascular identity",
    "dim_7 mural and smooth-muscle identity",
    "dim_8 fibroblast/stellate identity",
    "dim_9 activated extracellular-matrix stellate state",
    "dim_10 myeloid and antigen-presenting immune identity",
    "dim_11 T/NK lymphoid identity",
    "dim_12 B/plasma identity",
    "dim_13 mast-cell identity",
    "dim_14 Schwann/neural identity",
    "dim_15 cycling state",
    "dim_16 interferon/inflammatory state",
)


def _as_array(matrix):
    if hasattr(matrix, "toarray"):
        return matrix.toarray()
    return matrix


def _robust_z(values):
    import numpy as np

    values = np.asarray(values, dtype=np.float32)
    if values.size == 0:
        return values
    median = np.nanmedian(values)
    mad = np.nanmedian(np.abs(values - median))
    if not np.isfinite(mad) or mad <= 1e-6:
        std = np.nanstd(values)
        if not np.isfinite(std) or std <= 1e-6:
            return np.zeros(values.shape, dtype=np.float32)
        scaled = (values - np.nanmean(values)) / std
    else:
        scaled = (values - median) / (1.4826 * mad)
    return np.clip(scaled, -8.0, 8.0).astype(np.float32)


def _present_indices(var_names: Sequence[str], genes: Iterable[str]) -> List[int]:
    first_index_by_upper: Dict[str, int] = {}
    for i, gene in enumerate(var_names):
        first_index_by_upper.setdefault(str(gene).upper(), i)

    indices: List[int] = []
    seen = set()
    for gene in genes:
        key = gene.upper()
        if key in first_index_by_upper and key not in seen:
            indices.append(first_index_by_upper[key])
            seen.add(key)
    return indices


def _score_signature(matrix, var_names: Sequence[str], signature: Signature):
    import numpy as np

    n_cells = matrix.shape[0]
    indices = _present_indices(var_names, signature.genes)
    if len(indices) < signature.min_genes:
        return np.zeros(n_cells, dtype=np.float32)

    block = _as_array(matrix[:, indices])
    block = np.asarray(block, dtype=np.float32)
    if block.ndim == 1:
        block = block.reshape(-1, 1)

    detected = (block > 0).mean(axis=1)
    mean_expr = block.mean(axis=1)
    score = _robust_z(mean_expr) + 0.35 * _robust_z(detected)
    return _robust_z(score)


def _combine_dimensions(scores: Mapping[str, object]):
    import numpy as np

    zeros = np.zeros_like(next(iter(scores.values())), dtype=np.float32)

    def s(name: str):
        return np.asarray(scores.get(name, zeros), dtype=np.float32)

    endocrine = s("pan_endocrine")
    beta = s("beta_cell")
    alpha = s("alpha_cell")
    delta = s("delta_cell")
    gamma = s("gamma_pp_cell")
    epsilon = s("epsilon_cell")
    acinar = s("acinar")
    ductal = s("ductal")
    endothelial = s("endothelial")
    mural = s("mural_smooth_muscle")
    fibro = s("fibroblast_stellate")
    activated_ecm = s("activated_ecm_stellate")
    myeloid = s("myeloid_antigen_presenting")
    lymphoid = s("t_nk_lymphoid")
    b_plasma = s("b_plasma")
    mast = s("mast_cell")
    neural = s("schwann_neural")
    cycling = s("cell_cycle")
    interferon = s("interferon_inflammatory")

    dimensions = [
        endocrine - 0.20 * (acinar + ductal + endothelial + fibro),
        beta - alpha,
        delta + 0.60 * gamma + 0.60 * epsilon - 0.25 * (beta + alpha),
        acinar - 0.15 * ductal,
        ductal - 0.15 * acinar,
        acinar - ductal,
        endothelial - 0.20 * mural,
        mural - 0.20 * endothelial,
        fibro - 0.20 * activated_ecm,
        activated_ecm - 0.20 * fibro,
        myeloid - 0.20 * (lymphoid + b_plasma + mast),
        lymphoid - 0.20 * (myeloid + b_plasma + mast),
        b_plasma - 0.20 * (myeloid + lymphoid + mast),
        mast - 0.20 * (myeloid + lymphoid + b_plasma),
        neural,
        cycling,
        interferon,
    ]
    return np.vstack([_robust_z(dim) for dim in dimensions]).T.astype(np.float32)


def build_embedding(input_path: str):
    import anndata as ad
    import pandas as pd

    adata = ad.read_h5ad(input_path)
    matrix = adata.X
    if "counts" in adata.layers:
        # The contract says X is log1p-normalized and counts are available.
        # Scoring uses X so gene modules are comparable across library sizes.
        matrix = adata.X

    var_names = [str(x) for x in adata.var_names]
    scores = {
        signature.name: _score_signature(matrix, var_names, signature)
        for signature in SIGNATURES
    }
    embedding = _combine_dimensions(scores)

    columns = {f"dim_{i}": embedding[:, i].astype("float32") for i in range(embedding.shape[1])}
    frame = pd.DataFrame({"cell_id": adata.obs_names.astype(str).to_numpy(), **columns})
    return frame


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Build a literature-derived embedding for healthy human pancreas."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad AnnData file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", required=False, default=None, help="Run tag; ignored.")
    args = parser.parse_args(argv)

    frame = build_embedding(args.input)
    frame.to_parquet(args.output, index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
