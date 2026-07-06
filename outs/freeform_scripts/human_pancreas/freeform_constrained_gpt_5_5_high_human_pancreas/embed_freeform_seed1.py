#!/usr/bin/env python3
"""
Literature-grounded embedding for a healthy human pancreas single-cell matrix.

The scenario_config.jsonc for this run specifies human pancreas with no
restricted age, sex, condition, or cell-type subset. The embedding therefore
uses marker modules for the expected normal pancreatic epithelial, endocrine,
vascular, stromal, immune, neuroglial, and generic cycling/stress states.

No model is trained and no pretrained model is run. Coordinates are robustly
scaled literature module scores and a small set of signed contrasts.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Iterable

import anndata as ad
import numpy as np
import pandas as pd
from scipy import sparse


@dataclass(frozen=True)
class ModuleSpec:
    name: str
    genes: tuple[str, ...]


MODULES: tuple[ModuleSpec, ...] = (
    ModuleSpec(
        name="endocrine_pancreatic",
        # Baron et al., Cell Systems, 10.1016/j.cels.2016.08.011; Muraro et al.,
        # Cell Systems, 10.1016/j.cels.2016.09.002; Segerstolpe et al., Cell
        # Metabolism, 10.1016/j.cmet.2016.08.020.
        genes=("CHGA", "CHGB", "SCG2", "PCSK1N", "PAX6", "NEUROD1", "ISL1", "RFX6"),
    ),
    ModuleSpec(
        name="beta_cell",
        # Baron et al., Cell Systems, 10.1016/j.cels.2016.08.011; Muraro et al.,
        # Cell Systems, 10.1016/j.cels.2016.09.002; Segerstolpe et al., Cell
        # Metabolism, 10.1016/j.cmet.2016.08.020.
        genes=("INS", "IAPP", "MAFA", "NKX6-1", "PDX1", "SLC30A8", "PCSK1", "G6PC2", "UCN3"),
    ),
    ModuleSpec(
        name="alpha_cell",
        # Baron et al., Cell Systems, 10.1016/j.cels.2016.08.011; Muraro et al.,
        # Cell Systems, 10.1016/j.cels.2016.09.002; Segerstolpe et al., Cell
        # Metabolism, 10.1016/j.cmet.2016.08.020.
        genes=("GCG", "ARX", "IRX2", "MAFB", "PCSK2", "TTR", "ALDH1A1"),
    ),
    ModuleSpec(
        name="minor_endocrine_delta_pp_epsilon",
        # Baron et al., Cell Systems, 10.1016/j.cels.2016.08.011; Muraro et al.,
        # Cell Systems, 10.1016/j.cels.2016.09.002; Segerstolpe et al., Cell
        # Metabolism, 10.1016/j.cmet.2016.08.020.
        genes=("SST", "HHEX", "RBP4", "BCHE", "PPY", "PYY", "GHRL", "MLN"),
    ),
    ModuleSpec(
        name="acinar_exocrine",
        # Baron et al., Cell Systems, 10.1016/j.cels.2016.08.011; Muraro et al.,
        # Cell Systems, 10.1016/j.cels.2016.09.002; KEGG pancreatic secretion,
        # hsa04972.
        genes=(
            "PRSS1",
            "PRSS2",
            "CPA1",
            "CPB1",
            "CTRB1",
            "CTRB2",
            "CEL",
            "PNLIP",
            "CLPS",
            "CELA3A",
            "REG1A",
            "REG1B",
        ),
    ),
    ModuleSpec(
        name="ductal_centroacinar",
        # Baron et al., Cell Systems, 10.1016/j.cels.2016.08.011; Muraro et al.,
        # Cell Systems, 10.1016/j.cels.2016.09.002; KEGG pancreatic secretion,
        # hsa04972.
        genes=("KRT19", "KRT7", "KRT8", "KRT18", "SOX9", "CFTR", "SLC4A4", "SLC26A6", "MUC1", "EPCAM", "ANXA4"),
    ),
    ModuleSpec(
        name="endothelial_vascular_lymphatic",
        # Baron et al., Cell Systems, 10.1016/j.cels.2016.08.011; The Tabula
        # Sapiens Consortium, Science, 10.1126/science.abl4896; GO vasculature
        # development, GO:0001944.
        genes=("PECAM1", "VWF", "KDR", "CDH5", "EMCN", "PLVAP", "ESAM", "RAMP2", "ACKR1", "PROX1", "LYVE1", "PDPN"),
    ),
    ModuleSpec(
        name="stromal_stellate_pericyte_smooth_muscle",
        # Baron et al., Cell Systems, 10.1016/j.cels.2016.08.011; The Tabula
        # Sapiens Consortium, Science, 10.1126/science.abl4896; MSigDB
        # HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.
        genes=("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "SPARC", "PDGFRA", "PDGFRB", "RGS5", "MCAM", "ACTA2", "TAGLN", "DES"),
    ),
    ModuleSpec(
        name="immune_lymphoid",
        # The Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
        # Reactome adaptive immune system, R-HSA-1280215.
        genes=("PTPRC", "CD3D", "CD3E", "TRAC", "IL7R", "LTB", "NKG7", "GNLY", "MS4A1", "CD79A"),
    ),
    ModuleSpec(
        name="immune_myeloid",
        # The Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
        # Reactome innate immune system, R-HSA-168249.
        genes=("PTPRC", "LST1", "CSF1R", "AIF1", "FCER1G", "TYROBP", "C1QA", "C1QB", "LYZ", "CD14", "FCGR3A", "MS4A7"),
    ),
    ModuleSpec(
        name="neuroglial_schwann",
        # Baron et al., Cell Systems, 10.1016/j.cels.2016.08.011; The Tabula
        # Sapiens Consortium, Science, 10.1126/science.abl4896.
        genes=("SOX10", "S100B", "PLP1", "MPZ", "MBP", "NGFR", "GFRA3"),
    ),
    ModuleSpec(
        name="cycling",
        # The Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
        # MSigDB HALLMARK_E2F_TARGETS and HALLMARK_G2M_CHECKPOINT.
        genes=("MKI67", "TOP2A", "PCNA", "TYMS", "STMN1", "CENPF", "HMGB2", "UBE2C"),
    ),
    ModuleSpec(
        name="stress_interferon_unfolded_protein",
        # Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
        # MSigDB HALLMARK_UNFOLDED_PROTEIN_RESPONSE and
        # HALLMARK_INTERFERON_ALPHA_RESPONSE.
        genes=("HSPA1A", "HSPA1B", "DDIT3", "XBP1", "ATF3", "JUN", "FOS", "ISG15", "IFIT1", "IFIT3", "MX1"),
    ),
)


AXIS_NAMES: tuple[str, ...] = (
    "dim_0_pancreatic_endocrine",
    "dim_1_beta_cell",
    "dim_2_alpha_cell",
    "dim_3_delta_pp_epsilon",
    "dim_4_acinar_exocrine",
    "dim_5_ductal_centroacinar",
    "dim_6_endothelial",
    "dim_7_stromal_pericyte_smooth_muscle",
    "dim_8_lymphoid_immune",
    "dim_9_myeloid_immune",
    "dim_10_neuroglial_schwann",
    "dim_11_cycling",
    "dim_12_stress_interferon_upr",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build a literature-defined human pancreas cell-type/state embedding."
    )
    parser.add_argument("--input", required=True, help="Input AnnData .h5ad file.")
    parser.add_argument("--output", required=True, help="Output parquet path.")
    parser.add_argument("--seed", type=int, default=0, help="Run tag; not used by this deterministic scorer.")
    return parser.parse_args()


def canonical_symbol(symbol: object) -> str:
    text = str(symbol).strip().upper()
    if "-" in text:
        head, tail = text.rsplit("-", 1)
        if tail.isdigit():
            return head
    return text


def build_gene_index(var_names: Iterable[object]) -> dict[str, list[int]]:
    index: dict[str, list[int]] = {}
    for i, symbol in enumerate(var_names):
        index.setdefault(canonical_symbol(symbol), []).append(i)
    return index


def module_mean(x, indices: list[int], n_cells: int) -> np.ndarray:
    if not indices:
        return np.zeros(n_cells, dtype=np.float32)
    sub = x[:, indices]
    if sparse.issparse(sub):
        values = np.asarray(sub.mean(axis=1)).ravel()
    else:
        values = np.asarray(sub, dtype=np.float32).mean(axis=1)
    return values.astype(np.float32, copy=False)


def robust_z(values: np.ndarray) -> np.ndarray:
    values = values.astype(np.float32, copy=False)
    median = np.nanmedian(values)
    q25, q75 = np.nanpercentile(values, [25.0, 75.0])
    scale = (q75 - q25) / 1.349
    if not np.isfinite(scale) or scale < 1e-6:
        scale = np.nanstd(values)
    if not np.isfinite(scale) or scale < 1e-6:
        return np.zeros_like(values, dtype=np.float32)
    z = (values - median) / scale
    return np.clip(z, -8.0, 8.0).astype(np.float32)


def score_modules(adata: ad.AnnData) -> dict[str, np.ndarray]:
    x = adata.X
    if x is None:
        raise ValueError("AnnData.X must contain log1p-normalized expression.")
    gene_index = build_gene_index(adata.var_names)
    scores: dict[str, np.ndarray] = {}
    for module in MODULES:
        indices: list[int] = []
        for gene in module.genes:
            indices.extend(gene_index.get(canonical_symbol(gene), []))
        raw = module_mean(x, sorted(set(indices)), adata.n_obs)
        scores[module.name] = robust_z(raw)
    return scores


def build_embedding(scores: dict[str, np.ndarray]) -> np.ndarray:
    endocrine = scores["endocrine_pancreatic"]
    beta = scores["beta_cell"]
    alpha = scores["alpha_cell"]
    minor = scores["minor_endocrine_delta_pp_epsilon"]
    acinar = scores["acinar_exocrine"]
    ductal = scores["ductal_centroacinar"]
    endothelial = scores["endothelial_vascular_lymphatic"]
    stromal = scores["stromal_stellate_pericyte_smooth_muscle"]
    lymphoid = scores["immune_lymphoid"]
    myeloid = scores["immune_myeloid"]
    neuroglial = scores["neuroglial_schwann"]
    cycling = scores["cycling"]
    stress = scores["stress_interferon_unfolded_protein"]

    broad_non_endocrine = (acinar + ductal + endothelial + stromal + lymphoid + myeloid) / 6.0
    broad_non_acinar = (endocrine + ductal + endothelial + stromal + lymphoid + myeloid) / 6.0
    broad_non_ductal = (endocrine + acinar + endothelial + stromal + lymphoid + myeloid) / 6.0

    axes = (
        endocrine - 0.35 * broad_non_endocrine,
        beta - 0.30 * (alpha + minor) / 2.0,
        alpha - 0.30 * (beta + minor) / 2.0,
        minor - 0.25 * (beta + alpha) / 2.0,
        acinar - 0.25 * broad_non_acinar,
        ductal - 0.25 * broad_non_ductal,
        endothelial - 0.20 * (stromal + lymphoid + myeloid) / 3.0,
        stromal - 0.20 * (endothelial + acinar + ductal) / 3.0,
        lymphoid - 0.30 * myeloid,
        myeloid - 0.30 * lymphoid,
        neuroglial - 0.15 * (stromal + endothelial) / 2.0,
        cycling,
        stress,
    )
    return np.column_stack([robust_z(axis) for axis in axes]).astype(np.float32, copy=False)


def main() -> None:
    args = parse_args()
    adata = ad.read_h5ad(args.input)
    scores = score_modules(adata)
    embedding = build_embedding(scores)

    data: dict[str, object] = {"cell_id": pd.Index(adata.obs_names.astype(str), dtype="string")}
    for i in range(embedding.shape[1]):
        data[f"dim_{i}"] = embedding[:, i].astype(np.float32, copy=False)
    out = pd.DataFrame(data, columns=["cell_id", *[f"dim_{i}" for i in range(embedding.shape[1])]])
    out.to_parquet(args.output, index=False)


if __name__ == "__main__":
    main()
