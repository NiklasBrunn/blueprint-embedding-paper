#!/usr/bin/env python3
"""Literature-derived embedding for baseline human pancreas single-cell data."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


# Sources: Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004; Tabula Sapiens Consortium,
# Science, 10.1126/science.abl4896; KEGG hsa04972 PANCREATIC_SECRETION.
ACINAR_GENES = [
    "PRSS1", "PRSS2", "PRSS3", "CPA1", "CPA2", "CPB1", "CELA3A",
    "CELA3B", "CTRB1", "CTRB2", "CTRC", "PNLIP", "PNLIPRP1", "CLPS",
    "CEL", "PLA2G1B", "REG1A", "REG1B",
]

# Sources: Enge et al., Cell, 10.1016/j.cell.2017.09.004; Tabula Sapiens
# Consortium, Science, 10.1126/science.abl4896.
DUCTAL_GENES = [
    "KRT19", "KRT7", "KRT8", "KRT18", "EPCAM", "SOX9", "CFTR",
    "SLC4A4", "MUC1", "MUC6", "ANXA4", "TFF1", "TFF2", "S100A6",
    "SPP1", "MSLN",
]

# Sources: Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004; Tabula Sapiens Consortium,
# Science, 10.1126/science.abl4896; KEGG hsa04911 INSULIN_SECRETION.
BETA_GENES = [
    "INS", "IAPP", "MAFA", "PDX1", "NKX6-1", "SLC30A8", "PCSK1",
    "UCN3", "G6PC2", "HADH", "RBP4", "ADCYAP1", "DLK1", "NPTX2",
]

# Sources: Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004; Tabula Sapiens Consortium,
# Science, 10.1126/science.abl4896.
ALPHA_GENES = [
    "GCG", "ARX", "IRX2", "MAFB", "PCSK2", "TTR", "GC", "LOXL4",
    "DPP4", "TM4SF4", "SLC38A5", "CRYBA2",
]

# Sources: Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
# Enge et al., Cell, 10.1016/j.cell.2017.09.004; Tabula Sapiens Consortium,
# Science, 10.1126/science.abl4896.
DELTA_GENES = [
    "SST", "HHEX", "RBP4", "LEPR", "GHSR", "SLC18A2", "PCSK1",
    "RGS2", "LY6H",
]

# Sources: Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
GAMMA_PP_GENES = [
    "PPY", "PYY", "TSPAN8", "FAP", "SLITRK6", "TM4SF4", "PCSK1N",
]

# Sources: Segerstolpe et al., Cell Metabolism, 10.1016/j.cmet.2016.08.020;
# Tabula Sapiens Consortium, Science, 10.1126/science.abl4896.
EPSILON_GENES = [
    "GHRL", "MLN", "IAPP", "RBP4", "PAX6", "PCSK1",
]

# Sources: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# MSigDB HALLMARK_ANGIOGENESIS.
ENDOTHELIAL_GENES = [
    "PECAM1", "VWF", "KDR", "FLT1", "ENG", "CDH5", "EMCN", "RAMP2",
    "CLDN5", "ESAM", "PLVAP", "CLEC14A", "SOX17",
]

# Sources: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# GO:0006939 SMOOTH_MUSCLE_CONTRACTION.
PERICYTE_SMOOTH_GENES = [
    "RGS5", "PDGFRB", "CSPG4", "MCAM", "NOTCH3", "ACTA2", "TAGLN",
    "MYH11", "MYL9", "COL4A1", "COL4A2", "NDUFA4L2",
]

# Sources: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# MSigDB HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.
STELLATE_FIBROBLAST_GENES = [
    "COL1A1", "COL1A2", "COL3A1", "COL6A1", "COL6A2", "DCN", "LUM",
    "PDGFRA", "MMP2", "FAP", "VIM", "THY1", "CXCL12", "CFD",
]

# Sources: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# MSigDB HALLMARK_INFLAMMATORY_RESPONSE; MSigDB HALLMARK_COMPLEMENT.
MYELOID_GENES = [
    "PTPRC", "LST1", "TYROBP", "AIF1", "C1QA", "C1QB", "C1QC", "CD14",
    "CD68", "FCGR3A", "MS4A7", "LYZ", "CSF1R", "SPI1",
]

# Sources: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# MSigDB HALLMARK_IL2_STAT5_SIGNALING.
LYMPHOID_GENES = [
    "PTPRC", "CD3D", "CD3E", "CD2", "TRAC", "IL7R", "LTB", "NKG7",
    "GNLY", "GZMB", "MS4A1", "CD79A", "CD79B", "BANK1",
]

# Sources: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# MSigDB HALLMARK_INFLAMMATORY_RESPONSE.
MAST_GENES = [
    "TPSAB1", "TPSB2", "CPA3", "KIT", "MS4A2", "HDC", "GATA2", "HPGDS",
]

# Sources: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# GO:0007272 ENSHEATHMENT_OF_NEURONS.
SCHWANN_NEURAL_GENES = [
    "SOX10", "S100B", "PLP1", "MPZ", "GFRA3", "NGFR", "NRXN1", "SCN7A",
    "NCAM1", "GPM6B",
]

# Sources: Enge et al., Cell, 10.1016/j.cell.2017.09.004; Segerstolpe et al.,
# Cell Metabolism, 10.1016/j.cmet.2016.08.020; MSigDB
# HALLMARK_UNFOLDED_PROTEIN_RESPONSE.
SECRETORY_STRESS_GENES = [
    "HSPA5", "XBP1", "DDIT3", "ATF4", "HERPUD1", "DNAJB9", "PDIA4",
    "P4HB", "FOS", "FOSB", "JUN", "JUND", "HSPA1A", "HSPA1B", "ATF3",
    "CLU",
]

# Sources: Tabula Sapiens Consortium, Science, 10.1126/science.abl4896;
# MSigDB HALLMARK_G2M_CHECKPOINT; KEGG hsa04110 CELL_CYCLE.
CELL_CYCLE_GENES = [
    "MKI67", "TOP2A", "UBE2C", "CENPF", "NUSAP1", "PCNA", "MCM5",
    "TYMS", "STMN1", "HMGB2", "BIRC5", "AURKB",
]


SIGNATURES = {
    "acinar": ACINAR_GENES,
    "ductal": DUCTAL_GENES,
    "beta": BETA_GENES,
    "alpha": ALPHA_GENES,
    "delta": DELTA_GENES,
    "gamma_pp": GAMMA_PP_GENES,
    "epsilon": EPSILON_GENES,
    "endothelial": ENDOTHELIAL_GENES,
    "pericyte_smooth": PERICYTE_SMOOTH_GENES,
    "stellate_fibroblast": STELLATE_FIBROBLAST_GENES,
    "myeloid": MYELOID_GENES,
    "lymphoid": LYMPHOID_GENES,
    "mast": MAST_GENES,
    "schwann_neural": SCHWANN_NEURAL_GENES,
    "secretory_stress": SECRETORY_STRESS_GENES,
    "cell_cycle": CELL_CYCLE_GENES,
}


# Fourteen axes are used because adult human pancreas spans endocrine,
# exocrine, vascular, mesenchymal, immune, neural, and stress/cycling states;
# fewer axes would collapse rare endocrine and non-epithelial compartments.
AXES = [
    ("endocrine_islet_vs_other", {
        "beta": 0.8, "alpha": 0.8, "delta": 0.7, "gamma_pp": 0.7,
        "epsilon": 0.7, "acinar": -0.35, "ductal": -0.35,
        "endothelial": -0.2, "stellate_fibroblast": -0.2,
        "myeloid": -0.2, "lymphoid": -0.2,
    }),
    ("acinar_vs_ductal", {"acinar": 1.0, "ductal": -1.0}),
    ("beta_vs_alpha", {"beta": 1.0, "alpha": -1.0}),
    ("delta_identity", {"delta": 1.0, "beta": -0.2, "alpha": -0.2}),
    ("gamma_pp_vs_epsilon", {"gamma_pp": 1.0, "epsilon": -1.0}),
    ("endothelial_vascular", {"endothelial": 1.0}),
    ("pericyte_smooth_muscle", {"pericyte_smooth": 1.0}),
    ("stellate_fibroblast", {"stellate_fibroblast": 1.0}),
    ("myeloid_macrophage", {"myeloid": 1.0, "lymphoid": -0.25}),
    ("lymphoid", {"lymphoid": 1.0, "myeloid": -0.25}),
    ("mast_cell", {"mast": 1.0}),
    ("schwann_neural", {"schwann_neural": 1.0}),
    ("secretory_stress", {"secretory_stress": 1.0}),
    ("cell_cycle", {"cell_cycle": 1.0}),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create a literature-derived embedding for human pancreas cells."
    )
    parser.add_argument("--input", required=True, help="Input .h5ad file")
    parser.add_argument("--output", required=True, help="Output parquet path")
    parser.add_argument("--seed", type=int, required=True, help="Run tag; ignored")
    return parser.parse_args()


def dense_view(matrix):
    import numpy as np
    from scipy import sparse

    if sparse.issparse(matrix):
        matrix = matrix.toarray()
    arr = np.asarray(matrix, dtype=np.float32)
    if arr.ndim == 1:
        arr = arr.reshape(-1, 1)
    return arr


def gene_index(var_names) -> dict[str, int]:
    index: dict[str, int] = {}
    for i, gene in enumerate(var_names):
        key = str(gene).strip().upper()
        if key and key not in index:
            index[key] = i
    return index


def score_signature(x, index: dict[str, int], genes: list[str], n_cells: int, name: str):
    import numpy as np

    cols = [index[g.upper()] for g in genes if g.upper() in index]
    if not cols:
        print(f"warning: no genes matched for signature {name}", file=sys.stderr)
        return np.zeros(n_cells, dtype=np.float32)

    block = dense_view(x[:, cols])
    mean = block.mean(axis=0, dtype=np.float64)
    std = block.std(axis=0, dtype=np.float64)
    std = np.where(std < 1e-6, 1.0, std)
    z = (block - mean) / std
    z = np.clip(z, -5.0, 5.0)
    return z.mean(axis=1, dtype=np.float64).astype(np.float32)


def robust_unit(values):
    import numpy as np

    values = np.asarray(values, dtype=np.float32)
    median = float(np.nanmedian(values))
    q25, q75 = np.nanpercentile(values, [25.0, 75.0])
    scale = float((q75 - q25) / 1.349) if q75 > q25 else float(np.nanstd(values))
    if not np.isfinite(scale) or scale < 1e-6:
        return np.zeros(values.shape[0], dtype=np.float32)
    centered = (values - median) / (2.5 * scale)
    return np.tanh(centered).astype(np.float32)


def run_embedding(input_path: Path, output_path: Path) -> None:
    import anndata as ad
    import numpy as np
    import pandas as pd

    adata = ad.read_h5ad(input_path)
    n_cells = int(adata.n_obs)
    index = gene_index(adata.var_names)
    x = adata.X

    scores = {
        name: score_signature(x, index, genes, n_cells, name)
        for name, genes in SIGNATURES.items()
    }

    coords = []
    for _, weights in AXES:
        raw = np.zeros(n_cells, dtype=np.float32)
        for name, weight in weights.items():
            raw += np.float32(weight) * scores[name]
        coords.append(robust_unit(raw))

    data = {"cell_id": np.asarray(adata.obs_names, dtype=str)}
    for i, coord in enumerate(coords):
        data[f"dim_{i}"] = coord.astype(np.float32, copy=False)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(data).to_parquet(output_path, index=False)


def main() -> None:
    args = parse_args()
    run_embedding(Path(args.input), Path(args.output))


if __name__ == "__main__":
    main()
