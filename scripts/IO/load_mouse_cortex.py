"""
Load the integrated three-source adult-mouse-isocortex dataset.

Combines:

    - Tasic et al. (2016) "Adult mouse cortical cell taxonomy revealed by
      single cell transcriptomics", Nature Neuroscience,
      doi:10.1038/nn.4216 — VISp primary visual cortex, SMART-Seq
      (full-length), ~1.5k cells post-QC. GEO GSE71585.
    - Tasic et al. (2018) "Shared and distinct transcriptomic cell types
      across neocortical areas", Nature,
      doi:10.1038/s41586-018-0654-5 — VISp + ALM, SMART-Seq v4
      (full-length), ~24k cells, ~133 transcriptomic types. GEO GSE115746.
    - Zeisel et al. (2018) "Molecular Architecture of the Mouse Nervous
      System", Cell, doi:10.1016/j.cell.2018.06.021 — whole CNS+PNS atlas,
      10x Genomics 3' droplet, ~500k cells, filtered here to the
      isocortex slice (~70-160k depending on the non-neuronal allowlist).
      l5_all.loom from mousebrain.org.

The dataset deliberately mixes two distinct protocols (full-length
SMART-Seq* vs 3' droplet 10x) on overlapping biological scope so that
downstream metrics can assess cross-protocol integratability of
literature-grounded blueprint embeddings without explicit batch
correction.

obs columns harmonized across all three studies:
    study              — "tasic2016" / "tasic2018" / "zeisel2018"; BATCH_KEY
    protocol           — "smartseq" / "smartseq_v4" / "10x_v2"
    area               — "VISp" / "ALM" / "cortex" (Zeisel cortex slice has
                         no per-area resolution)
    celltype_tasic2016 — Tasic 2016 native subclass; NaN for other studies
    celltype_tasic2018 — Tasic 2018 native subclass; NaN for other studies
    celltype_zeisel    — Zeisel native ClusterName; NaN for other studies
    celltype_harmonized — Tasic-2018-subclass canonical vocabulary
                          (~24 labels: L2/3 IT, L4, L5 IT, L5 PT, L5 NP,
                          L6 IT, L6 CT, L6b, Pvalb, Sst, Sst Chodl, Vip,
                          Lamp5, Sncg, Meis2, Astro, Oligo, OPC, Microglia,
                          Endo, Peri, SMC, VLMC, Macrophage). The canonical
                          LABEL_KEY.
    class_harmonized   — coarse 3-bucket class derived from
                         celltype_harmonized: Glutamatergic / GABAergic /
                         Non-neuronal. Useful as the root level of a
                         hierarchical gating tree.

Label-cleaning policy: only cells with a usable harmonized cell-type label are
emitted. Cells with no annotation (nan / "") and Zeisel's generic "Neurons"
label (real cells but no class-resolvable subtype) are dropped at integration
time, so the LABEL_KEY ground truth carries no spurious "Unknown" cluster. This
mirrors the standalone Tasic loader (load_mouse_brain.py), which drops
"Unclassified" cells. Rare but genuine fine types (e.g. Cajal-Retzius "CR") are
kept and bucketed to a coarse class.

Outputs:
    - data/raw/raw_mouse_cortex.h5ad   (raw counts in adata.X and adata.layers["counts"])
    - data/raw/raw_mouse_cortex.json   (provenance + per-study counts + harmonization audit)

Cache (shared across reruns; cleared by --force-download):
    - data/raw/_mousebrain_cortex/tasic2018/{counts,metadata}.csv.gz
    - data/raw/_mousebrain_cortex/zeisel/l5_all.loom
    - data/raw/_mousebrain_cortex/zeisel/zeisel_cortex_filtered.h5ad  (filtered subset)

The Tasic 2016 sources share load_mouse_brain.py's data/raw/GSE71585_* cache;
--force-download re-fetches them.

Usage:
    python scripts/IO/load_mouse_cortex.py [--overwrite] [--force-download]
"""
from __future__ import annotations

import argparse
import json
import re
import urllib.request
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


REPO_ROOT = Path(__file__).resolve().parents[2]
RAW_DIR = REPO_ROOT / "data" / "raw"
RAW_OUTPUT = RAW_DIR / "raw_mouse_cortex.h5ad"
METADATA_OUTPUT = RAW_DIR / "raw_mouse_cortex.json"
CACHE_DIR = RAW_DIR / "_mousebrain_cortex"

DATASET_NAME = "mouse_cortex"
CITATION = (
    "Tasic et al. (2016) Nature Neuroscience, doi:10.1038/nn.4216; "
    "Tasic et al. (2018) Nature, doi:10.1038/s41586-018-0654-5; "
    "Zeisel et al. (2018) Cell, doi:10.1016/j.cell.2018.06.021"
)
SPECIES = "Mus musculus"
TISSUE = "isocortex"
PROTOCOL = "scRNA-seq cross-protocol (SMART-Seq + SMART-Seq v4 + 10x Genomics 3')"
BATCH_KEY = "study"
LABEL_KEY = "celltype_harmonized"

# --- Source URLs ---------------------------------------------------------

GEO_BASE = "https://www.ncbi.nlm.nih.gov/geo/download/?acc={acc}&format=file&file={file}"

TASIC2016_COUNTS_FILENAME = "GSE71585_RefSeq_counts.csv.gz"
TASIC2016_CLUSTER_FILENAME = "GSE71585_Clustering_Results.csv.gz"

TASIC2018_COUNTS_FILENAME = "GSE115746_cells_exon_counts.csv.gz"
TASIC2018_META_FILENAME = "GSE115746_complete_metadata_28706-cells.csv.gz"

ZEISEL_LOOM_URL = "https://storage.googleapis.com/linnarsson-lab-loom/l5_all.loom"
ZEISEL_LOOM_FILENAME = "l5_all.loom"

# --- Zeisel isocortex filter ---------------------------------------------
# Region=="Cortex" in the level-5 atlas tags Neurons only; non-neurons are
# tagged CNS-wide and cannot be cleanly attributed to a specific brain region
# from this file. We therefore restrict the Zeisel contribution to its cortex
# NEURONS and rely on Tasic 2018 for non-neuronal coverage (Tasic 2018 covers
# Astro/Oligo/OPC/Endo/Micro/Peri/SMC/VLMC in VISp+ALM).
#
# Region containing "cortex" (case-insensitive substring) is the broad filter,
# but Zeisel's "Cortex" Region also includes paleocortex (Piriform), allocortex
# (Entorhinal superficial, Anterior olfactory), and the hippocampal formation
# (Subiculum). These are NOT isocortex and are excluded by Probable_location.

ZEISEL_NONISOCORTEX_LOCATION_PATTERNS = (
    "piriform",
    "subiculum",
    "entorhinal",
    "anterior olfactory",
)

# --- Harmonized vocabulary -----------------------------------------------
# The Tasic-2018 subclass list (area suffixes stripped) is the canonical
# label set. Native labels that fall outside this vocabulary pass through
# unchanged and are surfaced in the harmonization audit in the JSON sidecar.

HARMONIZED_GLUTAMATERGIC = {
    "L2/3 IT", "L4", "L5 IT", "L5 PT", "L5 NP", "L6 IT", "L6 CT", "L6b",
    # Cajal-Retzius: a rare, real excitatory (reelin+) type kept with its own
    # fine `celltype_harmonized` label "CR"; bucketed here so it carries a class
    # rather than falling into "Unknown".
    "CR",
}
HARMONIZED_GABAERGIC = {
    "Pvalb", "Sst", "Sst Chodl", "Vip", "Lamp5", "Sncg", "Meis2",
}
HARMONIZED_NONNEURONAL = {
    "Astro", "Oligo", "OPC", "Microglia", "Endo", "Peri", "SMC", "VLMC",
    "Macrophage",
}
HARMONIZED_VOCAB = (
    HARMONIZED_GLUTAMATERGIC | HARMONIZED_GABAERGIC | HARMONIZED_NONNEURONAL
)

# celltype_harmonized values (lowercased) that are NOT a usable ground-truth
# label: missing annotations and generic / unresolvable coarse labels (notably
# Zeisel's "Neurons"). Cells carrying these are dropped at integration time.
UNUSABLE_LABELS = {"", "nan", "none", "na", "unknown", "unclassified", "neurons"}

# Tasic 2016 subclass (derived from primary_type leading token in
# load_mouse_brain.py) → harmonized. Lower-cased keys; unmapped → pass through.
_HARMONIZE_TASIC2016 = {
    "l2": "L2/3 IT", "l2/3": "L2/3 IT",
    "l4": "L4",
    "l5": "L5 IT", "l5a": "L5 IT",
    "l5b": "L5 PT",
    "l6a": "L6 CT",
    "l6b": "L6b",
    "vip": "Vip",
    "sst": "Sst",
    "pvalb": "Pvalb",
    "ndnf": "Lamp5",
    "igtp": "Lamp5",
    "smad3": "Sncg",
    "sncg": "Sncg",
    "lamp5": "Lamp5",
    "astro": "Astro",
    "opc": "OPC",
    "oligo": "Oligo",
    "micro": "Microglia",
    "endo": "Endo",
    "smc": "SMC",
}

# Tasic 2018 subclass (with area suffix stripped) → harmonized. Lower-cased keys.
_HARMONIZE_TASIC2018 = {
    "l2/3 it": "L2/3 IT",
    "l4": "L4",
    "l5 it": "L5 IT",
    "l5 pt": "L5 PT",
    "l5/6 np": "L5 NP", "l5 np": "L5 NP", "np": "L5 NP",
    "l6 it": "L6 IT",
    "l6 ct": "L6 CT",
    "l6b": "L6b",
    "pvalb": "Pvalb",
    "sst": "Sst",
    "sst chodl": "Sst Chodl",
    "vip": "Vip",
    "lamp5": "Lamp5",
    "sncg": "Sncg",
    "serpinf1": "Sncg",
    "meis2": "Meis2",
    "astro": "Astro",
    "oligo": "Oligo",
    "opc": "OPC",
    "endo": "Endo",
    "peri": "Peri",
    "smc": "SMC",
    "vlmc": "VLMC",
    "microglia": "Microglia",
    "macrophage": "Macrophage",
}

# Zeisel Class → harmonized for non-neurons. The level-5 atlas uses coarse
# Class names (Astrocytes / Oligos / Immune / Vascular / Ependymal /
# PeripheralGlia / Neurons). Cortex non-neurons aren't isolable from this
# file (see ZEISEL_NONISOCORTEX_LOCATION_PATTERNS comment), so this map is
# kept for completeness but the loader only emits cortex Neurons.
_HARMONIZE_ZEISEL_CLASS = {
    "astrocytes": "Astro",
    "oligos": "Oligo",
    "immune": "Microglia",
    "vascular": "Endo",
    "ependymal": "Other",
    "peripheralglia": "Other",
}

# Zeisel EXCITATORY-neuron Probable_location → harmonized subclass. The
# level-5 atlas tags excitatory cortex clusters as e.g.
# "Cortical pyramidal layer 4" / "Cortex pyramidal layer 6b" /
# "Cingulate/Retrosplenial area, layer 5". First case-insensitive substring
# match wins; ordered most-specific first so "layer 6b" beats "layer 6".
_HARMONIZE_ZEISEL_EXC_LOCATION = (
    ("layer 6b", "L6b"),
    ("layer 2/3", "L2/3 IT"),
    ("layer 2-3", "L2/3 IT"),
    ("layer 2", "L2/3 IT"),
    ("layer 3", "L2/3 IT"),
    ("layer 4", "L4"),
    # Lateral cortex L6 (gustatory/barrel-field/auditory) is IT, not CT.
    # Must come BEFORE the generic "layer 6" rule.
    ("lateral cortex", "L6 IT"),
    ("layer 5", "L5 IT"),
    ("layer 6", "L6 CT"),
)

# Zeisel INHIBITORY-interneuron Description → harmonized subclass. The
# level-5 atlas annotates interneurons by functional/morphological class
# (e.g. "Axo-axonic, cortex/hippocampus" → Pvalb chandelier;
# "Sleep-active, long-range projection interneurons" → Sst Chodl).
# Ordered most-specific first.
_HARMONIZE_ZEISEL_INH_DESCRIPTION = (
    ("sleep-active, long-range projection", "Sst Chodl"),
    ("hippocamposeptal projection", "Sst"),
    ("axo-axonic", "Pvalb"),
    ("basket and bistratified", "Pvalb"),
    ("ivy and mge-derived neurogliaform", "Lamp5"),
    ("cge-derived neurogliaform", "Lamp5"),
    ("interneuron-selective", "Vip"),
    ("border cck", "Vip"),
    ("cck interneurons", "Vip"),
    ("pvalb", "Pvalb"),
    ("parvalbumin", "Pvalb"),
    ("sst chodl", "Sst Chodl"),
    ("somatostatin", "Sst"),
    ("vip", "Vip"),
    ("lamp5", "Lamp5"),
    ("ndnf", "Lamp5"),
    ("sncg", "Sncg"),
    ("meis2", "Meis2"),
)


def _harmonize_subclass_tasic2016(label: object) -> str | float:
    if not isinstance(label, str) or not label.strip():
        return np.nan
    return _HARMONIZE_TASIC2016.get(label.strip().lower(), label.strip())


def _harmonize_subclass_tasic2018(label: object) -> str | float:
    if not isinstance(label, str) or not label.strip():
        return np.nan
    # Strip area suffix ("L2/3 IT VISp" → "L2/3 IT"); area is captured separately.
    s = label.strip()
    for suffix in (" VISp", " ALM"):
        if s.endswith(suffix):
            s = s[: -len(suffix)]
            break
    return _HARMONIZE_TASIC2018.get(s.lower(), s)


def _harmonize_zeisel(
    cls: object,
    taxonomy_rank4: object = "",
    probable_location: object = "",
    description: object = "",
) -> str | float:
    """Map a Zeisel-atlas row to the harmonized subclass.

    - Non-neurons (Class != "Neurons") map by Class via _HARMONIZE_ZEISEL_CLASS.
    - Excitatory neurons (TaxonomyRank4 contains "excitatory") map by
      Probable_location via _HARMONIZE_ZEISEL_EXC_LOCATION (layer-based).
    - Interneurons (TaxonomyRank4 contains "inhibitory" or "interneuron") map
      by Description via _HARMONIZE_ZEISEL_INH_DESCRIPTION (functional class).
    - Anything else returns the Class label unchanged, surfacing in the audit.
    """
    if not isinstance(cls, str):
        return np.nan
    cls_key = cls.strip().lower()
    tax4 = str(taxonomy_rank4 or "").lower()
    loc = str(probable_location or "").lower()
    desc = str(description or "").lower()

    if cls_key != "neurons":
        return _HARMONIZE_ZEISEL_CLASS.get(cls_key, cls.strip())

    if "excitatory" in tax4:
        for substring, target in _HARMONIZE_ZEISEL_EXC_LOCATION:
            if substring in loc:
                return target
        return cls.strip()  # unmapped excitatory — surfaces in audit

    if "inhibitory" in tax4 or "interneuron" in tax4:
        for substring, target in _HARMONIZE_ZEISEL_INH_DESCRIPTION:
            if substring in desc:
                return target
        return cls.strip()  # unmapped interneuron — surfaces in audit

    return cls.strip()


def _derive_class_harmonized(label: object) -> str:
    if not isinstance(label, str):
        return "Unknown"
    if label in HARMONIZED_GLUTAMATERGIC:
        return "Glutamatergic"
    if label in HARMONIZED_GABAERGIC:
        return "GABAergic"
    if label in HARMONIZED_NONNEURONAL:
        return "Non-neuronal"
    return "Unknown"


def _download(url: str, dest: Path, label: str = "") -> Path:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        print(f"[skip] {dest}")
        return dest
    print(f"[download] {label or dest.name} <- {url}")
    urllib.request.urlretrieve(url, dest)
    if dest.stat().st_size == 0:
        dest.unlink(missing_ok=True)
        raise RuntimeError(f"Downloaded file is empty: {url}")
    return dest


# --- Tasic 2016 ----------------------------------------------------------

def _load_tasic2016(force_download: bool = False) -> ad.AnnData:
    """Read the GEO GSE71585 CSVs the existing load_mouse_brain.py caches.

    Reuses the same parsing as load_mouse_brain.py (no shared helper to keep
    this loader self-contained).
    """
    counts_path = RAW_DIR / TASIC2016_COUNTS_FILENAME
    cluster_path = RAW_DIR / TASIC2016_CLUSTER_FILENAME

    counts_url = GEO_BASE.format(acc="GSE71585", file=TASIC2016_COUNTS_FILENAME)
    cluster_url = GEO_BASE.format(acc="GSE71585", file=TASIC2016_CLUSTER_FILENAME)

    if force_download or not counts_path.exists():
        _download(counts_url, counts_path, "tasic2016 counts")
    if force_download or not cluster_path.exists():
        _download(cluster_url, cluster_path, "tasic2016 metadata")

    counts_df = pd.read_csv(counts_path, compression="gzip")
    gene_names = pd.Index(counts_df.pop("gene").astype(str), name="gene")
    sample_ids = counts_df.columns.astype(str).tolist()

    meta = pd.read_csv(cluster_path, compression="gzip")
    meta = meta.rename(columns={"sample_title": "sample_id"}).set_index("sample_id")
    meta = meta.loc[sample_ids].copy()
    meta.index = meta.index.astype(str)
    meta.index.name = "cell_id"
    meta["sample_id"] = meta.index

    # Drop Unclassified cells (same filter as load_mouse_brain.py).
    keep = meta["broad_type"].astype(str) != "Unclassified"
    n_dropped = int((~keep).sum())
    meta = meta.loc[keep].copy()
    counts_df = counts_df.loc[:, meta.index]
    print(f"[tasic2016] dropped {n_dropped} Unclassified cells; kept {meta.shape[0]}")

    # Derive subclass from primary_type's leading token (same as load_mouse_brain).
    primary = meta["primary_type"].fillna("Unclassified").astype(str)
    meta["subclass_native"] = primary.str.split().str[0].replace({"nan": "Unclassified"})

    matrix = counts_df.to_numpy(dtype=np.float32, copy=False).T
    obs = pd.DataFrame(index=meta.index)
    obs["study"] = "tasic2016"
    obs["protocol"] = "smartseq"
    obs["area"] = "VISp"
    obs["celltype_tasic2016"] = meta["subclass_native"].astype(str).values
    obs["celltype_tasic2018"] = ""
    obs["celltype_zeisel"] = ""
    obs["celltype_harmonized"] = meta["subclass_native"].apply(
        _harmonize_subclass_tasic2016
    ).values
    obs["sample_id"] = meta["sample_id"].astype(str).values

    adata = ad.AnnData(
        X=sp.csr_matrix(matrix),
        obs=obs,
        var=pd.DataFrame(index=gene_names),
    )
    adata.layers["counts"] = adata.X.copy()
    return adata


# --- Tasic 2018 ----------------------------------------------------------

def _load_tasic2018(force_download: bool = False) -> ad.AnnData:
    """Download GSE115746 supplementary CSVs and parse into AnnData."""
    cache = CACHE_DIR / "tasic2018"
    cache.mkdir(parents=True, exist_ok=True)
    counts_path = cache / TASIC2018_COUNTS_FILENAME
    meta_path = cache / TASIC2018_META_FILENAME

    counts_url = GEO_BASE.format(acc="GSE115746", file=TASIC2018_COUNTS_FILENAME)
    meta_url = GEO_BASE.format(acc="GSE115746", file=TASIC2018_META_FILENAME)

    if force_download or not counts_path.exists():
        _download(counts_url, counts_path, "tasic2018 counts")
    if force_download or not meta_path.exists():
        _download(meta_url, meta_path, "tasic2018 metadata")

    meta = pd.read_csv(meta_path, compression="gzip")
    # GEO supplementary files use sample_name as the canonical cell ID.
    if "sample_name" in meta.columns:
        meta = meta.set_index("sample_name")
    elif "sample_id" in meta.columns:
        meta = meta.set_index("sample_id")
    else:
        meta = meta.set_index(meta.columns[0])
    meta.index = meta.index.astype(str)
    meta.index.name = "cell_id"

    # Orientation: Tasic 2018 ships the counts as genes-in-rows / cells-in-cols
    # with no header on the gene column (pandas reads it as "Unnamed: 0").
    # Other vendors may ship cells-in-rows. Detect by looking at the first
    # column's name and a value sample.
    head = pd.read_csv(counts_path, compression="gzip", nrows=2)
    first_col_name = str(head.columns[0])
    first_col_sample = str(head.iloc[0, 0]) if len(head) else ""
    second_col_name = str(head.columns[1]) if head.shape[1] > 1 else ""
    cell_id_re = re.compile(r"^[A-Z0-9]+_\d{6}_\d{3}_[A-Z]\d{2}$|^[A-Z]+\.[A-Z0-9.]+$")
    cells_first = first_col_name.lower() in ("sample_name", "cell", "cell_id") or bool(
        cell_id_re.match(first_col_sample)
    )

    if cells_first:
        # Cells in rows, genes in columns.
        counts_df = pd.read_csv(counts_path, compression="gzip")
        counts_df = counts_df.set_index(counts_df.columns[0])
        counts_df.index = counts_df.index.astype(str)
        gene_names = counts_df.columns.astype(str).tolist()
        cell_ids = counts_df.index.tolist()
        matrix = counts_df.to_numpy(dtype=np.float32, copy=False)
    else:
        # Genes in rows, cells in columns (Tasic 2018 + Tasic 2016 layout).
        # Header[0] is either a gene-label column ("gene", "gene_symbol") or
        # unnamed (pandas: "Unnamed: 0").
        counts_df = pd.read_csv(counts_path, compression="gzip")
        gene_names = counts_df.pop(counts_df.columns[0]).astype(str).tolist()
        cell_ids = counts_df.columns.astype(str).tolist()
        matrix = counts_df.to_numpy(dtype=np.float32, copy=False).T

    # Align metadata to count matrix.
    meta = meta.reindex(cell_ids).copy()

    # Drop technical-control / QC-fail cells. The Tasic 2018 metadata uses
    # `cell_class` (not `class`) and includes spike-in / no-template / batch
    # rows alongside the real Glutamatergic / GABAergic / Non-Neuronal /
    # Endothelial cells. Drop everything that isn't a real cell class.
    valid_cell_classes = {"glutamatergic", "gabaergic", "non-neuronal", "endothelial"}
    if "cell_class" in meta.columns:
        cls_lower = meta["cell_class"].astype(str).str.lower()
        keep_mask = cls_lower.isin(valid_cell_classes)
    else:
        keep_mask = pd.Series(True, index=meta.index)
    n_dropped = int((~keep_mask).sum())
    meta = meta.loc[keep_mask].copy()
    matrix = matrix[keep_mask.to_numpy(), :]
    cell_ids = [c for c, k in zip(cell_ids, keep_mask) if k]
    print(
        f"[tasic2018] dropped {n_dropped} controls / QC-fail / unmapped-class cells; "
        f"kept {len(cell_ids)}"
    )

    subclass_native = meta.get("cell_subclass", pd.Series([np.nan] * len(meta), index=meta.index))
    brain_region = meta.get("dissected_region", pd.Series([np.nan] * len(meta), index=meta.index))

    obs = pd.DataFrame(index=pd.Index(cell_ids, name="cell_id"))
    obs["study"] = "tasic2018"
    obs["protocol"] = "smartseq_v4"
    obs["area"] = brain_region.astype(str).values
    obs["celltype_tasic2016"] = ""
    obs["celltype_tasic2018"] = subclass_native.astype(str).values
    obs["celltype_zeisel"] = ""
    obs["celltype_harmonized"] = subclass_native.apply(_harmonize_subclass_tasic2018).values
    obs["sample_id"] = obs.index

    adata = ad.AnnData(
        X=sp.csr_matrix(matrix),
        obs=obs,
        var=pd.DataFrame(index=pd.Index(gene_names, name="gene")),
    )
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    return adata


# --- Zeisel 2018 ---------------------------------------------------------

def _zeisel_isocortex_mask(region: pd.Series, probable_location: pd.Series) -> np.ndarray:
    """Bool mask: Region contains "cortex" AND Probable_location is not paleo/archicortex."""
    region_lc = region.astype(str).str.lower()
    loc_lc = probable_location.astype(str).str.lower()
    in_cortex_region = region_lc.str.contains("cortex", regex=False, na=False)
    is_nonisocortex = loc_lc.apply(
        lambda s: any(p in s for p in ZEISEL_NONISOCORTEX_LOCATION_PATTERNS)
    )
    return (in_cortex_region & ~is_nonisocortex).to_numpy()


def _load_zeisel(force_download: bool = False) -> ad.AnnData:
    """Load the Zeisel atlas and filter to the cortex Neuron slice.

    The mousebrain.org canonical source is l5_all.loom (HDF5 with /matrix +
    /col_attrs/* + /row_attrs/*). This loader also accepts a pre-converted
    AnnData h5ad with the same field names in obs (Region, Class,
    TaxonomyRank4, Probable_location, Description, ClusterName); format
    is auto-detected from the file's HDF5 top-level keys.

    The filtered isocortex-Neurons subset is cached as a small h5ad next
    to the input; subsequent reruns skip the parse.
    """
    cache = CACHE_DIR / "zeisel"
    cache.mkdir(parents=True, exist_ok=True)
    src_path = cache / ZEISEL_LOOM_FILENAME
    filtered_path = cache / "zeisel_cortex_filtered.h5ad"

    if force_download and filtered_path.exists():
        filtered_path.unlink()
    if filtered_path.exists():
        print(f"[skip] {filtered_path}")
        return ad.read_h5ad(filtered_path)

    if not src_path.exists():
        if force_download:
            _download(ZEISEL_LOOM_URL, src_path, "zeisel l5_all.loom (~2.5 GB)")
        else:
            raise FileNotFoundError(
                f"{src_path} not found. The canonical mousebrain.org URL is "
                f"intermittent; either pass --force-download to retry the "
                f"canonical URL, or place the Zeisel atlas at this path "
                f"manually (l5_all.loom or a pre-converted AnnData h5ad with "
                f"obs.Region / Class / TaxonomyRank4 / Probable_location / "
                f"Description / ClusterName)."
            )

    import h5py
    with h5py.File(src_path, "r") as f:
        top_keys = set(f.keys())
    is_anndata = {"X", "obs", "var"}.issubset(top_keys)
    is_loom = {"matrix", "row_attrs", "col_attrs"}.issubset(top_keys)
    if not (is_anndata or is_loom):
        raise RuntimeError(
            f"Unrecognised HDF5 layout in {src_path}: top_keys={sorted(top_keys)}; "
            f"expected AnnData (X/obs/var) or Loom (matrix/row_attrs/col_attrs)."
        )

    if is_anndata:
        adata = _load_zeisel_from_anndata(src_path)
    else:
        adata = _load_zeisel_from_loom(src_path)

    print(f"[zeisel] writing filtered cache: {filtered_path}")
    adata.write_h5ad(filtered_path)
    return adata


def _load_zeisel_from_anndata(src_path: Path) -> ad.AnnData:
    print(f"[zeisel] reading {src_path.name} as AnnData h5ad")
    a = ad.read_h5ad(src_path)
    # Coerce expected obs columns to str (may come back as Categorical).
    obs_cols = ("Region", "Tissue", "Class", "Subclass", "ClusterName",
                "Description", "Probable_location", "TaxonomyRank1",
                "TaxonomyRank2", "TaxonomyRank3", "TaxonomyRank4", "SampleID")
    for col in obs_cols:
        if col in a.obs.columns:
            a.obs[col] = a.obs[col].astype(str)
    region = a.obs.get("Region", pd.Series("", index=a.obs.index))
    prob_loc = a.obs.get("Probable_location", pd.Series("", index=a.obs.index))
    keep = _zeisel_isocortex_mask(region, prob_loc)
    n_kept = int(keep.sum())
    print(f"[zeisel] keeping {n_kept}/{a.n_obs} isocortex cells (Region~cortex \\ paleo/archicortex)")
    a = a[keep].copy()

    # Use the dataset's raw counts. Prefer layers["counts"] if present (the
    # pre-converted h5ad ships it); otherwise treat X as counts.
    if "counts" in a.layers:
        X = a.layers["counts"]
    else:
        X = a.X

    return _zeisel_assemble(X, a.var_names.astype(str).tolist(), a.obs)


def _load_zeisel_from_loom(src_path: Path) -> ad.AnnData:
    print(f"[zeisel] reading column attrs from {src_path.name}")
    import h5py
    with h5py.File(src_path, "r") as f:
        col_attrs_grp = f["col_attrs"]
        col_attrs: dict[str, np.ndarray] = {}
        for name in col_attrs_grp:
            arr = col_attrs_grp[name][:]
            if arr.dtype.kind == "S":
                arr = arr.astype("U")
            col_attrs[name] = arr
        gene_names = f["row_attrs/Gene"][:]
        if gene_names.dtype.kind == "S":
            gene_names = gene_names.astype("U")
        # Build keep mask.
        region_arr = col_attrs.get("Region", np.array([""] * f["matrix"].shape[1], dtype="U"))
        prob_loc_arr = col_attrs.get("Probable_location", np.array([""] * len(region_arr), dtype="U"))
        keep = _zeisel_isocortex_mask(pd.Series(region_arr), pd.Series(prob_loc_arr))
        keep_idx = np.where(keep)[0]
        print(f"[zeisel] keeping {len(keep_idx)}/{len(region_arr)} isocortex cells")

        # Stream the matrix in column chunks (loom matrix is genes × cells).
        matrix_ds = f["matrix"]
        _, n_cells = matrix_ds.shape
        chunk_size = 5000
        kept_chunks: list[sp.csr_matrix] = []
        for start in range(0, n_cells, chunk_size):
            end = min(start + chunk_size, n_cells)
            chunk_keep = keep[start:end]
            if not chunk_keep.any():
                continue
            block = matrix_ds[:, start:end][:, chunk_keep]        # genes × kept_in_chunk
            kept_chunks.append(sp.csr_matrix(block.T, dtype=np.float32))
        if not kept_chunks:
            raise RuntimeError("[zeisel] empty isocortex slice")
        X = sp.vstack(kept_chunks, format="csr")

    obs = pd.DataFrame(index=pd.Index([f"zeisel_{i}" for i in keep_idx], name="cell_id"))
    for name in ("Region", "Tissue", "Class", "Subclass", "ClusterName",
                 "Description", "Probable_location", "TaxonomyRank1",
                 "TaxonomyRank2", "TaxonomyRank3", "TaxonomyRank4", "SampleID"):
        if name in col_attrs:
            obs[name] = col_attrs[name][keep_idx].astype(str)
    return _zeisel_assemble(X, gene_names.astype(str).tolist(), obs)


def _zeisel_assemble(X, gene_names: list[str], obs: pd.DataFrame) -> ad.AnnData:
    """Common tail: stamp study/protocol/area columns, harmonize labels,
    wrap into an AnnData with layers['counts']."""
    obs = obs.copy()
    obs["study"] = "zeisel2018"
    obs["protocol"] = "10x_v2"
    obs["area"] = "cortex"
    obs["celltype_tasic2016"] = ""
    obs["celltype_tasic2018"] = ""
    obs["celltype_zeisel"] = obs.get("ClusterName", pd.Series("", index=obs.index)).astype(str)
    obs["celltype_harmonized"] = [
        _harmonize_zeisel(c, t4, p, d)
        for c, t4, p, d in zip(
            obs.get("Class", pd.Series("", index=obs.index)).astype(str),
            obs.get("TaxonomyRank4", pd.Series("", index=obs.index)).astype(str),
            obs.get("Probable_location", pd.Series("", index=obs.index)).astype(str),
            obs.get("Description", pd.Series("", index=obs.index)).astype(str),
        )
    ]
    obs["sample_id"] = obs.get("SampleID", obs.index).astype(str)

    X_sp = sp.csr_matrix(X, dtype=np.float32) if not sp.issparse(X) else X.tocsr().astype(np.float32)
    adata = ad.AnnData(
        X=X_sp,
        obs=obs,
        var=pd.DataFrame(index=pd.Index(gene_names, name="gene")),
    )
    adata.var_names_make_unique()
    adata.layers["counts"] = adata.X.copy()
    return adata


# --- Gene-symbol harmonization + concatenation ---------------------------

def _intersect_genes(
    adatas: dict[str, ad.AnnData],
    *,
    min_sources: int = 2,
) -> tuple[dict[str, ad.AnnData], dict[str, int]]:
    """K-of-N gene intersection across the input AnnDatas (default K=2).

    A gene (compared case-insensitively) is kept if it appears in at least
    ``min_sources`` source datasets. Sources lacking the gene receive a
    zero column at the corresponding position so the per-source AnnData is
    aligned to the same var axis as the rest. This keeps cross-protocol
    markers that one older-annotation source (e.g. Tasic 2016) lacks but
    the other two carry — without forcing a costly reannotation.

    Casing precedence stays Tasic 2018 → Tasic 2016 → Zeisel (most-modern
    annotation wins).

    Returns ``(subsetted adatas, audit)`` where ``audit`` reports:
      - ``n_kept``: total genes in output (≥K=min_sources sources)
      - ``n_recovered_kof2``: genes that strict K=N would have dropped but
        K=min_sources keeps (i.e. genes in exactly N-1 sources, mostly)
      - ``n_dropped_singleton``: genes present in fewer than ``min_sources``
        sources and therefore dropped from output
    """
    from collections import Counter

    per_source_lower: dict[str, dict[str, str]] = {}
    for name, ad_ in adatas.items():
        canon = {}
        for g in ad_.var_names.astype(str):
            canon.setdefault(g.lower(), g)
        per_source_lower[name] = canon

    n_sources = len(per_source_lower)

    # Count distinct sources holding each lower-cased gene.
    source_counts: Counter[str] = Counter()
    for d in per_source_lower.values():
        source_counts.update(d.keys())

    common_lower = {g for g, n in source_counts.items() if n >= min_sources}
    full_intersection_count = sum(1 for n in source_counts.values() if n == n_sources)
    n_recovered = len(common_lower) - full_intersection_count
    n_dropped_singleton = sum(1 for n in source_counts.values() if n < min_sources)

    canon_priority = ["tasic2018", "tasic2016", "zeisel"]
    canon_for: dict[str, str] = {}
    for low in sorted(common_lower):
        for src in canon_priority:
            if src in per_source_lower and low in per_source_lower[src]:
                canon_for[low] = per_source_lower[src][low]
                break

    sorted_common = sorted(common_lower)
    canon_genes = [canon_for[low] for low in sorted_common]
    n_common = len(canon_genes)

    out: dict[str, ad.AnnData] = {}
    for name, ad_ in adatas.items():
        canon_low = {g.lower(): g for g in ad_.var_names.astype(str)}
        n_cells = ad_.n_obs

        # For each common gene, find this source's column index (-1 if absent).
        src_idx_arr = np.full(n_common, -1, dtype=np.int64)
        for j, low in enumerate(sorted_common):
            if low in canon_low:
                src_idx_arr[j] = ad_.var_names.get_loc(canon_low[low])

        has_mask = src_idx_arr >= 0
        existing_src_cols = src_idx_arr[has_mask]   # source-side column indices
        out_positions = np.where(has_mask)[0]       # destination positions

        # Build the source's contribution to the aligned matrix. Sparse-
        # friendly: extract the existing columns, then remap COO col indices
        # to the full common-width shape; missing-source columns stay zero
        # because they're simply never inserted.
        X_src = ad_.X
        if sp.issparse(X_src):
            X_existing = X_src[:, existing_src_cols].tocoo()
            new_col = out_positions[X_existing.col]
            new_X = sp.csr_matrix(
                (X_existing.data.astype(np.float32, copy=False),
                 (X_existing.row, new_col)),
                shape=(n_cells, n_common),
            )
        else:
            X_arr = np.asarray(X_src, dtype=np.float32)
            dense = np.zeros((n_cells, n_common), dtype=np.float32)
            dense[:, has_mask] = X_arr[:, existing_src_cols]
            new_X = sp.csr_matrix(dense)

        obs_copy = ad_.obs.copy()
        sub = ad.AnnData(
            X=new_X,
            obs=obs_copy,
            var=pd.DataFrame(index=pd.Index(canon_genes, name="gene")),
        )
        sub.layers["counts"] = new_X.copy()
        out[name] = sub

    audit = {
        "n_kept": n_common,
        "n_recovered_kof2": n_recovered,
        "n_dropped_singleton": n_dropped_singleton,
    }
    return out, audit


def _build_integrated(adatas: dict[str, ad.AnnData]) -> tuple[ad.AnnData, dict[str, int]]:
    aligned, audit = _intersect_genes(adatas, min_sources=2)
    n_per = {n: a.n_obs for n, a in aligned.items()}
    print(
        f"[concat] K-of-N gene intersection (>=2 of 3 sources): "
        f"kept {audit['n_kept']} genes; "
        f"recovered {audit['n_recovered_kof2']} that strict K=3 would have dropped; "
        f"dropped {audit['n_dropped_singleton']} singletons (in 1 source only)"
    )
    print(f"[concat] cells per study: {n_per}")

    adata = ad.concat(
        list(aligned.values()),
        axis=0,
        join="inner",
        label=None,           # study col already set in each obs
        index_unique=None,
    )
    # Make obs_names globally unique (Tasic-Tasic-Zeisel may collide on barcodes).
    adata.obs_names = pd.Index(
        [f"{study}::{cid}" for study, cid in zip(adata.obs["study"].astype(str), adata.obs_names.astype(str))],
        name="cell_id",
    )
    # Coarse class for the gating-tree root.
    adata.obs["class_harmonized"] = adata.obs["celltype_harmonized"].astype(str).apply(_derive_class_harmonized)

    # Label-cleaning policy (same principle as the standalone Tasic loader,
    # load_mouse_brain.py): restrict to cells carrying a USABLE harmonized
    # cell-type label, which is the metric ground truth (LABEL_KEY =
    # celltype_harmonized). This drops (i) cells with no annotation (nan / "")
    # and (ii) Zeisel's generic "Neurons" label — real cells that lack a
    # class-resolvable subtype and would otherwise enter the ground-truth
    # partition as a spurious "Unknown" cluster, confounding bio-conservation
    # metrics. Genuinely-labelled rare types (e.g. Cajal-Retzius "CR") are kept.
    ct = adata.obs["celltype_harmonized"].astype(str).str.strip().str.lower()
    keep = ~ct.isin(UNUSABLE_LABELS)
    n_drop = int((~keep).sum())
    breakdown = (adata.obs.loc[~keep, "celltype_harmonized"].astype(str)
                 .value_counts().to_dict())
    adata = adata[keep.to_numpy()].copy()
    print(f"[filter] dropped {n_drop} cells without a usable harmonized label "
          f"{breakdown}; kept {adata.n_obs}")

    adata.layers["counts"] = adata.X.copy()
    return adata, audit


# --- Metadata sidecar ----------------------------------------------------

def _harmonization_audit(adata: ad.AnnData) -> dict[str, list[str]]:
    """Per-study list of native labels that fell OUTSIDE the harmonized vocab.

    These are surfaced in the JSON sidecar so the next iteration of the
    harmonization tables can shrink the audit list.
    """
    audit: dict[str, list[str]] = {}
    for study in ("tasic2016", "tasic2018", "zeisel2018"):
        mask = adata.obs["study"].astype(str) == study
        if not mask.any():
            continue
        labels = adata.obs.loc[mask, "celltype_harmonized"].astype(str)
        unmapped = sorted(set(labels) - HARMONIZED_VOCAB - {"nan", ""})
        if unmapped:
            audit[study] = unmapped
    return audit


def build_metadata(adata: ad.AnnData, intersection_audit: dict[str, int] | None = None) -> dict:
    per_study = (
        adata.obs["study"].astype(str).value_counts().to_dict()
    )
    per_celltype_harm = (
        adata.obs["celltype_harmonized"].astype(str).value_counts().to_dict()
    )
    per_class = (
        adata.obs["class_harmonized"].astype(str).value_counts().to_dict()
    )
    per_area = (
        adata.obs["area"].astype(str).value_counts().to_dict()
    )
    audit = intersection_audit or {}
    return {
        "dataset": DATASET_NAME,
        "citation": CITATION,
        "species": SPECIES,
        "tissue": TISSUE,
        "protocol": PROTOCOL,
        "n_cells_raw": int(adata.n_obs),
        "n_genes_raw": int(adata.n_vars),
        "n_genes_kept_kof2": int(audit.get("n_kept", adata.n_vars)),
        "n_genes_recovered_kof2": int(audit.get("n_recovered_kof2", 0)),
        "n_genes_dropped_singleton": int(audit.get("n_dropped_singleton", 0)),
        "source_urls": {
            "tasic2016_counts": GEO_BASE.format(acc="GSE71585", file=TASIC2016_COUNTS_FILENAME),
            "tasic2016_metadata": GEO_BASE.format(acc="GSE71585", file=TASIC2016_CLUSTER_FILENAME),
            "tasic2018_counts": GEO_BASE.format(acc="GSE115746", file=TASIC2018_COUNTS_FILENAME),
            "tasic2018_metadata": GEO_BASE.format(acc="GSE115746", file=TASIC2018_META_FILENAME),
            "zeisel2018_loom": ZEISEL_LOOM_URL,
        },
        "download_date": datetime.today().strftime("%Y-%m-%d"),
        "counts_layer": "counts",
        "batch_key": BATCH_KEY,
        "label_key": LABEL_KEY,
        "per_study_counts": per_study,
        "per_celltype_harmonized_counts": per_celltype_harm,
        "per_class_harmonized_counts": per_class,
        "per_area_counts": per_area,
        "harmonization_audit": _harmonization_audit(adata),
    }


# --- CLI ----------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=f"Load {DATASET_NAME} (Tasic 2016 + Tasic 2018 + Zeisel 2018 isocortex).")
    parser.add_argument("--force-download", action="store_true",
                        help="Re-fetch source files even if cached.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Overwrite existing data/raw/raw_mouse_cortex.h5ad.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if RAW_OUTPUT.exists() and not args.overwrite:
        raise FileExistsError(f"{RAW_OUTPUT} exists. Pass --overwrite to replace it.")

    t16 = _load_tasic2016(force_download=args.force_download)
    t18 = _load_tasic2018(force_download=args.force_download)
    zei = _load_zeisel(force_download=args.force_download)
    adata, intersection_audit = _build_integrated(
        {"tasic2016": t16, "tasic2018": t18, "zeisel": zei}
    )

    adata.uns["source"] = {
        "dataset": DATASET_NAME,
        "citation": CITATION,
        "batch_key": BATCH_KEY,
        "label_key": LABEL_KEY,
        "download_date": datetime.today().strftime("%Y-%m-%d"),
    }

    RAW_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    METADATA_OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(RAW_OUTPUT)
    with METADATA_OUTPUT.open("w") as handle:
        json.dump(build_metadata(adata, intersection_audit=intersection_audit), handle, indent=2)

    print(f"[ok] wrote {RAW_OUTPUT} ({adata.n_obs} cells x {adata.n_vars} genes)")
    print(f"[ok] wrote {METADATA_OUTPUT}")


if __name__ == "__main__":
    main()
