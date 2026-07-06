"""Unit tests for scripts/IO/load_mouse_cortex.py — the pure helpers only
(no network, no h5py loom read: those live inside the per-source loaders).

Uses a lazy module load via importlib so this test file stays dependency-light
when the bigger optional deps aren't needed.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse as sp


REPO_ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def load_mouse_cortex():
    spec = importlib.util.spec_from_file_location(
        "load_mouse_cortex", REPO_ROOT / "scripts" / "IO" / "load_mouse_cortex.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["load_mouse_cortex"] = module
    spec.loader.exec_module(module)
    return module


# --- _harmonize_subclass_tasic2016 ----------------------------------------

@pytest.mark.parametrize("native,expected", [
    ("L2", "L2/3 IT"),
    ("L2/3", "L2/3 IT"),
    ("L4", "L4"),
    ("L5", "L5 IT"),
    ("L5a", "L5 IT"),
    ("L5b", "L5 PT"),
    ("L6a", "L6 CT"),
    ("L6b", "L6b"),
    ("Vip", "Vip"),
    ("Sst", "Sst"),
    ("Pvalb", "Pvalb"),
    ("Ndnf", "Lamp5"),
    ("Astro", "Astro"),
    ("OPC", "OPC"),
    ("Oligo", "Oligo"),
    ("Micro", "Microglia"),
    ("Endo", "Endo"),
    ("SMC", "SMC"),
    ("  Pvalb  ", "Pvalb"),                  # whitespace stripped
    ("UnknownNovelType", "UnknownNovelType"),  # unmapped passes through
])
def test_harmonize_tasic2016(load_mouse_cortex, native, expected):
    assert load_mouse_cortex._harmonize_subclass_tasic2016(native) == expected


def test_harmonize_tasic2016_nan_inputs(load_mouse_cortex):
    assert pd.isna(load_mouse_cortex._harmonize_subclass_tasic2016(None))
    assert pd.isna(load_mouse_cortex._harmonize_subclass_tasic2016(""))
    assert pd.isna(load_mouse_cortex._harmonize_subclass_tasic2016(np.nan))


# --- _harmonize_subclass_tasic2018 ----------------------------------------

@pytest.mark.parametrize("native,expected", [
    ("L2/3 IT VISp", "L2/3 IT"),       # area suffix stripped
    ("L2/3 IT ALM", "L2/3 IT"),        # area suffix stripped
    ("L4", "L4"),
    ("L5 IT VISp", "L5 IT"),
    ("L5 PT ALM", "L5 PT"),
    ("L5/6 NP", "L5 NP"),
    ("L5 NP", "L5 NP"),
    ("L6 IT", "L6 IT"),
    ("L6 CT VISp", "L6 CT"),
    ("L6b", "L6b"),
    ("Pvalb", "Pvalb"),
    ("Sst", "Sst"),
    ("Sst Chodl", "Sst Chodl"),
    ("Vip", "Vip"),
    ("Lamp5", "Lamp5"),
    ("Sncg", "Sncg"),
    ("Serpinf1", "Sncg"),
    ("Meis2", "Meis2"),
    ("Macrophage", "Macrophage"),
    ("Microglia", "Microglia"),
    ("VLMC", "VLMC"),
    ("Endo", "Endo"),
    ("Peri", "Peri"),
    ("SMC", "SMC"),
    ("UnknownType", "UnknownType"),    # unmapped passes through
])
def test_harmonize_tasic2018(load_mouse_cortex, native, expected):
    assert load_mouse_cortex._harmonize_subclass_tasic2018(native) == expected


# --- _harmonize_zeisel ----------------------------------------------------
# Signature: (cls, taxonomy_rank4, probable_location, description).
# Non-neurons resolve by Class; excitatory neurons by Probable_location;
# interneurons by Description. Mirrors the actual l5_all atlas layout.

@pytest.mark.parametrize("cls,tax4,prob_loc,desc,expected", [
    # Non-neurons: Class-based (coarse atlas vocab)
    ("Astrocytes",     "Astrocytes",       "",                              "",                                                 "Astro"),
    ("Oligos",         "Oligodendrocytes", "",                              "",                                                 "Oligo"),
    ("Immune",         "Microglia",        "",                              "",                                                 "Microglia"),
    ("Vascular",       "Vascular cells",   "",                              "",                                                 "Endo"),
    # Excitatory cortex Neurons resolve by Probable_location (layer string)
    ("Neurons", "Telencephalon projecting excitatory neurons", "Cortical pyramidal layer 2/3",   "Excitatory neurons, cerebral cortex", "L2/3 IT"),
    ("Neurons", "Telencephalon projecting excitatory neurons", "Cortical pyramidal layer 4",     "Excitatory neurons, cerebral cortex", "L4"),
    ("Neurons", "Telencephalon projecting excitatory neurons", "Cortical pyramidal layer 5",     "Excitatory neurons, cerebral cortex", "L5 IT"),
    ("Neurons", "Telencephalon projecting excitatory neurons", "Cortex pyramidal layer 6",       "Excitatory neurons, cerebral cortex", "L6 CT"),
    ("Neurons", "Telencephalon projecting excitatory neurons", "Cortex pyramidal layer 6b",      "Excitatory neurons, cerebral cortex", "L6b"),
    ("Neurons", "Telencephalon projecting excitatory neurons", "Cingulate/Retrosplenial area, layer 5", "Excitatory neurons, cerebral cortex", "L5 IT"),
    ("Neurons", "Telencephalon projecting excitatory neurons", "Lateral cortex layer 6: gustatory, barrel field, auditory", "Excitatory neurons, cerebral cortex", "L6 IT"),
    ("Neurons", "Telencephalon projecting excitatory neurons", "Cortical pyramidal (poor markers)", "Excitatory neurons, cerebral cortex", "Neurons"),  # unmapped
    # Inhibitory interneurons resolve by Description (functional class)
    ("Neurons", "Telencephalon inhibitory interneurons", "Hippocampus,Cortex",       "Axo-axonic, cortex/hippocampus",                          "Pvalb"),
    ("Neurons", "Telencephalon inhibitory interneurons", "Hippocampus,Cortex",       "Basket and bistratified cells, cortex/hippocampus",       "Pvalb"),
    ("Neurons", "Telencephalon inhibitory interneurons", "Hippocampus,Cortex",       "CGE-derived neurogliaform cells, cortex/hippocampus",      "Lamp5"),
    ("Neurons", "Telencephalon inhibitory interneurons", "Hippocampus,Cortex",       "Ivy and MGE-derived neurogliaform cells, cortex/hippocampus", "Lamp5"),
    ("Neurons", "Telencephalon inhibitory interneurons", "Hippocampus,Cortex",       "Interneuron-selective interneurons, cortex/hippocampus",  "Vip"),
    ("Neurons", "Telencephalon inhibitory interneurons", "Hippocampus (sr/slm),Cortex", "R-LM border Cck interneurons, cortex/hippocampus",    "Vip"),
    ("Neurons", "Telencephalon inhibitory interneurons", "Hippocampus,Cortex",       "Non-border Cck interneurons, cortex/hippocampus",         "Vip"),
    ("Neurons", "Telencephalon inhibitory interneurons", "Hippocampus (so),Cortex",  "Hippocamposeptal projection, cortex/hippocampus",         "Sst"),
    ("Neurons", "Telencephalon inhibitory interneurons", "Hippocampus,Cortex",       "Sleep-active, long-range projection interneurons, cortex/hippocampus", "Sst Chodl"),
    ("Neurons", "Telencephalon inhibitory interneurons", "Hippocampus,Cortex",       "Some novel interneuron description with no keyword",       "Neurons"),  # unmapped
])
def test_harmonize_zeisel(load_mouse_cortex, cls, tax4, prob_loc, desc, expected):
    assert load_mouse_cortex._harmonize_zeisel(cls, tax4, prob_loc, desc) == expected


def test_harmonize_zeisel_nan_inputs(load_mouse_cortex):
    assert pd.isna(load_mouse_cortex._harmonize_zeisel(None, "", "", ""))


# --- _derive_class_harmonized --------------------------------------------

@pytest.mark.parametrize("subclass,expected_class", [
    ("L2/3 IT", "Glutamatergic"),
    ("L4", "Glutamatergic"),
    ("L5 IT", "Glutamatergic"),
    ("L5 PT", "Glutamatergic"),
    ("L5 NP", "Glutamatergic"),
    ("L6 IT", "Glutamatergic"),
    ("L6 CT", "Glutamatergic"),
    ("L6b", "Glutamatergic"),
    ("Pvalb", "GABAergic"),
    ("Sst", "GABAergic"),
    ("Sst Chodl", "GABAergic"),
    ("Vip", "GABAergic"),
    ("Lamp5", "GABAergic"),
    ("Sncg", "GABAergic"),
    ("Meis2", "GABAergic"),
    ("Astro", "Non-neuronal"),
    ("Oligo", "Non-neuronal"),
    ("OPC", "Non-neuronal"),
    ("Microglia", "Non-neuronal"),
    ("Endo", "Non-neuronal"),
    ("Peri", "Non-neuronal"),
    ("SMC", "Non-neuronal"),
    ("VLMC", "Non-neuronal"),
    ("Macrophage", "Non-neuronal"),
    ("UnknownNovelLabel", "Unknown"),
    ("", "Unknown"),
])
def test_derive_class_harmonized(load_mouse_cortex, subclass, expected_class):
    assert load_mouse_cortex._derive_class_harmonized(subclass) == expected_class


# --- _intersect_genes (case-insensitive intersection) --------------------

def _make_dummy(genes: list[str], n_cells: int, study: str) -> ad.AnnData:
    X = sp.csr_matrix(np.zeros((n_cells, len(genes)), dtype=np.float32))
    obs = pd.DataFrame(
        {
            "study": [study] * n_cells,
            "protocol": [study] * n_cells,
            "area": ["x"] * n_cells,
            "celltype_tasic2016": pd.NA,
            "celltype_tasic2018": pd.NA,
            "celltype_zeisel": pd.NA,
            "celltype_harmonized": ["Pvalb"] * n_cells,
            "sample_id": [f"{study}_{i}" for i in range(n_cells)],
        },
        index=[f"{study}_cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    return ad.AnnData(X=X, obs=obs, var=var)


def _make_dummy_with_values(genes: list[str], values: np.ndarray, study: str) -> ad.AnnData:
    """Like _make_dummy but with caller-supplied dense X (so K-of-N zero-fill
    is visually distinguishable from a real zero in the source data)."""
    assert values.shape == (values.shape[0], len(genes))
    X = sp.csr_matrix(values.astype(np.float32))
    n_cells = values.shape[0]
    obs = pd.DataFrame(
        {
            "study": [study] * n_cells,
            "protocol": [study] * n_cells,
            "area": ["x"] * n_cells,
            "celltype_tasic2016": pd.NA,
            "celltype_tasic2018": pd.NA,
            "celltype_zeisel": pd.NA,
            "celltype_harmonized": ["Pvalb"] * n_cells,
            "sample_id": [f"{study}_{i}" for i in range(n_cells)],
        },
        index=[f"{study}_cell_{i}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=pd.Index(genes, name="gene"))
    return ad.AnnData(X=X, obs=obs, var=var)


def test_intersect_genes_kof2_case_insensitive_and_dropped_singletons(load_mouse_cortex):
    # Three sources share Snap25 (case-mixed), Pvalb, Gad1; each has a
    # singleton-only gene. Under K=2 the shared three pass, the three
    # singletons drop. Canonical casing comes from tasic2018.
    t16 = _make_dummy(["Snap25", "Pvalb", "Gad1", "OnlyIn16"], n_cells=2, study="tasic2016")
    t18 = _make_dummy(["snap25", "Pvalb", "Gad1", "OnlyIn18"], n_cells=3, study="tasic2018")
    zei = _make_dummy(["SNAP25", "Pvalb", "Gad1", "OnlyInZeisel"], n_cells=4, study="zeisel")

    aligned, audit = load_mouse_cortex._intersect_genes(
        {"tasic2016": t16, "tasic2018": t18, "zeisel": zei}
    )
    expected_genes = {"snap25", "pvalb", "gad1"}
    for name, a in aligned.items():
        assert {g.lower() for g in a.var_names} == expected_genes
        # Canonical casing from tasic2018 (snap25 lowercased there).
        assert "snap25" in a.var_names

    assert audit["n_kept"] == 3
    assert audit["n_dropped_singleton"] == 3        # OnlyIn16, OnlyIn18, OnlyInZeisel
    assert audit["n_recovered_kof2"] == 0           # nothing in exactly 2 of 3 here


def test_intersect_genes_kof2_preserves_two_of_three_with_zero_fill(load_mouse_cortex):
    # Gene 'Lamp5' is in tasic2018 + zeisel but NOT tasic2016 (the actual
    # Lamp5/Ndnf situation in the real mouse_cortex pipeline). Under K=2 it
    # must survive, with the tasic2016 column zero-filled.
    t16 = _make_dummy_with_values(
        ["Snap25", "Pvalb"],
        np.array([[1.0, 2.0], [3.0, 4.0]]),  # no Lamp5 column at all
        study="tasic2016",
    )
    t18 = _make_dummy_with_values(
        ["Snap25", "Pvalb", "Lamp5"],
        np.array([[10.0, 20.0, 30.0], [40.0, 50.0, 60.0], [70.0, 80.0, 90.0]]),
        study="tasic2018",
    )
    zei = _make_dummy_with_values(
        ["Snap25", "Pvalb", "Lamp5"],
        np.array([[100.0, 200.0, 300.0], [400.0, 500.0, 600.0]]),
        study="zeisel",
    )

    aligned, audit = load_mouse_cortex._intersect_genes(
        {"tasic2016": t16, "tasic2018": t18, "zeisel": zei}
    )
    expected_genes = {"snap25", "pvalb", "lamp5"}
    assert audit["n_kept"] == 3
    assert audit["n_recovered_kof2"] == 1           # Lamp5 (in 2 of 3, not in t16)
    assert audit["n_dropped_singleton"] == 0

    # Every aligned adata has all three genes, in canonical Tasic-2018 case.
    for name, a in aligned.items():
        assert list(a.var_names) == ["Lamp5", "Pvalb", "Snap25"] or \
               {g.lower() for g in a.var_names} == expected_genes
        assert a.n_vars == 3

    # Tasic 2016: Lamp5 column is all zero (zero-filled), Snap25 + Pvalb
    # carry the real values 1-4. Order alphabetic in canonical: Lamp5, Pvalb, Snap25.
    t16_aligned = aligned["tasic2016"]
    t16_dense = t16_aligned.X.toarray()
    lamp5_col = list(t16_aligned.var_names).index("Lamp5")
    snap25_col = list(t16_aligned.var_names).index("Snap25")
    pvalb_col = list(t16_aligned.var_names).index("Pvalb")
    assert (t16_dense[:, lamp5_col] == 0).all(), "Lamp5 should be zero-filled in Tasic 2016"
    assert t16_dense[0, snap25_col] == 1.0
    assert t16_dense[1, snap25_col] == 3.0
    assert t16_dense[0, pvalb_col] == 2.0
    assert t16_dense[1, pvalb_col] == 4.0

    # Tasic 2018: Lamp5 should carry its real values 30, 60, 90.
    t18_aligned = aligned["tasic2018"]
    t18_dense = t18_aligned.X.toarray()
    lamp5_col_18 = list(t18_aligned.var_names).index("Lamp5")
    assert t18_dense[0, lamp5_col_18] == 30.0
    assert t18_dense[1, lamp5_col_18] == 60.0
    assert t18_dense[2, lamp5_col_18] == 90.0


def test_intersect_genes_min_sources_3_reproduces_strict_intersection(load_mouse_cortex):
    # Sanity: setting min_sources=3 falls back to the old strict-intersection
    # behavior — the Lamp5-in-2 example should now drop Lamp5.
    t16 = _make_dummy(["Snap25", "Pvalb"],          n_cells=2, study="tasic2016")
    t18 = _make_dummy(["Snap25", "Pvalb", "Lamp5"], n_cells=3, study="tasic2018")
    zei = _make_dummy(["Snap25", "Pvalb", "Lamp5"], n_cells=4, study="zeisel")

    aligned, audit = load_mouse_cortex._intersect_genes(
        {"tasic2016": t16, "tasic2018": t18, "zeisel": zei},
        min_sources=3,
    )
    assert audit["n_kept"] == 2
    for a in aligned.values():
        assert {g.lower() for g in a.var_names} == {"snap25", "pvalb"}


# --- _build_integrated (concat preserves obs columns) --------------------

def test_build_integrated_preserves_obs_columns(load_mouse_cortex):
    genes = ["Snap25", "Pvalb", "Gad1"]
    t16 = _make_dummy(genes, n_cells=2, study="tasic2016")
    t18 = _make_dummy(genes, n_cells=3, study="tasic2018")
    zei = _make_dummy(genes, n_cells=4, study="zeisel")

    # Set per-source celltype_harmonized so class_harmonized derives correctly.
    t16.obs["celltype_harmonized"] = ["L2/3 IT", "Pvalb"]
    t18.obs["celltype_harmonized"] = ["L5 PT", "Vip", "Astro"]
    zei.obs["celltype_harmonized"] = ["L6 CT", "Microglia", "Oligo", "Sst"]

    integrated, audit = load_mouse_cortex._build_integrated(
        {"tasic2016": t16, "tasic2018": t18, "zeisel": zei}
    )
    assert audit["n_kept"] == 3

    assert integrated.n_obs == 9
    assert integrated.n_vars == 3
    for col in (
        "study", "protocol", "area",
        "celltype_tasic2016", "celltype_tasic2018", "celltype_zeisel",
        "celltype_harmonized", "class_harmonized", "sample_id",
    ):
        assert col in integrated.obs.columns, f"missing obs column: {col}"

    # class_harmonized is derived correctly per row.
    expected_classes = [
        "Glutamatergic", "GABAergic",                # from t16
        "Glutamatergic", "GABAergic", "Non-neuronal",  # from t18
        "Glutamatergic", "Non-neuronal", "Non-neuronal", "GABAergic",  # from zei
    ]
    assert integrated.obs["class_harmonized"].tolist() == expected_classes

    # obs_names are globally unique with study:: prefix.
    assert integrated.obs_names.is_unique
    assert all("::" in n for n in integrated.obs_names)

    # counts layer present.
    assert "counts" in integrated.layers


# --- _harmonization_audit ------------------------------------------------

def test_harmonization_audit_lists_unmapped_labels(load_mouse_cortex):
    # Build a minimal integrated adata with some unmapped labels.
    genes = ["Snap25", "Pvalb"]
    t16 = _make_dummy(genes, n_cells=2, study="tasic2016")
    t18 = _make_dummy(genes, n_cells=2, study="tasic2018")
    zei = _make_dummy(genes, n_cells=2, study="zeisel2018")
    t16.obs["celltype_harmonized"] = ["L2/3 IT", "WeirdT16Label"]
    t18.obs["celltype_harmonized"] = ["L5 PT", "Vip"]
    zei.obs["celltype_harmonized"] = ["Oligo", "WeirdZeiselLabel"]

    integrated, _ = load_mouse_cortex._build_integrated(
        {"tasic2016": t16, "tasic2018": t18, "zeisel": zei}
    )
    integrated.obs["study"] = integrated.obs["study"].replace(
        {"zeisel": "zeisel2018"}
    )

    audit = load_mouse_cortex._harmonization_audit(integrated)
    assert "tasic2016" in audit and "WeirdT16Label" in audit["tasic2016"]
    assert "zeisel2018" in audit and "WeirdZeiselLabel" in audit["zeisel2018"]
    # Tasic 2018 was fully mapped.
    assert "tasic2018" not in audit or len(audit.get("tasic2018", [])) == 0
