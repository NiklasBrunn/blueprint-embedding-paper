"""Unit tests for scripts/IO/load_human_immune.py — the pure helpers only
(no network download; the figshare fetch is exercised separately)."""
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
def load_human_immune():
    spec = importlib.util.spec_from_file_location(
        "load_human_immune", REPO_ROOT / "scripts" / "IO" / "load_human_immune.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["load_human_immune"] = module
    spec.loader.exec_module(module)
    return module


# --- as_csr_float32 --------------------------------------------------------

def test_as_csr_float32_from_dense(load_human_immune):
    out = load_human_immune.as_csr_float32(np.array([[0, 2], [3, 0]], dtype=np.int64))
    assert sp.isspmatrix_csr(out)
    assert out.dtype == np.float32
    assert np.array_equal(out.toarray(), np.array([[0, 2], [3, 0]], dtype=np.float32))


def test_as_csr_float32_from_sparse_preserves_values(load_human_immune):
    src = sp.coo_matrix(np.array([[1, 0, 0], [0, 0, 4]], dtype=np.float64))
    out = load_human_immune.as_csr_float32(src)
    assert sp.isspmatrix_csr(out)
    assert out.dtype == np.float32
    assert np.array_equal(out.toarray(), src.toarray().astype(np.float32))


# --- build_metadata --------------------------------------------------------

def _synthetic_adata(load_human_immune) -> ad.AnnData:
    X = sp.csr_matrix(np.array([[1, 0], [0, 2], [3, 0], [0, 4]], dtype=np.float32))
    obs = pd.DataFrame(
        {
            load_human_immune.BATCH_KEY: ["b1", "b1", "b2", "b2"],
            load_human_immune.LABEL_KEY: ["T cell", "B cell", "T cell", "NK cell"],
        },
        index=[f"cell_{i}" for i in range(4)],
    )
    var = pd.DataFrame(index=["GENE1", "GENE2"])
    return ad.AnnData(X=X, obs=obs, var=var)


def test_build_metadata_keys_and_counts(load_human_immune):
    adata = _synthetic_adata(load_human_immune)
    meta = load_human_immune.build_metadata(adata)

    assert meta["dataset"] == "human_immune"
    assert meta["batch_key"] == load_human_immune.BATCH_KEY == "batch"
    assert meta["label_key"] == load_human_immune.LABEL_KEY == "final_annotation"
    assert meta["n_cells_raw"] == 4
    assert meta["n_genes_raw"] == 2
    assert meta["batch_counts"] == {"b1": 2, "b2": 2}
    assert meta["celltype_counts"]["T cell"] == 2
