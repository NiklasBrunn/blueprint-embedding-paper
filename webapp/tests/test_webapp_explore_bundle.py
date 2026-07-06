"""Tests for the client-side explore bundle assembly (explore_bundle.build)."""
from __future__ import annotations

from pathlib import Path

import pytest

from webapp import explorer
from webapp.mapping import REPO_ROOT, map_h5ad, validate_blueprint
from webapp.tests.fixtures.make_tiny_h5ad import write_tiny

REPO = REPO_ROOT


@pytest.fixture(scope="module")
def bundle_inputs(tmp_path_factory):
    """Build a real MapResult + summary + adata once for the module."""
    import anndata as ad

    td = tmp_path_factory.mktemp("explore_bundle")
    h5 = write_tiny(td / "tiny.h5ad")
    bp = next((REPO / "outs" / "blueprints").rglob("embedder_blueprint_seed*.py"))
    ok, summary, errors = validate_blueprint(bp)
    assert ok, errors
    res = map_h5ad(blueprint=bp, h5ad=h5, outdir=td)
    adata = ad.read_h5ad(h5)
    return res, adata, summary


def test_bundle_shape(bundle_inputs):
    from webapp.explore_bundle import build

    res, adata, summary = bundle_inputs
    b = build(res, adata, summary)
    for key in ("coords", "axisOrder", "values", "branches", "axisMeta",
                "metadata", "modules", "geneList", "shown", "total"):
        assert key in b, key
    shown = b["shown"]
    assert len(b["coords"]["x"]) == shown
    assert len(b["coords"]["y"]) == shown
    assert b["axisOrder"] == [a["name"] for a in summary["blueprint"]["axes"]]
    for a in b["axisOrder"]:
        assert len(b["values"][a]) == shown
        for v in b["values"][a][:20]:
            assert round(v, 3) == v
    for v in b["coords"]["x"][:20]:
        assert round(v, 3) == v


def test_axisMeta_roles(bundle_inputs):
    from webapp.explore_bundle import build

    res, adata, summary = bundle_inputs
    b = build(res, adata, summary)
    for name in b["axisOrder"]:
        m = b["axisMeta"][name]
        assert set(m) == {"kind", "positive", "gate", "negatives"}
        assert isinstance(m["negatives"], list)


def test_geneList_orders_module_genes_first(bundle_inputs):
    from webapp.explore_bundle import build

    res, adata, summary = bundle_inputs
    b = build(res, adata, summary)
    module_genes = []
    seen = set()
    for genes in b["modules"].values():
        for g in genes:
            if g not in seen:
                seen.add(g)
                module_genes.append(g)
    assert b["geneList"][: len(module_genes)] == module_genes
    # no duplicates overall
    assert len(b["geneList"]) == len(set(b["geneList"]))
    # every var name is present
    assert set(str(v) for v in adata.var_names).issubset(set(b["geneList"]))


def test_branches_match_helper(bundle_inputs):
    from webapp.explore_bundle import build

    res, adata, summary = bundle_inputs
    b = build(res, adata, summary)
    assert b["branches"] == explorer.branches_from_summary(summary)


def test_metadata_categorical_only(bundle_inputs):
    from webapp.explore_bundle import build

    res, adata, summary = bundle_inputs
    b = build(res, adata, summary)
    # cell_type is a known low-card categorical -> included
    assert "cell_type" in b["metadata"]
    meta = b["metadata"]["cell_type"]
    assert meta["kind"] == "label"
    assert len(meta["labels"]) == len(meta["colors"])
    assert len(meta["cellIdx"]) == b["shown"]
    assert all(0 <= i < len(meta["labels"]) for i in meta["cellIdx"])
    # an all-unique id column is excluded
    import numpy as np

    adata2 = adata.copy()
    adata2.obs["uid"] = [f"u{i}" for i in range(adata2.n_obs)]
    b2 = build(res, adata2, summary)
    assert "uid" not in b2["metadata"]


def test_metadata_handles_categorical_dtype(bundle_inputs):
    """Real h5ads store obs (e.g. cell_type) as pandas Categorical; build() must
    not crash on it. Regression for: 'Cannot setitem on a Categorical with a new
    category (unmatched)'."""
    import pandas as pd
    from webapp.explore_bundle import build

    res, adata, summary = bundle_inputs
    adata2 = adata.copy()
    adata2.obs["cell_type"] = pd.Categorical(adata2.obs["cell_type"].astype(str))
    b = build(res, adata2, summary)  # must NOT raise
    assert "cell_type" in b["metadata"]
    assert len(b["metadata"]["cell_type"]["cellIdx"]) == b["shown"]


def test_metadata_categorical_with_missing_cells():
    """A cell absent from obs (reindex -> NaN in a Categorical column) must map
    to 'unmatched' rather than raising."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    from webapp import explore_bundle

    obs = pd.DataFrame({"grp": pd.Categorical(["A", "B", "A", "B"])},
                       index=["c0", "c1", "c2", "c3"])
    a = ad.AnnData(X=np.zeros((4, 3), dtype=float), obs=obs)
    md = explore_bundle._build_metadata(
        a, np.array(["c0", "c1", "MISSING", "c3"]), np.array([0, 1, 2, 3]))
    assert "grp" in md
    assert "unmatched" in md["grp"]["labels"]
    assert len(md["grp"]["cellIdx"]) == 4


def test_cell_cap_downsamples(bundle_inputs):
    from webapp.explore_bundle import build

    res, adata, summary = bundle_inputs
    n = adata.n_obs
    cap = n // 3
    b = build(res, adata, summary, cell_cap=cap)
    assert b["total"] == n
    assert b["shown"] == cap
    assert len(b["coords"]["x"]) == cap
    for a in b["axisOrder"]:
        assert len(b["values"][a]) == cap
    for meta in b["metadata"].values():
        assert len(meta["cellIdx"]) == cap
