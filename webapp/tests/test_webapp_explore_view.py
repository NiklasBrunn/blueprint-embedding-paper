from webapp import explore_view
from webapp.tests.fixtures.make_tiny_h5ad import write_tiny


def _collect_ids(component, out=None):
    """Recursively collect every component id in a Dash component tree."""
    if out is None:
        out = set()
    cid = getattr(component, "id", None)
    if isinstance(cid, str):
        out.add(cid)
    children = getattr(component, "children", None)
    if isinstance(children, (list, tuple)):
        for c in children:
            _collect_ids(c, out)
    elif children is not None and hasattr(children, "_traverse_ids"):
        _collect_ids(children, out)
    elif children is not None and hasattr(children, "children"):
        _collect_ids(children, out)
    return out


def test_panel_has_required_ids():
    ids = _collect_ids(explore_view.panel())
    required = {
        "exp-graph", "exp-sidebar", "exp-mode", "exp-branch", "exp-palette",
        "exp-time-host", "explore-bundle", "gene-vector", "exp-paths",
        "exp-gene-select",
    }
    missing = required - ids
    assert not missing, f"missing ids: {missing}"


def _real_gene(h5ad):
    import anndata as ad
    a = ad.read_h5ad(h5ad)
    return str(a.var_names[0])


def test_gene_vector_valid(tmp_path):
    h5 = write_tiny(tmp_path / "tiny.h5ad")
    gene = _real_gene(h5)
    out = explore_view.gene_vector(str(h5), gene)
    assert isinstance(out, dict)
    assert out["gene"] == gene
    assert isinstance(out["values"], list) and out["values"]


def test_gene_vector_missing(tmp_path):
    h5 = write_tiny(tmp_path / "tiny.h5ad")
    assert explore_view.gene_vector(str(h5), "NOT_A_GENE") is None


def test_gene_vector_bad_paths():
    assert explore_view.gene_vector(None, "x") is None
    assert explore_view.gene_vector("/no/such.h5ad", "x") is None


def test_gene_vector_aligns_to_cell_ids(tmp_path):
    """Values must be reindexed to the requested cell order (not h5ad order),
    and cells absent from the h5ad must come back as 0.0."""
    import anndata as ad
    import numpy as np

    h5 = write_tiny(tmp_path / "tiny.h5ad")
    gene = _real_gene(h5)
    a = ad.read_h5ad(str(h5))
    ids = [str(c) for c in a.obs_names]

    # baseline: full h5ad order
    base = explore_view.gene_vector(str(h5), gene)["values"]
    assert len(base) == len(ids)

    # reversed order -> values must reverse one-to-one
    rev = explore_view.gene_vector(str(h5), gene, cell_ids=ids[::-1])["values"]
    assert np.allclose(rev, base[::-1])

    # a missing cell id -> 0.0 in that slot, real value otherwise
    mixed = explore_view.gene_vector(str(h5), gene,
                                     cell_ids=[ids[0], "NOPE_CELL"])["values"]
    assert len(mixed) == 2
    assert np.isclose(mixed[0], base[0])
    assert mixed[1] == 0.0


def test_displayed_cell_ids_matches_bundle_shown(tmp_path):
    """displayed_cell_ids must reproduce build()'s plot order/length so gene
    vectors line up with the bundle coords."""
    import pandas as pd
    from webapp import explore_bundle

    pq = tmp_path / "emb.parquet"
    pd.DataFrame({"cell_id": [f"c{i}" for i in range(10)],
                  "dim_0": range(10)}).to_parquet(pq)
    full = explore_bundle.displayed_cell_ids(pq)
    assert full == [f"c{i}" for i in range(10)]
    capped = explore_bundle.displayed_cell_ids(pq, cell_cap=4)
    assert len(capped) == 4
    assert set(capped) <= set(full)
