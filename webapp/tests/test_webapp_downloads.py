def test_augment_adata_roundtrip(tmp_path):
    import anndata as ad, numpy as np, pandas as pd
    from webapp.downloads import augment_adata
    from webapp.tests.fixtures.make_tiny_h5ad import write_tiny
    h5 = write_tiny(tmp_path/"tiny.h5ad"); a = ad.read_h5ad(h5)
    emb = pd.DataFrame({"cell_id": a.obs_names, "dim_0": np.zeros(a.n_obs),
                        "dim_1": np.ones(a.n_obs)})
    umap = pd.DataFrame({"cell_id": a.obs_names, "umap_0": np.zeros(a.n_obs),
                         "umap_1": np.zeros(a.n_obs)})
    out = augment_adata(h5, emb, umap, axis_names=["ax0", "ax1"],
                        out_path=tmp_path/"aug.h5ad")
    b = ad.read_h5ad(out)
    assert "X_blueprint" in b.obsm and "X_umap" in b.obsm
    assert "ax0" in b.obs and "ax1" in b.obs
    assert b.obsm["X_blueprint"].shape == (a.n_obs, 2)


def test_export_figure_png_and_html():
    import plotly.graph_objects as go
    from webapp.downloads import export_figure
    fig = go.Figure(go.Scatter(x=[0, 1], y=[0, 1]))
    png = export_figure(fig, "png"); assert isinstance(png, (bytes, bytearray)) and len(png) > 0
    html = export_figure(fig, "html"); assert isinstance(html, str) and "<html" in html.lower()
