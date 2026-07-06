"""Assemble the augmented AnnData download and export figures."""
from __future__ import annotations

from pathlib import Path


def augment_adata(h5ad, emb_df, umap_df, axis_names, out_path):
    """Write a copy of the h5ad with embedding coords + axis scores attached.

    obsm['X_blueprint'] = dim_* matrix, obsm['X_umap'] = umap coords,
    obs[<axis_name>] = each axis score, obs['blueprint_axis_assignment'] = the
    argmax axis name per cell. Aligns rows by cell_id == obs_names.
    Returns out_path (Path).
    """
    import anndata as ad
    import pandas as pd

    from webapp.mapping import axis_assignment

    a = ad.read_h5ad(h5ad)
    ids = a.obs_names.astype(str)
    emb = emb_df.set_index(emb_df["cell_id"].astype(str)).reindex(ids)
    umap = umap_df.set_index(umap_df["cell_id"].astype(str)).reindex(ids)
    dim_cols = [c for c in emb_df.columns if c.startswith("dim_")]
    a.obsm["X_blueprint"] = emb[dim_cols].to_numpy(dtype=float)
    a.obsm["X_umap"] = umap[["umap_0", "umap_1"]].to_numpy(dtype=float)
    for i, name in enumerate(axis_names):
        if i < len(dim_cols):
            a.obs[name] = emb[dim_cols[i]].to_numpy(dtype=float)
    a.obs["blueprint_axis_assignment"] = pd.Categorical(
        axis_assignment(emb_df, axis_names).reindex(ids).to_numpy())
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    a.write_h5ad(out_path)
    return out_path


def export_figure(fig, fmt):
    """Export a Plotly figure. fmt in {'png','svg','html'}.

    png/svg -> bytes (kaleido); html -> str.
    """
    fmt = fmt.lower()
    if fmt == "html":
        return fig.to_html(full_html=True, include_plotlyjs="cdn")
    if fmt in ("png", "svg"):
        return fig.to_image(format=fmt)  # bytes via kaleido
    raise ValueError(f"unsupported format: {fmt}")
