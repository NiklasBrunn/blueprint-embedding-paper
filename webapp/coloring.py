"""Obs/gene color schemes and discrete palette for the explorer."""
from __future__ import annotations

import numpy as np


def scheme_for_obs(values):
    """Categorical (object/str/low-card) -> discrete; numeric -> sequential."""
    arr = list(values)
    is_numeric = all(
        isinstance(v, (int, float, np.integer, np.floating)) and not isinstance(v, bool)
        for v in arr
    ) and len(arr) > 0
    if is_numeric:
        return {"type": "sequential", "colorscale": "Viridis", "symmetric": False}
    return {"type": "categorical"}


def gene_log1p(X, var_names, gene):
    """log1p expression of one gene (case-insensitive var lookup)."""
    lower = [str(v).lower() for v in var_names]
    idx = lower.index(str(gene).lower())
    col = np.asarray(X)[:, idx]
    return np.log1p(col)


def discrete_colors(labels):
    """Return a {label: color} dict for categorical labels.

    Tries scripts.plots.umap_plots.build_categorical_colors, which returns a
    (color_map, point_colors, ordered_labels) tuple; the first element is the
    {label: color} mapping. Falls back to a Plotly qualitative palette.
    """
    cats = list(dict.fromkeys(labels))
    try:
        from scripts.plots.umap_plots import build_categorical_colors

        result = build_categorical_colors(labels)
        color_map = None
        if isinstance(result, dict):
            color_map = result
        elif isinstance(result, (tuple, list)) and result:
            first = result[0]
            if isinstance(first, dict):
                color_map = first
        if color_map:
            return dict(color_map)
    except Exception:
        pass
    import plotly.express as px

    palette = px.colors.qualitative.Dark24
    return {c: palette[i % len(palette)] for i, c in enumerate(cats)}
