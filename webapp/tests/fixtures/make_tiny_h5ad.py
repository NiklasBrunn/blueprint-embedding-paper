"""Generate a tiny AnnData fixture whose genes overlap an archived blueprint.

The var_names deliberately reuse gene names from the first archived
embedder_blueprint, so running that blueprint's embedder over the fixture
produces non-degenerate dim_* values (modules actually find their genes).
"""
from __future__ import annotations

from pathlib import Path

from webapp.mapping import REPO_ROOT, validate_blueprint

REPO = REPO_ROOT


def _blueprint_genes() -> list[str]:
    bp = next((REPO / "outs" / "blueprints").rglob("embedder_blueprint_seed*.py"))
    ok, summary, errors = validate_blueprint(bp)
    if not ok or summary is None:
        raise RuntimeError(f"could not validate blueprint {bp}: {errors}")
    genes: list[str] = []
    seen: set[str] = set()
    for module in summary["blueprint"]["modules"]:
        for g in module["genes"]:
            if g not in seen:
                seen.add(g)
                genes.append(g)
    return genes


def write_tiny(path) -> Path:
    """Write a small AnnData h5ad whose genes overlap the first archived
    blueprint, and return its Path."""
    import anndata as ad
    import numpy as np
    import pandas as pd

    path = Path(path)

    var_names = _blueprint_genes()[:60]
    i = 0
    while len(var_names) < 30:
        var_names.append(f"GENEpad{i}")
        i += 1

    n_cells = 120
    n_genes = len(var_names)
    rng = np.random.default_rng(0)
    X = rng.poisson(lam=2.0, size=(n_cells, n_genes)).astype("float32")

    obs_names = [f"cell{i}" for i in range(n_cells)]
    labels = ["A", "B", "C", "D"]
    cell_type = [labels[i % len(labels)] for i in range(n_cells)]
    obs = pd.DataFrame(
        {
            "cell_type": pd.Categorical(cell_type),
            "n_counts": X.sum(axis=1).astype("float64"),
        },
        index=obs_names,
    )
    var = pd.DataFrame(index=pd.Index(var_names, name="gene"))

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names = obs_names

    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
    return path
