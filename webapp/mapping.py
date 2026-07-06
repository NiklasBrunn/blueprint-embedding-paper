"""Map uploaded h5ad data through a blueprint and prep it for the explorer."""
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
# Prefer the vendored copy (standalone); fall back to the parent repo in-tree.
_VENDOR = Path(__file__).resolve().parent / "vendor"
_SRC_ROOT = _VENDOR if _VENDOR.exists() else REPO_ROOT
EXTRACT_RUN = _SRC_ROOT / "scripts" / "agentic_embedding" / "extract_run.py"


def validate_blueprint(path):
    """Static-validate a populated embedder_blueprint.py via extract_run.py.

    Returns (ok: bool, summary: dict | None, errors: list[str])."""
    path = Path(path)
    errors = []
    with tempfile.TemporaryDirectory() as td:
        wd = Path(td)
        (wd / "embedder_blueprint.py").write_bytes(path.read_bytes())
        # extract_run.py also needs a scenario_config.jsonc present.
        (wd / "scenario_config.jsonc").write_text(
            '{"species":["mouse"],"tissue":null,"celltype":null}')
        proc = subprocess.run(
            [sys.executable, str(EXTRACT_RUN), str(wd)],
            capture_output=True, text=True, timeout=120)
        if proc.returncode != 0:
            errors.append(proc.stderr.strip() or "extract_run.py failed")
            return False, None, errors
        summary_path = wd / "summary.json"
        if not summary_path.exists():
            errors.append("no summary.json produced")
            return False, None, errors
        summary = json.loads(summary_path.read_text())
    sanity = summary.get("sanity", {})
    if sanity.get("placeholder_strings_present"):
        errors.append("placeholder strings present")
    if sanity.get("orphaned_modules"):
        errors.append(f"orphaned modules: {sanity['orphaned_modules']}")
    if not summary.get("blueprint", {}).get("n_axes"):
        errors.append("no axes registered — this looks like an unpopulated "
                      "template (generation may not have finished, or the "
                      "wrong .py was selected)")
    ok = not errors
    return ok, summary, errors


@dataclass
class MapResult:
    parquet: Path
    umap_coords: Path
    coverage_json: Path | None


def axis_assignment(emb_df, axis_names):
    """Argmax assignment of each cell to its single strongest axis.

    The embedding parquet's dim_* columns are already z-scored + clipped by the
    embedder, so an argmax over them is the same axis-ownership notion used in
    scripts/temporal_metrics/metrics.py (``owner = X.argmax(axis=1)``). dim_i
    maps positionally to axis_names[i]. Returns a pandas Series of axis names
    indexed by cell_id (str), named 'blueprint_axis_assignment'."""
    import numpy as np
    import pandas as pd

    dim_cols = sorted((c for c in emb_df.columns if c.startswith("dim_")),
                      key=lambda c: int(c.split("_")[1]))
    owner = emb_df[dim_cols].to_numpy().argmax(axis=1)
    names = np.asarray(axis_names, dtype=object)
    return pd.Series(names[owner], index=emb_df["cell_id"].astype(str),
                     name="blueprint_axis_assignment")


def map_h5ad(*, blueprint, h5ad, outdir):
    """Run the embedder on an h5ad through a blueprint, then compute + cache UMAP.

    Returns MapResult(parquet, umap_coords, coverage_json)."""
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    parquet = outdir / "embedding.parquet"
    proc = subprocess.run(
        [sys.executable, str(blueprint), "--input", str(h5ad),
         "--output", str(parquet)],
        capture_output=True, text=True, timeout=1200)
    if proc.returncode != 0:
        raise RuntimeError(f"embedder failed: {proc.stderr.strip()}")

    import pandas as pd
    import umap

    df = pd.read_parquet(parquet)
    dim_cols = [c for c in df.columns if c.startswith("dim_")]
    X = df[dim_cols].to_numpy()
    n_neighbors = min(30, max(2, len(df) - 1))
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=0.15,
                        metric="euclidean", random_state=42, n_components=2)
    xy = reducer.fit_transform(X)
    coords = pd.DataFrame({"cell_id": df["cell_id"].to_numpy(),
                           "umap_0": xy[:, 0], "umap_1": xy[:, 1]})
    umap_coords = outdir / "umap_coords.parquet"
    coords.to_parquet(umap_coords)

    cov = parquet.with_name(parquet.stem + "_coverage.json")
    return MapResult(parquet=parquet, umap_coords=umap_coords,
                     coverage_json=cov if cov.exists() else None)
