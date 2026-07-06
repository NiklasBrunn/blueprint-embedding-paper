"""
Smoke + behavioural tests for scripts/plots/plot_metrics.py.

We don't visually validate the figures here — that's a manual check on real data.
These tests pin the data-prep contract (seed aggregation, composite formula,
NaN handling) and confirm draw_* helpers produce non-empty PNG/PDF artifacts on
a synthetic CSV.
"""
from __future__ import annotations

import importlib.util
import os
import sys
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import numpy as np
import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, rel_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / rel_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def plot_metrics():
    return _load("plot_metrics", "scripts/plots/plot_metrics.py")


@pytest.fixture
def synthetic_csv(plot_metrics, tmp_path: Path):
    """Three models × two seeds × celltype label, batch metrics populated.

    Deliberately omits the canonical bio_score/batch_score_scib/composite_scib
    columns so aggregate_seeds exercises its scIB-set fallback derivation.
    """
    rows = []
    rng = np.random.default_rng(0)
    for model, base_bio, base_batch in [("pca", 0.40, 0.55), ("scvi", 0.70, 0.80),
                                        ("blueprint", 0.85, 0.75)]:
        for seed in (0, 1):
            rows.append({
                "model": model, "dataset": "human_pancreas", "seed": seed,
                "label_col": "celltype", "n_cells": 1000, "embedding_dim": 10,
                "walltime_s": 1.0, "participation_ratio": 5.0,
                "ari":   base_bio + rng.normal(scale=0.02),
                "nmi":   base_bio + rng.normal(scale=0.02),
                "asw_label": base_bio + rng.normal(scale=0.02),
                "clisi":     base_bio + rng.normal(scale=0.02),
                "isolated_label_f1": base_bio + rng.normal(scale=0.02),
                "asw_batch":         base_batch + rng.normal(scale=0.02),
                "ilisi":             base_batch + rng.normal(scale=0.02),
                "kbet_acceptance":   base_batch + rng.normal(scale=0.02),
                "graph_connectivity": base_batch + rng.normal(scale=0.02),
                "pcr":               base_batch + rng.normal(scale=0.02),
            })
    csv_path = tmp_path / "metrics.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    return csv_path


@pytest.fixture
def synthetic_no_batch_csv(tmp_path: Path):
    """mouse_brain-shaped fixture: no batch metric columns finite."""
    rows = []
    rng = np.random.default_rng(0)
    for model, base in [("pca", 0.4), ("expimap", 0.7), ("blueprint", 0.85)]:
        for seed in (0, 1, 2):
            rows.append({
                "model": model, "dataset": "mouse_brain", "seed": seed,
                "label_col": "subclass", "n_cells": 1000, "embedding_dim": 10,
                "walltime_s": 1.0, "participation_ratio": 5.0,
                "ari": base + rng.normal(scale=0.02),
                "nmi": base + rng.normal(scale=0.02),
                "asw_label": base + rng.normal(scale=0.02),
                "clisi": base + rng.normal(scale=0.02),
                "isolated_label_f1":  base + rng.normal(scale=0.02),
                "asw_batch": np.nan, "ilisi": np.nan, "kbet_acceptance": np.nan,
                "graph_connectivity": np.nan, "pcr": np.nan,
            })
    csv_path = tmp_path / "metrics.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    return csv_path


def test_aggregate_seeds_collapses_per_model(plot_metrics, synthetic_csv):
    df = pd.read_csv(synthetic_csv)
    agg = plot_metrics.aggregate_seeds(df)
    assert set(agg["model"]) == {"pca", "scvi", "blueprint"}
    assert (agg["n_seeds"] == 2).all()
    # blueprint had highest base_bio AND base_batch → should win composite.
    top = agg.sort_values("composite_scib", ascending=False).iloc[0]
    assert top["model"] == "blueprint"


def test_aggregate_seeds_composite_uses_60_40_when_batch_present(plot_metrics, synthetic_csv):
    df = pd.read_csv(synthetic_csv)
    agg = plot_metrics.aggregate_seeds(df)
    row = agg.iloc[0]
    expected = 0.6 * row["bio_score"] + 0.4 * row["batch_score_scib"]
    assert row["composite_scib"] == pytest.approx(expected, rel=1e-9)


def test_aggregate_seeds_composite_falls_back_to_bio_alone_without_batch(
    plot_metrics, synthetic_no_batch_csv
):
    df = pd.read_csv(synthetic_no_batch_csv)
    agg = plot_metrics.aggregate_seeds(df)
    assert agg["batch_score_scib"].isna().all()
    # composite should equal bio_score in that case.
    np.testing.assert_allclose(agg["composite_scib"], agg["bio_score"], rtol=1e-9)


def test_resolve_label_cols_rejects_missing(plot_metrics, synthetic_csv):
    df = pd.read_csv(synthetic_csv)
    agg = plot_metrics.aggregate_seeds(df)
    args = type("A", (), {"label_cols": "celltype,nonsense", "dataset": "human_pancreas"})()
    with pytest.raises(ValueError, match="not present"):
        plot_metrics.resolve_label_cols(agg, args)


def test_resolve_label_cols_default_uses_dataset_params(plot_metrics, synthetic_csv):
    df = pd.read_csv(synthetic_csv)
    agg = plot_metrics.aggregate_seeds(df)
    args = type("A", (), {"label_cols": None, "dataset": "human_pancreas"})()
    assert plot_metrics.resolve_label_cols(agg, args) == ["celltype"]


def test_draw_heatmap_writes_png_and_pdf(plot_metrics, synthetic_csv, tmp_path):
    df = pd.read_csv(synthetic_csv)
    agg = plot_metrics.aggregate_seeds(df)
    agg_label = agg[agg["label_col"] == "celltype"]
    plot_metrics.draw_heatmap(agg_label, "human_pancreas", "celltype", tmp_path / "h")
    assert (tmp_path / "h.png").stat().st_size > 1000
    assert (tmp_path / "h.pdf").stat().st_size > 500


def test_draw_bars_writes_png_and_pdf(plot_metrics, synthetic_csv, tmp_path):
    df = pd.read_csv(synthetic_csv)
    agg = plot_metrics.aggregate_seeds(df)
    agg_label = agg[agg["label_col"] == "celltype"]
    plot_metrics.draw_bars(agg_label, "human_pancreas", "celltype", tmp_path / "b")
    assert (tmp_path / "b.png").stat().st_size > 1000


def test_draw_tradeoff_skips_when_no_batch(plot_metrics, synthetic_no_batch_csv, tmp_path):
    df = pd.read_csv(synthetic_no_batch_csv)
    agg = plot_metrics.aggregate_seeds(df)
    agg_label = agg[agg["label_col"] == "subclass"]
    plot_metrics.draw_tradeoff(agg_label, "mouse_brain", "subclass", tmp_path / "t")
    assert not (tmp_path / "t.png").exists()
    assert not (tmp_path / "t.pdf").exists()


def test_draw_tradeoff_writes_when_batch_present(plot_metrics, synthetic_csv, tmp_path):
    df = pd.read_csv(synthetic_csv)
    agg = plot_metrics.aggregate_seeds(df)
    agg_label = agg[agg["label_col"] == "celltype"]
    plot_metrics.draw_tradeoff(agg_label, "human_pancreas", "celltype", tmp_path / "t")
    assert (tmp_path / "t.png").stat().st_size > 1000
