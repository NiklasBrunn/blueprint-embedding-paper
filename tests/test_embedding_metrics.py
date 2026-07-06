"""
Smoke + behavioural tests for scripts/embedding_metrics/metrics.py.

The pure-Python helpers (parse_stem, prettify_model, participation_ratio,
read_sibling_report, merge_csv, the CSV schema) are tested without heavy deps.
The end-to-end scib-metrics Benchmarker path (build_obsm -> run_benchmarker ->
results_to_rows) is exercised on a tiny synthetic dataset behind @needs_scib,
so the file stays a useful CI sanity check even where scib_metrics is absent.
"""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]
SCIB_AVAILABLE = importlib.util.find_spec("scib_metrics") is not None


def _load(name: str, rel_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / rel_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def metrics():
    return _load("embedding_metrics", "scripts/embedding_metrics/metrics.py")


needs_scib = pytest.mark.skipif(not SCIB_AVAILABLE, reason="scib_metrics not installed")


# --------------------------------------------------------------------
# Pure-Python helpers
# --------------------------------------------------------------------
def test_parse_stem_seeded_and_unseeded(metrics):
    assert metrics.parse_stem("pca_mouse_brain_seed3", "mouse_brain") == ("pca", 3)
    assert metrics.parse_stem("ucell_human_pancreas", "human_pancreas") == ("ucell", None)
    assert metrics.parse_stem("scvi_census_mouse_brain", "mouse_brain") == ("scvi_census", None)


def test_prettify_model_baselines(metrics):
    assert metrics.prettify_model("pca") == "PCA"
    assert metrics.prettify_model("scvi") == "scVI"
    assert metrics.prettify_model("scvi_census") == "scVI-census"
    assert metrics.prettify_model("expimap") == "expiMap"


def test_prettify_model_agentic_abbreviations(metrics):
    # noref->N, unrefined->U, refine->R
    assert metrics.prettify_model("blueprint_gpt_5_5_noref_unrefined_high") == "BP GPT-5.5 N U"
    assert metrics.prettify_model("blueprint_gpt_5_5_noref_refine_high") == "BP GPT-5.5 N R"
    assert metrics.prettify_model("blueprint_gpt_5_5_unrefined_high") == "BP GPT-5.5 U"
    assert metrics.prettify_model("blueprint_gpt_5_5_refine_high") == "BP GPT-5.5 R"
    assert metrics.prettify_model("freeform_naive_gpt_5_5_high") == "FF naïve GPT-5.5"
    assert metrics.prettify_model("freeform_noref_gpt_5_5_high") == "FF GPT-5.5 N"


def test_prettify_model_unknown_passthrough(metrics):
    assert metrics.prettify_model("some_other_model") == "some_other_model"


def test_csv_columns_schema(metrics):
    cols = set(metrics.CSV_COLUMNS)
    # New scIB metrics + canonical scores present.
    assert {"clisi", "pcr", "bio_score", "batch_score_scib", "composite_scib"} <= cols
    # Dropped non-scIB extras and _ext variants are gone.
    assert not ({"knn_acc", "knn_batch_entropy", "knn_mixing_rate",
                 "batch_score_ext", "composite_ext"} & cols)


def test_discover_embeddings_filters_umap_coords(metrics, tmp_path: Path):
    (tmp_path / "a_mouse_brain.parquet").write_bytes(b"")
    (tmp_path / "a_mouse_brain_umap_coords.parquet").write_bytes(b"")
    found = metrics.discover_embeddings(tmp_path, override=None)
    assert [p.name for p in found] == ["a_mouse_brain.parquet"]


def test_participation_ratio_isotropic_close_to_d(metrics):
    rng = np.random.default_rng(1)
    X = rng.normal(size=(2000, 7)).astype(np.float32)
    pr = metrics.participation_ratio(X)
    assert 5.5 < pr < 8.0


def test_participation_ratio_rank_one_close_to_one(metrics):
    rng = np.random.default_rng(2)
    direction = rng.normal(size=(8,))
    direction /= np.linalg.norm(direction)
    scale = rng.normal(size=(500, 1))
    X = (scale @ direction[None, :]).astype(np.float32)
    pr = metrics.participation_ratio(X)
    assert pr < 1.5


def test_read_sibling_report_missing_returns_nan(metrics, tmp_path: Path):
    out = metrics.read_sibling_report(tmp_path / "not_real.parquet")
    assert np.isnan(out["walltime_s"])
    assert out["report_path"] is None


def test_read_sibling_report_parses_walltime(metrics, tmp_path: Path):
    parquet = tmp_path / "scvi_mouse_brain_seed0.parquet"
    parquet.write_bytes(b"")
    md = tmp_path / "scvi_mouse_brain_seed0_report.md"
    md.write_text("# scVI embedding report\n- walltime: 12.34 s\n")
    out = metrics.read_sibling_report(parquet)
    assert out["walltime_s"] == pytest.approx(12.34)


def test_load_embedding_returns_matrix_and_ids(metrics, tmp_path: Path):
    p = tmp_path / "emb_testds_seed0.parquet"
    cols = {f"dim_{i}": np.arange(3, dtype=float) + i for i in range(4)}
    pq.write_table(pa.Table.from_pydict({"cell_id": ["a", "b", "c"], **cols}), p)
    X, ids = metrics.load_embedding(p)
    assert X.shape == (3, 4)
    assert ids == ["a", "b", "c"]


def test_merge_csv_replaces_matching_rows(metrics, tmp_path: Path):
    csv = tmp_path / "metrics.csv"
    pd.DataFrame([
        {c: None for c in metrics.CSV_COLUMNS} | {
            "model": "pca", "dataset": "mouse_brain", "seed": 0, "label_col": "subclass", "ari": 0.1,
        },
        {c: None for c in metrics.CSV_COLUMNS} | {
            "model": "scvi", "dataset": "mouse_brain", "seed": 0, "label_col": "subclass", "ari": 0.2,
        },
    ]).to_csv(csv, index=False)
    new_rows = [{c: None for c in metrics.CSV_COLUMNS} | {
        "model": "pca", "dataset": "mouse_brain", "seed": 0, "label_col": "subclass", "ari": 0.9,
    }]
    df = metrics.merge_csv(csv, new_rows)
    pca_row = df[(df["model"] == "pca") & (df["seed"] == 0) & (df["label_col"] == "subclass")]
    assert len(pca_row) == 1
    assert pca_row["ari"].iloc[0] == 0.9
    # untouched scvi row preserved
    assert ((df["model"] == "scvi") & (df["ari"] == 0.2)).any()


# --------------------------------------------------------------------
# End-to-end Benchmarker pipeline (synthetic, behind scib_metrics)
# --------------------------------------------------------------------
@pytest.fixture
def dataset_fixture(tmp_path: Path):
    """Three well-separated 8-dim blob embeddings (3 cell types, 2 batches) plus a
    PCA 'Unintegrated' reference, written as parquets + a pp h5ad — the layout
    build_obsm/run_benchmarker expect. Dataset token is 'testds'."""
    rng = np.random.default_rng(0)
    per, k, d = 40, 3, 8
    centres = np.zeros((k, d))
    centres[0, 0] = centres[1, 1] = centres[2, 2] = 8.0
    blocks, labels = [], []
    for ci in range(k):
        blocks.append(rng.normal(loc=centres[ci], scale=0.3, size=(per, d)))
        labels += [["A", "B", "C"][ci]] * per
    base = np.vstack(blocks).astype(np.float32)
    n = base.shape[0]
    cell_ids = [f"c{i}" for i in range(n)]
    labels = np.array(labels)
    batch = np.array((["t1", "t2"] * (n // 2 + 1))[:n])

    def _write(stem, X):
        p = tmp_path / f"{stem}.parquet"
        cols = {f"dim_{i}": X[:, i] for i in range(X.shape[1])}
        pq.write_table(pa.Table.from_pydict({"cell_id": cell_ids, **cols}), p)
        return p

    parquets = [
        _write("emba_testds_seed0", base + rng.normal(scale=0.05, size=base.shape).astype(np.float32)),
        _write("emba_testds_seed1", base + rng.normal(scale=0.05, size=base.shape).astype(np.float32)),
    ]
    pca = _write("pca_testds", base + rng.normal(scale=0.5, size=base.shape).astype(np.float32))

    h5ad = tmp_path / "testds_pp.h5ad"
    obs = pd.DataFrame({"celltype": labels, "batch": batch},
                       index=pd.Index(cell_ids, name="cell_id"))
    ad.AnnData(X=base.copy(), obs=obs).write_h5ad(h5ad)
    return {"parquets": parquets, "pca": pca, "h5ad": h5ad}


@needs_scib
def test_benchmarker_pipeline_finite_and_composite(metrics, dataset_fixture):
    import warnings
    obs_df = metrics.load_obs(dataset_fixture["h5ad"], ["celltype"], "batch")
    obsm, info, common, unint, _order = metrics.build_obsm(
        dataset_fixture["parquets"], obs_df, dataset_fixture["pca"], trajectory_col=None)
    assert len(common) == obs_df.shape[0]
    assert unint is not None  # PCA reference loaded

    obs_sub = obs_df[["celltype", "batch"]]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bm = metrics.run_benchmarker(obsm, common, obs_sub, "celltype", "batch",
                                     unint, n_neighbors=15, n_jobs=1)
        results = bm.get_results(min_max_scale=False, clean_names=True)
    rows = metrics.results_to_rows(results, "testds", "celltype", info, has_batch=True)

    assert len(rows) == len(dataset_fixture["parquets"])
    for r in rows:
        for k in ("ari", "nmi", "asw_label", "clisi", "isolated_label_f1"):
            assert np.isfinite(r[k]), f"{k} not finite"
        # At least the core batch metrics + PCR are present with a reference.
        for k in ("asw_batch", "ilisi", "kbet_acceptance", "graph_connectivity", "pcr"):
            assert np.isfinite(r[k]), f"{k} not finite"
        # Composite is the package Total = 0.6*bio + 0.4*batch.
        assert r["composite_scib"] == pytest.approx(
            0.6 * r["bio_score"] + 0.4 * r["batch_score_scib"], abs=1e-6)
        assert r["ari"] > 0.5  # well-separated blobs


@needs_scib
def test_benchmarker_batchless_composite_equals_bio(metrics, dataset_fixture):
    import warnings
    obs_df = metrics.load_obs(dataset_fixture["h5ad"], ["celltype"], batch_key=None)
    obsm, info, common, unint, _order = metrics.build_obsm(
        dataset_fixture["parquets"], obs_df, dataset_fixture["pca"], trajectory_col=None)
    obs_sub = obs_df[["celltype"]]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        bm = metrics.run_benchmarker(obsm, common, obs_sub, "celltype", None,
                                     unint, n_neighbors=15, n_jobs=1)
        results = bm.get_results(min_max_scale=False, clean_names=True)
    rows = metrics.results_to_rows(results, "testds", "celltype", info, has_batch=False)
    for r in rows:
        assert np.isnan(r["asw_batch"]) and np.isnan(r["ilisi"])
        assert r["composite_scib"] == pytest.approx(r["bio_score"], abs=1e-9)
