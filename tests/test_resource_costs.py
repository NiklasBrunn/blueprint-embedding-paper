"""Smoke + behavioural tests for the resource-cost pipeline.

Covers scripts/agentic_embedding/extract_cost.py (token / timing / footprint
introspection of an agent run dir), scripts/resource_costs/collect.py
(report-regex parsing — both the new `## Resources` block and the legacy
`## Stats` fallback — + agent-run aggregation into the long CSV), and
scripts/plots/plot_resource_costs.py (figures render and emit PNG + PDF) on
small synthetic fixtures — no heavy single-cell deps.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

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
def extract_cost():
    return _load("extract_cost", "scripts/agentic_embedding/extract_cost.py")


@pytest.fixture(scope="module")
def collect():
    return _load("rc_collect", "scripts/resource_costs/collect.py")


@pytest.fixture(scope="module")
def plotter():
    return _load("rc_plot", "scripts/plots/plot_resource_costs.py")


# --------------------------------------------------------------- extract_cost
_NEW_TIMING = (
    "[timing] walltime_s=2.5100 load_data_s=0.4000 obtain_model_s=0.0010 "
    "build_embedding_s=2.0600 other_s=0.0490 model_source=authored model_download_s=0.0 "
    "peak_rss_mb=205.4 rss_after_load_mb=120.0 rss_after_obtain_mb=121.0 rss_after_build_mb=205.4 "
    "model_storage_bytes=1100 embedding_bytes=99999 device=cpu host=testhost\n"
)


def _make_agent_run(tmp_path: Path, *, timing_line: str = _NEW_TIMING) -> Path:
    cohort = tmp_path / "work" / "gpt_5_5_high_mouse_brain"
    run = cohort / "260512_blueprint_run01"
    run.mkdir(parents=True)
    (run / "embedder_blueprint.py").write_text("# populated blueprint\n" * 50)  # ~1 KB
    (run / "codex_edit.log").write_text(
        "OpenAI Codex v0.1\n...transcript...\ntokens used\n12.345\nEdited build_blueprint().\n"
    )
    (run / ".stage_a_rc").write_text("0\n")
    (run / ".stage_a_secs").write_text("180\n")
    (run / ".embed_secs").write_text("4\n")
    (run / ".embed_status").write_text("OK\n")
    (run / ".embed_detail").write_text("outs/embeddings/mouse_brain/blueprint_x_seed1.parquet\n")
    (run / "embed.log").write_text(
        "[coverage] module 'a': 3/3 genes present\nWrote out.parquet\n" + timing_line
    )
    return run


def test_extract_cost_basic(extract_cost, tmp_path):
    run = _make_agent_run(tmp_path)
    cost = extract_cost._build_cost(run)
    assert cost["variant"] == "blueprint"
    assert cost["agent_cli"] == "codex"
    assert cost["agent_tokens_stage_a"] == 12345          # "12.345" -> 12345
    # Stage A only — no refinement transcripts, no Stage B fields.
    assert cost["agent_tokens_stage_b_refine"] is None
    assert cost["agent_tokens_total"] == 12345
    assert cost["agent_tokens_note"] is None
    assert cost["stage_a_s"] == 180
    assert cost["stage_b_s"] is None
    assert cost["embed_s"] == 4
    assert cost["stage_b_status"] is None
    assert cost["stage_b_iters"] is None
    assert cost["stage_b_final_rating"] is None
    # New 4-bucket [timing] line:
    assert cost["embed_walltime_s"] == 2.51
    assert cost["embed_load_data_s"] == 0.4
    assert cost["embed_obtain_model_s"] == 0.001
    assert cost["embed_build_embedding_s"] == 2.06
    assert cost["embed_other_s"] == 0.049
    assert cost["embed_model_source"] == "authored"
    assert cost["embed_model_download_s"] == 0.0
    assert cost["embed_peak_rss_mb"] == 205.4
    assert cost["embed_rss_after_load_mb"] == 120.0 and cost["embed_rss_after_build_mb"] == 205.4
    assert cost["embed_device"] == "cpu" and cost["embed_host"] == "testhost"
    assert cost["agent_compute_device"] == "remote-llm"
    assert cost["embedding_bytes"] == 99999               # from the [timing] line (parquet path absent)
    assert cost["model_storage_bytes"] == (run / "embedder_blueprint.py").stat().st_size

    # The CLI entrypoint writes cost.json next to the run dir.
    import subprocess
    subprocess.run([sys.executable, str(REPO_ROOT / "scripts/agentic_embedding/extract_cost.py"), str(run)],
                   check=True, capture_output=True)
    out = run / "cost.json"
    assert out.exists()
    assert json.loads(out.read_text())["agent_tokens_total"] == 12345


def test_extract_cost_legacy_timing_line(extract_cost, tmp_path):
    """A pre-taxonomy embed.log (fit_s/transform_s) still maps onto the new keys."""
    legacy = ("[timing] fit_s=0.01 transform_s=2.50 walltime_s=2.51 peak_rss_mb=205.4 "
              "model_storage_bytes=1100 embedding_bytes=99999 device=cpu host=testhost\n")
    run = _make_agent_run(tmp_path, timing_line=legacy)
    cost = extract_cost._build_cost(run)
    assert cost["embed_obtain_model_s"] == 0.01      # <- fit_s
    assert cost["embed_build_embedding_s"] == 2.5    # <- transform_s
    assert cost["embed_walltime_s"] == 2.51
    assert cost["embed_load_data_s"] is None         # not in the legacy line
    assert cost["embed_model_source"] == "authored"  # default for an agent embed.log


def test_extract_cost_no_tokens(extract_cost, tmp_path):
    cohort = tmp_path / "work" / "claude_opus_4_7_high_mouse_brain"
    run = cohort / "260512_blueprint_run02"
    run.mkdir(parents=True)
    (run / "embedder_blueprint.py").write_text("# bp\n")
    (run / "agent_edit.log").write_text("...claude transcript with no token footer...\n")
    cost = extract_cost._build_cost(run)
    assert cost["agent_tokens_total"] is None
    assert "no 'tokens used'" in cost["agent_tokens_note"]
    assert cost["agent_cli"] == "claude"


# Stage A sibling embed: a refine workdir snapshotted embedder_blueprint_unrefined.py
# at refine-loop entry and Phase 6 produced a second parquet from it. Distinct
# timing values let us verify each branch is parsed from its own log file.
_UNREFINED_TIMING = (
    "[timing] walltime_s=1.7700 load_data_s=0.3500 obtain_model_s=0.0009 "
    "build_embedding_s=1.3700 other_s=0.0491 model_source=authored model_download_s=0.0 "
    "peak_rss_mb=180.2 rss_after_load_mb=110.0 rss_after_obtain_mb=110.5 rss_after_build_mb=180.2 "
    "model_storage_bytes=900 embedding_bytes=88888 device=cpu host=testhost\n"
)


def _make_refine_run(tmp_path: Path) -> Path:
    """Refine workdir under the new convention (.stage_b_iters present), with
    both refined and Stage A artifacts. Cohort slug `_refine` is what triggers
    the `_unrefined` sibling row in collect.agent_rows."""
    cohort = tmp_path / "work" / "gpt_5_5_refine_high_mouse_brain"
    run = cohort / "260520_blueprint_run01"
    logs = run / "logs"
    logs.mkdir(parents=True)
    (run / "embedder_blueprint.py").write_text("# refined blueprint\n" * 80)
    (run / "embedder_blueprint_unrefined.py").write_text("# stage A snapshot\n" * 60)
    (run / "scenario_config.jsonc").write_text("{}\n")
    (logs / "codex_edit.log").write_text(
        "OpenAI Codex\n...\ntokens used\n50.000\nDone.\n"
    )
    (logs / "codex_discriminate_01.log").write_text("...\ntokens used\n10.000\nDone.\n")
    (logs / "codex_generate_01.log").write_text("...\ntokens used\n8.000\nDone.\n")
    (logs / "embed.log").write_text("[coverage] ok\n" + _NEW_TIMING)
    (logs / "embed_unrefined.log").write_text("[coverage] ok\n" + _UNREFINED_TIMING)
    (run / ".stage_a_rc").write_text("0\n")
    (run / ".stage_a_secs").write_text("210\n")
    (run / ".stage_b_secs").write_text("75\n")
    (run / ".stage_b_status").write_text("OK 1 fine\n")
    (run / ".stage_b_iters").write_text("1\n")  # triggers new-convention parsing
    (run / ".stage_b_final_rating").write_text("fine\n")
    (run / ".embed_secs").write_text("3\n")
    (run / ".embed_status").write_text("OK\n")
    (run / ".embed_detail").write_text(
        "outs/embeddings/mouse_brain/"
        "blueprint_gpt_5_5_refine_high_mouse_brain_seed1.parquet\n"
    )
    (run / ".embed_unrefined_secs").write_text("2\n")
    (run / ".embed_unrefined_status").write_text("OK\n")
    (run / ".embed_unrefined_detail").write_text(
        "outs/embeddings/mouse_brain/"
        "blueprint_gpt_5_5_unrefined_high_mouse_brain_seed1.parquet\n"
    )
    return run


def test_extract_cost_emits_unrefined_keys(extract_cost, tmp_path):
    run = _make_refine_run(tmp_path)
    cost = extract_cost._build_cost(run)
    # Refined branch is unaffected by the Stage A addition.
    assert cost["embed_walltime_s"] == 2.51
    assert cost["embed_build_embedding_s"] == 2.06
    # Stage A sibling parsed from logs/embed_unrefined.log.
    assert cost["embed_unrefined_walltime_s"] == 1.77
    assert cost["embed_unrefined_load_data_s"] == 0.35
    assert cost["embed_unrefined_build_embedding_s"] == 1.37
    assert cost["embed_unrefined_peak_rss_mb"] == 180.2
    assert cost["embed_unrefined_device"] == "cpu"
    assert cost["embed_unrefined_s"] == 2
    assert cost["embed_unrefined_status"] == "OK"
    assert cost["unrefined_model_file"] == "embedder_blueprint_unrefined.py"
    assert (
        cost["unrefined_model_storage_bytes"]
        == (run / "embedder_blueprint_unrefined.py").stat().st_size
    )


def test_collect_emits_dual_rows_for_refine_workdir(collect, tmp_path):
    repo = tmp_path
    _make_refine_run(repo)
    collect.REPO_ROOT = repo
    collect.EMB_ROOT = repo / "outs" / "embeddings"
    collect.WORK_ROOT = repo / "work"
    collect.METRICS_ROOT = repo / "outs" / "metrics"

    rows = collect.agent_rows("mouse_brain")
    by_method = {r["method"]: r for r in rows}
    assert set(by_method) == {
        "blueprint_gpt_5_5_refine_high",
        "blueprint_gpt_5_5_unrefined_high",
    }, "refine workdir should emit one refined + one Stage A row"

    refined = by_method["blueprint_gpt_5_5_refine_high"]
    unrefined = by_method["blueprint_gpt_5_5_unrefined_high"]

    # Refined row keeps the full Stage A + Stage B picture.
    assert refined["family"] == "blueprint"
    assert refined["agent_stage_a_s"] == 210
    assert refined["agent_stage_b_refine_s"] == 75
    assert refined["agent_stage_b_iters"] == 1
    assert refined["agent_tokens_stage_a"] == 50000
    assert refined["agent_tokens_stage_b_refine"] == 18000  # 10000 + 8000
    assert refined["agent_tokens_total"] == 50000 + 18000
    assert refined["walltime_s"] == 2.51  # from embed.log
    assert refined["cohort"] == "gpt_5_5_refine_high_mouse_brain"

    # Stage A row shares Stage A authoring cost but nulls refinement fields and
    # carries the Stage A embed timings (from embed_unrefined.log).
    assert unrefined["family"] == "blueprint"
    assert unrefined["variant"] == "blueprint"
    assert unrefined["seed"] == refined["seed"]
    assert unrefined["cohort"] == "gpt_5_5_refine_high_mouse_brain"
    assert unrefined["agent_model"] == "gpt_5_5_unrefined"
    assert unrefined["agent_stage_a_s"] == 210
    assert unrefined["agent_tokens_stage_a"] == 50000
    assert unrefined["agent_stage_b_refine_s"] is None
    assert unrefined["agent_stage_b_iters"] is None
    assert unrefined["agent_stage_b_final_rating"] is None
    assert unrefined["agent_tokens_stage_b_refine"] is None
    assert unrefined["agent_tokens_total"] == 50000  # Stage A only
    assert unrefined["walltime_s"] == 1.77            # from embed_unrefined.log
    assert unrefined["build_embedding_s"] == 1.37
    assert unrefined["embedding_bytes"] == 88888
    assert unrefined["model_source"] == "authored"
    assert unrefined["status"].startswith("refine=n/a")


# ------------------------------------------------------------------- collect
def _make_repo_fixture(tmp_path: Path) -> Path:
    """A tiny stand-in repo: model reports (new `## Resources` block for pca/scvi,
    legacy `## Stats` for scvi_census) + one agent cohort."""
    ds = "mouse_brain"
    emb = tmp_path / "outs" / "embeddings" / ds
    emb.mkdir(parents=True)
    # pca — NEW format. walltime = 0.30 + 0.50 + 0.40 + 2.01 = 3.21
    (emb / f"pca_{ds}_report.md").write_text(
        "# PCA embedding report — mouse_brain\n\n## Stats\n"
        "- cells: 1809\n- input genes (after HVG mask): 2000\n- embedding dim: 50\n\n"
        "## Resources\n"
        "- walltime: 3.21 s   (= load 0.30 s + obtain-model 0.50 s + build-embedding 0.40 s + other 2.010 s)\n"
        "- load (read input data from disk): 0.30 s\n"
        "- obtain model — compute PC loadings + save model [trained]: 0.50 s\n"
        "- build embedding (project cells onto the loadings): 0.40 s\n"
        "- other (data prep + parquet write + glue): 2.010 s\n"
        "- peak RSS: 512.5 MB   (RSS trajectory: after-load 200.0 → after-model 300.0 → after-embedding 400.0 MB)\n"
        "- model storage (pca_model.npz): 80000 bytes\n"
        "- model download: 0.00 s (0 — model already cached on disk / nothing to download)\n"
        "- embedding parquet: 12345 bytes\n- device: cpu\n- host: testhost\n"
    )
    # scvi — NEW format. walltime = 1.00 + 92.10 + 3.30 + 0.00 = 96.40. (The
    # "what" labels here deliberately contain nested parens — regression for the
    # collect.py regex that used to choke on "...scores (rank scoring)...".)
    (emb / f"scvi_{ds}_seed0_report.md").write_text(
        "# scVI embedding report — mouse_brain\n\n## Stats\n"
        "- cells: 1809\n- input genes (HVGs): 2000\n- embedding dim: 10\n\n"
        "## Resources\n"
        "- walltime: 96.40 s   (= load 1.00 s + obtain-model 92.10 s + build-embedding 3.30 s + other 0.000 s)\n"
        "- load (read input data from disk): 1.00 s\n"
        "- obtain model — setup_anndata + construct (heuristic max_epochs) + train + save [trained]: 92.10 s\n"
        "- build embedding (encoder forward pass (mean of the latent posterior)): 3.30 s\n"
        "- other (data prep + parquet write + glue): 0.000 s\n"
        "- peak RSS: 1801.0 MB   (RSS trajectory: after-load 300.0 → after-model 1700.0 → after-embedding 1801.0 MB)\n"
        "- model storage (scVI checkpoint dir): 2500000 bytes\n"
        "- model download: 0.00 s (0 — model already cached on disk / nothing to download)\n"
        "- embedding parquet: 222000 bytes\n- device: cpu\n- host: testhost\n"
    )
    # scvi_census — LEGACY format (exercises the fallback + model_source inference).
    (emb / f"scvi_census_{ds}_report.md").write_text(
        "# scVI Census embedding report — mouse_brain\n\n## Stats\n"
        "- cells: 1809\n- embedding dim: 50\n- walltime: 14.7 s\n"
        "- fit (download + load pretrained checkpoint + query adaptation): 12.50 s\n"
        "- transform (encoder forward, data → embedding): 0.15 s\n"
        "- peak RSS: 2900.0 MB\n- model storage (pretrained model.pt): 196379229 bytes\n"
        "- embedding parquet: 478210 bytes\n- device: cpu\n- host: testhost\n"
    )
    # blueprint_* per-stem report rendered by extract_run.py must be ignored here.
    (emb / f"blueprint_gpt_5_5_high_{ds}_seed1_report.md").write_text("# Blueprint Report\n")
    _make_agent_run(tmp_path)  # work/gpt_5_5_high_mouse_brain/260512_blueprint_run01/
    return tmp_path


def test_collect_build_df(collect, tmp_path):
    repo = _make_repo_fixture(tmp_path)
    # Re-point the module's data roots at the fixture repo (REPO_ROOT is left
    # alone so _load_extract_cost still finds the real extract_cost.py).
    collect.REPO_ROOT = repo
    collect.EMB_ROOT = repo / "outs" / "embeddings"
    collect.WORK_ROOT = repo / "work"
    collect.OUT_ROOT = repo / "outs" / "resource_costs"
    collect.METRICS_ROOT = repo / "outs" / "metrics"   # absent -> no embedding_dim backfill
    df = collect.build_dataset_df("mouse_brain")

    methods = set(df["method"])
    assert {"pca", "scvi", "scvi_census", "blueprint_gpt_5_5_high"} <= methods
    assert "load_s" not in df.columns and "fit_s" not in df.columns and "transform_s" not in df.columns

    pca = df[df["method"] == "pca"].iloc[0]
    assert pca["walltime_s"] == 3.21 and pca["peak_rss_mb"] == 512.5
    assert pca["load_data_s"] == 0.30 and pca["obtain_model_s"] == 0.50
    assert pca["build_embedding_s"] == 0.40 and pca["other_s"] == 2.01
    assert pca["model_source"] == "trained"
    assert pca["rss_after_load_mb"] == 200.0 and pca["rss_after_build_mb"] == 400.0
    assert pca["model_storage_bytes"] == 80000 and pca["embedding_bytes"] == 12345
    assert pca["device"] == "cpu" and pca["host"] == "testhost"
    # the four buckets sum to walltime
    assert abs((pca["load_data_s"] + pca["obtain_model_s"] + pca["build_embedding_s"]
                + pca["other_s"]) - pca["walltime_s"]) < 1e-6

    scvi = df[df["method"] == "scvi"].iloc[0]
    assert scvi["obtain_model_s"] == 92.1 and scvi["build_embedding_s"] == 3.3
    assert scvi["model_source"] == "trained" and scvi["model_storage_bytes"] == 2500000

    # Legacy report: old "fit (...)" -> obtain_model_s; model_source inferred as "loaded".
    census = df[df["method"] == "scvi_census"].iloc[0]
    assert census["obtain_model_s"] == 12.5 and census["build_embedding_s"] == 0.15
    assert census["model_source"] == "loaded"
    assert pd.isna(census["load_data_s"])  # not present in the legacy block

    bp = df[df["method"] == "blueprint_gpt_5_5_high"].iloc[0]
    assert bp["family"] == "blueprint" and bp["model_source"] == "authored"
    # The fixture run has Stage A only; agent_tokens_total = Stage A tokens
    # and agent_stage_b_refine_s is NaN.
    assert bp["agent_tokens_total"] == 12345
    assert bp["agent_stage_a_s"] == 180 and pd.isna(bp["agent_stage_b_refine_s"])
    assert bp["walltime_s"] == 2.51 and bp["load_data_s"] == 0.4 and bp["build_embedding_s"] == 2.06
    assert bp["device"] == "cpu"  # the local Phase-5 device

    # Summary + markdown render without error.
    summ = collect.summarize(df)
    assert {"family", "method", "n", "walltime_s_median", "obtain_model_s_median"} <= set(summ.columns)
    md_path = tmp_path / "report.md"
    collect.write_markdown(summ, md_path, "test")
    txt = md_path.read_text()
    assert "Cost / wall-clock / storage" not in txt  # we passed a custom title
    assert "| method |" in txt and "obtain-model" in txt


# -------------------------------------------------------------------- plotter
def _row(**kw):
    base = dict(seed=None, run_id="r", cohort=None, agent_model=None, reasoning=None,
                variant="conventional", device="cpu", host="h1", load_data_s=np.nan, obtain_model_s=np.nan,
                model_source="", build_embedding_s=np.nan, other_s=np.nan, model_download_s=0.0,
                peak_rss_mb=np.nan, rss_after_load_mb=np.nan, rss_after_obtain_mb=np.nan,
                rss_after_build_mb=np.nan, model_storage_bytes=np.nan, embedding_bytes=np.nan,
                embedding_dim=np.nan, agent_stage_a_s=np.nan, agent_stage_b_refine_s=np.nan,
                agent_stage_b_iters=np.nan, agent_stage_b_final_rating="",
                agent_tokens_stage_a=np.nan, agent_tokens_stage_b_refine=np.nan, agent_tokens_total=np.nan,
                status="ok", source="report")
    base.update(kw)
    return base


def test_plot_runs(plotter, tmp_path):
    ds = "mouse_brain"
    costs_dir = tmp_path / "outs" / "resource_costs" / ds
    costs_dir.mkdir(parents=True)
    rows = [
        _row(method="pca", family="conventional", dataset=ds, model_source="trained",
             walltime_s=3.0, load_data_s=0.3, obtain_model_s=0.5, build_embedding_s=0.4, other_s=1.8,
             peak_rss_mb=500.0, rss_after_load_mb=200.0, rss_after_obtain_mb=300.0, rss_after_build_mb=400.0,
             model_storage_bytes=80000, embedding_bytes=12345, embedding_dim=50),
        _row(method="scvi", family="conventional", dataset=ds, seed=0, model_source="trained",
             walltime_s=96.0, load_data_s=1.0, obtain_model_s=92.0, build_embedding_s=3.0, other_s=0.0,
             peak_rss_mb=1800.0, rss_after_load_mb=300.0, rss_after_obtain_mb=1700.0, rss_after_build_mb=1800.0,
             model_storage_bytes=2_500_000, embedding_bytes=222000, embedding_dim=10),
        _row(method="blueprint_gpt_5_5_high", family="blueprint", dataset=ds, seed=1, source="cost",
             variant="blueprint", agent_model="gpt_5_5", reasoning="high",
             model_source="authored", walltime_s=2.5, load_data_s=0.4, obtain_model_s=0.001,
             build_embedding_s=2.0, other_s=0.099, peak_rss_mb=200.0, rss_after_load_mb=120.0,
             rss_after_obtain_mb=121.0, rss_after_build_mb=200.0, model_storage_bytes=30000,
             embedding_bytes=400000, embedding_dim=23, agent_stage_a_s=300, agent_stage_b_refine_s=140,
             agent_tokens_stage_a=70000, agent_tokens_stage_b_refine=40000, agent_tokens_total=110000),
        _row(method="blueprint_gpt_5_5_high", family="blueprint", dataset=ds, seed=2, source="cost",
             variant="blueprint", agent_model="gpt_5_5", reasoning="high",
             model_source="authored", walltime_s=2.6, load_data_s=0.4, obtain_model_s=0.001,
             build_embedding_s=2.1, other_s=0.099, peak_rss_mb=210.0, rss_after_load_mb=120.0,
             rss_after_obtain_mb=121.0, rss_after_build_mb=210.0, model_storage_bytes=31000,
             embedding_bytes=410000, embedding_dim=24, agent_stage_a_s=320, agent_stage_b_refine_s=150,
             agent_tokens_stage_a=72000, agent_tokens_stage_b_refine=42000, agent_tokens_total=114000),
    ]
    pd.DataFrame(rows).to_csv(costs_dir / "costs.csv", index=False)
    plotter.COSTS_ROOT = tmp_path / "outs" / "resource_costs"
    plotter.METRICS_ROOT = tmp_path / "outs" / "metrics"   # absent -> quality axis falls back to n
    plotter.OUT_ROOT = tmp_path / "outs" / "plots" / "resource_costs"
    plotter.matplotlib.rcParams.update(plotter.PLOT_RC)
    plotter.run_dataset(ds, overwrite=True)
    out = plotter.OUT_ROOT / ds
    for stem in ("tradeoff_total_vs_quality", "tradeoff_obtain_model_vs_quality",
                 "tradeoff_build_embedding_vs_quality", "time_breakdown", "memory_storage", "agent_cost"):
        assert (out / f"{stem}.png").exists() and (out / f"{stem}.pdf").exists()

    # Derived analysis columns. All time axes EXCLUDE the disk read (load_data_s).
    # Default: the optional Stage B refinement loop is EXCLUDED.
    df = plotter._load_costs(ds)
    assert {"total_s", "obtain_s_plot", "build_s_plot", "walltime_xl_s"} <= set(df.columns)
    scvi_row = df[df["method"] == "scvi"].iloc[0]
    assert scvi_row["walltime_xl_s"] == 95.0           # 96.0 walltime − 1.0 disk read
    assert scvi_row["obtain_s_plot"] == 92.0           # de-novo: fit/train time (no disk read in it)
    assert scvi_row["build_s_plot"] == 3.0             # trained -> forward pass only
    assert scvi_row["total_s"] == 95.0                 # = walltime − load_data
    bp_row = df[df["family"] == "blueprint"].iloc[0]
    assert bp_row["obtain_s_plot"] == 300.0            # Stage A only (Stage B refinement excluded by default)
    # total = Stage A + (Phase-5 walltime − its disk read) = 300 + (2.5 − 0.4) = 302.1
    assert abs(bp_row["total_s"] - (300.0 + bp_row["walltime_s"] - bp_row["load_data_s"])) < 1e-6
    # ...included with --include-refine:
    df_b = plotter._load_costs(ds, include_refine=True)
    bp_row_b = df_b[df_b["family"] == "blueprint"].iloc[0]
    assert bp_row_b["obtain_s_plot"] == 440.0          # Stage A + Stage B refinement  (300 + 140)
