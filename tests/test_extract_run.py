"""Tests for scripts/agentic_embedding/extract_run.py — summary.json schema."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load_extract_module():
    spec = importlib.util.spec_from_file_location(
        "extract_run", REPO_ROOT / "scripts" / "agentic_embedding" / "extract_run.py"
    )
    module = importlib.util.module_from_spec(spec)
    sys.modules["extract_run"] = module
    spec.loader.exec_module(module)
    return module


def _make_run_dir(tmp_path: Path, body: str) -> Path:
    """Create a synthetic run dir with embedder_blueprint.py + scenario_config.jsonc.

    `body` is the Python source that goes inside build_blueprint() between
    the placeholder marker and `return blueprint`. The canonical script
    is copied into the run dir, then its TODO(AGENT) stub is replaced
    with the test-supplied body.
    """
    rd = tmp_path / "260505_blueprint_run01"
    rd.mkdir()
    (rd / "scenario_config.jsonc").write_text(
        '{\n  "species": ["mouse"],\n  "tissue": ["cortex"],\n'
        '  "celltype": null, "age": "adult",\n  "sex": null, "condition": null\n}\n'
    )
    canonical = (REPO_ROOT / "scripts" / "agentic_embedding" / "embedder_blueprint.py").read_text()
    marker = "    # TODO(AGENT): register your resources, gene modules, and axes here."
    assert marker in canonical
    head, _, tail = canonical.partition(marker)
    return_marker = "    return blueprint"
    tail_after_return = tail.partition(return_marker)[2]
    new_script = head + body + "\n" + return_marker + tail_after_return
    (rd / "embedder_blueprint.py").write_text(new_script)
    return rd


def test_summary_emits_resources_block(tmp_path):
    body = """
    blueprint.add_resource("k1", kind="paper",
                            citation="Author1 Year1 Journal1",
                            doi="10.1000/aaa")
    blueprint.add_resource("k2", kind="database",
                            citation="Author2 Year2 Journal2",
                            doi="10.1000/bbb")
    blueprint.add_gene_module("m1", ["A", "B"], resources=["k1", "k2"])
    blueprint.add_gene_module("m2", ["C"], resources=["k1"])
    blueprint.add_axis("ax1", module="m1")
    blueprint.add_axis("ax2", module="m2")
"""
    rd = _make_run_dir(tmp_path, body)
    em = _load_extract_module()
    summary = em._build_summary(rd)
    bp = summary["blueprint"]

    assert bp["n_resources"] == 2
    assert bp["n_resources_by_kind"] == {"paper": 1, "atlas": 0, "database": 1}
    assert bp["n_dois_unique"] == 2
    keys = {r["key"] for r in bp["resources"]}
    assert keys == {"k1", "k2"}
    k1 = next(r for r in bp["resources"] if r["key"] == "k1")
    assert k1["kind"] == "paper"
    assert k1["doi"] == "10.1000/aaa"
    assert set(k1["used_in_modules"]) == {"m1", "m2"}
    k2 = next(r for r in bp["resources"] if r["key"] == "k2")
    assert k2["used_in_modules"] == ["m1"]


def test_summary_module_records_use_resource_keys(tmp_path):
    body = """
    blueprint.add_resource("k1", kind="paper", citation="C1", doi="10.1000/aaa")
    blueprint.add_gene_module("m1", ["A"], resources=["k1"])
    blueprint.add_axis("ax1", module="m1")
"""
    rd = _make_run_dir(tmp_path, body)
    em = _load_extract_module()
    summary = em._build_summary(rd)
    m1 = next(m for m in summary["blueprint"]["modules"] if m["name"] == "m1")
    assert m1["resource_keys"] == ["k1"]
    assert m1["n_resource_keys"] == 1
    # Old fields must not be present
    assert "references" not in m1
    assert "n_references" not in m1


def test_summary_flags_unused_resources(tmp_path):
    body = """
    blueprint.add_resource("k_used", kind="paper", citation="C", doi="10.1000/aaa")
    blueprint.add_resource("k_unused", kind="paper", citation="C", doi="10.1000/bbb")
    blueprint.add_gene_module("m", ["A"], resources=["k_used"])
    blueprint.add_axis("ax", module="m")
"""
    rd = _make_run_dir(tmp_path, body)
    em = _load_extract_module()
    summary = em._build_summary(rd)
    assert summary["sanity"]["unused_resources"] == ["k_unused"]
    # Replaces modules_missing_refs:
    assert "modules_missing_refs" not in summary["sanity"]
    assert summary["sanity"]["modules_missing_resources"] == []


def _write_coverage(rd: Path, payload: dict) -> Path:
    cov = rd / "emb_coverage.json"
    cov.write_text(json.dumps(payload))
    return cov


def test_summary_prunes_excluded_axes_and_reflects_clip(tmp_path):
    """A coverage sidecar that excludes a MIDDLE axis must prune summary.json's
    axes list to the surviving axes in dim order, fix the counts, and adopt the
    sidecar's actual clip."""
    body = """
    blueprint.add_resource("k1", kind="paper", citation="C1", doi="10.1000/aaa")
    blueprint.add_gene_module("m_a", ["A"], resources=["k1"])
    blueprint.add_gene_module("m_z", ["Z"], resources=["k1"])
    blueprint.add_gene_module("m_c", ["C"], resources=["k1"])
    blueprint.add_axis("ax_a", module="m_a")
    blueprint.add_contrast_axis("ax_z", positive_module="m_z", negative_module="m_a")
    blueprint.add_axis("ax_c", module="m_c")
"""
    rd = _make_run_dir(tmp_path, body)
    cov = _write_coverage(rd, {
        "h5ad_path": "data.h5ad", "layer": None, "n_cells": 3, "n_vars": 3,
        "n_min_genes": 1, "clip": 2.0,
        "modules": [],
        "axes": [
            {"name": "ax_a", "positive_module": "m_a", "n_present": 1, "kept": True, "dim": 0},
            {"name": "ax_z", "positive_module": "m_z", "n_present": 0, "kept": False, "dim": None},
            {"name": "ax_c", "positive_module": "m_c", "n_present": 1, "kept": True, "dim": 1},
        ],
    })
    em = _load_extract_module()
    summary = em._build_summary(rd, coverage_path=cov)
    bp = summary["blueprint"]

    assert [a["name"] for a in bp["axes"]] == ["ax_a", "ax_c"]   # pruned, in dim order
    assert bp["n_axes"] == 2
    assert bp["n_axes_declared"] == 3
    assert bp["excluded_axes"] == ["ax_z"]
    assert bp["n_axes_by_kind"] == {"plain": 2, "contrast": 0, "gated": 0}
    assert bp["n_min_genes"] == 1
    assert bp["clip"] == 2.0
    assert summary["sanity"]["set_clip_overridden"] is True


def test_summary_no_exclusion_leaves_axes_untouched(tmp_path):
    body = """
    blueprint.add_resource("k1", kind="paper", citation="C1", doi="10.1000/aaa")
    blueprint.add_gene_module("m_a", ["A"], resources=["k1"])
    blueprint.add_gene_module("m_c", ["C"], resources=["k1"])
    blueprint.add_axis("ax_a", module="m_a")
    blueprint.add_axis("ax_c", module="m_c")
"""
    rd = _make_run_dir(tmp_path, body)
    cov = _write_coverage(rd, {
        "h5ad_path": "data.h5ad", "layer": None, "n_cells": 3, "n_vars": 3,
        "n_min_genes": 0, "clip": 4.0,
        "modules": [],
        "axes": [
            {"name": "ax_a", "positive_module": "m_a", "n_present": 1, "kept": True, "dim": 0},
            {"name": "ax_c", "positive_module": "m_c", "n_present": 1, "kept": True, "dim": 1},
        ],
    })
    em = _load_extract_module()
    summary = em._build_summary(rd, coverage_path=cov)
    bp = summary["blueprint"]
    assert [a["name"] for a in bp["axes"]] == ["ax_a", "ax_c"]
    assert bp["n_axes"] == 2
    assert "n_axes_declared" not in bp          # no rewrite when nothing excluded
    assert "excluded_axes" not in bp
    assert bp["n_min_genes"] == 0
    assert summary["sanity"]["set_clip_overridden"] is False


def test_summary_n_dois_unique_dedupes_shared_dois(tmp_path):
    """Two resources sharing a DOI count as one DOI (resources stay separate)."""
    body = """
    blueprint.add_resource("k1", kind="paper",
                            citation="Author 2024 Cell",
                            doi="10.1000/shared")
    blueprint.add_resource("k2", kind="atlas",
                            citation="Author 2024 Cell (atlas entry)",
                            doi="10.1000/shared")
    blueprint.add_resource("k3", kind="database",
                            citation="Curated 2025 NAR",
                            doi="10.1000/distinct")
    blueprint.add_gene_module("m", ["A"], resources=["k1", "k2", "k3"])
    blueprint.add_axis("ax", module="m")
"""
    rd = _make_run_dir(tmp_path, body)
    em = _load_extract_module()
    summary = em._build_summary(rd)
    bp = summary["blueprint"]
    assert bp["n_resources"] == 3
    assert bp["n_dois_unique"] == 2
