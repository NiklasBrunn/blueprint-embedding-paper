"""Tests for scripts/blueprint_variability_metrics/metrics.py — scenario-field
normalisation. `species` is array-typed in current scenario configs but older
`summary.json` files on disk still carry it as a bare string; the helpers below
must accept both forms (and keep `scenario_key` hashable, since it is used as a
dict key when grouping summaries by scenario)."""
from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[1]


def _load(name: str, rel_path: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / rel_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def metrics():
    return _load("blueprint_variability_metrics",
                 "scripts/blueprint_variability_metrics/metrics.py")


def test_as_seq_accepts_str_list_and_none(metrics):
    assert metrics._as_seq(None) == []
    assert metrics._as_seq("mouse") == ["mouse"]
    assert metrics._as_seq(["mouse", "human"]) == ["mouse", "human"]
    assert metrics._as_seq([]) == []


@pytest.mark.parametrize("species", ["mouse", ["mouse"], ["mouse", "human"], None])
def test_scenario_key_is_hashable_for_any_species_form(metrics, species):
    sc = {"species": species, "tissue": ["cortex"], "celltype": None,
          "age": "adult", "sex": None, "condition": None}
    key = metrics.scenario_key(sc)
    # Must be usable as a dict key (this is how summaries get grouped).
    hash(key)
    {key: 1}


def test_scenario_key_groups_string_and_one_element_list_species_together(metrics):
    a = metrics.scenario_key({"species": "mouse", "tissue": ["cortex"]})
    b = metrics.scenario_key({"species": ["mouse"], "tissue": ["cortex"]})
    assert a == b


def test_build_run_df_joins_species_for_both_forms(metrics):
    summaries = [
        {"run_id": "r1", "scenario": {"species": "mouse", "tissue": ["cortex"]},
         "blueprint": {}, "sanity": {}, "report": {}},
        {"run_id": "r2", "scenario": {"species": ["mouse", "human"], "tissue": ["cortex"]},
         "blueprint": {}, "sanity": {}, "report": {}},
    ]
    df = metrics.build_run_df(summaries)
    assert list(df["species"]) == ["mouse", "mouse|human"]
