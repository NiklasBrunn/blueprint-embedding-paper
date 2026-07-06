"""Tests for scripts/agentic_embedding/embedder_blueprint.py."""
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
import anndata as ad

import argparse

from scripts.agentic_embedding.embedder_blueprint import (
    AxisSpec,
    EmbedderBlueprint,
    ModuleSpec,
    Resource,
    _check_placeholders,
    _clip_value,
    build_blueprint,
)


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "agentic_embedding" / "embedder_blueprint.py"

_TEST_DOI = "10.1000/test"
_VALID_DOI = "10.1038/s41586-018-0654-5"
_VALID_DOI_2 = "10.1016/j.cell.2021.04.021"


def _add_default_resource(bp, key="ref"):
    bp.add_resource(key, kind="paper", citation=key, doi=_TEST_DOI)
    return key


def _blueprint_with_modules(*names: str) -> EmbedderBlueprint:
    blueprint = EmbedderBlueprint("test")
    _add_default_resource(blueprint)
    for n in names:
        blueprint.add_gene_module(n, ["DUMMY"], resources=["ref"])
    return blueprint


def _make_h5ad(tmp_path, X=None, var_names=None, obs_names=None):
    """Build a tiny synthetic AnnData and write it to tmp_path/data.h5ad."""
    if X is None:
        X = np.array([
            [1.0, 0.0, 2.0],
            [3.0, 0.0, 1.0],
            [2.0, 0.0, 5.0],
        ])
    if var_names is None:
        var_names = ["A", "B", "C"]
    if obs_names is None:
        obs_names = [f"c{i+1}" for i in range(X.shape[0])]
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(index=obs_names),
        var=pd.DataFrame(index=var_names),
    )
    path = tmp_path / "data.h5ad"
    adata.write_h5ad(path)
    return path


# ---------------------------------------------------------------------------
# CLI smoke
# ---------------------------------------------------------------------------

def test_help_exits_cleanly_and_lists_flags():
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--help"],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0
    assert "--input" in result.stdout
    assert "--output" in result.stdout
    assert "--describe" in result.stdout
    assert "--layer" in result.stdout
    assert "--clip" in result.stdout
    assert "--n-min-genes" in result.stdout


def test_cli_describe_runs_without_input_file(tmp_path):
    # Even with no --input path existing, --describe should succeed.
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--describe",
         "--input", str(tmp_path / "nonexistent.h5ad")],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0


# ---------------------------------------------------------------------------
# Dataclass API
# ---------------------------------------------------------------------------

def test_module_spec_is_frozen_and_stores_genes_as_tuple():
    module = ModuleSpec(name="inhibitory", genes=("Gad1", "Gad2"))
    assert module.name == "inhibitory"
    assert module.genes == ("Gad1", "Gad2")
    with pytest.raises(Exception):
        module.name = "other"  # frozen


def test_module_spec_stores_resource_keys_as_tuple():
    module = ModuleSpec(name="m", genes=("A",), resource_keys=("tasic2018",))
    assert module.resource_keys == ("tasic2018",)


def test_module_spec_default_resource_keys_empty():
    module = ModuleSpec(name="m", genes=("A",))
    assert module.resource_keys == ()


def test_axis_spec_defaults_for_plain_axis():
    axis = AxisSpec(name="vip", positive="vip_cge")
    assert axis.negatives == ()
    assert axis.negative_weight == 1.0
    assert axis.gate is None


def test_resource_dataclass_is_frozen_with_expected_fields():
    r = Resource(
        key="tasic2018",
        kind="paper",
        citation="Tasic 2018 Nature",
        doi="10.1038/s41586-018-0654-5",
    )
    assert r.key == "tasic2018"
    assert r.kind == "paper"
    assert r.citation == "Tasic 2018 Nature"
    assert r.doi == "10.1038/s41586-018-0654-5"
    with pytest.raises(Exception):
        r.kind = "atlas"  # frozen


# ---------------------------------------------------------------------------
# EmbedderBlueprint construction + add_resource / add_gene_module
# ---------------------------------------------------------------------------

def test_blueprint_starts_empty():
    blueprint = EmbedderBlueprint(name="test")
    assert blueprint.name == "test"
    assert blueprint.modules == {}
    assert blueprint.axes == ()
    assert blueprint.resources == {}


def test_build_blueprint_returns_empty_blueprint():
    bp = build_blueprint("default")
    assert isinstance(bp, EmbedderBlueprint)
    assert bp.axes == ()


def test_add_resource_registers_a_resource():
    bp = EmbedderBlueprint(name="t")
    bp.add_resource(
        "tasic2018", kind="paper",
        citation="Tasic 2018 Nature", doi=_VALID_DOI,
    )
    assert "tasic2018" in bp.resources
    r = bp.resources["tasic2018"]
    assert r.kind == "paper"
    assert r.citation == "Tasic 2018 Nature"
    assert r.doi == _VALID_DOI


def test_add_resource_rejects_duplicate_key():
    bp = EmbedderBlueprint(name="t")
    bp.add_resource("k", kind="paper", citation="C", doi=_VALID_DOI)
    with pytest.raises(ValueError, match="already registered"):
        bp.add_resource("k", kind="paper", citation="C", doi=_VALID_DOI)


@pytest.mark.parametrize("bad_kind", ["", "Paper", "book", "preprint", None])
def test_add_resource_rejects_invalid_kind(bad_kind):
    bp = EmbedderBlueprint(name="t")
    with pytest.raises(ValueError, match="invalid kind"):
        bp.add_resource("k", kind=bad_kind, citation="C", doi=_VALID_DOI)


@pytest.mark.parametrize("field", ["key", "citation"])
def test_add_resource_rejects_empty_required_string(field):
    bp = EmbedderBlueprint(name="t")
    args = {"key": "k", "kind": "paper", "citation": "C", "doi": _VALID_DOI}
    args[field] = ""
    with pytest.raises(ValueError, match="must be non-empty"):
        bp.add_resource(**args)


@pytest.mark.parametrize("kind", ["paper", "atlas"])
def test_add_resource_requires_doi_for_paper_and_atlas(kind):
    bp = EmbedderBlueprint(name="t")
    with pytest.raises(ValueError, match="doi is required"):
        bp.add_resource("k", kind=kind, citation="C", doi=None)
    bp2 = EmbedderBlueprint(name="t")
    with pytest.raises(ValueError, match="doi is required"):
        bp2.add_resource("k", kind=kind, citation="C", doi="")


def test_add_resource_allows_none_doi_for_database():
    bp = EmbedderBlueprint(name="t")
    bp.add_resource(
        "msigdb_hallmark_glycolysis",
        kind="database",
        citation="MSigDB HALLMARK_GLYCOLYSIS",
        doi=None,
    )
    assert bp.resources["msigdb_hallmark_glycolysis"].doi is None


def test_add_resource_database_with_doi_still_format_validated():
    # If a database resource provides a DOI, it still has to be format-valid.
    bp = EmbedderBlueprint(name="t")
    bp.add_resource(
        "zenodo_deposit",
        kind="database",
        citation="DB GO:0007049",
        doi="10.5281/zenodo.1234567",
    )
    bp2 = EmbedderBlueprint(name="t")
    with pytest.raises(ValueError, match="invalid DOI"):
        bp2.add_resource(
            "bad_db", kind="database", citation="DB GO:0007049",
            doi="not a doi",
        )


@pytest.mark.parametrize("bad_doi", [
    "doi:10.1038/s41586-018-0654-5",
    "https://doi.org/10.1038/x",
    "10.x/abc",
    "10/abc",
    "not a doi",
    "10.1038",
])
def test_add_resource_rejects_malformed_doi(bad_doi):
    bp = EmbedderBlueprint(name="t")
    with pytest.raises(ValueError, match="invalid DOI"):
        bp.add_resource("k", kind="paper", citation="C", doi=bad_doi)


@pytest.mark.parametrize("good_doi", [
    "10.1038/s41586-018-0654-5",
    "10.1016/j.cell.2021.04.021",
    "10.5281/zenodo.1234567",
    "10.1093/nar/gkx1132",
])
def test_add_resource_accepts_valid_dois(good_doi):
    bp = EmbedderBlueprint(name="t")
    bp.add_resource("k", kind="paper", citation="C", doi=good_doi)
    assert bp.resources["k"].doi == good_doi


def test_add_gene_module_stores_genes_as_tuple():
    blueprint = EmbedderBlueprint(name="test")
    _add_default_resource(blueprint)
    blueprint.add_gene_module("inhibitory", ["Gad1", "Gad2"], resources=["ref"])
    assert blueprint.modules["inhibitory"].genes == ("Gad1", "Gad2")


def test_add_gene_module_rejects_duplicates():
    blueprint = EmbedderBlueprint(name="test")
    _add_default_resource(blueprint)
    blueprint.add_gene_module("inhibitory", ["Gad1"], resources=["ref"])
    with pytest.raises(ValueError, match="already registered"):
        blueprint.add_gene_module("inhibitory", ["Gad2"], resources=["ref"])


def test_add_gene_module_requires_non_empty_resources():
    bp = EmbedderBlueprint(name="t")
    bp.add_resource("k", kind="paper", citation="C", doi=_VALID_DOI)
    with pytest.raises(ValueError, match="resources must be non-empty"):
        bp.add_gene_module("m", ["A"], resources=[])


def test_add_gene_module_stores_resource_keys():
    bp = EmbedderBlueprint(name="t")
    bp.add_resource("k1", kind="paper", citation="C1", doi=_VALID_DOI)
    bp.add_resource("k2", kind="atlas", citation="C2", doi=_VALID_DOI_2)
    bp.add_gene_module("m", ["A"], resources=["k1", "k2"])
    assert bp.modules["m"].resource_keys == ("k1", "k2")


def test_add_gene_module_rejects_unknown_resource_key():
    bp = EmbedderBlueprint(name="t")
    bp.add_resource("k1", kind="paper", citation="C", doi=_VALID_DOI)
    with pytest.raises(ValueError, match="unknown resource"):
        bp.add_gene_module("m", ["A"], resources=["k1", "missing_key"])


def test_add_gene_module_rejects_when_no_resources_registered():
    bp = EmbedderBlueprint(name="t")
    with pytest.raises(ValueError, match="unknown resource"):
        bp.add_gene_module("m", ["A"], resources=["k"])


# ---------------------------------------------------------------------------
# add_axis / add_contrast_axis / add_gated_axis
# ---------------------------------------------------------------------------

def test_add_axis_plain():
    b = _blueprint_with_modules("vip_cge")
    b.add_axis("vip", module="vip_cge")
    (axis,) = b.axes
    assert axis.name == "vip"
    assert axis.positive == "vip_cge"
    assert axis.negatives == ()
    assert axis.gate is None


def test_add_axis_rejects_unknown_module():
    b = EmbedderBlueprint("test")
    with pytest.raises(ValueError, match="unknown module"):
        b.add_axis("vip", module="vip_cge")


def test_add_contrast_axis():
    b = _blueprint_with_modules("pan_neuron", "non_neuron")
    b.add_contrast_axis("neuron", positive_module="pan_neuron", negative_module="non_neuron")
    (axis,) = b.axes
    assert axis.positive == "pan_neuron"
    assert axis.negatives == ("non_neuron",)
    assert axis.negative_weight == 1.0
    assert axis.gate is None


def test_add_gated_axis_defaults():
    b = _blueprint_with_modules("inhibitory", "vip_cge", "sst", "pvalb")
    b.add_gated_axis(
        "vip",
        gate_module="inhibitory",
        positive_module="vip_cge",
        negative_modules=["sst", "pvalb"],
    )
    (axis,) = b.axes
    assert axis.gate == "inhibitory"
    assert axis.positive == "vip_cge"
    assert axis.negatives == ("sst", "pvalb")
    assert axis.negative_weight == 0.5


def test_axis_name_must_be_unique():
    b = _blueprint_with_modules("m")
    b.add_axis("dup", module="m")
    with pytest.raises(ValueError, match="already registered"):
        b.add_axis("dup", module="m")


def test_set_clip_can_be_disabled():
    b = EmbedderBlueprint("t")
    b.set_clip(None)
    assert b._clip is None
    b.set_clip(3.0)
    assert b._clip == 3.0


# ---------------------------------------------------------------------------
# --clip CLI value parser (_clip_value)
# ---------------------------------------------------------------------------

def test_clip_value_parses_positive_float():
    assert _clip_value("2.5") == 2.5


@pytest.mark.parametrize("s", ["none", "NONE", "off", "disable", " none "])
def test_clip_value_returns_none_sentinel(s):
    # The sentinel is the string "none", distinct from the flag default None.
    assert _clip_value(s) == "none"


@pytest.mark.parametrize("s", ["0", "-1", "abc", ""])
def test_clip_value_rejects_invalid(s):
    with pytest.raises(argparse.ArgumentTypeError):
        _clip_value(s)


def test_normalize_disable_clip_keeps_outlier():
    b = EmbedderBlueprint("t")
    b.set_clip(None)
    raw = np.array([[0.0], [0.0], [0.0], [0.0], [100.0]])
    out = b._normalize(raw)
    # With clipping off, the z-scored outlier stays well above any clip bound.
    assert out.max() > 1.5


# ---------------------------------------------------------------------------
# _normalize
# ---------------------------------------------------------------------------

def test_normalize_zscore_with_clip():
    b = EmbedderBlueprint("t")
    b.set_clip(1.0)
    raw = np.array([[0.0], [1.0], [10.0]])
    out = b._normalize(raw)
    # outlier clipped to +1; bottom value goes below zero
    assert out[-1, 0] == pytest.approx(1.0)
    assert out[0, 0] < 0


# ---------------------------------------------------------------------------
# End-to-end build_embedding using h5ad fixture
# ---------------------------------------------------------------------------

def test_build_embedding_output_dim_matches_axis_count(tmp_path):
    """One dim column per registered axis, in declaration order."""
    h5ad_path = _make_h5ad(tmp_path)  # default 3-cell fixture
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("m_a", ["A"], resources=["ref"])
    bp.add_gene_module("m_c", ["C"], resources=["ref"])
    bp.add_axis("dim_a", module="m_a")
    bp.add_axis("dim_c", module="m_c")
    bp.add_contrast_axis("dim_a_vs_c", positive_module="m_a", negative_module="m_c")
    out = tmp_path / "embedding.parquet"
    bp.build_embedding(h5ad_path, out)

    df = pd.read_parquet(out)
    assert list(df.columns) == ["cell_id", "dim_0", "dim_1", "dim_2"]
    for col in ("dim_0", "dim_1", "dim_2"):
        assert df[col].std() > 0, f"{col} is constant; expected nonzero variance"


def test_build_embedding_keeps_constant_axes_as_zero_columns(tmp_path):
    h5ad_path = _make_h5ad(tmp_path)
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("m_a", ["A"], resources=["ref"])
    bp.add_gene_module("m_missing", ["Z"], resources=["ref"])  # Z not in var_names
    bp.add_axis("dim_a", module="m_a")
    bp.add_axis("dim_missing", module="m_missing")
    out = tmp_path / "embedding.parquet"
    bp.build_embedding(h5ad_path, out)

    df = pd.read_parquet(out)
    assert list(df.columns) == ["cell_id", "dim_0", "dim_1"]
    assert (df["dim_1"] == 0.0).all()
    assert df["dim_0"].std() > 0


def test_build_embedding_matches_genes_case_insensitively(tmp_path):
    """The engine matches gene names case-insensitively. Same data, same
    expected output regardless of whether the blueprint says 'A' (matches
    var_names exactly), 'a' (lowercase), or 'a' / 'A' mixed across modules.
    """
    # Dataset uses uppercase var_names (the default fixture).
    h5ad_path = _make_h5ad(tmp_path)

    out_exact = tmp_path / "exact.parquet"
    bp1 = EmbedderBlueprint("exact")
    _add_default_resource(bp1)
    bp1.add_gene_module("m_a", ["A"], resources=["ref"])
    bp1.add_gene_module("m_c", ["C"], resources=["ref"])
    bp1.add_axis("dim_a", module="m_a")
    bp1.add_axis("dim_c", module="m_c")
    bp1.build_embedding(h5ad_path, out_exact)
    df_exact = pd.read_parquet(out_exact)

    out_lower = tmp_path / "lower.parquet"
    bp2 = EmbedderBlueprint("lower")
    _add_default_resource(bp2)
    bp2.add_gene_module("m_a", ["a"], resources=["ref"])     # lowercase
    bp2.add_gene_module("m_c", ["c"], resources=["ref"])     # lowercase
    bp2.add_axis("dim_a", module="m_a")
    bp2.add_axis("dim_c", module="m_c")
    bp2.build_embedding(h5ad_path, out_lower)
    df_lower = pd.read_parquet(out_lower)

    # Identical columns, identical values: case-insensitive matching produced
    # the same output as the exact-case blueprint.
    assert list(df_lower.columns) == list(df_exact.columns)
    for col in ("dim_0", "dim_1"):
        assert df_lower[col].std() > 0, f"{col} should have nonzero variance under case-blind match"
        assert (df_lower[col].to_numpy() == df_exact[col].to_numpy()).all(), \
            f"{col} differs between exact-case and lowercase blueprint — case-insensitive match failed"


def test_build_embedding_case_insensitive_coverage_sidecar(tmp_path):
    """When a gene is registered in different case from var_names, it must
    NOT show up as missing in the coverage sidecar."""
    h5ad_path = _make_h5ad(tmp_path)  # var_names = ["A", "B", "C"]
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("m_mixed", ["a", "C"], resources=["ref"])  # a vs A
    bp.add_axis("dim_m", module="m_mixed")
    out = tmp_path / "embedding.parquet"
    bp.build_embedding(h5ad_path, out)

    import json
    sidecar = tmp_path / "embedding_coverage.json"
    record = json.loads(sidecar.read_text())
    mod = next(m for m in record["modules"] if m["name"] == "m_mixed")
    assert mod["n_genes"] == 2
    assert mod["n_genes_present"] == 2, "Lowercase 'a' should match 'A' in var_names case-blind"
    assert mod["n_genes_missing"] == 0
    assert mod["missing_genes"] == []


def test_build_embedding_rejects_missing_file(tmp_path):
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("m", ["A"], resources=["ref"])
    bp.add_axis("dim_m", module="m")
    with pytest.raises(FileNotFoundError):
        bp.build_embedding(tmp_path / "missing.h5ad", tmp_path / "out.parquet")


def test_build_embedding_rejects_empty_axes(tmp_path):
    h5ad_path = _make_h5ad(tmp_path)
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("m", ["A"], resources=["ref"])  # no axis registered
    with pytest.raises(ValueError, match="no axes registered"):
        bp.build_embedding(h5ad_path, tmp_path / "out.parquet")


def test_build_embedding_writes_all_zero_parquet_when_all_axes_constant(tmp_path):
    """No more RuntimeError; constant axes are zero columns."""
    h5ad_path = _make_h5ad(tmp_path)
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("m_z", ["Z"], resources=["ref"])  # Z absent
    bp.add_axis("dim_z", module="m_z")
    out = tmp_path / "embedding.parquet"
    bp.build_embedding(h5ad_path, out)

    df = pd.read_parquet(out)
    assert list(df.columns) == ["cell_id", "dim_0"]
    assert (df["dim_0"] == 0.0).all()


def test_build_embedding_prints_coverage_line_with_did_you_mean(tmp_path, capsys):
    """Coverage telemetry reports missing genes and suggests close matches from var_names."""
    X = np.array([
        [1.0, 2.0, 3.0],
        [2.0, 1.0, 4.0],
    ])
    h5ad_path = _make_h5ad(tmp_path, X=X, var_names=["A", "Gad1", "B"])
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("m_a", ["A"], resources=["ref"])             # all present
    bp.add_gene_module("m_partial", ["Gad2", "ZZZ"], resources=["ref"])  # 0/2 present, one close to Gad1
    bp.add_axis("dim_a", module="m_a")
    bp.add_axis("dim_partial", module="m_partial")
    out = tmp_path / "embedding.parquet"
    bp.build_embedding(h5ad_path, out)

    captured = capsys.readouterr().out
    assert "module 'm_a'" not in captured
    assert "module 'm_partial'" in captured
    assert "0/2 genes present" in captured
    assert "Gad2 (did you mean: 'Gad1')" in captured
    assert "ZZZ" in captured


# ---------------------------------------------------------------------------
# describe()
# ---------------------------------------------------------------------------

def test_describe_lists_axes_with_kind_and_composition():
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp, "r1")
    _add_default_resource(bp, "r2")
    bp.add_gene_module("pos", ["A"], resources=["r1"])
    bp.add_gene_module("neg", ["B"], resources=["r1", "r2"])
    bp.add_axis("dim_pos", module="pos")
    bp.add_contrast_axis("dim_pn", positive_module="pos", negative_module="neg")
    out = bp.describe()
    assert "Blueprint: t" in out
    assert "modules: 2 registered" in out
    assert "axes: 2 registered" in out
    assert "[plain]" in out
    assert "[contrast]" in out
    assert "resource_keys" in out


def test_describe_lists_module_resource_key_counts():
    b = EmbedderBlueprint("t")
    _add_default_resource(b, "r1")
    _add_default_resource(b, "r2")
    _add_default_resource(b, "r3")
    b.add_gene_module("m1", ["A"], resources=["r1", "r2"])
    b.add_gene_module("m2", ["B"], resources=["r3"])
    b.add_axis("plain_m1", module="m1")
    b.add_axis("plain_m2", module="m2")

    text = b.describe()
    assert "2 resource_keys" in text
    assert "1 resource_keys" in text


# ---------------------------------------------------------------------------
# Placeholder fail-fast gate (called on every `--describe` invocation)
# ---------------------------------------------------------------------------

def test_check_placeholders_passes_on_clean_blueprint():
    bp = _blueprint_with_modules("real_module")
    bp.add_axis("real_axis", module="real_module")
    _check_placeholders(bp)  # must not raise


@pytest.mark.parametrize("marker", ["<<EXAMPLE_MOD>>", "MOD_DELETE_ME", "PLACEHOLDER_X"])
def test_check_placeholders_raises_on_module_name(marker):
    bp = _blueprint_with_modules(marker)
    with pytest.raises(RuntimeError, match="Placeholder text still present"):
        _check_placeholders(bp)


def test_check_placeholders_raises_on_gene_name():
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("real_module", ["<<EXAMPLE_GENE_A>>", "REAL"], resources=["ref"])
    with pytest.raises(RuntimeError, match="gene in module"):
        _check_placeholders(bp)


def test_check_placeholders_raises_on_axis_name():
    bp = _blueprint_with_modules("real_module")
    bp.add_axis("<<EXAMPLE_AXIS>>", module="real_module")
    with pytest.raises(RuntimeError, match="axis name"):
        _check_placeholders(bp)


def test_check_placeholders_raises_on_resource_citation():
    bp = EmbedderBlueprint("t")
    bp.add_resource("ref", kind="paper",
                    citation="<<EXAMPLE_CITATION>>", doi=_TEST_DOI)
    bp.add_gene_module("real_module", ["GENE"], resources=["ref"])
    with pytest.raises(RuntimeError, match="citation"):
        _check_placeholders(bp)


def test_check_placeholders_lists_every_offender():
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("<<EXAMPLE_MOD>>", ["<<EXAMPLE_GENE>>"], resources=["ref"])
    bp.add_axis("<<EXAMPLE_AXIS>>", module="<<EXAMPLE_MOD>>")
    with pytest.raises(RuntimeError) as exc:
        _check_placeholders(bp)
    msg = str(exc.value)
    assert "module name" in msg
    assert "gene in module" in msg
    assert "axis name" in msg


# ---------------------------------------------------------------------------
# n_min_genes axis exclusion
# ---------------------------------------------------------------------------

def test_n_min_genes_zero_keeps_all_axes(tmp_path):
    h5ad_path = _make_h5ad(tmp_path)  # var_names A, B, C
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("m_a", ["A"], resources=["ref"])   # present
    bp.add_gene_module("m_z", ["Z"], resources=["ref"])   # absent
    bp.add_axis("ax_a", module="m_a")
    bp.add_axis("ax_z", module="m_z")
    out = tmp_path / "e.parquet"
    bp.build_embedding(h5ad_path, out, n_min_genes=0)
    df = pd.read_parquet(out)
    assert list(df.columns) == ["cell_id", "dim_0", "dim_1"]


def test_n_min_genes_excludes_absent_positive_module(tmp_path):
    h5ad_path = _make_h5ad(tmp_path)
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("m_a", ["A"], resources=["ref"])   # 1 present
    bp.add_gene_module("m_z", ["Z"], resources=["ref"])   # 0 present
    bp.add_axis("ax_a", module="m_a")
    bp.add_axis("ax_z", module="m_z")
    out = tmp_path / "e.parquet"
    bp.build_embedding(h5ad_path, out, n_min_genes=1)
    df = pd.read_parquet(out)
    assert list(df.columns) == ["cell_id", "dim_0"]   # only the kept axis
    assert df["dim_0"].std() > 0


def test_n_min_genes_compact_renumbers_after_middle_drop(tmp_path):
    """Dropping a MIDDLE axis must compactly renumber surviving dims, and the
    surviving dim must carry the right axis's score (no positional slip)."""
    h5ad_path = _make_h5ad(tmp_path)  # A, B, C present
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("m_a", ["A"], resources=["ref"])   # present
    bp.add_gene_module("m_z", ["Z"], resources=["ref"])   # absent (middle)
    bp.add_gene_module("m_c", ["C"], resources=["ref"])   # present
    bp.add_axis("ax_a", module="m_a")
    bp.add_axis("ax_z", module="m_z")     # dropped
    bp.add_axis("ax_c", module="m_c")
    out = tmp_path / "e.parquet"
    bp.build_embedding(h5ad_path, out, n_min_genes=1)
    df = pd.read_parquet(out)
    assert list(df.columns) == ["cell_id", "dim_0", "dim_1"]

    # dim_1 must equal the C-axis score from an A+C-only blueprint (not B/Z).
    ref = tmp_path / "ref.parquet"
    bp2 = EmbedderBlueprint("t2")
    _add_default_resource(bp2)
    bp2.add_gene_module("m_a", ["A"], resources=["ref"])
    bp2.add_gene_module("m_c", ["C"], resources=["ref"])
    bp2.add_axis("ax_a", module="m_a")
    bp2.add_axis("ax_c", module="m_c")
    bp2.build_embedding(h5ad_path, ref)
    df_ref = pd.read_parquet(ref)
    assert (df["dim_1"].to_numpy() == df_ref["dim_1"].to_numpy()).all()


def test_n_min_genes_records_axes_in_sidecar(tmp_path):
    h5ad_path = _make_h5ad(tmp_path)
    bp = EmbedderBlueprint("t")
    bp.set_clip(2.0)
    _add_default_resource(bp)
    bp.add_gene_module("m_a", ["A"], resources=["ref"])
    bp.add_gene_module("m_z", ["Z"], resources=["ref"])
    bp.add_gene_module("m_c", ["C"], resources=["ref"])
    bp.add_axis("ax_a", module="m_a")
    bp.add_axis("ax_z", module="m_z")
    bp.add_axis("ax_c", module="m_c")
    out = tmp_path / "e.parquet"
    bp.build_embedding(h5ad_path, out, n_min_genes=1)

    import json
    rec = json.loads((tmp_path / "e_coverage.json").read_text())
    assert rec["n_min_genes"] == 1
    assert rec["clip"] == 2.0
    axmap = {a["name"]: a for a in rec["axes"]}
    assert axmap["ax_a"]["kept"] is True and axmap["ax_a"]["dim"] == 0
    assert axmap["ax_z"]["kept"] is False and axmap["ax_z"]["dim"] is None
    assert axmap["ax_c"]["kept"] is True and axmap["ax_c"]["dim"] == 1


def test_n_min_genes_only_considers_positive_module(tmp_path):
    """Gate/negative module presence is irrelevant; only the positive counts."""
    h5ad_path = _make_h5ad(tmp_path)  # A, B, C present
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("gate", ["A"], resources=["ref"])         # present
    bp.add_gene_module("pos_absent", ["Z"], resources=["ref"])   # absent
    bp.add_gene_module("pos_present", ["C"], resources=["ref"])  # present
    bp.add_gated_axis("g_absent", gate_module="gate", positive_module="pos_absent")
    bp.add_gated_axis("g_present", gate_module="gate", positive_module="pos_present")
    out = tmp_path / "e.parquet"
    bp.build_embedding(h5ad_path, out, n_min_genes=1)
    df = pd.read_parquet(out)
    assert list(df.columns) == ["cell_id", "dim_0"]

    import json
    rec = json.loads((tmp_path / "e_coverage.json").read_text())
    kept = {a["name"]: a["kept"] for a in rec["axes"]}
    assert kept == {"g_absent": False, "g_present": True}


def test_n_min_genes_all_excluded_raises(tmp_path):
    h5ad_path = _make_h5ad(tmp_path)
    bp = EmbedderBlueprint("t")
    _add_default_resource(bp)
    bp.add_gene_module("m_a", ["A"], resources=["ref"])  # only 1 present
    bp.add_axis("ax_a", module="m_a")
    with pytest.raises(ValueError, match="all 1 axes excluded"):
        bp.build_embedding(h5ad_path, tmp_path / "o.parquet", n_min_genes=2)


# ---------------------------------------------------------------------------
# CLI integration for --clip / --n-min-genes (populated skeleton)
# ---------------------------------------------------------------------------

def _populated_skeleton(tmp_path):
    """Write a stub that monkeypatches build_blueprint to a populated one and
    calls main(), so the full CLI path (flags → set_clip/n_min_genes) runs."""
    skeleton = tmp_path / "stub_blueprint.py"
    skeleton.write_text(
        "import sys; sys.path.insert(0, "
        + repr(str(SCRIPT.parent.parent.parent)) + ")\n"
        "from scripts.agentic_embedding import embedder_blueprint as eb\n"
        "def _patched(name):\n"
        "    bp = eb.EmbedderBlueprint(name)\n"
        f"    bp.add_resource('ref', kind='paper', citation='ok', doi={_TEST_DOI!r})\n"
        "    bp.add_gene_module('m_a', ['A'], resources=['ref'])\n"
        "    bp.add_gene_module('m_z', ['Z'], resources=['ref'])\n"
        "    bp.add_gene_module('m_c', ['C'], resources=['ref'])\n"
        "    bp.add_axis('ax_a', module='m_a')\n"
        "    bp.add_axis('ax_z', module='m_z')\n"
        "    bp.add_axis('ax_c', module='m_c')\n"
        "    return bp\n"
        "eb.build_blueprint = _patched\n"
        "eb.main()\n"
    )
    return skeleton


def test_cli_n_min_genes_excludes_axis(tmp_path):
    h5ad_path = _make_h5ad(tmp_path)
    skeleton = _populated_skeleton(tmp_path)
    out = tmp_path / "o.parquet"
    r = subprocess.run(
        [sys.executable, str(skeleton), "--input", str(h5ad_path),
         "--output", str(out), "--n-min-genes", "1"],
        capture_output=True, text=True, check=False,
    )
    assert r.returncode == 0, r.stderr
    df = pd.read_parquet(out)
    assert list(df.columns) == ["cell_id", "dim_0", "dim_1"]
    assert "[n_min_genes] excluded 1/3 axes" in r.stdout


def test_cli_clip_override_bounds_output(tmp_path):
    h5ad_path = _make_h5ad(tmp_path)
    skeleton = _populated_skeleton(tmp_path)
    out_default = tmp_path / "def.parquet"
    out_clip = tmp_path / "clip.parquet"
    subprocess.run(
        [sys.executable, str(skeleton), "--input", str(h5ad_path),
         "--output", str(out_default)],
        check=True, capture_output=True, text=True,
    )
    subprocess.run(
        [sys.executable, str(skeleton), "--input", str(h5ad_path),
         "--output", str(out_clip), "--clip", "0.5"],
        check=True, capture_output=True, text=True,
    )
    a = pd.read_parquet(out_default)
    b = pd.read_parquet(out_clip)
    dim_cols = [c for c in b.columns if c.startswith("dim_")]
    # Default (clip 4.0) leaves |z| > 0.5 on this fixture; --clip 0.5 clamps it.
    assert a[dim_cols].abs().to_numpy().max() > 0.5
    assert b[dim_cols].abs().to_numpy().max() == pytest.approx(0.5)


def test_cli_clip_none_runs(tmp_path):
    h5ad_path = _make_h5ad(tmp_path)
    skeleton = _populated_skeleton(tmp_path)
    out = tmp_path / "o.parquet"
    r = subprocess.run(
        [sys.executable, str(skeleton), "--input", str(h5ad_path),
         "--output", str(out), "--clip", "none"],
        capture_output=True, text=True, check=False,
    )
    assert r.returncode == 0, r.stderr
    assert out.exists()


def test_cli_describe_fails_when_skeleton_has_placeholders(tmp_path, monkeypatch):
    # Drop a tiny skeleton with an un-replaced sentinel and run --describe; the
    # gate fires before the describe branch, so the CLI must exit non-zero.
    skeleton = tmp_path / "stub_blueprint.py"
    skeleton.write_text(
        "import sys; sys.path.insert(0, "
        + repr(str(SCRIPT.parent.parent.parent)) + ")\n"
        "from scripts.agentic_embedding import embedder_blueprint as eb\n"
        "def _patched(name):\n"
        "    bp = eb.EmbedderBlueprint(name)\n"
        f"    bp.add_resource('ref', kind='paper', citation='ok', doi={_TEST_DOI!r})\n"
        "    bp.add_gene_module('<<EXAMPLE_MOD>>', ['G'], resources=['ref'])\n"
        "    return bp\n"
        "eb.build_blueprint = _patched\n"
        "eb.main()\n"
    )
    result = subprocess.run(
        [sys.executable, str(skeleton), "--describe",
         "--input", str(tmp_path / "nonexistent.h5ad")],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode != 0
    assert "Placeholder text still present" in result.stderr
