import json
import os
import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]


def _has_archived_bp():
    return next(
        (REPO / "outs" / "blueprints").rglob("embedder_blueprint_seed*.py"), None
    ) is not None


pytestmark = pytest.mark.skipif(
    not _has_archived_bp(), reason="no archived blueprint fixture"
)


@pytest.fixture(scope="module")
def mapped(tmp_path_factory):
    """Run the real token-free app flow once and share the mapped artifacts.

    Mirrors the real app flow exactly: prepare_run_workdir creates run01/ with
    the bare template, the driver writes the POPULATED blueprint to its own
    cohort-named run dir under --work-root, and the app must resolve that real
    dir from the captured log (regression: the app used to read run01/ and hit
    "no axes registered"). Module-scoped so the expensive embed + UMAP map runs
    only once for the whole file.
    """
    import webapp.pipeline as P

    tmp_path = tmp_path_factory.mktemp("e2e")
    orig_runtime = P.RUNTIME
    P.RUNTIME = tmp_path / "runtime"
    try:
        session = "sess1"
        app_workdir, scen = P.prepare_run_workdir(
            session=session, run=1, species=["mouse"], tissue=["cortex"],
            temporal=False, noref=False)
        env = dict(os.environ)
        env["WEBAPP_MOCK_CODEX"] = "1"
        cmd = P.build_driver_command(
            scenario_file=scen, input_h5ad=scen, work_root=P.RUNTIME / session,
            model="gpt-5.5", reasoning="high", with_refine=False, max_iterations=3,
            noref=False, no_embed=True, dataset="webapp", n_parallel=1)
        log_path = app_workdir / "logs" / "codex_edit.log"
        with open(log_path, "w") as fh:  # same capture the JobRegistry does
            r = subprocess.run(cmd, cwd=str(REPO), env=env, stdout=fh,
                               stderr=subprocess.STDOUT, text=True, timeout=600)
        assert r.returncode == 0, log_path.read_text()

        # The app resolves the driver's real output dir, NOT the bare template.
        run_dir = P.resolve_run_dir(app_workdir)
        assert run_dir != app_workdir, "resolved to the bare-template dir"
        bp = run_dir / "embedder_blueprint.py"
        assert bp.exists()
        summary = json.loads((run_dir / "summary.json").read_text())
        assert summary["blueprint"]["n_axes"] >= 1

        # No re-upload: that in-session blueprint maps directly.
        from webapp.mapping import map_h5ad
        from webapp.tests.fixtures.make_tiny_h5ad import write_tiny

        h5 = write_tiny(tmp_path / "tiny.h5ad")
        res = map_h5ad(blueprint=bp, h5ad=h5, outdir=tmp_path / "map")

        import anndata as ad

        adata = ad.read_h5ad(h5)
        return {
            "summary": summary,
            "map_result": res,
            "h5ad": h5,
            "adata": adata,
        }
    finally:
        P.RUNTIME = orig_runtime


def test_generation_e2e_with_mock(mapped):
    import pandas as pd

    df = pd.read_parquet(mapped["map_result"].parquet)
    assert any(c.startswith("dim_") for c in df.columns)
    assert mapped["map_result"].umap_coords.exists()


def test_explore_bundle_shape_and_json(mapped):
    from webapp import explore_bundle

    bundle = explore_bundle.build(
        mapped["map_result"], mapped["adata"], mapped["summary"])

    expected = {"coords", "axisOrder", "values", "branches", "axisMeta",
                "metadata", "modules", "geneList", "shown", "total"}
    assert expected <= set(bundle)

    # Tiny dataset: nothing downsampled, all coords present.
    assert len(bundle["coords"]["x"]) == bundle["shown"]
    assert len(bundle["coords"]["y"]) == bundle["shown"]
    assert bundle["shown"] == bundle["total"]
    assert bundle["total"] == mapped["adata"].n_obs

    # One value array per axis, each as long as the shown cells.
    assert len(bundle["values"]) == len(bundle["axisOrder"])
    for arr in bundle["values"].values():
        assert len(arr) == bundle["shown"]

    # Must survive a dcc.Store round-trip (JSON-serializable).
    assert json.loads(json.dumps(bundle)) == bundle


def test_explore_gene_vector_path(mapped):
    from webapp import explore_view

    bundle = explore_bundle_build(mapped)
    gene = bundle["geneList"][0]
    n_cells = mapped["adata"].n_obs

    vec = explore_view.gene_vector(str(mapped["h5ad"]), gene)
    assert vec is not None
    assert vec["gene"] == gene
    assert len(vec["values"]) == n_cells

    assert explore_view.gene_vector(str(mapped["h5ad"]), "NOT_A_GENE") is None


def explore_bundle_build(mapped):
    from webapp import explore_bundle

    return explore_bundle.build(
        mapped["map_result"], mapped["adata"], mapped["summary"])


def test_branches_from_summary_empty_axes():
    from webapp.explorer import branches_from_summary

    empty = {"blueprint": {"axes": [], "modules": []}}
    assert branches_from_summary(empty) == []
