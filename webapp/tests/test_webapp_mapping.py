import shutil
from pathlib import Path

from webapp.mapping import validate_blueprint

REPO = Path(__file__).resolve().parents[2]


def test_validate_clean_blueprint(tmp_path):
    src = next((REPO / "outs" / "blueprints").rglob("embedder_blueprint_seed*.py"))
    bp = tmp_path / "embedder_blueprint.py"
    shutil.copy(src, bp)
    ok, summary, errors = validate_blueprint(bp)
    assert ok, errors
    assert summary["blueprint"]["n_axes"] >= 1


def test_validate_rejects_placeholder(tmp_path):
    bp = tmp_path / "embedder_blueprint.py"
    bp.write_text("# template with <<EXAMPLE_MODULE>> placeholder, no build_blueprint\n")
    ok, summary, errors = validate_blueprint(bp)
    assert not ok


def test_map_h5ad_produces_dims_and_umap(tmp_path):
    from webapp.mapping import map_h5ad
    from webapp.tests.fixtures.make_tiny_h5ad import write_tiny
    h5 = write_tiny(tmp_path/"tiny.h5ad")
    bp = next((REPO/"outs"/"blueprints").rglob("embedder_blueprint_seed*.py"))
    res = map_h5ad(blueprint=bp, h5ad=h5, outdir=tmp_path)
    import pandas as pd
    df = pd.read_parquet(res.parquet)
    assert any(c.startswith("dim_") for c in df.columns)
    assert "cell_id" in df.columns
    assert res.umap_coords.exists()
    u = pd.read_parquet(res.umap_coords)
    assert {"cell_id","umap_0","umap_1"}.issubset(u.columns)
    assert len(u) == len(df)


def test_axis_assignment_argmax():
    import pandas as pd

    from webapp.mapping import axis_assignment
    # dim_0..dim_2 with an obvious argmax per row; dim_cols sorted by int suffix.
    emb = pd.DataFrame({
        "cell_id": ["c1", "c2", "c3"],
        "dim_0": [2.0, 0.0, 0.1],
        "dim_2": [0.5, 0.5, 3.0],
        "dim_1": [1.0, 9.0, 0.2],
    })
    s = axis_assignment(emb, ["root", "mid", "leaf"])
    assert list(s.index) == ["c1", "c2", "c3"]
    assert list(s) == ["root", "mid", "leaf"]
    assert s.name == "blueprint_axis_assignment"


def test_validate_rejects_unpopulated_template():
    # The vendored embedder_blueprint.py is the bare template (no axes); it must
    # be rejected so mapping never runs it and hits "no axes registered".
    from webapp.mapping import validate_blueprint, _SRC_ROOT
    tmpl = _SRC_ROOT / "scripts" / "agentic_embedding" / "embedder_blueprint.py"
    ok, _summary, errors = validate_blueprint(tmpl)
    assert not ok
    assert any("no axes" in e for e in errors)
