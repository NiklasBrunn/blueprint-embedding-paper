from pathlib import Path

import webapp.editing as E
from webapp.mapping import REPO_ROOT


def _archived_source() -> str:
    bp = next((REPO_ROOT / "outs" / "blueprints").rglob("embedder_blueprint_seed*.py"))
    return bp.read_text()


def test_edited_path_under_token(tmp_path, monkeypatch):
    monkeypatch.setattr(E, "RUNTIME", tmp_path)
    assert E.edited_blueprint_path("tok") == \
        tmp_path / "tok" / "edited" / "embedder_blueprint.py"


def test_save_edited_valid_writes_and_validates(tmp_path, monkeypatch):
    monkeypatch.setattr(E, "RUNTIME", tmp_path)
    ok, errors = E.save_edited("tok", _archived_source())
    assert ok, errors
    assert E.edited_blueprint_path("tok").exists()


def test_save_edited_invalid_still_written_but_flagged(tmp_path, monkeypatch):
    monkeypatch.setattr(E, "RUNTIME", tmp_path)
    bad = ("def build_blueprint(name):\n"
           "    return None  # <<EXAMPLE_DELETE_ME>> placeholder left in\n")
    ok, errors = E.save_edited("tok", bad)
    assert not ok and errors
    assert E.edited_blueprint_path("tok").exists()  # work-in-progress persists


def test_resolve_prefers_edited_then_generated(tmp_path, monkeypatch):
    monkeypatch.setattr(E, "RUNTIME", tmp_path)
    gen_dir = tmp_path / "sess" / "run01"
    gen_dir.mkdir(parents=True)
    (gen_dir / "embedder_blueprint.py").write_text("# generated\n")
    active = {"workdir": str(gen_dir)}
    # No edited copy -> the generated file.
    assert E.resolve_session_blueprint(active, "tok") == \
        gen_dir / "embedder_blueprint.py"
    # Edited copy present -> it wins.
    ep = E.edited_blueprint_path("tok")
    ep.parent.mkdir(parents=True, exist_ok=True)
    ep.write_text("# edited\n")
    assert E.resolve_session_blueprint(active, "tok") == ep


def test_revert_removes_edited(tmp_path, monkeypatch):
    monkeypatch.setattr(E, "RUNTIME", tmp_path)
    p = E.edited_blueprint_path("tok")
    p.parent.mkdir(parents=True)
    p.write_text("x")
    E.revert_edited("tok")
    assert not p.exists()
    E.revert_edited("tok")  # idempotent / no error when absent
