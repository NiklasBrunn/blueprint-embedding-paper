def test_build_app_has_tabs():
    from webapp.server import build_app
    html = str(build_app().layout)
    assert "Generate" in html and "Explore" in html

def test_phase_tracker_parses_log():
    from webapp.server import parse_phase
    assert parse_phase(["Stage A (populator) for 1 workdir(s)."]) == "Stage A"
    assert parse_phase(["stage_b iter=01 discriminator starting (isolated)"]) == "Stage B"
    assert parse_phase(["Embed: embedding accepted blueprints"]) == "Embedding"

def test_phase_tracker_done_and_idle():
    from webapp.server import parse_phase
    assert parse_phase(["Session log: /x/y.log"]) == "Done"
    assert parse_phase([]) == "Idle"


def test_layout_has_branded_header():
    from webapp.server import build_app
    app = build_app()
    html_str = str(app.layout)
    assert "ukf-header" in html_str
    assert "Blueprint embeddings" in html_str


def test_bp_upload_visibility_toggle():
    from webapp.server import _bp_upload_style
    assert _bp_upload_style("upload") == {"display": "block"}
    assert _bp_upload_style("session") == {"display": "none"}


def test_bp_ready_chip_after_generate():
    from webapp.server import _bp_ready_chip
    chip = _bp_ready_chip({"workdir": "/tmp/run01"})
    assert "ready" in str(chip).lower()
    assert _bp_ready_chip(None) == ""


def test_stepper_marks_active_and_done():
    from webapp.server import _stepper
    kids = _stepper("Stage B")
    classes = [k.className for k in kids]
    # Stage A done, Stage B active, later steps plain.
    assert any("done" in c for c in classes)
    assert any("active" in c for c in classes)


def test_stepper_done_marks_all_done():
    from webapp.server import _stepper
    kids = _stepper("Done")
    assert all("done" in k.className for k in kids)


def test_autofill_notes_toggle_on_off():
    from webapp.server import _autofill_notes_value
    from webapp.pipeline import TEMPORAL_NOTE
    # Toggle ON with empty notes -> notes becomes the temporal note.
    assert _autofill_notes_value(["on"], [], "") == TEMPORAL_NOTE
    # Toggle ON preserves custom text and appends the note once.
    v = _autofill_notes_value(["on"], [], "focus on myeloid")
    assert "focus on myeloid" in v and TEMPORAL_NOTE in v
    # Toggle ON again when already present -> no change (no duplicate).
    assert _autofill_notes_value(["on"], [], v) is None
    # Toggle OFF removes the temporal block, keeps custom text.
    assert _autofill_notes_value([], [], v) == "focus on myeloid"
    # Toggle OFF with only the temporal note -> empty.
    assert _autofill_notes_value([], [], TEMPORAL_NOTE) == ""
    # Toggle OFF when no note present -> no change.
    assert _autofill_notes_value([], [], "plain") is None


def test_autofill_notes_celltype_toggle():
    from webapp.server import _autofill_notes_value
    from webapp.pipeline import CELLTYPE_NOTE, TEMPORAL_NOTE
    # Cell-type toggle ON with empty notes -> notes becomes the cell-type note.
    assert _autofill_notes_value([], ["on"], "") == CELLTYPE_NOTE
    # OFF removes only the cell-type block, keeps custom text.
    v = _autofill_notes_value([], ["on"], "keep me")
    assert "keep me" in v and CELLTYPE_NOTE in v
    assert _autofill_notes_value([], [], v) == "keep me"
    # Both toggles ON -> both blocks present; the two notes are independent.
    both = _autofill_notes_value(["on"], ["on"], "")
    assert TEMPORAL_NOTE in both and CELLTYPE_NOTE in both
    # Turning temporal OFF leaves the cell-type block intact.
    only_ct = _autofill_notes_value([], ["on"], both)
    assert TEMPORAL_NOTE not in only_ct and CELLTYPE_NOTE in only_ct


def test_layout_has_blueprint_editor():
    from webapp.server import build_app
    app = build_app()
    s = str(app.layout)
    assert "gen-editor" in s and "gen-done" in s


def test_resolve_h5ad_prefers_local_path(tmp_path):
    from webapp.server import _resolve_h5ad_source
    f = tmp_path / "x.h5ad"; f.write_text("d")
    p, err = _resolve_h5ad_source("tok", None, None, str(f))
    assert err is None and p == f


def test_resolve_h5ad_path_missing():
    from webapp.server import _resolve_h5ad_source
    p, err = _resolve_h5ad_source("tok", None, None, "/no/such/x.h5ad")
    assert p is None and "not found" in err.lower()


def test_resolve_h5ad_path_wrong_ext(tmp_path):
    from webapp.server import _resolve_h5ad_source
    f = tmp_path / "x.txt"; f.write_text("d")
    p, err = _resolve_h5ad_source("tok", None, None, str(f))
    assert p is None and "h5ad" in err.lower()


def test_resolve_h5ad_upload_branch(tmp_path, monkeypatch):
    import base64
    import webapp.server as S
    monkeypatch.setattr(S, "RUNTIME", tmp_path)
    data = b"hello-bytes"
    uri = "data:application/octet-stream;base64," + base64.b64encode(data).decode()
    p, err = S._resolve_h5ad_source("tok", uri, "data.h5ad", None)
    assert err is None and p.read_bytes() == data


def test_resolve_h5ad_nothing_provided():
    from webapp.server import _resolve_h5ad_source
    p, err = _resolve_h5ad_source("tok", None, None, "")
    assert p is None and err


def test_layout_has_browse_button():
    from webapp.server import build_app
    assert "exp-h5ad-browse" in str(build_app().layout)


def test_h5ad_shape_reads_dims(tmp_path):
    from webapp.server import _h5ad_shape
    from webapp.tests.fixtures.make_tiny_h5ad import write_tiny
    h5 = write_tiny(tmp_path / "t.h5ad")
    s = _h5ad_shape(h5)
    assert "120 cells" in s and "genes" in s


def test_resolve_bp_prefers_local_path(tmp_path):
    from webapp.server import _resolve_bp_source
    f = tmp_path / "embedder_blueprint.py"; f.write_text("x")
    p, err = _resolve_bp_source("tok", None, str(f))
    assert err is None and p == f


def test_resolve_bp_path_missing():
    from webapp.server import _resolve_bp_source
    p, err = _resolve_bp_source("tok", None, "/no/such/bp.py")
    assert p is None and "not found" in err.lower()


def test_resolve_bp_wrong_ext(tmp_path):
    from webapp.server import _resolve_bp_source
    f = tmp_path / "bp.txt"; f.write_text("x")
    p, err = _resolve_bp_source("tok", None, str(f))
    assert p is None and ".py" in err


def test_resolve_bp_upload_branch(tmp_path, monkeypatch):
    import base64
    import webapp.server as S
    monkeypatch.setattr(S, "RUNTIME", tmp_path)
    uri = "data:text/x-python;base64," + base64.b64encode(b"# bp").decode()
    p, err = S._resolve_bp_source("tok", uri, None)
    assert err is None and p.read_bytes() == b"# bp"


def test_resolve_bp_nothing_provided():
    from webapp.server import _resolve_bp_source
    p, err = _resolve_bp_source("tok", None, "")
    assert p is None and err


def test_layout_has_bp_browse():
    from webapp.server import build_app
    assert "exp-bp-browse" in str(build_app().layout)


def test_completion_status_success(tmp_path):
    from webapp.server import _completion_status
    app_wd = tmp_path / "run01"; (app_wd / "logs").mkdir(parents=True)
    run_dir = tmp_path / "cohort" / "260622_blueprint_run01"
    run_dir.mkdir(parents=True)
    (run_dir / "embedder_blueprint.py").write_text("# bp")
    (run_dir / ".stage_a_rc").write_text("0")
    (app_wd / "logs" / "codex_edit.log").write_text(
        f"Prepared working directory: {run_dir}\n")
    label, tail = _completion_status({"workdir": str(app_wd)})
    assert label == "Done" and tail == ""


def test_completion_status_failure_surfaces_codex_log(tmp_path):
    from webapp.server import _completion_status
    app_wd = tmp_path / "run01"; (app_wd / "logs").mkdir(parents=True)
    run_dir = tmp_path / "cohort" / "260622_blueprint_run01"
    (run_dir / "logs").mkdir(parents=True)
    (run_dir / "embedder_blueprint.py").write_text("# bp")
    (run_dir / ".stage_a_rc").write_text("1")
    (run_dir / "logs" / "codex_edit.log").write_text(
        "stream error: Not in a trusted directory\n")
    (app_wd / "logs" / "codex_edit.log").write_text(
        f"Prepared working directory: {run_dir}\n")
    label, tail = _completion_status({"workdir": str(app_wd)})
    assert "failed" in label.lower() and "codex exit 1" in label
    assert "trusted directory" in tail
