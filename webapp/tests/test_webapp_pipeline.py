import json
from webapp.pipeline import write_scenario_config, TEMPORAL_NOTE

def test_write_scenario_minimal(tmp_path):
    p = tmp_path / "scenario_config.jsonc"
    write_scenario_config(p, species=["mouse"], tissue=["cortex"], temporal=False)
    cfg = json.loads(p.read_text())
    assert cfg["species"] == ["mouse"] and cfg["tissue"] == ["cortex"]
    assert cfg["celltype"] is None and cfg["notes"] is None

def test_write_scenario_temporal_injects_note(tmp_path):
    p = tmp_path / "scenario_config.jsonc"
    write_scenario_config(p, species=["mouse"], tissue=["embryo"], temporal=True)
    assert json.loads(p.read_text())["notes"] == TEMPORAL_NOTE

def test_write_scenario_custom_note_appended(tmp_path):
    p = tmp_path / "scenario_config.jsonc"
    write_scenario_config(p, species=["human"], tissue=None, temporal=True,
                          notes="emphasize myeloid")
    notes = json.loads(p.read_text())["notes"]
    assert "emphasize myeloid" in notes and TEMPORAL_NOTE in notes

def test_build_driver_command(tmp_path):
    from webapp.pipeline import build_driver_command
    cmd = build_driver_command(
        scenario_file=tmp_path/"s.jsonc", input_h5ad=tmp_path/"d.h5ad",
        work_root=tmp_path/"runs", model="gpt-5.5", reasoning="high",
        with_refine=True, max_iterations=3, noref=False,
        no_embed=True, dataset="webapp", n_parallel=1)
    assert "--with-refine" in cmd and "--non-interactive" in cmd
    assert "--no-embed" in cmd and "--work-root" in cmd
    assert not any("sk-" in str(c) for c in cmd)  # key never in argv

def test_driver_env_carries_key():
    from webapp.pipeline import driver_env
    assert driver_env("sk-TESTKEY")["OPENAI_API_KEY"] == "sk-TESTKEY"

def test_driver_env_no_key_when_none():
    from webapp.pipeline import driver_env
    env = driver_env(None)
    # should not raise; OPENAI_API_KEY only set if provided
    assert isinstance(env, dict)

def test_prepare_run_workdir(tmp_path, monkeypatch):
    import webapp.pipeline as P
    monkeypatch.setattr(P, "RUNTIME", tmp_path/"runtime")
    wd, scen = P.prepare_run_workdir(session="sess1", run=1,
                                     species=["mouse"], tissue=["cortex"],
                                     temporal=False, noref=False)
    assert (wd/"embedder_blueprint.py").exists()
    assert (wd/"scenario_config.jsonc").exists() and scen == wd/"scenario_config.jsonc"
    assert (wd/"agent_prompt_blueprint_edit.md").exists()
    assert (wd/"logs").is_dir()


def test_write_scenario_no_duplicate_when_note_already_present(tmp_path):
    # The UI may auto-fill the temporal note into the notes field; writing with
    # temporal=True must not append a second copy.
    p = tmp_path / "scenario_config.jsonc"
    preset = f"custom rule\n\n{TEMPORAL_NOTE}"
    write_scenario_config(p, species=["mouse"], tissue=None, temporal=True,
                          notes=preset)
    notes = json.loads(p.read_text())["notes"]
    assert notes.count(TEMPORAL_NOTE) == 1
    assert "custom rule" in notes


def test_resolve_run_dir_from_log(tmp_path):
    from webapp.pipeline import resolve_run_dir
    app_wd = tmp_path / "run01"
    (app_wd / "logs").mkdir(parents=True)
    driver_dir = tmp_path / "cohort" / "260620_blueprint_run01"
    driver_dir.mkdir(parents=True)
    (driver_dir / "embedder_blueprint.py").write_text("# populated\n")
    (app_wd / "logs" / "codex_edit.log").write_text(
        f"Prepared working directory: {driver_dir}\nStage A: ...\n")
    assert resolve_run_dir(app_wd) == driver_dir


def test_resolve_run_dir_fallback_rglob(tmp_path):
    from webapp.pipeline import resolve_run_dir
    app_wd = tmp_path / "run01"
    (app_wd / "logs").mkdir(parents=True)
    (app_wd / "logs" / "codex_edit.log").write_text("no marker line here\n")
    driver_dir = tmp_path / "cohort" / "260620_blueprint_run01"
    driver_dir.mkdir(parents=True)
    (driver_dir / "embedder_blueprint.py").write_text("# populated\n")
    (driver_dir / "summary.json").write_text("{}")
    assert resolve_run_dir(app_wd) == driver_dir


def test_resolve_run_dir_falls_back_to_self(tmp_path):
    from webapp.pipeline import resolve_run_dir
    app_wd = tmp_path / "sess" / "run01"
    app_wd.mkdir(parents=True)
    assert resolve_run_dir(app_wd) == app_wd
