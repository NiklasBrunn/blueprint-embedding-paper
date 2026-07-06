"""Write scenario configs and (later) launch the vendored blueprint driver."""
from __future__ import annotations
import json, os, shutil, subprocess
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
WEBAPP_DIR = Path(__file__).resolve().parent
DRIVER = WEBAPP_DIR / "pipeline_driver.sh"
RUNTIME = WEBAPP_DIR / "runtime"
# Canonical pipeline sources: prefer the vendored copies (standalone), falling
# back to the parent repo when running in-tree without a vendor/ dir.
VENDOR = WEBAPP_DIR / "vendor"
SRC_ROOT = VENDOR if VENDOR.exists() else REPO_ROOT

# Verbatim notes value from configs/scenario_config_mouse_gastrulation.jsonc.
TEMPORAL_NOTE = "IMPORTANT AGENT CONSTRAINTS: organize axes as a time-forward hierarchical tree matching the scenario; tree depth = time (root = earliest state, leaves = latest). Every non-root axis MUST declare its parent via gate — the parent's module becomes the gate, which both anchors tree parentage and scopes the signal to cells expressing the parent's program (the biologically correct behavior for both sub-states and sequential transitions). For sibling-vs-sibling discrimination at a branch point, pass the sibling modules as negatives on the same gated axis. Use a plain axis (no gate, no parent) ONLY for roots: scenario-given starting states (one or several, as the scenario allows), populations that do not progress within the window, terminal populations whose precursor isn't represented, and cycling/proliferating populations. Naming: <chainPrefix>_<descriptiveTail>; <chainPrefix> is one camelCase token (no underscores) grouping axes from the same lineage. When you gate across chains (parent has a different <chainPrefix>), the child MUST use a fresh, unique <chainPrefix> — reusing a prefix that already exists in a different subtree creates ambiguity. Cover the full cell landscape implied by the scenario; preserve continuous maturation gradients rather than treating states as discrete clusters."  # noqa: E501

# Verbatim celltype-tree notes directive for the webapp's Generate mode.
CELLTYPE_NOTE = "IMPORTANT AGENT CONSTRAINTS: organize axes as a hierarchical cell-type tree; tree depth = annotation granularity (root = broadest class, leaves = finest grounded subtype). Roots are plain axes (no gate) — the tissue's top-level cell classes. Every non-root axis MUST gate on its parent's module (anchoring parentage and scoping its signal to the parent class) and pass its sibling modules as negatives for within-parent discrimination. Refine broad→fine, splitting only where subtypes are molecularly distinct; stop a branch when no further grounded subdivision exists. Naming: <chainPrefix>_<descriptiveTail>, one camelCase token per lineage; a child gated under a parent of a different prefix MUST take a fresh, unique prefix. Cover every level from broad classes to terminal subtypes."  # noqa: E501

def write_scenario_config(path, *, species, tissue, temporal,
                          celltype=None, age=None, sex=None,
                          condition=None, notes=None):
    # Idempotent: if the temporal note is already embedded in `notes` (e.g. the
    # UI auto-filled it into the notes field), pass it through unchanged instead
    # of appending a duplicate.
    if temporal and notes and TEMPORAL_NOTE not in notes:
        note = f"{notes}\n\n{TEMPORAL_NOTE}"
    elif temporal and not notes:
        note = TEMPORAL_NOTE
    else:
        note = notes
    cfg = {"species": species, "tissue": tissue, "celltype": celltype,
           "age": age, "sex": sex, "condition": condition, "notes": note}
    Path(path).write_text(json.dumps(cfg, indent=2))
    return cfg

def prepare_run_workdir(*, session, run, species, tissue, temporal,
                        noref=False, celltype=None, age=None, sex=None,
                        condition=None, notes=None):
    """Create RUNTIME/<session>/run<NN>/ (+ logs/), copy the canonical
    blueprint source + edit prompt in, write scenario_config.jsonc.
    Returns (workdir: Path, scenario_path: Path)."""
    workdir = RUNTIME / session / f"run{run:02d}"
    workdir.mkdir(parents=True, exist_ok=True)
    (workdir / "logs").mkdir(exist_ok=True)
    blueprint_src = (SRC_ROOT / "scripts" / "agentic_embedding" /
                     ("embedder_blueprint_noref.py" if noref
                      else "embedder_blueprint.py"))
    prompt_src = (SRC_ROOT / "prompt_templates" /
                  ("agent_prompt_blueprint_edit_noref.md" if noref
                   else "agent_prompt_blueprint_edit.md"))
    shutil.copy(blueprint_src, workdir / "embedder_blueprint.py")
    shutil.copy(prompt_src, workdir / "agent_prompt_blueprint_edit.md")
    scen = workdir / "scenario_config.jsonc"
    write_scenario_config(scen, species=species, tissue=tissue,
                          temporal=temporal, celltype=celltype, age=age,
                          sex=sex, condition=condition, notes=notes)
    return workdir, scen

def resolve_run_dir(app_workdir):
    """Return the driver's actual run dir (where the populated blueprint and
    summary.json live).

    The vendored driver, given ``--work-root``, creates its OWN cohort-named run
    dir (e.g. ``<work_root>/<cohort>/<ts>_blueprint_run01/``) and edits the
    blueprint there — it never touches the app's prepared ``run<NN>/`` dir, which
    only holds the bare template + scenario + the captured driver log. The driver
    prints its run dir in that log as ``Prepared working directory: <path>``;
    parse it, falling back to the newest populated run dir under the session
    root, and finally to ``app_workdir`` itself.
    """
    app_workdir = Path(app_workdir)
    log = app_workdir / "logs" / "codex_edit.log"
    if log.exists():
        for line in log.read_text().splitlines():
            if line.startswith("Prepared working directory:"):
                cand = Path(line.split(":", 1)[1].strip())
                if (cand / "embedder_blueprint.py").exists():
                    return cand
    session_root = app_workdir.parent
    populated = [p.parent for p in session_root.rglob("summary.json")
                 if (p.parent / "embedder_blueprint.py").exists()]
    if populated:
        return max(populated, key=lambda d: d.stat().st_mtime)
    return app_workdir


def build_driver_command(*, scenario_file, input_h5ad, work_root, model,
                         reasoning, with_refine, max_iterations, noref,
                         no_embed, dataset="webapp", n_parallel=1):
    # The OpenAI API key is NEVER added to this command. It is passed to the
    # subprocess exclusively through driver_env() (the OPENAI_API_KEY env var).
    cmd = ["bash", str(DRIVER), "--dataset", dataset, "--model", model,
           "--reasoning-effort", reasoning,
           "--scenario-file", str(scenario_file), "--input", str(input_h5ad),
           "--work-root", str(work_root),
           "--non-interactive", "--n-parallel", str(n_parallel)]
    if with_refine:
        cmd += ["--with-refine", "--max-iterations", str(max_iterations)]
    if noref:
        cmd += ["--noref"]
    if no_embed:
        cmd += ["--no-embed"]
    return cmd

def driver_env(api_key):
    # The key lives ONLY in the returned env dict; never on argv, never on disk,
    # never logged.
    env = dict(os.environ)
    if api_key:
        env["OPENAI_API_KEY"] = api_key
    return env

def launch_driver(*, log_path, api_key, **kw):
    """Start the vendored driver as a subprocess; stdout+stderr -> log_path.
    Returns the Popen handle. The API key is injected via env only (never on
    the command line, never written into the log)."""
    cmd = build_driver_command(**kw)
    log = open(log_path, "w")
    return subprocess.Popen(cmd, cwd=str(REPO_ROOT), env=driver_env(api_key),
                            stdout=log, stderr=subprocess.STDOUT, text=True)
