"""Plotly Dash app: Generate tab (launch the blueprint driver) + Explore stub.

The OpenAI API key never touches argv or disk: the Generate callback stashes it
in the module-level ``_KEYS`` dict keyed by a per-session token, passes it to the
subprocess via ``driver_env`` (env-only), clears it from the browser field on
launch, and pops it from ``_KEYS`` the moment the job stops running.
"""
from __future__ import annotations

import base64
import json
import uuid
from pathlib import Path

import dash_ace
import dash_bootstrap_components as dbc
from dash import Dash, dcc, html, Input, Output, State, no_update

from webapp import editing, explore_view, filepicker, jobs
from webapp.explore_bundle import build as build_bundle
from webapp.mapping import axis_assignment, map_h5ad, validate_blueprint
from webapp.pipeline import (
    CELLTYPE_NOTE, REPO_ROOT, RUNTIME, TEMPORAL_NOTE, build_driver_command,
    driver_env, prepare_run_workdir, resolve_run_dir,
)

# Single shared registry for all background driver subprocesses.
REGISTRY = jobs.JobRegistry()

# Server-side, in-memory key store: {session_token: api_key}. Never persisted.
_KEYS: dict[str, str] = {}


def clear_key(token: str) -> None:
    """Drop the stored API key for a session (called when a job stops running)."""
    _KEYS.pop(token, None)


def parse_phase(log_lines: list[str]) -> str:
    """Return the most advanced pipeline phase reached in the driver log."""
    joined = "\n".join(log_lines)
    if "Session log:" in joined:
        return "Done"
    if "Embed" in joined:
        return "Embedding"
    if "stage_b iter=" in joined:
        return "Stage B"
    if "Stage A:" in joined or "Stage A" in joined:
        return "Stage A"
    return "Idle"


_PHASE_ORDER = ["Stage A", "Stage B", "Embedding", "Done"]


def _stepper(phase):
    """Render the phase stepper: steps before the current phase are 'done',
    the current one is 'active' (unless 'Done', which marks all done)."""
    try:
        cur = _PHASE_ORDER.index(phase)
    except ValueError:
        cur = -1
    steps = []
    for i, label in enumerate(_PHASE_ORDER):
        cls = "ukf-step"
        if phase == "Done" or i < cur:
            cls += " done"
        elif i == cur:
            cls += " active"
        steps.append(html.Div(label, className=cls))
    return steps


def _autofill_notes_value(temporal, celltype, notes):
    """Keep the notes field in sync with the structural-tree toggles.

    Toggling a mode ON appends its rules block (once); toggling OFF removes that
    block, preserving any custom text. The temporal-tree and cell-type-tree
    blocks are handled independently (either, both, or neither may be present).
    Returns the new notes string, or None when nothing should change (avoids
    clobbering edits)."""
    notes = notes or ""
    changed = False
    for on, note in ((bool(temporal), TEMPORAL_NOTE),
                     (bool(celltype), CELLTYPE_NOTE)):
        has = note in notes
        if on and not has:
            notes = f"{notes}\n\n{note}".strip() if notes.strip() else note
            changed = True
        elif not on and has:
            notes = notes.replace(note, "").strip()
            changed = True
    return notes if changed else None


def _completion_status(active):
    """After the driver exits, inspect the real run dir so a failed Stage A
    surfaces codex's ACTUAL error (from the cohort workdir's codex log) instead
    of a misleading 'Done'. Returns ``(label, codex_log_tail)``; label == 'Done'
    only on success."""
    try:
        run_dir = resolve_run_dir(active["workdir"])
    except Exception:
        return "Done", ""
    rc_file = run_dir / ".stage_a_rc"
    rc = rc_file.read_text().strip() if rc_file.exists() else None
    bp = run_dir / "embedder_blueprint.py"
    if rc == "0" and bp.exists():
        return "Done", ""
    clog = run_dir / "logs" / "codex_edit.log"
    tail = ""
    if clog.exists():
        tail = "\n".join(clog.read_text().splitlines()[-40:]).strip()
    detail = f" (codex exit {rc})" if rc not in (None, "0") else ""
    return f"Generation failed{detail}", tail


def _bp_upload_style(source):
    """Show the .py upload only when 'uploaded blueprint' is selected."""
    return {"display": "block"} if source == "upload" else {"display": "none"}


def _bp_ready_chip(active):
    """Green chip signalling an in-session blueprint is ready to map."""
    if not active or not active.get("workdir"):
        return ""
    return html.Span("✓ Blueprint ready from this session",
                     className="ukf-chip")


def _csv_list(value):
    """Comma-separated string -> list of trimmed items, or None if empty."""
    if not value:
        return None
    items = [v.strip() for v in str(value).split(",") if v.strip()]
    return items or None


def _next_run(session: str) -> int:
    """Pick the next run number for a session by counting existing run dirs."""
    base = RUNTIME / session
    if not base.exists():
        return 1
    existing = [p for p in base.iterdir() if p.is_dir() and p.name.startswith("run")]
    return len(existing) + 1


def _generate_form() -> html.Div:
    sex_opts = [{"label": s, "value": s} for s in ("male", "female", "mixed")]
    sex_opts.append({"label": "(unspecified)", "value": ""})
    scenario_card = dbc.Card([
        dbc.CardHeader("Scenario"),
        dbc.CardBody([
            dbc.Label("Species (comma-separated, e.g. mouse or mouse,human)"),
            dbc.Input(id="gen-species", type="text", value="mouse"),
            dbc.Label("Tissue (comma-separated; blank = none)",
                      class_name="mt-2"),
            dbc.Input(id="gen-tissue", type="text", value=""),
            dbc.Label("Age (optional)", class_name="mt-2"),
            dbc.Input(id="gen-age", type="text", value=""),
            dbc.Label("Sex (optional)", class_name="mt-2"),
            dcc.Dropdown(id="gen-sex", options=sex_opts, value=""),
            dbc.Label("Condition (comma-separated, optional)", class_name="mt-2"),
            dbc.Input(id="gen-condition", type="text", value=""),
            dbc.Checklist(id="gen-temporal", class_name="mt-2",
                          options=[{"label": " Temporal analysis",
                                    "value": "on"}], value=[]),
            dbc.Checklist(id="gen-celltype-tree", class_name="mt-2",
                          options=[{"label": " Cell-type tree",
                                    "value": "on"}], value=[]),
            dbc.Label("Notes / additional rules (optional)", class_name="mt-2"),
            dcc.Textarea(id="gen-notes", value="",
                         placeholder="Extra instructions for the agent. "
                                     "Toggling 'Temporal analysis' or "
                                     "'Cell-type tree' fills this with the "
                                     "corresponding tree rules, which you can "
                                     "then edit.",
                         style={"width": "100%", "minHeight": "90px"}),
        ]),
    ])
    model_card = dbc.Card([
        dbc.CardHeader("Model"),
        dbc.CardBody([
            dbc.Label("Model"),
            dcc.Dropdown(id="gen-model",
                         options=[{"label": m, "value": m} for m in
                                  ("gpt-5.5", "gpt-5.4", "gpt-5.4-mini",
                                   "gpt-5.3-codex")], value="gpt-5.5"),
            dbc.Label("Reasoning effort", class_name="mt-2"),
            dcc.Dropdown(id="gen-reasoning",
                         options=[{"label": r, "value": r} for r in
                                  ("low", "medium", "high", "xhigh")],
                         value="high"),
            dbc.Checklist(id="gen-refine", class_name="mt-2",
                          options=[{"label": " Stage B refinement",
                                    "value": "on"}], value=[]),
            dbc.Label("Max iterations", class_name="mt-2"),
            dbc.Input(id="gen-maxiter", type="number", value=3),
            dbc.Checklist(id="gen-noref", class_name="mt-2",
                          options=[{"label": " Resource-unconstrained (noref)",
                                    "value": "on"}], value=[]),
        ]),
    ])
    creds_card = dbc.Card([
        dbc.CardHeader("Credentials"),
        dbc.CardBody([
            dbc.Input(id="gen-apikey", type="password",
                      placeholder="OpenAI API key (kept in memory only)",
                      value=""),
            dbc.Button("Generate", id="gen-run-btn", n_clicks=0,
                       color="success", class_name="mt-3 w-100"),
        ]),
    ])
    status_card = dbc.Card([
        dbc.CardHeader("Status"),
        dbc.CardBody([
            html.Div(id="gen-stepper", className="ukf-stepper"),
            html.Div(id="gen-eta", className="text-muted mb-2"),
            html.Div(id="gen-phase", children="Idle", className="fw-bold mb-2"),
            html.Pre(id="gen-log", className="ukf-log"),
            html.Hr(),
            html.Div(className="d-flex flex-wrap gap-2", children=[
                dbc.Button("Download blueprint", id="gen-dl-bp-btn",
                           n_clicks=0, color="primary", size="sm"),
                dcc.Download(id="gen-dl-bp"),
                dbc.Button("Download summary.json", id="gen-dl-summary-btn",
                           n_clicks=0, color="primary", size="sm"),
                dcc.Download(id="gen-dl-summary"),
                dbc.Button("Download report.md", id="gen-dl-report-btn",
                           n_clicks=0, color="primary", size="sm"),
                dcc.Download(id="gen-dl-report"),
            ]),
        ]),
    ])
    editor_card = dbc.Card(class_name="mt-3", children=[
        dbc.CardHeader("Edit blueprint"),
        dbc.CardBody([
            html.Div("After generation the blueprint loads here. Edit it, then "
                     "Save & validate — “Use in-session blueprint” on "
                     "the Explore tab will map your edited version.",
                     className="text-muted mb-2"),
            dash_ace.DashAceEditor(
                id="gen-editor", value="", mode="python", theme="github",
                tabSize=4, fontSize=13, showGutter=True, showPrintMargin=False,
                width="100%", height="460px",
                placeholder="Generate a blueprint first…"),
            html.Div(className="d-flex gap-2 mt-2 align-items-center", children=[
                dbc.Button("Save & validate", id="gen-editor-save-btn",
                           n_clicks=0, color="success", size="sm"),
                dbc.Button("Revert to generated", id="gen-editor-revert-btn",
                           n_clicks=0, color="primary", outline=True, size="sm"),
                html.Span(id="gen-editor-status", className="ms-2"),
            ]),
        ]),
    ])
    return html.Div([
        dcc.Store(id="gen-done"),
        dbc.Row([
            dbc.Col([scenario_card, html.Div(className="mt-3"), model_card,
                     html.Div(className="mt-3"), creds_card], md=5),
            dbc.Col(status_card, md=7),
        ]),
        editor_card,
    ])


def _explore_panel() -> html.Div:
    return explore_view.panel()


def _save_upload(token: str, contents: str, filename: str) -> Path:
    """Decode a dcc.Upload data-URI and write it under RUNTIME/<token>/uploads/."""
    raw = base64.b64decode(contents.split(",", 1)[1])
    updir = RUNTIME / token / "uploads"
    updir.mkdir(parents=True, exist_ok=True)
    dest = updir / filename
    dest.write_bytes(raw)
    return dest


def _resolve_bp_source(token, contents, path):
    """Resolve the uploaded-blueprint source: a typed/picked local .py path
    (used directly from disk) takes precedence over a dragged-in upload.
    Returns ``(Path | None, error | None)``."""
    if path and path.strip():
        p = Path(path.strip()).expanduser()
        if not p.exists():
            return None, f"File not found: {p}"
        if p.suffix != ".py":
            return None, "Blueprint must be a .py file."
        return p, None
    if contents:
        return _save_upload(token, contents, "embedder_blueprint.py"), None
    return None, "Select a blueprint .py (Browse or drag-and-drop) first."


def _h5ad_shape(path) -> str:
    """' — N cells × M genes' read cheaply from an h5ad header, or '' if the
    shape can't be read quickly (e.g. a .gz)."""
    try:
        import anndata as ad
        a = ad.read_h5ad(str(path), backed="r")
        n, m = a.shape
        try:
            a.file.close()
        except Exception:
            pass
        return f" — {n:,} cells × {m:,} genes"
    except Exception:
        return ""


def _resolve_h5ad_source(token, contents, filename, path):
    """Resolve the h5ad to map. A typed local path (read straight from disk, no
    size limit — best for large files) takes precedence over a dragged-in
    upload. Returns ``(Path | None, error | None)``."""
    if path and path.strip():
        p = Path(path.strip()).expanduser()
        if not p.exists():
            return None, f"File not found: {p}"
        if not p.name.endswith((".h5ad", ".h5ad.gz")):
            return None, "Local file must be a .h5ad or .h5ad.gz."
        return p, None
    if contents:
        return _save_upload(token, contents, filename or "data.h5ad"), None
    return None, "Drop an h5ad file or paste a local path first."


def _header() -> html.Div:
    return html.Div(className="ukf-header", children=[
        html.Div("Blueprint embeddings", className="eyebrow"),
        html.H1("Interpretable single-cell embeddings",
                className="ukf-title"),
        html.P("Generate a literature-grounded blueprint using coding agents, "
               "then explore the embedding.", className="ukf-subtitle"),
    ])


def build_app() -> Dash:
    app = Dash(__name__, suppress_callback_exceptions=True,
               external_stylesheets=[dbc.themes.FLATLY])
    token = uuid.uuid4().hex

    app.layout = html.Div([
        _header(),
        dcc.Store(id="session-token", data=token, storage_type="session"),
        dcc.Store(id="active-job"),
        dcc.Interval(id="gen-poll", interval=1500, disabled=True),
        dbc.Container(fluid=True, className="py-3", children=[
            dcc.Tabs(id="main-tabs", value="generate",
                     parent_className="ukf-tabs",
                     className="ukf-tabs-bar", children=[
                dcc.Tab(label="Generate", value="generate",
                        className="ukf-tab",
                        selected_className="ukf-tab--active",
                        children=[_generate_form()]),
                dcc.Tab(label="Explore", value="explore",
                        className="ukf-tab",
                        selected_className="ukf-tab--active", children=[
                    html.Div(id="explore-root", children=[_explore_panel()]),
                ]),
            ]),
        ]),
    ])

    _register_callbacks(app)
    return app


def _register_callbacks(app: Dash) -> None:
    explore_view.register(app)

    @app.callback(
        Output("active-job", "data"),
        Output("gen-poll", "disabled"),
        Output("gen-apikey", "value"),
        Input("gen-run-btn", "n_clicks"),
        State("gen-species", "value"), State("gen-tissue", "value"),
        State("gen-age", "value"), State("gen-sex", "value"),
        State("gen-condition", "value"), State("gen-temporal", "value"),
        State("gen-model", "value"), State("gen-reasoning", "value"),
        State("gen-refine", "value"), State("gen-maxiter", "value"),
        State("gen-noref", "value"),
        State("gen-notes", "value"),
        State("gen-apikey", "value"), State("session-token", "data"),
        prevent_initial_call=True,
    )
    def _launch(n_clicks, species, tissue, age, sex, condition, temporal,
                model, reasoning, refine, maxiter, noref, notes, apikey,
                token):
        if not n_clicks:
            return no_update, no_update, no_update
        token = token or uuid.uuid4().hex
        if apikey:
            _KEYS[token] = apikey
        # A new run supersedes any earlier in-browser edits for this session.
        editing.revert_edited(token)

        session = token
        run = _next_run(session)
        workdir, scen = prepare_run_workdir(
            session=session, run=run,
            species=_csv_list(species) or ["mouse"],
            tissue=_csv_list(tissue),
            temporal=bool(temporal),
            noref=bool(noref),
            age=age or None,
            sex=sex or None,
            condition=_csv_list(condition),
            notes=(notes or None),
        )
        # Generation runs with --no-embed, so no real h5ad is consumed. We pass
        # the scenario file itself as --input: it is a guaranteed-existing path,
        # so the driver's required-file check passes, and --no-embed means it is
        # never actually read as data.
        cmd = build_driver_command(
            scenario_file=scen, input_h5ad=scen,
            work_root=RUNTIME / session, model=model, reasoning=reasoning,
            with_refine=bool(refine), max_iterations=int(maxiter or 3),
            noref=bool(noref), no_embed=True,
            dataset="webapp", n_parallel=1,
        )
        log_path = workdir / "logs" / "codex_edit.log"
        jid = REGISTRY.start_subprocess(
            cmd, cwd=REPO_ROOT, env=driver_env(_KEYS.get(token)),
            log_path=log_path,
        )
        active = {"jid": jid, "workdir": str(workdir)}
        # Clear the key from the browser field after launching.
        return active, False, ""

    @app.callback(
        Output("gen-stepper", "children"),
        Output("gen-phase", "children"),
        Output("gen-log", "children"),
        Output("gen-poll", "disabled", allow_duplicate=True),
        Output("gen-done", "data"),
        Input("gen-poll", "n_intervals"),
        State("active-job", "data"), State("session-token", "data"),
        prevent_initial_call=True,
    )
    def _poll(_n, active, token):
        if not active:
            return (no_update,) * 5
        jid = active["jid"]
        lines = REGISTRY.tail(jid, 60)
        status = REGISTRY.get(jid).status
        phase = parse_phase(lines)
        if status != "running":
            if token:
                clear_key(token)
            label, codex_tail = _completion_status(active)
            ok = label == "Done"
            display = "\n".join(lines)
            if codex_tail:
                display += "\n\n--- codex log (Stage A error) ---\n" + codex_tail
            phase_text = "Done" if ok else f"⚠ {label} — see log below"
            return (_stepper("Done" if ok else phase), phase_text,
                    display, True, ok)
        return _stepper(phase), phase, "\n".join(lines), no_update, False

    def _make_download(out_id, btn_id, filename):
        @app.callback(
            Output(out_id, "data"),
            Input(btn_id, "n_clicks"),
            State("active-job", "data"),
            prevent_initial_call=True,
        )
        def _dl(n_clicks, active):
            if not n_clicks or not active:
                return no_update
            path = resolve_run_dir(active["workdir"]) / filename
            if not path.exists():
                return no_update
            return dcc.send_file(str(path))
        return _dl

    _make_download("gen-dl-bp", "gen-dl-bp-btn", "embedder_blueprint.py")
    _make_download("gen-dl-summary", "gen-dl-summary-btn", "summary.json")
    _make_download("gen-dl-report", "gen-dl-report-btn", "reports/report.md")

    # --- In-browser blueprint editor -------------------------------------
    @app.callback(
        Output("gen-editor", "value"),
        Output("gen-editor-status", "children"),
        Input("gen-done", "data"),
        State("active-job", "data"),
        prevent_initial_call=True,
    )
    def _load_editor(done, active):
        if not done:
            return no_update, no_update
        src = editing.load_generated_source(active)
        if not src:
            return no_update, no_update
        return src, html.Span("Loaded generated blueprint — edit, then "
                              "Save & validate.", className="text-muted")

    @app.callback(
        Output("gen-editor-status", "children", allow_duplicate=True),
        Input("gen-editor-save-btn", "n_clicks"),
        State("gen-editor", "value"),
        State("session-token", "data"),
        prevent_initial_call=True,
    )
    def _save_editor(n_clicks, source, token):
        if not n_clicks:
            return no_update
        token = token or uuid.uuid4().hex
        if not source or not source.strip():
            return html.Span("Nothing to save yet.", className="text-danger")
        ok, errors = editing.save_edited(token, source)
        if ok:
            return html.Span("✓ Saved & valid — “Use in-session blueprint” "
                             "will map your edits.", className="text-success")
        return html.Span("✗ Saved, but invalid: " + "; ".join(errors),
                         className="text-danger")

    @app.callback(
        Output("gen-editor", "value", allow_duplicate=True),
        Output("gen-editor-status", "children", allow_duplicate=True),
        Input("gen-editor-revert-btn", "n_clicks"),
        State("active-job", "data"),
        State("session-token", "data"),
        prevent_initial_call=True,
    )
    def _revert_editor(n_clicks, active, token):
        if not n_clicks:
            return no_update, no_update
        token = token or uuid.uuid4().hex
        editing.revert_edited(token)
        src = editing.load_generated_source(active)
        return src, html.Span("Reverted to generated blueprint.",
                              className="text-muted")

    # --- Explore tab -------------------------------------------------------
    @app.callback(
        Output("exp-h5ad-name", "children"),
        Input("exp-upload-h5ad", "filename"),
        prevent_initial_call=True,
    )
    def _h5ad_name(filename):
        return f"Selected: {filename}" if filename else no_update

    @app.callback(
        Output("exp-h5ad-path", "value"),
        Input("exp-h5ad-browse", "n_clicks"),
        prevent_initial_call=True,
    )
    def _browse_h5ad(n_clicks):
        if not n_clicks:
            return no_update
        path = filepicker.pick_file_native("Select an .h5ad / .h5ad.gz file")
        return path or no_update  # cancelled -> leave the field unchanged

    @app.callback(
        Output("exp-h5ad-info", "children"),
        Input("exp-h5ad-path", "value"),
        prevent_initial_call=True,
    )
    def _h5ad_path_info(path):
        if not path or not path.strip():
            return ""
        p = Path(path.strip()).expanduser()
        if not p.exists():
            return html.Span(f"File not found: {p}", className="text-danger")
        if not p.name.endswith((".h5ad", ".h5ad.gz")):
            return html.Span("Not a .h5ad / .h5ad.gz file",
                             className="text-warning")
        return html.Span(f"✓ Selected {p.name}{_h5ad_shape(p)} — click Map "
                         "to load and embed.", className="text-success")

    @app.callback(
        Output("exp-bp-name", "children"),
        Input("exp-upload-bp", "filename"),
        prevent_initial_call=True,
    )
    def _bp_name(filename):
        return f"Selected: {filename}" if filename else no_update

    @app.callback(
        Output("exp-bp-path", "value"),
        Input("exp-bp-browse", "n_clicks"),
        prevent_initial_call=True,
    )
    def _browse_bp(n_clicks):
        if not n_clicks:
            return no_update
        path = filepicker.pick_file_native("Select a blueprint .py file")
        return path or no_update

    @app.callback(
        Output("exp-bp-info", "children"),
        Input("exp-bp-path", "value"),
        prevent_initial_call=True,
    )
    def _bp_path_info(path):
        if not path or not path.strip():
            return ""
        p = Path(path.strip()).expanduser()
        if not p.exists():
            return html.Span(f"File not found: {p}", className="text-danger")
        if p.suffix != ".py":
            return html.Span("Blueprint must be a .py file",
                             className="text-warning")
        return html.Span(f"✓ Selected {p.name} — click Map to validate & embed.",
                         className="text-success")

    @app.callback(
        Output("gen-notes", "value"),
        Input("gen-temporal", "value"),
        Input("gen-celltype-tree", "value"),
        State("gen-notes", "value"),
        prevent_initial_call=True,
    )
    def _sync_notes(temporal, celltype, notes):
        new = _autofill_notes_value(temporal, celltype, notes)
        return new if new is not None else no_update

    @app.callback(
        Output("exp-bp-upload-wrap", "style"),
        Input("exp-bp-source", "value"),
        prevent_initial_call=True,
    )
    def _toggle_bp_upload(source):
        return _bp_upload_style(source)

    @app.callback(
        Output("exp-bp-ready", "children"),
        Input("active-job", "data"),
        Input("main-tabs", "value"),
    )
    def _bp_ready(active, _tab):
        return _bp_ready_chip(active)

    @app.callback(
        Output("gen-eta", "children"),
        Input("gen-model", "value"),
        Input("gen-reasoning", "value"),
        Input("gen-refine", "value"),
    )
    def _eta(model, effort, refine):
        from webapp.eta import estimate, format_eta
        secs, n = estimate(model or "gpt-5.5", effort or "high",
                           refine=bool(refine))
        return format_eta(secs, n)

    @app.callback(
        Output("explore-bundle", "data"),
        Output("exp-paths", "data"),
        Output("exp-coverage", "children"),
        Output("exp-current", "children"),
        Input("exp-map-btn", "n_clicks"),
        State("exp-upload-h5ad", "contents"),
        State("exp-upload-h5ad", "filename"),
        State("exp-h5ad-path", "value"),
        State("exp-upload-bp", "contents"),
        State("exp-bp-path", "value"),
        State("exp-bp-source", "value"),
        State("active-job", "data"),
        State("session-token", "data"),
        prevent_initial_call=True,
    )
    def _map(n_clicks, h5ad_contents, h5ad_filename, h5ad_path, bp_contents,
             bp_path, bp_source, active, token):
        if not n_clicks:
            return (no_update,) * 4
        token = token or uuid.uuid4().hex
        err = lambda msg: (no_update, no_update,
                           html.Div(msg, style={"color": "crimson"}), no_update)
        h5ad, h5err = _resolve_h5ad_source(token, h5ad_contents,
                                           h5ad_filename, h5ad_path)
        if h5err:
            return err(h5err)

        # Resolve the blueprint path.
        if bp_source == "upload":
            bp, bperr = _resolve_bp_source(token, bp_contents, bp_path)
            if bperr:
                return err(bperr)
        else:
            if not active or not active.get("workdir"):
                return err("No in-session blueprint yet — run Generate first.")
            bp = editing.resolve_session_blueprint(active, token)
            if not bp or not bp.exists():
                return err("In-session blueprint not found — run Generate first.")

        ok, summary, errors = validate_blueprint(bp)
        if not ok:
            return err("Blueprint invalid: " + "; ".join(errors))

        try:
            result = map_h5ad(blueprint=bp, h5ad=h5ad,
                              outdir=RUNTIME / token / "map")
        except Exception as exc:  # surface mapping failures in the banner
            return err(f"Mapping failed: {exc}")

        import anndata as ad
        import pandas as pd
        a = ad.read_h5ad(h5ad)
        # Argmax axis assignment: label each cell with the axis it scores most
        # strongly on, stored in obs so it can color the Explore plot (and ride
        # along into the downloaded h5ad). Aligned by cell_id == obs_names.
        axis_names = [ax["name"] for ax in summary["blueprint"]["axes"]]
        assign = axis_assignment(pd.read_parquet(result.parquet), axis_names)
        a.obs["blueprint_axis_assignment"] = pd.Categorical(
            assign.reindex(a.obs_names.astype(str)).to_numpy())
        bundle = build_bundle(result, a, summary)
        # Nonce so the explore-bundle Store value changes on EVERY Map click —
        # otherwise Dash deep-compares and, when the same data is re-mapped
        # (UMAP is deterministic), the value is identical and the clientside
        # re-render callback never fires (the UMAP appears not to refresh).
        bundle["_nonce"] = n_clicks
        paths = {"parquet": str(result.parquet),
                 "umap": str(result.umap_coords),
                 "h5ad": str(h5ad),
                 "axis_names": [ax["name"]
                                for ax in summary["blueprint"]["axes"]]}

        banner = _coverage_banner(result.coverage_json)
        bp_label = ("in-session blueprint" if bp_source != "upload"
                    else Path(bp).name)
        current = html.Span(["Currently showing: ",
                             html.Strong(Path(h5ad).name), "  ·  ",
                             html.Strong(bp_label)])
        return bundle, paths, banner, current

    @app.callback(
        Output("exp-dl-h5ad", "data"),
        Input("exp-dl-h5ad-btn", "n_clicks"),
        State("exp-paths", "data"),
        State("session-token", "data"),
        prevent_initial_call=True,
    )
    def _dl_h5ad(n_clicks, paths, token):
        if not n_clicks or not paths:
            return no_update
        import pandas as pd
        from webapp.downloads import augment_adata

        token = token or uuid.uuid4().hex
        emb_df = pd.read_parquet(paths["parquet"])
        umap_df = pd.read_parquet(paths["umap"])
        out = RUNTIME / token / "download" / "augmented.h5ad"
        augment_adata(paths["h5ad"], emb_df, umap_df, paths["axis_names"],
                      out_path=out)
        return dcc.send_file(str(out))

    @app.callback(
        Output("exp-tree-kind", "value"),
        Input("explore-bundle", "data"),
        prevent_initial_call=True,
    )
    def _sync_tree_kind(bundle):
        # Pre-select the tree mode from the mapped blueprint (user can override).
        return (bundle or {}).get("treeKind", "celltype")

    @app.callback(
        Output("exp-dl-tree", "data"),
        Input("exp-dl-tree-btn", "n_clicks"),
        State("exp-tree-kind", "value"),
        State("exp-tree-fmt", "value"),
        State("explore-bundle", "data"),
        prevent_initial_call=True,
    )
    def _dl_tree(n_clicks, kind, fmt, bundle):
        if not n_clicks or not bundle or not bundle.get("axisOrder"):
            return no_update
        import matplotlib.pyplot as plt

        from webapp.tree_plot import figure_from_bundle

        kind = kind or "celltype"
        fmt = (fmt or "png").lower()
        fig = figure_from_bundle(bundle, kind)
        try:
            return dcc.send_bytes(
                lambda buf: fig.savefig(buf, format=fmt, dpi=300,
                                        bbox_inches="tight"),
                f"blueprint_tree_{kind}.{fmt}")
        finally:
            plt.close(fig)


def _coverage_banner(coverage_json):
    """Non-blocking warning banner for modules with low non-zero coverage."""
    if not coverage_json:
        return no_update
    try:
        cov = json.loads(Path(coverage_json).read_text())
    except Exception:
        return no_update
    modules = cov.get("modules", cov) if isinstance(cov, dict) else cov
    low = []
    items = modules.items() if isinstance(modules, dict) else enumerate(modules)
    for key, entry in items:
        frac = (entry or {}).get("fraction_nonzero_cells")
        if frac is not None and frac < 0.1:
            name = entry.get("name", key) if isinstance(entry, dict) else key
            low.append(f"{name} ({frac:.1%})")
    if not low:
        return no_update
    return html.Div("Low coverage modules: " + ", ".join(low),
                    style={"color": "#a60", "marginTop": "6px"})
