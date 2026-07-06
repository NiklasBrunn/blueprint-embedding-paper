#!/usr/bin/env python3
"""Interactive HTML scrubber: walk along a developmental branch and watch the
UMAP coloring slide between per-axis activations via convex combination.

The HTML embeds the UMAP coords + per-axis activation arrays + every root-to-
leaf branch found in the blueprint's tree. A dropdown selects a branch, a
time slider scrubs t in [0, K-1] across the branch's K axes. At slider
value t with i = floor(t), alpha = t - i, every cell's color is the convex
combination:

    color(cell) = (1 - alpha) * values[axis_i] + alpha * values[axis_{i+1}]

so moving the slider visually morphs the activation pattern of axis_i into
that of axis_{i+1}. The plot library is Plotly (scattergl), loaded from CDN;
all interactivity is client-side JS embedded in the same file. No Python
plotly dependency.

Inputs (everything but the parquet is optional / auto-derived):
    outs/embeddings/<dataset>/<stem>.parquet                       REQUIRED  cell_id + dim_*
    outs/embeddings/<dataset>/<stem>_report.md                     REQUIRED  axis table (tree)
    outs/plots/umap/<dataset>/<model>/<stem>_umap_coords.parquet   2D layout cache — COMPUTED
        in-script from the parquet if absent (same UMAP recipe as umap_plots.py), and reused
        if present.
    data/preprocessed/<dataset>_pp.h5ad                            optional: obs columns for the
        "Metadata" coloring mode (best-effort; skipped if missing).
    outs/blueprints/<dataset>/<model>/embedder_blueprint_seed<N>.py  optional: module->gene chips
        in the sidebar (auto-discovered; skipped if missing).

Output:
    outs/plots/umap_branch_interactive/<dataset>/<model>/<stem>_interactive.html

Plotly: loaded from CDN by default (needs internet when the HTML is opened). Pass
--inline-plotly to embed plotly.min.js so the file is fully self-contained / offline-ready.

Usage:
    # From pipeline outputs + dataset alone; computes coords on the fly:
    python scripts/plots/umap_branch_interactive.py --dataset mouse_gastrulation \
        --input outs/embeddings/mouse_gastrulation/blueprint_<...>_seed1.parquet
    # all cells, no downsample, fully self-contained offline HTML:
    python scripts/plots/umap_branch_interactive.py --dataset mouse_gastrulation \
        --input <parquet> --max-cells -1 --inline-plotly
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import re
import sys
import tempfile
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]

# Reuse the static-UMAP palette builders so categorical colors match panel
# figures one-to-one. The import has to come after MPLCONFIGDIR is set because
# umap_plots itself pulls in matplotlib at import time.
os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import matplotlib.colors as mcolors  # noqa: E402
from umap_plots import (  # noqa: E402
    DEFAULT_PARAMS as UMAP_DEFAULT_PARAMS,
    build_categorical_colors,
    compute_or_load_umap,
    load_obs_columns,
)

DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation")
_TABLE_HEADER = "## Embedding dimensions"
_KIND_PRIORITY = {"plain": 0, "gated": 1, "contrast": 2}
_PLOTLY_CDN = "https://cdn.plot.ly/plotly-2.35.2.min.js"
_UNMATCHED_HEX = "#d0d0d0"


# ============================================================ axis-table parsing


def parse_axis_table(report_path: Path) -> list[dict]:
    if not report_path.exists():
        raise FileNotFoundError(report_path)
    lines = report_path.read_text().splitlines()
    try:
        start = next(i for i, ln in enumerate(lines) if ln.strip() == _TABLE_HEADER)
    except StopIteration as exc:
        raise SystemExit(f"No '{_TABLE_HEADER}' section in {report_path}") from exc
    rows: list[dict] = []
    seen_separator = False
    for ln in lines[start + 1:]:
        s = ln.strip()
        if not s.startswith("|"):
            if rows or seen_separator:
                break
            continue
        cells = [c.strip() for c in s.strip("|").split("|")]
        if not cells:
            continue
        if set(cells[0]) <= {"-", ":"}:
            seen_separator = True
            continue
        if cells[0].lower() == "axis":
            continue
        name = cells[0]
        kind = cells[1] if len(cells) > 1 else ""
        mods_str = cells[2] if len(cells) > 2 else ""
        mods = [m.strip() for m in mods_str.split(",") if m.strip()]
        positive = mods[0] if mods else None
        gate: str | None = None
        negatives: list[str] = []
        for m in mods[1:]:
            if m.startswith("gate="):
                gate = m[len("gate="):].strip()
            else:
                negatives.append(m)
        chain, ord_n = _parse_axis_name(name)
        rows.append({
            "name": name, "kind": kind,
            "positive": positive, "gate": gate, "negatives": negatives,
            "chain": chain, "ord": ord_n,
        })
    return rows


def _parse_axis_name(name: str) -> tuple[str, int | None]:
    """Extract (chain_prefix, ordinal) from an axis name.

    Current convention is `<chainPrefix>_<descriptiveTail>`; ordinal is None
    and parentage is carried by `gate=`. Legacy `<prefix>_T<n>_<tail>` and
    `<prefix>_<NN>_<tail>` formats are still parsed so archived runs keep
    plotting correctly (ordinal then drives same-prefix stitching when no
    gate is set).
    """
    m = re.match(r"^([A-Za-z][A-Za-z0-9]*)_[Tt]?(\d+)(?:_|$)", name)
    if m:
        return m.group(1), int(m.group(2))
    return name.split("_", 1)[0], None


def build_parent_of(axes: list[dict], stitch: bool) -> dict[str, str]:
    by_positive: dict[str, str] = {}
    for ax in sorted(axes, key=lambda a: _KIND_PRIORITY.get(a["kind"], 9)):
        if ax["positive"] and ax["positive"] not in by_positive:
            by_positive[ax["positive"]] = ax["name"]
    parent_of: dict[str, str] = {}
    for ax in axes:
        gate_mod = ax["gate"]
        if not gate_mod:
            continue
        parent = by_positive.get(gate_mod)
        if parent and parent != ax["name"]:
            parent_of.setdefault(ax["name"], parent)
    if stitch:
        by_chain: dict[str, list[dict]] = {}
        for ax in axes:
            if ax["ord"] is None:
                continue
            by_chain.setdefault(ax["chain"], []).append(ax)
        for members in by_chain.values():
            members.sort(key=lambda a: (a["ord"], a["name"]))
            for ax in members:
                if ax["name"] in parent_of:
                    continue
                pred = None
                for other in members:
                    if other["name"] == ax["name"] or other["ord"] is None:
                        continue
                    if other["ord"] < ax["ord"] and (pred is None or other["ord"] > pred["ord"]):
                        pred = other
                if pred is not None:
                    parent_of[ax["name"]] = pred["name"]
    return parent_of


def enumerate_branches(
    axes: list[dict], parent_of: dict[str, str]
) -> list[dict]:
    """Every leaf -> root path. Skips singletons (paths of length 1)."""
    by_name = {ax["name"]: ax for ax in axes}
    # leaf = node with no incoming child-of-itself in parent_of (i.e. not a parent of anything).
    parents = set(parent_of.values())
    leaves = [a["name"] for a in axes if a["name"] not in parents]
    branches: list[dict] = []
    for leaf in leaves:
        path = [leaf]
        cur = leaf
        while cur in parent_of:
            cur = parent_of[cur]
            if cur in path:
                break
            path.append(cur)
        path = list(reversed(path))
        if len(path) < 2:
            continue
        label = f"{path[0]}  →  {path[-1]}  ({len(path)} steps)"
        branches.append({
            "label": label, "leaf": leaf, "axes": path,
            "chains": [by_name[n]["chain"] for n in path],
        })
    branches.sort(key=lambda b: (-len(b["axes"]), b["label"]))
    return branches


# ============================================================ data loading


def load_coords_and_values(
    parquet: Path, dataset: str, axes: list[dict], umap_params: dict,
) -> tuple[pd.DataFrame, np.ndarray]:
    emb = pd.read_parquet(parquet)
    dim_cols = [c for c in emb.columns if c.startswith("dim_")]
    if len(dim_cols) != len(axes):
        raise SystemExit(
            f"axis table has {len(axes)} rows but parquet has {len(dim_cols)} dim_* columns"
        )
    model_slug = _model_slug_from_stem(parquet.stem, dataset)
    coords_path = (
        REPO_ROOT / "outs" / "plots" / "umap" / dataset / model_slug
        / f"{parquet.stem}_umap_coords.parquet"
    )
    # Reuse umap_plots.py's coords cache if it exists; otherwise compute the 2D
    # layout in-script via the *same* UMAP recipe and write the *same* cache
    # path, so (a) the script needs only the embedding parquet — no prior
    # umap_plots.py run — and (b) the coords match the static UMAP panels
    # one-to-one and are reused by a later umap_plots sweep. compute_or_load_umap
    # returns cell_id-sorted coords; the merge below re-aligns them to the
    # embedding's row order.
    if not coords_path.exists():
        print(f"[umap] coords cache miss -> computing in-script: {coords_path}")
    cell_ids_sorted, xy, _, cached = compute_or_load_umap(
        parquet, coords_path, umap_params, overwrite=False
    )
    if not cached:
        print(f"[umap] computed + cached UMAP coords for {len(cell_ids_sorted)} cells")
    coords = pd.DataFrame(
        {"cell_id": cell_ids_sorted.to_numpy(), "x": xy[:, 0], "y": xy[:, 1]}
    )
    merged = emb[["cell_id"]].merge(coords, on="cell_id", how="left")
    if merged[["x", "y"]].isna().any().any():
        raise SystemExit("Some cells lack UMAP coords after computing the layout.")
    values = emb[dim_cols].to_numpy(dtype=float)
    return merged, values


def _model_slug_from_stem(stem: str, dataset: str) -> str:
    suffix = f"_{dataset}_seed"
    if suffix in stem:
        return stem.split(suffix, 1)[0]
    return stem


# ============================================================ obs-metadata loading


def resolve_metadata_cols(dataset: str, override: str | None) -> list[str]:
    if override is not None:
        cols = [c.strip() for c in override.split(",") if c.strip()]
        return cols
    defaults = list(UMAP_DEFAULT_PARAMS[dataset]["color_keys"])
    # Skip composite "a+b" keys in v1; the UI mode list stays short.
    return [c for c in defaults if "+" not in c]


def load_metadata_payload(
    dataset: str,
    cols: list[str],
    cell_ids: pd.Series,
) -> dict[str, dict]:
    """Return {col: {kind, labels, colors, cellIdx}} aligned to cell_ids.

    Palettes come from umap_plots.build_categorical_colors so the interactive
    coloring matches the static UMAP panels. cellIdx is a per-cell int index
    into the labels list (with the sentinel len(labels) meaning unmatched).
    """
    if not cols:
        return {}
    h5ad_path = REPO_ROOT / "data" / "preprocessed" / f"{dataset}_pp.h5ad"
    obs_by_col = load_obs_columns(h5ad_path, cols)
    column_kinds = UMAP_DEFAULT_PARAMS[dataset].get("column_kinds", {})

    cell_id_idx = pd.Index(cell_ids.astype(str))
    payload: dict[str, dict] = {}
    for col in cols:
        series = obs_by_col[col]
        kind = column_kinds.get(col)
        # Align series values to coords ordering. Missing cell_ids → unmatched.
        aligned = series.reindex(cell_id_idx).fillna("unmatched").astype(str)
        color_map, _point_colors, ordered_labels = build_categorical_colors(
            aligned, dataset=dataset, column=col, kind=kind,
        )
        labels = list(ordered_labels)
        hex_colors = [mcolors.to_hex(color_map[lbl]) for lbl in labels]
        # Append unmatched bucket only if any cell's value isn't in color_map
        # (build_categorical_colors normally includes "unmatched" already when
        # the fillna sentinel appears, but guard against palette functions that
        # might drop it).
        label_to_idx = {lbl: i for i, lbl in enumerate(labels)}
        unmatched_idx = len(labels)
        needs_unmatched = False
        cell_idx: list[int] = []
        for v in aligned.tolist():
            i = label_to_idx.get(v, unmatched_idx)
            cell_idx.append(i)
            if i == unmatched_idx:
                needs_unmatched = True
        if needs_unmatched:
            labels.append("unmatched")
            hex_colors.append(_UNMATCHED_HEX)
        payload[col] = {
            "kind": kind or "label",
            "labels": labels,
            "colors": hex_colors,
            "cellIdx": cell_idx,
        }
    return payload


# ============================================================ blueprint discovery (for module → genes)


def find_blueprint_for_parquet(parquet: Path, dataset: str) -> Path | None:
    """Locate the populated blueprint script that produced this parquet.

    Canonical archive layout (written by Finalize of the orchestrators):
        outs/blueprints/<dataset>/<model_slug>_<reasoning>_<mode>/embedder_blueprint_seed<N>.py
    Falls back to globbing the run's workdir for legacy runs predating the
    archive convention.
    """
    sm = re.match(rf"^blueprint_(.+?)_{re.escape(dataset)}_seed(\d+)$", parquet.stem)
    if sm:
        model_dir, seed = sm.group(1), int(sm.group(2))
        archived = (REPO_ROOT / "outs" / "blueprints" / dataset / model_dir
                    / f"embedder_blueprint_seed{seed}.py")
        if archived.exists():
            return archived
    # Fallback: glob the workdir.
    m = re.match(r"^blueprint_(.*)_seed(\d+)$", parquet.stem)
    if not m:
        return None
    parent_name, seed = m.group(1), int(m.group(2))
    parent_dir = REPO_ROOT / "work" / parent_name
    if not parent_dir.is_dir():
        return None
    candidates = sorted(parent_dir.glob(f"*_blueprint_run{seed:02d}/embedder_blueprint.py"))
    return candidates[-1] if candidates else None


def load_modules_from_script(script: Path) -> dict[str, list[str]]:
    """Import a populated blueprint script and return module_name -> [genes].
    Same importlib trick as scripts/agentic_embedding/extract_run.py.

    Sets sys.dont_write_bytecode for the duration of the import so no
    __pycache__/ directories are left inside outs/blueprints/.
    """
    if not script.exists():
        return {}
    prev = sys.dont_write_bytecode
    sys.dont_write_bytecode = True
    try:
        spec = importlib.util.spec_from_file_location(
            f"_bp_{script.parent.name}_{script.stem}", script
        )
        if spec is None or spec.loader is None:
            return {}
        mod = importlib.util.module_from_spec(spec)
        sys.modules[spec.name] = mod
        spec.loader.exec_module(mod)
        bp = mod.build_blueprint("default")
        return {name: list(m.genes) for name, m in bp.modules.items()}
    finally:
        sys.dont_write_bytecode = prev


# ============================================================ HTML emission


# Round payload numbers to 3 decimal places to keep the embedded JSON compact
# without visibly degrading the activation gradient.
def _round_list(arr: np.ndarray, decimals: int = 3) -> list[float]:
    return np.round(np.asarray(arr, dtype=float), decimals).tolist()


_HTML_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8" />
<title>{title}</title>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, sans-serif; margin: 0; padding: 12px 16px; background: #fafafa; color: #222; }}
  h1 {{ font-size: 14px; margin: 0 0 6px 0; font-weight: 600; }}
  .controls {{ display: flex; gap: 18px; align-items: center; flex-wrap: wrap; margin-bottom: 8px; font-size: 12px; }}
  .controls label {{ display: flex; align-items: center; gap: 6px; }}
  .controls select, .controls input[type=range] {{ font-family: ui-monospace, monospace; font-size: 11px; }}
  #mode {{ min-width: 220px; padding: 3px 6px; }}
  #branch {{ min-width: 380px; padding: 3px 6px; }}
  #palette {{ padding: 3px 6px; }}
  #time {{ width: 380px; }}
  .control-group {{ display: flex; align-items: center; gap: 8px; }}
  .control-group.hidden {{ display: none; }}
  #ticks {{ font-family: ui-monospace, monospace; font-size: 11px; color: #555; margin-top: 4px; display: flex; flex-wrap: wrap; align-items: center; gap: 4px; }}
  #ticks.hidden {{ display: none; }}
  .tick {{ padding: 2px 8px; border-radius: 10px; background: #eeeeee; color: #444; cursor: pointer; transition: background-color 50ms linear, color 50ms linear; user-select: none; border: 1px solid transparent; }}
  .tick:hover {{ border-color: #888; }}
  .tick-arrow {{ color: #999; padding: 0 1px; }}
  #stateLine {{ font-family: ui-monospace, monospace; font-size: 11px; color: #333; min-width: 360px; }}
  #catLegend {{ display: none; flex-wrap: wrap; gap: 6px 10px; margin: 6px 0 4px 0; font-size: 11px; font-family: ui-monospace, monospace; }}
  #catLegend.visible {{ display: flex; }}
  .cat-item {{ display: inline-flex; align-items: center; gap: 5px; padding: 2px 8px; border-radius: 10px; background: #f4f4f4; }}
  .cat-swatch {{ width: 10px; height: 10px; border-radius: 50%; display: inline-block; border: 1px solid rgba(0,0,0,0.15); }}
  #catTitle {{ font-size: 11px; color: #555; margin: 8px 0 2px 0; font-family: ui-monospace, monospace; display: none; }}
  #catTitle.visible {{ display: block; }}
  #plot {{ width: 100%; height: calc(100vh - 170px); background: white; border: 1px solid #ddd; }}
  .note {{ font-size: 11px; color: #666; }}
  .app {{ display: flex; gap: 14px; align-items: stretch; }}
  .main {{ flex: 1; min-width: 0; }}
  #sidebar {{
    flex: 0 0 320px; background: white; border: 1px solid #ddd;
    padding: 12px 14px; max-height: calc(100vh - 170px);
    overflow-y: auto; font-size: 11.5px; line-height: 1.45;
  }}
  #sidebar.hidden {{ display: none; }}
  #sidebar h2 {{ font-size: 12.5px; margin: 0 0 4px 0; font-family: ui-monospace, monospace; word-break: break-all; }}
  #sidebar .kind-badge {{ display: inline-block; font-size: 10px; padding: 1px 6px; border-radius: 8px; background: #eee; color: #555; text-transform: uppercase; letter-spacing: 0.4px; margin-bottom: 12px; }}
  .picker-title {{ font-size: 10px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.6px; color: #666; margin: 0 0 6px 0; }}
  .picker {{ display: flex; flex-direction: column; gap: 2px; margin-bottom: 12px; }}
  .picker-item {{
    padding: 4px 8px; border-radius: 4px; cursor: pointer;
    font-family: ui-monospace, monospace; font-size: 11.5px; color: #222;
    border-left: 3px solid transparent; background: #fafafa;
    transition: background-color 60ms linear, border-color 60ms linear;
    word-break: break-all;
  }}
  .picker-item:hover {{ background: #eef3f9; }}
  .picker-item.selected {{ background: #e6f0fb; border-left-color: #0072B2; font-weight: 600; }}
  .picker-divider {{ border-top: 1px solid #ddd; margin: 8px 0 10px 0; }}
  .role {{ margin-bottom: 14px; }}
  .role-heading {{ font-size: 10px; font-weight: 700; text-transform: uppercase; letter-spacing: 0.6px; margin-bottom: 4px; }}
  .role-mod {{ font-size: 11.5px; font-family: ui-monospace, monospace; margin: 4px 0 5px 0; color: #222; }}
  .gene-chips {{ display: flex; flex-wrap: wrap; gap: 3px; margin-bottom: 8px; }}
  .gene-chip {{
    padding: 1px 6px; border-radius: 4px; background: #f4f4f4;
    font-family: ui-monospace, monospace; font-size: 10.5px; color: #222;
    border-left: 2px solid;
  }}
  .role.positive .role-heading {{ color: #0072B2; }}
  .role.positive .gene-chip    {{ border-left-color: #0072B2; }}
  .role.gate     .role-heading {{ color: #CC79A7; }}
  .role.gate     .gene-chip    {{ border-left-color: #CC79A7; }}
  .role.negative .role-heading {{ color: #E69F00; }}
  .role.negative .gene-chip    {{ border-left-color: #E69F00; }}
  .no-genes {{ font-size: 10.5px; color: #999; font-style: italic; padding: 2px 0; }}
  .sidebar-note {{ font-size: 10px; color: #999; margin-top: 18px; line-height: 1.4; }}
</style>
{plotly_load}
</head>
<body>
<h1>{title}</h1>
<div class="app">
  <aside id="sidebar"></aside>
  <main class="main">
    <div class="controls">
      <label>Color by: <select id="mode"></select></label>
      <div class="control-group" id="branchGroup">
        <label>Branch: <select id="branch"></select></label>
      </div>
      <div class="control-group" id="paletteGroup">
        <label>Colors: <select id="palette"><option value="viridis">Viridis (sequential)</option><option value="divbo">Blue / orange (diverging)</option></select></label>
      </div>
      <div class="control-group" id="timeGroup">
        <label>t: <input type="range" id="time" min="0" max="1" step="0.01" value="0" /></label>
      </div>
      <span id="stateLine"></span>
    </div>
    <div id="ticks"></div>
    <div id="catTitle"></div>
    <div id="catLegend"></div>
    <div id="plot"></div>
    <p class="note" id="footNote">Color = convex combination (1-&alpha;)&middot;axis<sub>i</sub> + &alpha;&middot;axis<sub>i+1</sub> for i = floor(t), &alpha; = t-i. Drag the time slider to scrub; sidebar snaps to the dominant axis.</p>
  </main>
</div>
<script>
const DATA = {data_json};

const plot = document.getElementById("plot");
const modeSel = document.getElementById("mode");
const branchSel = document.getElementById("branch");
const paletteSel = document.getElementById("palette");
const timeSlider = document.getElementById("time");
const stateLine = document.getElementById("stateLine");
const ticksDiv = document.getElementById("ticks");
const sidebar = document.getElementById("sidebar");
const branchGroup = document.getElementById("branchGroup");
const paletteGroup = document.getElementById("paletteGroup");
const timeGroup = document.getElementById("timeGroup");
const catLegend = document.getElementById("catLegend");
const catTitle = document.getElementById("catTitle");
const footNote = document.getElementById("footNote");

function escapeHtml(s) {{
  return String(s).replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}}
function geneChipsHtml(moduleName) {{
  const genes = DATA.modules[moduleName] || [];
  if (genes.length === 0) return '<div class="no-genes">(no gene list available)</div>';
  return '<div class="gene-chips">' +
    genes.map(g => '<span class="gene-chip">' + escapeHtml(g) + '</span>').join("") +
    '</div>';
}}
function roleBlockHtml(role, moduleName) {{
  if (!moduleName) return "";
  return '<div class="role ' + role + '">' +
    '<div class="role-heading">' + role + '</div>' +
    '<div class="role-mod">' + escapeHtml(moduleName) + '</div>' +
    geneChipsHtml(moduleName) +
    '</div>';
}}
function negativeBlockHtml(modules) {{
  if (!modules || modules.length === 0) return "";
  const items = modules.map(nm =>
    '<div class="role-mod">' + escapeHtml(nm) + '</div>' + geneChipsHtml(nm)
  ).join("");
  return '<div class="role negative">' +
    '<div class="role-heading">negative</div>' + items +
    '</div>';
}}
function axisInfoHtml(axisName) {{
  const meta = DATA.axisMeta[axisName] || {{}};
  const parts = [];
  parts.push('<h2>' + escapeHtml(axisName) + '</h2>');
  parts.push('<div class="kind-badge">' + escapeHtml(meta.kind || "axis") + '</div>');
  parts.push(roleBlockHtml("positive", meta.positive));
  if (meta.gate) parts.push(roleBlockHtml("gate", meta.gate));
  if (meta.negatives && meta.negatives.length) parts.push(negativeBlockHtml(meta.negatives));
  if (Object.keys(DATA.modules || {{}}).length === 0) {{
    parts.push('<div class="sidebar-note">Gene lists not embedded — re-run with --workdir &lt;path&gt; or place the populated <code>embedder_blueprint.py</code> where auto-discovery can find it.</div>');
  }}
  return parts.join("");
}}

function renderBranchSidebar(axisName) {{
  sidebar.innerHTML = axisInfoHtml(axisName);
}}

function renderAxisSidebar(selectedAxis) {{
  const items = DATA.axisOrder.map(axn =>
    '<div class="picker-item' + (axn === selectedAxis ? ' selected' : '') + '" data-axis="' + escapeHtml(axn) + '">' +
    escapeHtml(axn) + '</div>'
  ).join("");
  sidebar.innerHTML =
    '<div class="picker-title">Color by axis (click to switch)</div>' +
    '<div class="picker" id="axisPicker">' + items + '</div>' +
    '<div class="picker-divider"></div>' +
    axisInfoHtml(selectedAxis);
  sidebar.querySelectorAll("#axisPicker .picker-item").forEach(el => {{
    el.addEventListener("click", () => {{
      currentAxis = el.getAttribute("data-axis");
      renderAxisSidebar(currentAxis);
      recolorAxis();
      const sel = sidebar.querySelector('.picker-item[data-axis="' + currentAxis + '"]');
      if (sel) sel.scrollIntoView({{ block: "nearest" }});
    }});
  }});
}}

function renderMetaSidebar(selectedCol) {{
  const cols = Object.keys(DATA.metadata || {{}});
  if (cols.length === 0) {{
    sidebar.innerHTML = '<div class="sidebar-note">No metadata columns were embedded for this dataset. Re-run with <code>--metadata-cols col1,col2</code>.</div>';
    return;
  }}
  const items = cols.map(col => {{
    const meta = DATA.metadata[col];
    const sub = '<div style="font-size:10px;color:#777;margin-top:2px;">' +
      escapeHtml(meta.kind) + ' · ' + meta.labels.length + ' categories</div>';
    return '<div class="picker-item' + (col === selectedCol ? ' selected' : '') + '" data-col="' + escapeHtml(col) + '">' +
      escapeHtml(col) + sub + '</div>';
  }}).join("");
  sidebar.innerHTML =
    '<div class="picker-title">Color by metadata (click to switch)</div>' +
    '<div class="picker" id="metaPicker">' + items + '</div>';
  sidebar.querySelectorAll("#metaPicker .picker-item").forEach(el => {{
    el.addEventListener("click", () => {{
      const col = el.getAttribute("data-col");
      currentMetaCol = col;
      renderMetaSidebar(col);
      recolorMeta(col);
    }});
  }});
}}

const PALETTES = {{
  viridis: {{
    colorscale: "Viridis",
    cmin: -1.0, cmax: 2.5,
    barTitle: "activation (z)",
    tickAccent: "68, 1, 84",
  }},
  divbo: {{
    colorscale: [
      [0.0, "#2166ac"], [0.25, "#67a9cf"], [0.5, "#f7f7f7"],
      [0.75, "#fdb863"], [1.0, "#b35806"]
    ],
    cmin: -2.5, cmax: 2.5,
    barTitle: "activation (z, signed)",
    tickAccent: "179, 88, 6",
  }},
}};
let currentPaletteKey = "viridis";

function applyPalette(key) {{
  currentPaletteKey = key;
  const p = PALETTES[key];
  Plotly.restyle(plot, {{
    "marker.colorscale": [p.colorscale],
    "marker.cmin": [p.cmin],
    "marker.cmax": [p.cmax],
    "marker.colorbar.title.text": [p.barTitle],
  }}, [0]);
}}

// --- mode setup ----------------------------------------------------------
const MODES = [
  {{key: "branch", label: "Branch (scrubber)"}},
  {{key: "axis",   label: "Single axis activation"}},
  {{key: "meta",   label: "Metadata"}},
];
for (const m of MODES) {{
  const opt = document.createElement("option");
  opt.value = m.key; opt.textContent = m.label;
  modeSel.appendChild(opt);
}}

for (let i = 0; i < DATA.branches.length; i++) {{
  const b = DATA.branches[i];
  const opt = document.createElement("option");
  opt.value = String(i); opt.textContent = b.label;
  branchSel.appendChild(opt);
}}

let currentMode = "branch";
let currentBranch = DATA.branches[0] || null;
let currentAxis = DATA.axisOrder[0];
const _metaCols = Object.keys(DATA.metadata || {{}});
let currentMetaCol = _metaCols.length ? _metaCols[0] : null;

function shortName(axn) {{
  const m = axn.match(/^([A-Za-z][A-Za-z0-9]*)_(.*)$/);
  return m ? m[2] : axn;
}}

function highlightTicks(i, alpha) {{
  const accent = PALETTES[currentPaletteKey].tickAccent;
  const chips = ticksDiv.querySelectorAll(".tick");
  chips.forEach((chip, k) => {{
    let w = 0;
    if (k === i) w = 1 - alpha;
    else if (k === i + 1) w = alpha;
    chip.style.backgroundColor = w === 0
      ? "#eeeeee"
      : "rgba(" + accent + ", " + (0.20 + w * 0.80).toFixed(3) + ")";
    chip.style.color = w > 0.5 ? "white" : "#333";
    chip.style.fontWeight = w > 0.5 ? "600" : "400";
  }});
}}

function renderCategoricalLegend(col, meta) {{
  catTitle.textContent = col + "  (" + meta.kind + ", " + meta.labels.length + " categories)";
  catTitle.classList.add("visible");
  catLegend.innerHTML = "";
  meta.labels.forEach((lbl, i) => {{
    const item = document.createElement("span");
    item.className = "cat-item";
    item.innerHTML = '<span class="cat-swatch" style="background:' + meta.colors[i] + '"></span>' + escapeHtml(lbl);
    catLegend.appendChild(item);
  }});
  catLegend.classList.add("visible");
}}
function hideCategoricalLegend() {{
  catLegend.classList.remove("visible");
  catTitle.classList.remove("visible");
}}

function recolorBranch() {{
  if (!currentBranch) return;
  const t = parseFloat(timeSlider.value);
  const K = currentBranch.axes.length;
  const tClamp = Math.max(0, Math.min(K - 1, t));
  const i = Math.min(K - 2, Math.floor(tClamp));
  const alpha = tClamp - i;
  const lo = DATA.values[currentBranch.axes[i]];
  const hi = DATA.values[currentBranch.axes[i + 1]];
  const n = lo.length;
  const interp = new Float32Array(n);
  for (let k = 0; k < n; k++) interp[k] = (1 - alpha) * lo[k] + alpha * hi[k];
  Plotly.restyle(plot, {{ "marker.color": [Array.from(interp)] }}, [0]);
  highlightTicks(i, alpha);
  const dominant = alpha < 0.5 ? i : i + 1;
  renderBranchSidebar(currentBranch.axes[dominant]);
  stateLine.textContent = "t=" + t.toFixed(2) +
    "   axis_i = " + currentBranch.axes[i] +
    "   axis_{{i+1}} = " + currentBranch.axes[i + 1] +
    "   α = " + alpha.toFixed(2);
}}

function recolorAxis() {{
  const arr = DATA.values[currentAxis];
  if (!arr) return;
  Plotly.restyle(plot, {{ "marker.color": [arr] }}, [0]);
  stateLine.textContent = "axis = " + currentAxis;
}}

function recolorMeta(col) {{
  const meta = DATA.metadata[col];
  if (!meta) return;
  const idx = meta.cellIdx;
  const colors = meta.colors;
  const n = idx.length;
  const out = new Array(n);
  for (let k = 0; k < n; k++) out[k] = colors[idx[k]] || "#d0d0d0";
  Plotly.restyle(plot, {{ "marker.color": [out], "marker.showscale": [false] }}, [0]);
  renderCategoricalLegend(col, meta);
  stateLine.textContent = "metadata = " + col;
}}

function renderTicks() {{
  if (!currentBranch) return;
  const axes = currentBranch.axes;
  ticksDiv.innerHTML = "";
  axes.forEach((axn, k) => {{
    if (k > 0) {{
      const arrow = document.createElement("span");
      arrow.className = "tick-arrow"; arrow.textContent = "→";
      ticksDiv.appendChild(arrow);
    }}
    const chip = document.createElement("span");
    chip.className = "tick";
    chip.textContent = "t=" + k + " · " + shortName(axn);
    chip.title = axn;
    chip.addEventListener("click", () => {{
      timeSlider.value = String(k); recolorBranch();
    }});
    ticksDiv.appendChild(chip);
  }});
}}

function selectBranch(idx) {{
  currentBranch = DATA.branches[idx];
  timeSlider.min = "0";
  timeSlider.max = String(currentBranch.axes.length - 1);
  timeSlider.step = "0.01";
  timeSlider.value = "0";
  renderTicks();
  recolorBranch();
}}

function applyMode(modeKey) {{
  currentMode = modeKey;
  const isBranch = modeKey === "branch";
  const isAxis = modeKey === "axis";
  const isMeta = modeKey === "meta";
  branchGroup.classList.toggle("hidden", !isBranch);
  timeGroup.classList.toggle("hidden", !isBranch);
  ticksDiv.classList.toggle("hidden", !isBranch);
  paletteGroup.classList.toggle("hidden", isMeta);
  if (isMeta) {{
    footNote.textContent = "Discrete coloring from preprocessed h5ad obs columns. Click a column on the left to switch.";
    renderMetaSidebar(currentMetaCol);
    if (currentMetaCol) recolorMeta(currentMetaCol);
  }} else {{
    hideCategoricalLegend();
    Plotly.restyle(plot, {{ "marker.showscale": [true] }}, [0]);
    applyPalette(currentPaletteKey);
    if (isBranch) {{
      footNote.textContent = "Color = convex combination (1-α)·axis_i + α·axis_{{i+1}}. Drag the time slider; sidebar snaps to the dominant axis.";
      if (currentBranch === null && DATA.branches.length) selectBranch(0);
      else {{ renderTicks(); recolorBranch(); }}
    }} else {{
      footNote.textContent = "Single-axis view: click an axis on the left to switch.";
      renderAxisSidebar(currentAxis);
      recolorAxis();
    }}
  }}
}}

modeSel.addEventListener("change", (e) => applyMode(e.target.value));
branchSel.addEventListener("change", (e) => selectBranch(parseInt(e.target.value, 10)));
timeSlider.addEventListener("input", recolorBranch);
paletteSel.addEventListener("change", (e) => {{
  applyPalette(e.target.value);
  if (currentMode === "branch") recolorBranch();
  else if (currentMode === "axis") recolorAxis();
}});

Plotly.newPlot(plot, [{{
  type: "scattergl",
  mode: "markers",
  x: DATA.coords.x,
  y: DATA.coords.y,
  marker: {{
    color: new Array(DATA.coords.x.length).fill(0),
    colorscale: "Viridis",
    reversescale: false,
    cmin: -1.0,
    cmax: 2.5,
    size: 2.5,
    showscale: true,
    colorbar: {{ title: "activation (z)", thickness: 12, len: 0.6 }}
  }}
}}], {{
  hovermode: false,
  showlegend: false,
  margin: {{ l: 24, r: 24, t: 8, b: 24 }},
  xaxis: {{ title: "UMAP 1", visible: true, showgrid: false, zeroline: false, showticklabels: false }},
  yaxis: {{ title: "UMAP 2", visible: true, showgrid: false, zeroline: false, showticklabels: false, scaleanchor: "x" }},
  plot_bgcolor: "white", paper_bgcolor: "white"
}}, {{ responsive: true, displaylogo: false }}).then(() => {{
  applyPalette("viridis");
  if (DATA.branches.length) {{ selectBranch(0); applyMode("branch"); }}
  else {{ modeSel.value = "axis"; applyMode("axis"); }}
}});
</script>
</body>
</html>
"""


def render_html(
    parquet: Path,
    dataset: str,
    axes: list[dict],
    branches: list[dict],
    coords: pd.DataFrame,
    values: np.ndarray,
    modules: dict[str, list[str]],
    metadata: dict[str, dict],
    inline_plotly: bool,
    out_html: Path,
) -> None:
    # Per-axis activation arrays keyed by axis name. Same row order as coords.
    # Embed *every* axis (not just those in branches) so the "Single axis" UI
    # mode can pick singletons and contrast axes too.
    values_payload: dict[str, list[float]] = {}
    for i, ax in enumerate(axes):
        values_payload[ax["name"]] = _round_list(values[:, i], decimals=3)

    branches_payload = [
        {"label": br["label"], "leaf": br["leaf"], "axes": br["axes"]}
        for br in branches
    ]
    # Per-axis role meta + module->genes table, for the left sidebar.
    axis_meta_payload = {
        ax["name"]: {
            "kind": ax["kind"],
            "positive": ax["positive"],
            "gate": ax["gate"],
            "negatives": ax["negatives"],
        }
        for ax in axes
    }
    data = {
        "coords": {
            "x": _round_list(coords["x"].to_numpy(), decimals=3),
            "y": _round_list(coords["y"].to_numpy(), decimals=3),
        },
        "values": values_payload,
        "branches": branches_payload,
        "axisMeta": axis_meta_payload,
        "modules": modules,
        "metadata": metadata,
        "axisOrder": [ax["name"] for ax in axes],
        "stem": parquet.stem,
    }

    if inline_plotly:
        plotly_js = urllib.request.urlopen(_PLOTLY_CDN).read().decode("utf-8")
        plotly_load = f"<script>\n{plotly_js}\n</script>"
    else:
        plotly_load = f'<script src="{_PLOTLY_CDN}"></script>'

    title = f"{parquet.stem} — interactive branch scrubber"
    html = _HTML_TEMPLATE.format(
        title=title,
        plotly_load=plotly_load,
        data_json=json.dumps(data, separators=(",", ":")),
    )
    out_html.parent.mkdir(parents=True, exist_ok=True)
    out_html.write_text(html)


# ============================================================ CLI


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Interactive branch scrubber on UMAP (HTML).")
    p.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--max-cells", type=int, default=25_000,
                   help="Random-subset to this many cells (seeded). -1 = use all.")
    p.add_argument("--seed", type=int, default=0)
    # UMAP recipe used only when no coords cache exists yet (defaults match
    # scripts/plots/umap_plots.py so in-script coords equal the static panels).
    p.add_argument("--n-neighbors", type=int, default=30,
                   help="UMAP n_neighbors for in-script coord computation (ignored if cache exists).")
    p.add_argument("--min-dist", type=float, default=0.15, help="UMAP min_dist (in-script compute).")
    p.add_argument("--metric", type=str, default="euclidean", help="UMAP metric (in-script compute).")
    p.add_argument("--random-state", type=int, default=42, help="UMAP random_state (in-script compute).")
    p.add_argument("--no-stitch", action="store_true")
    p.add_argument("--blueprint-script", type=Path, default=None,
                   help="Populated embedder_blueprint.py whose modules define the genes. "
                        "Default: auto-discover via outs/blueprints/<dataset>/<stem>.py, "
                        "falling back to work/<parent>/*_blueprint_run<NN>/.")
    p.add_argument("--metadata-cols", type=str, default=None,
                   help="Comma-separated obs columns to embed as discrete coloring options. "
                        "Defaults to umap_plots.DEFAULT_PARAMS[dataset]['color_keys'] "
                        "minus composite 'a+b' keys.")
    p.add_argument("--inline-plotly", action="store_true",
                   help="Embed plotly.min.js inline (offline-capable, +~3 MB HTML).")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    parquet: Path = args.input
    if not parquet.exists():
        raise SystemExit(f"--input not found: {parquet}")
    report_path = parquet.with_name(parquet.stem + "_report.md")
    axes = parse_axis_table(report_path)
    for idx, ax in enumerate(axes):
        ax["dim_index"] = idx
    parent_of = build_parent_of(axes, stitch=not args.no_stitch)
    branches = enumerate_branches(axes, parent_of)
    if not branches:
        print("[warn] no multi-axis branches found; 'Branch' mode will be unavailable.")
    umap_params = {
        "n_neighbors": args.n_neighbors, "min_dist": args.min_dist,
        "metric": args.metric, "random_state": args.random_state,
    }
    coords, values = load_coords_and_values(parquet, args.dataset, axes, umap_params)

    if args.max_cells > 0 and len(coords) > args.max_cells:
        n_orig = len(coords)
        rng = np.random.default_rng(args.seed)
        idx = np.sort(rng.choice(n_orig, size=args.max_cells, replace=False))
        coords = coords.iloc[idx].reset_index(drop=True)
        values = values[idx]
        print(f"[downsample] {args.max_cells} cells (from {n_orig} originals, seed={args.seed})")

    metadata_cols = resolve_metadata_cols(args.dataset, args.metadata_cols)
    metadata: dict[str, dict] = {}
    if metadata_cols:
        try:
            metadata = load_metadata_payload(args.dataset, metadata_cols, coords["cell_id"])
            print(f"[metadata] embedded {len(metadata)} columns: {list(metadata)}")
        except Exception as exc:
            print(f"[metadata] failed to load {metadata_cols}: {exc!r}")

    # Module -> genes table for the left sidebar. Canonical source =
    # outs/blueprints/<dataset>/<stem>.py written by Finalize; fallback =
    # the run's workdir embedder_blueprint.py.
    blueprint_script = args.blueprint_script or find_blueprint_for_parquet(parquet, args.dataset)
    modules: dict[str, list[str]] = {}
    if blueprint_script and blueprint_script.exists():
        try:
            modules = load_modules_from_script(blueprint_script)
            print(f"[modules] loaded {len(modules)} modules from {blueprint_script}")
        except Exception as exc:  # pragma: no cover (best-effort sidebar enrichment)
            print(f"[modules] failed to load from {blueprint_script}: {exc!r}")
    if not modules:
        print("[modules] no blueprint script located; sidebar will show module names without genes")

    out_dir = (
        REPO_ROOT / "outs" / "plots" / "umap_branch_interactive" / args.dataset
        / _model_slug_from_stem(parquet.stem, args.dataset)
    )
    out_html = out_dir / f"{parquet.stem}_interactive.html"
    if out_html.exists() and not args.overwrite:
        raise SystemExit(f"already exists (pass --overwrite): {out_html}")

    render_html(
        parquet=parquet, dataset=args.dataset,
        axes=axes, branches=branches,
        coords=coords, values=values, modules=modules,
        metadata=metadata,
        inline_plotly=args.inline_plotly, out_html=out_html,
    )
    size_mb = out_html.stat().st_size / 1e6
    print(f"[ok] {len(branches)} branches, {len(coords)} cells -> {out_html}  ({size_mb:.1f} MB)")
    for br in branches[:8]:
        print(f"  - {br['label']}")
    if len(branches) > 8:
        print(f"  ... and {len(branches) - 8} more")


if __name__ == "__main__":
    main()
