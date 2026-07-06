"""Render the blueprint gate-tree figure for download from the Explore tab.

Self-contained port of the pure tree logic in scripts/plots/blueprint_tree.py
(parentage from gate=module references, Reingold-Tilford layout) plus a renderer
that returns a matplotlib Figure with a mode-dependent x-axis label. The webapp
already reimplements gate-tree parentage in webapp/explorer.py; this module
follows the same "port compact, don't import scripts/plots" pattern so the
webapp stays standalone.

Input is the in-memory explore bundle (axisOrder + axisMeta), so no _report.md
on disk is needed. `detect_tree_kind` reads the blueprint's scenario notes to
pick the x-axis semantics (the Generate-mode toggle stamps the marker).
"""
from __future__ import annotations

import os
import re
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.lines as mlines
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt

from webapp.pipeline import TEMPORAL_NOTE

# x-axis semantics: what "tree depth" means for each tree kind.
TREE_KIND_XLABEL = {
    "celltype": "Tree depth (root = broad class → leaf = fine subtype)",
    "time": "Tree depth (root = earliest state → leaf = latest)",
}

# When several axes share the same positive module, prefer the one that most
# "owns" it: plain (state on its own) > gated (state within parent) > contrast.
_KIND_PRIORITY = {"plain": 0, "gated": 1, "contrast": 2}
_KIND_MARKER = {"plain": "o", "gated": "s", "contrast": "D"}


def detect_tree_kind(notes) -> str:
    """Pick the tree kind from a blueprint's scenario notes.

    The Generate-mode toggle injects TEMPORAL_NOTE (time tree) or CELLTYPE_NOTE
    (cell-type tree) into the notes. We only need to distinguish the time tree;
    everything else (cell-type marker, no marker, uploaded blueprint) defaults to
    "celltype" and the user can override the selector to "time"."""
    if notes and TEMPORAL_NOTE in notes:
        return "time"
    return "celltype"


def _parse_axis_name(name: str) -> tuple[str, int | None]:
    """Extract (chain_prefix, ordinal) from an axis name.

    Current convention is `<chainPrefix>_<descriptiveTail>` (no ordinal);
    parentage is carried by `gate=`, time progression by tree depth. Legacy
    `<chain>_T<n>_...` / `<chain>_<NN>_...` names still parse for archived runs.
    Returns (chain, None) when no numeric ordinal leads the second token."""
    m = re.match(r"^([A-Za-z][A-Za-z0-9]*)_[Tt]?(\d+)(?:_|$)", name)
    if m:
        return m.group(1), int(m.group(2))
    return name.split("_", 1)[0], None


def axes_from_bundle(bundle: dict) -> list[dict]:
    """Reshape the explore bundle's axisOrder + axisMeta into the axis dicts the
    tree builder expects: {name, kind, positive, gate, negatives, chain, ord}."""
    meta = bundle.get("axisMeta", {})
    order = bundle.get("axisOrder", []) or list(meta.keys())
    axes: list[dict] = []
    for name in order:
        m = meta.get(name, {})
        chain, ordn = _parse_axis_name(name)
        axes.append({
            "name": name,
            "kind": m.get("kind"),
            "positive": m.get("positive"),
            "gate": m.get("gate"),
            "negatives": list(m.get("negatives", [])),
            "chain": chain,
            "ord": ordn,
        })
    return axes


def build_edges(axes: list[dict], stitch: bool = True):
    """Return (edges, parent_of).

    edges: list of (parent, child, kind) with kind in {"gate", "stitch"}.
    parent_of: child -> parent (first edge in canonical order wins)."""
    by_positive: dict[str, str] = {}
    for ax in sorted(axes, key=lambda a: (_KIND_PRIORITY.get(a["kind"], 9),)):
        if ax["positive"] and ax["positive"] not in by_positive:
            by_positive[ax["positive"]] = ax["name"]

    edges: list[tuple[str, str, str]] = []
    parent_of: dict[str, str] = {}

    # 1) Gate edges (canonical): every gated axis's gate module -> its owner axis.
    for ax in axes:
        gate_mod = ax["gate"]
        if not gate_mod:
            continue
        parent_name = by_positive.get(gate_mod)
        if parent_name and parent_name != ax["name"]:
            edges.append((parent_name, ax["name"], "gate"))
            parent_of.setdefault(ax["name"], parent_name)

    # 2) Stitch edges (heuristic): same chain, predecessor-by-ordinal, no
    # gate-derived parent yet. A no-op under the current no-ordinal naming.
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
                    edges.append((pred["name"], ax["name"], "stitch"))
                    parent_of[ax["name"]] = pred["name"]

    return edges, parent_of


def roots_and_children(axes: list[dict], parent_of: dict[str, str]):
    children: dict[str, list[str]] = {ax["name"]: [] for ax in axes}
    for child, parent in parent_of.items():
        children[parent].append(child)
    roots = [ax["name"] for ax in axes if ax["name"] not in parent_of]
    return roots, children


def compute_positions(axes: list[dict], parent_of: dict[str, str]):
    """Tidy-tree layout: x = depth (root = earliest tier), y = Reingold-Tilford
    stacking (each leaf one slot, internal nodes centred on their children)."""
    roots, children = roots_and_children(axes, parent_of)

    depth: dict[str, int] = {}

    def _depth_of(name: str) -> int:
        if name in depth:
            return depth[name]
        depth[name] = 0 if name not in parent_of else _depth_of(parent_of[name]) + 1
        return depth[name]

    for ax in axes:
        _depth_of(ax["name"])

    by_name = {ax["name"]: ax for ax in axes}
    reg_index = {ax["name"]: i for i, ax in enumerate(axes)}

    def child_key(n: str) -> tuple[int, int]:
        ordn = by_name[n]["ord"]
        return (ordn if ordn is not None else 10 ** 6, reg_index[n])

    pos: dict[str, tuple[float, float]] = {}
    cursor = [0.0]

    def assign_y(name: str) -> float:
        kids = sorted(children.get(name, []), key=child_key)
        if not kids:
            y = cursor[0]
            cursor[0] += 1.0
        else:
            ys = [assign_y(k) for k in kids]
            y = sum(ys) / len(ys)
        pos[name] = (float(depth[name]), y)
        return y

    for root in sorted(roots, key=lambda n: reg_index[n]):
        assign_y(root)
        cursor[0] += 1.0  # gap between independent root subtrees
    return pos


def render_tree(axes, edges, positions, xlabel: str):
    """Render the gate tree and return the matplotlib Figure (no file writes).

    Ported from scripts/plots/blueprint_tree.render_figure with the x-axis label
    parameterized; the caller savefig()s the returned figure to a buffer."""
    chains: list[str] = []
    for ax in axes:
        if ax["chain"] not in chains:
            chains.append(ax["chain"])
    palette = (list(plt.get_cmap("tab20").colors)
               + list(plt.get_cmap("tab20b").colors)
               + list(plt.get_cmap("tab20c").colors))
    chain_color = {c: palette[i % len(palette)] for i, c in enumerate(chains)}

    xs = [p[0] for p in positions.values()]
    ys = [p[1] for p in positions.values()]
    xrange = max(xs) - min(xs) + 1
    yrange = max(ys) - min(ys) + 1
    fig, ax_fig = plt.subplots(
        figsize=(max(13.0, xrange * 3.4), max(5.5, yrange * 0.30)))

    for parent, child, kind in edges:
        if parent not in positions or child not in positions:
            continue
        x0, y0 = positions[parent]
        x1, y1 = positions[child]
        linestyle = "-" if kind == "gate" else "--"
        ax_fig.annotate(
            "", xy=(x1, y1), xytext=(x0, y0),
            arrowprops=dict(arrowstyle="-|>", linestyle=linestyle,
                            color="#666666", lw=1.0, shrinkA=12, shrinkB=12),
            zorder=2,
        )

    for ax_data in axes:
        name = ax_data["name"]
        if name not in positions:
            continue
        x, y = positions[name]
        marker = _KIND_MARKER.get(ax_data["kind"], "o")
        color = chain_color[ax_data["chain"]]
        ax_fig.scatter([x], [y], s=190, marker=marker, color=color,
                       edgecolor="black", linewidth=0.7, zorder=3)
        short = name
        prefix = f"{ax_data['chain']}_"
        if short.startswith(prefix):
            short = short[len(prefix):]
        ax_fig.text(x + 0.07, y, short, fontsize=10.5, ha="left", va="center",
                    color="#222222", zorder=4,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.75,
                              pad=0.4))

    ax_fig.invert_yaxis()
    ax_fig.set_xlabel(xlabel, fontsize=13)
    ax_fig.set_yticks([])
    xticks = sorted(set(int(p[0]) for p in positions.values()))
    ax_fig.set_xticks(xticks)
    ax_fig.set_xlim(min(xticks) - 0.3, max(xticks) + 1.6)
    ax_fig.tick_params(axis="x", labelsize=11)
    for spine in ("top", "right", "left"):
        ax_fig.spines[spine].set_visible(False)
    ax_fig.grid(axis="x", linestyle=":", color="#cccccc", alpha=0.6)

    present_kinds = {ax["kind"] for ax in axes}
    _kind_spec = [("plain", "o", 13), ("gated", "s", 13), ("contrast", "D", 12)]
    kind_handles = [
        mlines.Line2D([], [], marker=mk, color="w", markerfacecolor="white",
                      markeredgecolor="black", markersize=ms, label=lbl)
        for lbl, mk, ms in _kind_spec if lbl in present_kinds
    ]
    edge_handles = []
    if any(k == "gate" for _, _, k in edges):
        edge_handles.append(
            mlines.Line2D([], [], color="#666666", linestyle="-", label="gate edge"))
    if any(k == "stitch" for _, _, k in edges):
        edge_handles.append(
            mlines.Line2D([], [], color="#666666", linestyle="--", label="stitched edge"))
    chain_handles = [
        mpatches.Patch(facecolor=chain_color[c], edgecolor="black", label=c)
        for c in chains
    ]
    leg_kind = ax_fig.legend(
        handles=kind_handles + edge_handles,
        loc="upper left", bbox_to_anchor=(1.01, 1.0),
        fontsize=11, frameon=False,
        numpoints=1, handlelength=2.0, handletextpad=0.7, labelspacing=0.7,
    )
    ax_fig.add_artist(leg_kind)
    n_chain_cols = max(1, (len(chain_handles) + 14) // 15)
    ax_fig.legend(
        handles=chain_handles, title="lineage",
        loc="upper left", bbox_to_anchor=(1.01, 0.78),
        fontsize=10.5, title_fontsize=12, frameon=False,
        ncol=n_chain_cols, handlelength=1.6, handletextpad=0.6,
        columnspacing=1.4, labelspacing=0.7,
    )
    fig.tight_layout()
    return fig


def figure_from_bundle(bundle: dict, kind: str):
    """Build the gate-tree figure for a mapped blueprint's explore bundle."""
    axes = axes_from_bundle(bundle)
    if not axes:
        raise ValueError("bundle has no axes to render")
    edges, parent_of = build_edges(axes, stitch=True)
    positions = compute_positions(axes, parent_of)
    xlabel = TREE_KIND_XLABEL.get(kind, TREE_KIND_XLABEL["celltype"])
    return render_tree(axes, edges, positions, xlabel)
