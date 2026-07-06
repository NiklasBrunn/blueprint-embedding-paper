#!/usr/bin/env python3
"""Build a hierarchical tree from a populated blueprint's axis structure.

For one or more blueprint parquets under outs/embeddings/<dataset>/, reads the
`## Embedding dimensions` table from the run's <stem>_report.md companion,
recovers parent->child edges from the `gate=` module references, and
optionally adds "stitched" edges between same-chain plain axes whose ordinal
(T-tag or numeric infix) is one step apart. Renders:

    outs/plots/trees/<dataset>/<model>/<stem>_tree.txt   ASCII tree
    outs/plots/trees/<dataset>/<model>/<stem>_tree.png   matplotlib figure
    outs/plots/trees/<dataset>/<model>/<stem>_tree.pdf

Gate edges (solid arrows) are canonical: drawn whenever a gated axis has
`gate=<module>` and another axis owns that module as its positive program.
Stitched edges (dashed arrows, "(stitched)" tag in ASCII) bridge plain axes
when their chain prefix matches and their ordinal is consecutive — they
capture the developmental "next state" relationship the agent encoded in
axis names but not in `gate`. Disable with --no-stitch for a strictly
gate-faithful view.

Usage:
    python scripts/plots/blueprint_tree.py --dataset mouse_gastrulation
    python scripts/plots/blueprint_tree.py --dataset mouse_gastrulation \
        --input outs/embeddings/mouse_gastrulation/blueprint_<...>_seed1.parquet
    python scripts/plots/blueprint_tree.py --dataset mouse_gastrulation --no-stitch
"""
from __future__ import annotations

import argparse
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


REPO_ROOT = Path(__file__).resolve().parents[2]

DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation")
_TABLE_HEADER = "## Embedding dimensions"

# When several axes share the same positive module, prefer the one that most
# "owns" it: plain (state on its own) > gated (state within parent) > contrast.
_KIND_PRIORITY = {"plain": 0, "gated": 1, "contrast": 2}


# ============================================================ axis-table parsing


def _axes_from_summary(summary_path: Path) -> list[dict]:
    """Build the same axis dicts from a <stem>_summary.json when no report.md exists.

    The summary's blueprint.axes carry name/kind/positive/gate/negatives directly,
    so the tree can be recovered without the rendered markdown report (which is not
    shipped alongside the archived blueprints).
    """
    import json
    axes = json.loads(summary_path.read_text())["blueprint"]["axes"]
    rows: list[dict] = []
    for a in axes:
        chain, ord_n = _parse_axis_name(a["name"])
        rows.append({
            "name": a["name"], "kind": a.get("kind", ""),
            "positive": a.get("positive"), "gate": a.get("gate"),
            "negatives": list(a.get("negatives") or []),
            "chain": chain, "ord": ord_n,
        })
    return rows


def parse_axis_table(report_path: Path) -> list[dict]:
    """Parse '## Embedding dimensions' rows -> list of axis dicts.

    Each dict has: name, kind, positive, gate, negatives, chain, ord.
    The table is the same one umap_dimension_plots.py reads. Falls back to the
    sibling <stem>_summary.json when the report.md is absent (e.g. when embeddings
    are staged from the archived blueprints rather than freshly embedded).
    """
    if not report_path.exists():
        summary_path = report_path.with_name(
            report_path.name.replace("_report.md", "_summary.json"))
        if summary_path.exists():
            return _axes_from_summary(summary_path)
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

    Current convention is `<chainPrefix>_<descriptiveTail>` (no T-ordinal);
    parentage is carried by `gate=`, time progression by tree depth.
    Legacy formats `<chain>_T<n>_...` and `<chain>_<NN>_...` are still
    parsed for backwards compatibility with archived runs. Returns
    (chain, None) if no numeric ordinal is found in the name's first two
    underscore-separated tokens.
    """
    m = re.match(r"^([A-Za-z][A-Za-z0-9]*)_[Tt]?(\d+)(?:_|$)", name)
    if m:
        return m.group(1), int(m.group(2))
    return name.split("_", 1)[0], None


# ============================================================ tree edges


def build_edges(
    axes: list[dict], stitch: bool,
) -> tuple[list[tuple[str, str, str]], dict[str, str]]:
    """Return (edges, parent_of).

    edges: list of (parent_name, child_name, kind) with kind in {"gate", "stitch"}.
    parent_of: child_name -> parent_name (first edge in canonical order wins).
    """
    # Canonical "module -> axis" mapping: a module's owner is the highest-priority
    # axis that uses it as positive. Plain > gated > contrast; ties broken by
    # registration order.
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

    # 2) Stitch edges (heuristic): same chain, predecessor-by-ordinal, child has
    # no gate-derived parent yet. Captures plain-axis chains where the agent
    # encoded order in the name only (e.g. T0 -> T1 -> T2).
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


def roots_and_children(
    axes: list[dict], parent_of: dict[str, str],
) -> tuple[list[str], dict[str, list[str]]]:
    children: dict[str, list[str]] = {ax["name"]: [] for ax in axes}
    for child, parent in parent_of.items():
        children[parent].append(child)
    roots = [ax["name"] for ax in axes if ax["name"] not in parent_of]
    return roots, children


# ============================================================ ASCII rendering


def render_ascii_tree(
    axes: list[dict], edges: list[tuple[str, str, str]], parent_of: dict[str, str],
) -> str:
    by_name = {ax["name"]: ax for ax in axes}
    roots, children = roots_and_children(axes, parent_of)
    edge_kind = {(p, c): k for p, c, k in edges}

    def label(name: str) -> str:
        return f"{name}  [{by_name[name]['kind']}]"

    def child_sort_key(n: str) -> tuple[int, str]:
        ax = by_name[n]
        return (ax["ord"] if ax["ord"] is not None else 999, n)

    lines: list[str] = []

    def walk(name: str, prefix: str, is_last: bool, is_root: bool) -> None:
        if is_root:
            lines.append(label(name))
            new_prefix = ""
        else:
            parent = parent_of[name]
            connector = "└── " if is_last else "├── "
            tag = "" if edge_kind.get((parent, name)) == "gate" else "  (stitched)"
            lines.append(f"{prefix}{connector}{label(name)}{tag}")
            new_prefix = prefix + ("    " if is_last else "│   ")
        kids = sorted(children[name], key=child_sort_key)
        for i, kid in enumerate(kids):
            walk(kid, new_prefix, i == len(kids) - 1, is_root=False)

    for root in sorted(roots, key=child_sort_key):
        walk(root, "", True, is_root=True)
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


# ============================================================ figure layout


def compute_positions(
    axes: list[dict], parent_of: dict[str, str]
) -> dict[str, tuple[float, float]]:
    """Tidy-tree layout: x = tree depth (root = earliest = developmental time),
    y = a Reingold-Tilford-style stacking in which every leaf gets one slot and
    each internal node is centred on its children.

    This packs siblings compactly under their shared parent and removes the long
    edge crossings of a band-by-chain layout, so the gate tree reads as a clean
    dendrogram. Separate root subtrees are stacked with a one-slot gap.
    """
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
    cursor = [0.0]  # next free leaf slot

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


# ============================================================ figure rendering


_KIND_MARKER = {"plain": "o", "gated": "s", "contrast": "D"}


def render_figure(
    axes: list[dict],
    edges: list[tuple[str, str, str]],
    positions: dict[str, tuple[float, float]],
    title: str,
    out_png: Path,
    out_pdf: Path,
) -> None:
    chains: list[str] = []
    for ax in axes:
        if ax["chain"] not in chains:
            chains.append(ax["chain"])
    # Many lineages (~20+): concatenate qualitative palettes for distinct hues.
    palette = (list(plt.get_cmap("tab20").colors)
               + list(plt.get_cmap("tab20b").colors)
               + list(plt.get_cmap("tab20c").colors))
    chain_color = {c: palette[i % len(palette)] for i, c in enumerate(chains)}

    xs = [p[0] for p in positions.values()]
    ys = [p[1] for p in positions.values()]
    xrange = max(xs) - min(xs) + 1
    yrange = max(ys) - min(ys) + 1
    # Tidy layout: one slot per leaf in y, one depth tier per unit in x. Labels
    # sit to the right of each node, so width carries the depth tiers + label gutter.
    fig, ax_fig = plt.subplots(
        figsize=(max(13.0, xrange * 3.4), max(5.5, yrange * 0.30)))

    # Edges first so nodes render on top.
    for parent, child, kind in edges:
        if parent not in positions or child not in positions:
            continue
        x0, y0 = positions[parent]
        x1, y1 = positions[child]
        linestyle = "-" if kind == "gate" else "--"
        ax_fig.annotate(
            "", xy=(x1, y1), xytext=(x0, y0),
            arrowprops=dict(
                arrowstyle="-|>", linestyle=linestyle, color="#666666",
                lw=1.0, shrinkA=12, shrinkB=12,
            ),
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
        # Label to the right of the node (horizontal-dendrogram style), so labels
        # never collide vertically; the descriptive tail reads, color = lineage.
        short = name
        prefix = f"{ax_data['chain']}_"
        if short.startswith(prefix):
            short = short[len(prefix):]
        ax_fig.text(x + 0.07, y, short, fontsize=10.5, ha="left", va="center",
                    color="#222222", zorder=4,
                    bbox=dict(facecolor="white", edgecolor="none", alpha=0.75,
                              pad=0.4))

    ax_fig.invert_yaxis()
    ax_fig.set_xlabel("Tree depth (root = earliest state)", fontsize=13)
    ax_fig.set_yticks([])
    xticks = sorted(set(int(p[0]) for p in positions.values()))
    ax_fig.set_xticks(xticks)
    ax_fig.set_xlim(min(xticks) - 0.3, max(xticks) + 1.6)  # gutter for right labels
    ax_fig.tick_params(axis="x", labelsize=11)
    for spine in ("top", "right", "left"):
        ax_fig.spines[spine].set_visible(False)
    ax_fig.grid(axis="x", linestyle=":", color="#cccccc", alpha=0.6)

    # Only show legend entries for axis kinds / edge kinds actually present, so a
    # gate-only tree (no name-ordinal stitching, no contrast axes) does not carry a
    # dead "stitched edge" or "contrast" key.
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
    # Two stacked legends so the kind/edge markers don't get squeezed against
    # a long list of chain swatches, and the chain swatches can flow into
    # multiple columns instead of one tall stack.
    leg_kind = ax_fig.legend(
        handles=kind_handles + edge_handles,
        loc="upper left", bbox_to_anchor=(1.01, 1.0),
        fontsize=11, frameon=False,
        numpoints=1, handlelength=2.0, handletextpad=0.7, labelspacing=0.7,
    )
    ax_fig.add_artist(leg_kind)
    n_chain_cols = max(1, (len(chain_handles) + 14) // 15)  # ~15 rows per column
    ax_fig.legend(
        handles=chain_handles, title="lineage",
        loc="upper left", bbox_to_anchor=(1.01, 0.78),
        fontsize=10.5, title_fontsize=12, frameon=False,
        ncol=n_chain_cols,
        handlelength=1.6, handletextpad=0.6,
        columnspacing=1.4, labelspacing=0.7,
    )

    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    fig.savefig(out_pdf, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ============================================================ CLI


def discover_blueprints(embeddings_dir: Path, single: Path | None) -> list[Path]:
    if single is not None:
        # Tolerate a summary-only stem (no parquet) so the tree can be drawn from
        # staged blueprint metadata alone.
        if not single.exists() and not single.with_name(
                single.stem + "_summary.json").exists():
            raise SystemExit(f"--input not found: {single}")
        return [single]
    # Sweep stems from both parquets and staged summaries, so the tree renders whether
    # embeddings were freshly built or only the archived metadata was staged.
    stems = {p.name[: -len(".parquet")] for p in embeddings_dir.glob("blueprint_*.parquet")}
    stems |= {p.name[: -len("_summary.json")]
              for p in embeddings_dir.glob("blueprint_*_summary.json")}
    return [embeddings_dir / f"{stem}.parquet" for stem in sorted(stems)]


def model_slug_from_stem(stem: str, dataset: str) -> str:
    """`blueprint_<model>_<reasoning>_<dataset>_seed<N>` ->
    `blueprint_<model>_<reasoning>`."""
    suffix = f"_{dataset}_seed"
    if suffix in stem:
        return stem.split(suffix, 1)[0]
    return stem


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Hierarchical tree from blueprint axes.")
    p.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    p.add_argument("--input", type=Path, default=None,
                   help="Single parquet (default: sweep outs/embeddings/<dataset>/).")
    p.add_argument("--no-stitch", action="store_true",
                   help="Disable name-based stitching of consecutive plain axes.")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    stitch = not args.no_stitch
    embeddings_dir = REPO_ROOT / "outs" / "embeddings" / args.dataset
    out_root = REPO_ROOT / "outs" / "plots" / "trees" / args.dataset

    for parquet in discover_blueprints(embeddings_dir, args.input):
        stem = parquet.stem
        report_path = parquet.with_name(stem + "_report.md")
        summary_path = parquet.with_name(stem + "_summary.json")
        if not report_path.exists() and not summary_path.exists():
            print(f"[skip] {stem}: no report or summary companion")
            continue
        try:
            axes = parse_axis_table(report_path)
        except SystemExit as exc:
            print(f"[skip] {stem}: {exc}")
            continue
        if not axes:
            print(f"[skip] {stem}: empty axis table")
            continue

        edges, parent_of = build_edges(axes, stitch=stitch)
        model_dir = out_root / model_slug_from_stem(stem, args.dataset)
        out_txt = model_dir / f"{stem}_tree.txt"
        out_png = model_dir / f"{stem}_tree.png"
        out_pdf = model_dir / f"{stem}_tree.pdf"
        if out_png.exists() and not args.overwrite:
            print(f"[skip] {stem}: tree exists (use --overwrite)")
            continue
        model_dir.mkdir(parents=True, exist_ok=True)
        out_txt.write_text(render_ascii_tree(axes, edges, parent_of))
        positions = compute_positions(axes, parent_of)
        render_figure(axes, edges, positions, title=stem,
                      out_png=out_png, out_pdf=out_pdf)
        n_gate = sum(1 for _, _, k in edges if k == "gate")
        n_stitch = sum(1 for _, _, k in edges if k == "stitch")
        print(f"[ok] {stem}: {len(axes)} axes, "
              f"{n_gate} gate + {n_stitch} stitched edges -> {model_dir}")


if __name__ == "__main__":
    main()
