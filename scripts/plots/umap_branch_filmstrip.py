#!/usr/bin/env python3
"""Static branch filmstrip: a row of UMAP panels along one root→leaf lineage of
the gated blueprint tree, each coloured by that axis's activation.

This is the static, print-ready counterpart of the interactive branch explorer
(umap_branch_interactive.py). Reading left→right, it walks down one developmental
branch (e.g. epiblast → primitive streak → nascent mesoderm → haemato-endothelial
→ blood progenitor → erythroid) and shows how the axis activation sweeps across
the embedding as tree depth (= developmental time) increases. The point cloud is
the same expression-similarity UMAP throughout; only which axis is highlighted
changes — i.e. the temporal pattern is carried by the axis structure, not by a
separate trajectory layout.

Reuses _draw_dim / style_axes from umap_dimension_plots.py and the tree parser
from blueprint_tree.py. The 2D layout is the cached coords parquet that
umap_dimension_plots.py / umap_plots.py already write.

Inputs:
    outs/embeddings/<dataset>/<stem>.parquet                       (cell_id + dim_*)
    outs/embeddings/<dataset>/<stem>_report.md                     (axis table)
    outs/plots/umap/<dataset>/<model>/<stem>_umap_coords.parquet   (cached layout)

Output:
    outs/plots/umap_branch_filmstrip/<dataset>/<model>/<stem>_<leaf>_filmstrip.{png,pdf}

Usage:
    python scripts/plots/umap_branch_filmstrip.py --dataset mouse_gastrulation \
        --input outs/embeddings/mouse_gastrulation/blueprint_<...>_seed3.parquet
    # choose the branch by its leaf axis name:
    python scripts/plots/umap_branch_filmstrip.py --dataset mouse_gastrulation \
        --input <parquet> --leaf blood_erythroid
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import os
import tempfile
import textwrap
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
POINT_SIZE = {"mouse_gastrulation": 1.5}


def _load(name: str, relpath: str):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / relpath)
    mod = importlib.util.module_from_spec(spec)
    assert spec and spec.loader
    spec.loader.exec_module(mod)
    return mod


_UDP = _load("umap_dim", "scripts/plots/umap_dimension_plots.py")
_BT = _load("bp_tree", "scripts/plots/blueprint_tree.py")


def root_to_leaf_paths(axes, parent_of):
    _, children = _BT.roots_and_children(axes, parent_of)
    roots = [a["name"] for a in axes if a["name"] not in parent_of]
    paths = []

    def walk(n, acc):
        acc = acc + [n]
        ks = children.get(n, [])
        if not ks:
            paths.append(acc)
        for k in ks:
            walk(k, acc)

    for r in roots:
        walk(r, [])
    return paths


def pick_branch(paths, leaf: str | None) -> list[str]:
    if leaf:
        for p in paths:
            if p[-1] == leaf:
                return p
        raise SystemExit(f"--leaf {leaf} is not a leaf; leaves: "
                         + ", ".join(sorted(p[-1] for p in paths)))
    # default: the longest branch (deepest lineage), prefer a blood/erythroid leaf
    paths_sorted = sorted(paths, key=len, reverse=True)
    for p in paths_sorted:
        if "erythroid" in p[-1].lower():
            return p
    return paths_sorted[0]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", default="mouse_gastrulation")
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--leaf", type=str, default=None,
                   help="Leaf axis name selecting the branch (default: deepest, prefer blood).")
    p.add_argument("--n-genes", type=int, default=6,
                   help="How many marker genes of each axis's program to print under its panel.")
    args = p.parse_args()

    stem = args.input.stem
    model = _BT.model_slug_from_stem(stem, args.dataset)
    report = args.input.with_name(stem + "_report.md")
    axes = _BT.parse_axis_table(report)
    _, parent_of = _BT.build_edges(axes, stitch=True)
    branch = pick_branch(root_to_leaf_paths(axes, parent_of), args.leaf)
    axis_index = {ax["name"]: i for i, ax in enumerate(axes)}

    # literature-cited marker genes per program, from the run's summary sidecar,
    # so each panel shows the gene set the agent chose for that developmental state.
    summary_path = args.input.with_name(stem + "_summary.json")
    module_genes: dict[str, list[str]] = {}
    has_neg_by_name: dict[str, bool] = {}
    if summary_path.exists():
        bp = json.loads(summary_path.read_text()).get("blueprint", {})
        module_genes = {m["name"]: m.get("genes", []) for m in bp.get("modules", [])}
        # Axes with no negative modules (plain roots / gate-only) have no signed
        # pole -> one-sided colormap; sibling-contrasted gated axes stay diverging.
        has_neg_by_name = {a["name"]: len(a.get("negatives", [])) > 0
                           for a in bp.get("axes", [])}
    positive_of = {ax["name"]: ax["positive"] for ax in axes}

    coords_path = (REPO_ROOT / "outs" / "plots" / "umap" / args.dataset / model
                   / f"{stem}_umap_coords.parquet")
    if not coords_path.exists():
        raise SystemExit(f"Cached UMAP coords not found: {coords_path}\n"
                         f"Run umap_plots.py / umap_dimension_plots.py first.")

    emb = pd.read_parquet(args.input)
    coords_df = pd.read_parquet(coords_path)
    merged = emb.merge(coords_df, on="cell_id", how="inner")
    coords = merged[["umap_0", "umap_1"]].to_numpy()

    params = {"cmap": "RdBu_r", "seq_cmap": _UDP._upper_half_cmap("RdBu_r"),
              "point_size": POINT_SIZE.get(args.dataset, 2.0), "alpha": 0.6, "vmax": None}

    k = len(branch)
    fig, axs = plt.subplots(1, k, figsize=(3.0 * k, 3.5))
    axs = np.atleast_1d(axs)
    for j, name in enumerate(branch):
        values = merged[f"dim_{axis_index[name]}"].to_numpy()
        sc = _UDP._draw_dim(axs[j], coords, values, params, has_neg_by_name.get(name, True))
        short = name.split("_", 1)[1] if "_" in name else name
        _UDP.style_axes(axs[j], short, title_fontsize=9.5)
        # marker genes of this state's program, italicised under the panel
        genes = module_genes.get(positive_of.get(name, ""), [])[:args.n_genes]
        gene_txt = textwrap.fill(", ".join(genes), width=22) if genes else ""
        axs[j].set_xlabel(gene_txt, fontsize=8.5, style="italic", labelpad=4)
        axs[j].set_ylabel("")
    # a single shared "developmental time" arrow under the row
    fig.subplots_adjust(bottom=0.22)
    fig.text(0.5, 0.02, "developmental time  (tree depth, root → leaf)  →",
             ha="center", fontsize=10)
    cbar = fig.colorbar(sc, ax=axs.tolist(), fraction=0.012, pad=0.01, shrink=0.7)
    cbar.set_label("axis activation", fontsize=9)
    cbar.ax.tick_params(labelsize=8)

    out_dir = (REPO_ROOT / "outs" / "plots" / "umap_branch_filmstrip"
               / args.dataset / model)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_stem = out_dir / f"{stem}_{branch[-1]}_filmstrip"
    fig.savefig(f"{out_stem}.png", dpi=300, bbox_inches="tight", facecolor="white")
    fig.savefig(f"{out_stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"[ok] branch: {' → '.join(branch)}")
    print(f"[ok] wrote {out_stem}.png / .pdf")



if __name__ == "__main__":
    main()
