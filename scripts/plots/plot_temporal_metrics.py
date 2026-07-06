#!/usr/bin/env python3
"""Quantitative panels: the blueprint encodes developmental time in its
*axis topology* (gated-tree depth), not in the point-cloud geometry.

Each blueprint axis is a similarity-based readout (mean marker-gene expression),
so the embedding is expression-grounded like any other method; what differs is
that the agent arranges those axes into a time-forward gated tree whose depth
tracks embryonic stage. Shown two ways:

  (a) per-axis tree depth vs mean embryonic stage of the cells that axis owns
      (representative run), colored by lineage --- depth predicts time;
  (b) Spearman vs stage: the blueprint's *structural* depth ordering (by
      construction, no inference) against diffusion-pseudotime (DPT) run on each
      embedding's geometry (the route a black-box embedding needs to get time at
      all) --- the gated blueprint is deliberately a poor diffusion substrate
      because its temporal signal is in the tree, not the geometry.

Inputs (written by scripts/temporal_metrics/*.py):
    outs/temporal_metrics/<dataset>/temporal_metrics.csv
    outs/temporal_metrics/<dataset>/temporal_metrics_axes.csv
    outs/temporal_metrics/<dataset>/pseudotime_substrate.csv

Outputs:
    outs/plots/temporal_metrics/<dataset>/temporal_metrics.{png,pdf}

Usage:
    python scripts/plots/plot_temporal_metrics.py --dataset mouse_gastrulation
"""
from __future__ import annotations

import argparse
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]

HERO = "blueprint_gpt_5_5_refine_high"
BP_RED = "#d62728"
GREY = "#9e9e9e"

# pretty labels + bar colors for the substrate comparison
METHOD_LABEL = {
    "blueprint_gpt_5_5_refine_high": "Blueprint\n(refined)",
    "blueprint_gpt_5_5_unrefined_high": "Blueprint",
    "blueprint_gpt_5_5_noref_refine_high": "Blueprint\n(noref, ref.)",
    "blueprint_gpt_5_5_noref_unrefined_high": "Blueprint\n(noref)",
    "pca": "PCA",
    "scvi_nobatch": "scVI",
    "scvi_census": "scVI-census",
    "scvi": "scVI (batch)",
}
DPT_ORDER = ["pca", "scvi_nobatch", "scvi_census", "blueprint"]


def lineage_palette(chains: list[str]) -> dict[str, tuple]:
    palette = (list(plt.get_cmap("tab20").colors)
               + list(plt.get_cmap("tab20b").colors)
               + list(plt.get_cmap("tab20c").colors))
    return {c: palette[i % len(palette)] for i, c in enumerate(chains)}


def _hero_branches(axes_df: pd.DataFrame):
    """Reconstruct root→leaf branches of the representative run from the parent column."""
    from collections import defaultdict
    hero = axes_df[axes_df["model"] == HERO].copy()
    names = set(hero["axis"])
    depth_by = dict(zip(hero["axis"], hero["depth"]))
    stage_by = dict(zip(hero["axis"], hero["mean_stage"]))
    parent_by = {a: (p if isinstance(p, str) and p in names else "")
                 for a, p in zip(hero["axis"], hero["parent"])}
    children: dict[str, list[str]] = defaultdict(list)
    roots = []
    for a in hero["axis"]:
        if parent_by[a]:
            children[parent_by[a]].append(a)
        else:
            roots.append(a)
    paths: list[list[str]] = []

    def walk(n, acc):
        acc = acc + [n]
        ks = children.get(n, [])
        if not ks:
            paths.append(acc)
        for k in ks:
            walk(k, acc)

    for r in roots:
        walk(r, [])
    return paths, depth_by, stage_by


def panel_depth_stage(ax, axes_df: pd.DataFrame, summ: pd.DataFrame) -> None:
    """(a) per-lineage depth→stage trajectories (representative run): within every
    root→leaf branch, deeper axes own later-staged cells."""
    paths, depth_by, stage_by = _hero_branches(axes_df)
    chains = list(dict.fromkeys(axes_df["axis"].str.split("_").str[0]))
    cmap = lineage_palette(chains)
    for path in paths:
        pts = [(depth_by[a], stage_by[a]) for a in path
               if not np.isnan(stage_by[a])]
        if len(pts) < 2:
            continue
        xs, ys = zip(*pts)
        leaf_chain = path[-1].split("_")[0]
        ax.plot(xs, ys, "-o", color=cmap[leaf_chain], lw=1.2, ms=3.5,
                alpha=0.7, zorder=3)
    hero_row = summ[summ["model"] == HERO].iloc[0]
    ax.set_xlabel("Tree depth (root → leaf)")
    ax.set_ylabel("Mean embryonic stage of axis (E)")
    ax.set_xticks(sorted({int(d) for d in depth_by.values()}))
    # show the full developmental window (E6.5–E8.5), not just the range of means
    ax.set_ylim(6.45, 8.55)
    ax.set_yticks(np.arange(6.5, 8.51, 0.5))
    ax.annotate(f"per-lineage Spearman $\\bar\\rho$ = {hero_row['mean_branch_spearman']:.2f}\n"
                f"monotonic steps = {hero_row['monotonic_frac']:.0%}",
                xy=(0.04, 0.96), xycoords="axes fraction", va="top", fontsize=8.5,
                bbox=dict(boxstyle="round", fc="white", ec="#cccccc"))


def panel_structural_vs_dpt(ax, summ: pd.DataFrame, subs: pd.DataFrame,
                            refined_only: bool = False) -> None:
    """(b) blueprint structural depth ordering vs DPT-inferred ordering.

    refined_only=True uses only the refined run for the blueprint bars (no error
    bars); refined_only=False aggregates the blueprint over all re-authoring runs
    with s.d. error bars.
    """
    bp = summ[summ["model"] == HERO] if refined_only \
        else summ[summ["model"].str.startswith("blueprint")]
    struct_branch = bp["mean_branch_spearman"].mean()
    struct_branch_sd = 0.0 if refined_only else bp["mean_branch_spearman"].std(ddof=0)
    struct_pool = bp["pooled_depth_stage_spearman"].mean()
    struct_pool_sd = 0.0 if refined_only else bp["pooled_depth_stage_spearman"].std(ddof=0)

    # DPT global Spearman(pt, stage). Blueprint: refined run only (main) or
    # aggregated over runs with s.d.; baselines are single runs.
    def dpt_val(model_key):
        if model_key == "blueprint":
            sel = subs[subs["model"] == HERO] if refined_only \
                else subs[subs["model"].str.startswith("blueprint")]
            g = sel["spearman_stage"]
            return g.mean(), (0.0 if refined_only else g.std(ddof=0))
        g = subs[subs["model"] == model_key]["spearman_stage"]
        return (g.iloc[0], 0.0) if len(g) else (np.nan, 0.0)

    labels = ["Blueprint\nper-lineage\n(structural)", "Blueprint\nglobal\n(structural)"]
    vals = [struct_branch, struct_pool]
    errs = [struct_branch_sd, struct_pool_sd]
    colors = [BP_RED, "#f0a6a6"]
    for mk in DPT_ORDER:
        v, e = dpt_val(mk)
        labels.append(METHOD_LABEL.get(mk, mk).replace("\n", " ") + "\n+DPT")
        vals.append(v)
        errs.append(e)
        colors.append(BP_RED if mk == "blueprint" else GREY)

    x = np.arange(len(vals))
    ax.set_ylim(0, 1.22)
    ax.bar(x, vals, yerr=errs, color=colors, edgecolor="black", linewidth=0.6,
           capsize=3, width=0.7, zorder=3)
    ax.axvline(1.5, color="#cccccc", lw=1.0, ls=":")
    ax.text(0.5, 1.19, "structural (per-axis)\ndepth ↔ stage", ha="center",
            va="top", fontsize=8.5, style="italic", color=BP_RED, linespacing=0.95)
    ax.text(4.0, 1.19, "inferred (per-cell, +DPT)\npseudotime ↔ stage",
            ha="center", va="top", fontsize=8.5, style="italic", color="#555555",
            linespacing=0.95)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_ylabel("Spearman correlation with embryonic stage")
    ax.axhline(0, color="black", lw=0.6)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", default="mouse_gastrulation")
    args = p.parse_args()

    base = REPO_ROOT / "outs" / "temporal_metrics" / args.dataset
    summ = pd.read_csv(base / "temporal_metrics.csv")
    axes_df = pd.read_csv(base / "temporal_metrics_axes.csv")
    subs = pd.read_csv(base / "pseudotime_substrate.csv")

    out_dir = REPO_ROOT / "outs" / "plots" / "temporal_metrics" / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)

    # Refined-run panels, saved separately so each can be lettered on its own
    # (depth-stage trajectories + structural-vs-DPT bars).
    figc, axc = plt.subplots(figsize=(6.8, 4.6))
    panel_depth_stage(axc, axes_df, summ)
    figc.tight_layout()
    stemc = out_dir / "temporal_depthstage"
    figc.savefig(f"{stemc}.png", dpi=300, bbox_inches="tight")
    figc.savefig(f"{stemc}.pdf", bbox_inches="tight")
    plt.close(figc)

    figd, axd = plt.subplots(figsize=(6.8, 4.6))
    panel_structural_vs_dpt(axd, summ, subs, refined_only=True)
    figd.tight_layout()
    stemd = out_dir / "temporal_bars"
    figd.savefig(f"{stemd}.png", dpi=300, bbox_inches="tight")
    figd.savefig(f"{stemd}.pdf", bbox_inches="tight")
    plt.close(figd)
    print(f"[ok] wrote {stemc}.png and {stemd}.png  (refined run)")

    # Same contrast aggregated across all re-authoring runs (s.d. bars).
    figa, axa = plt.subplots(figsize=(6.4, 4.6))
    panel_structural_vs_dpt(axa, summ, subs, refined_only=False)
    figa.tight_layout()
    stema = out_dir / "temporal_allseeds"
    figa.savefig(f"{stema}.png", dpi=300, bbox_inches="tight")
    figa.savefig(f"{stema}.pdf", bbox_inches="tight")
    plt.close(figa)
    print(f"[ok] wrote {stema}.png / .pdf  (all runs)")


if __name__ == "__main__":
    main()
