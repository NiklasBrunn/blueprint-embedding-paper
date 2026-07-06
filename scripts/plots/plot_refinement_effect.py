#!/usr/bin/env python3
"""What the optional Stage-B refinement adds to a blueprint.

The refinement loop (discriminator → generator + regression guard) takes a
Stage-A blueprint and revises it. Because a refined run archives its pre-refine
Stage-A snapshot under the same run_id, the two pair exactly and show the
delta. On mouse_gastrulation the paired run (run03) is the representative one:
Stage A (unrefined seed3) → Stage A+B (refined seed3).

One panel: structural counts before vs after refinement (axes, modules,
resources, unique genes, gated axes). The named lineage states refinement
*added* (grouped by germ layer) and the "nothing removed" (regression guard)
and preserved-metric notes are written into the figure caption, not a panel.

Inputs:  outs/embeddings/<dataset>/blueprint_*_summary.json + *_report.md
         outs/temporal_metrics/<dataset>/temporal_metrics.csv   (fidelity)
         outs/metrics/<dataset>/metrics.csv                      (bio score)
Output:  outs/plots/temporal_metrics/<dataset>/refinement_effect.{png,pdf}

Usage:   python scripts/plots/plot_refinement_effect.py --dataset mouse_gastrulation
"""
from __future__ import annotations

import argparse
import importlib.util
import json
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
BP_RED = "#d62728"
A_GREY = "#9e9e9e"


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, REPO_ROOT / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


parse_axis_table = _load("bp_tree", "scripts/plots/blueprint_tree.py").parse_axis_table

# rough germ-layer grouping of chain prefixes, for the "added states" panel
GERM_LAYER = {
    "endo": "Endoderm", "gut": "Endoderm", "visc": "Endoderm", "parietal": "Endoderm",
    "meso": "Mesoderm", "cardiac": "Mesoderm", "hemo": "Mesoderm", "blood": "Mesoderm",
    "somite": "Mesoderm", "extra": "Mesoderm", "axial": "Mesoderm",
    "ecto": "Ectoderm", "neural": "Ectoderm", "crest": "Ectoderm", "surface": "Ectoderm",
    "pgc": "Germline", "streak": "Primitive streak", "caudal": "Mesoderm",
    "troph": "Extra-embryonic", "epi": "Epiblast", "cycle": "Cycling",
}


def find_pair(dataset: str):
    """Return (stageA_stem, stageAB_stem) sharing a run_id (refined run + its snapshot)."""
    emb = REPO_ROOT / "outs" / "embeddings" / dataset
    by_run = {}
    for f in emb.glob("blueprint_*_summary.json"):
        d = json.loads(f.read_text())
        stem = f.name.replace("_summary.json", "")
        kind = "AB" if "_refine_" in stem else "A"
        by_run.setdefault(d.get("run_id"), {})[kind] = (stem, d["blueprint"])
    for run_id, pair in by_run.items():
        if "A" in pair and "AB" in pair:
            return run_id, pair["A"], pair["AB"]
    raise SystemExit("No paired Stage-A / Stage-A+B run found (need a --with-refine run).")


def axis_names(stem: str, dataset: str) -> list[str]:
    rep = REPO_ROOT / "outs" / "embeddings" / dataset / f"{stem}_report.md"
    return [a["name"] for a in parse_axis_table(rep)]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", default="mouse_gastrulation")
    args = p.parse_args()

    run_id, (a_stem, a_bp), (b_stem, b_bp) = find_pair(args.dataset)
    print(f"[pair] run_id={run_id}\n  Stage A   : {a_stem}\n  Stage A+B : {b_stem}")

    metrics = [
        ("axes", a_bp["n_axes"], b_bp["n_axes"]),
        ("modules", a_bp["n_modules"], b_bp["n_modules"]),
        ("resources", a_bp["n_resources"], b_bp["n_resources"]),
        ("unique\ngenes", a_bp["n_unique_genes"], b_bp["n_unique_genes"]),
        ("gated\naxes", a_bp["n_axes_by_kind"]["gated"], b_bp["n_axes_by_kind"]["gated"]),
    ]
    a_axes, b_axes = axis_names(a_stem, args.dataset), axis_names(b_stem, args.dataset)
    added = [x for x in b_axes if x not in set(a_axes)]
    removed = [x for x in a_axes if x not in set(b_axes)]

    # fidelity + bio score before/after, matched to the *exact* paired run (same
    # seed) — not just the cohort, which would otherwise grab a different seed.
    def seed_of(stem):
        return int(stem.split("seed")[-1])

    def fidelity(stem):  # temporal_metrics.csv has a per-run `run` column == stem
        try:
            s = pd.read_csv(REPO_ROOT / "outs" / "temporal_metrics" / args.dataset
                            / "temporal_metrics.csv")
            return float(s[s["run"] == stem]["mean_branch_spearman"].iloc[0])
        except Exception:
            return float("nan")

    def bio(stem):  # metrics.csv keyed by (model, seed, label_col)
        try:
            model = stem.split(f"_{args.dataset}_seed")[0]
            m = pd.read_csv(REPO_ROOT / "outs" / "metrics" / args.dataset / "metrics.csv")
            m = m[(m["model"] == model) & (m["label_col"] == "celltype")
                  & (m["seed"].astype(float) == float(seed_of(stem)))]
            return float(m["bio_score"].iloc[0])
        except Exception:
            return float("nan")

    rho_a, rho_b = fidelity(a_stem), fidelity(b_stem)
    bio_a, bio_b = bio(a_stem), bio(b_stem)

    fig, ax0 = plt.subplots(figsize=(7.5, 4.4))

    # (a) grouped counts
    labels = [m[0] for m in metrics]
    av = [m[1] for m in metrics]
    bv = [m[2] for m in metrics]
    x = np.arange(len(labels))
    w = 0.38
    ax0.bar(x - w / 2, av, w, label="Stage A (unrefined)", color=A_GREY, edgecolor="black")
    ax0.bar(x + w / 2, bv, w, label="Stage A+B (refined)", color=BP_RED, edgecolor="black")
    for xi, (a_, b_) in zip(x, zip(av, bv)):
        ax0.text(xi - w / 2, a_ + 0.5, str(a_), ha="center", va="bottom", fontsize=8)
        ax0.text(xi + w / 2, b_ + 0.5, f"{b_}", ha="center", va="bottom", fontsize=8,
                 color=BP_RED, fontweight="bold")
    ax0.set_xticks(x)
    ax0.set_xticklabels(labels, fontsize=9)
    ax0.set_ylabel("count")
    ax0.legend(fontsize=8, frameon=False, loc="upper right")
    ax0.set_ylim(0, max(bv) * 1.18)

    # The list of added states and the preserved-metric numbers are reported in the
    # figure caption (text), not drawn as a second text panel.

    fig.tight_layout()
    out_dir = REPO_ROOT / "outs" / "plots" / "temporal_metrics" / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / "refinement_effect"
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] wrote {stem}.png / .pdf  (added {len(added)}, removed {len(removed)}; "
          f"rho {rho_a:.2f}->{rho_b:.2f}; bio {bio_a:.3f}->{bio_b:.3f}; "
          f"states: {', '.join(added)})")



if __name__ == "__main__":
    main()
