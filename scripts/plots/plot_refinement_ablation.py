#!/usr/bin/env python3
"""Stage-B refinement is coverage/audit completion, not a quality knob.

Across the cohort's paired re-authoring seeds, Stage-B refinement moves the
blueprint's *epistemic* properties (specification coverage and provenance:
cited resources, DOIs, gene modules, axes, unique marker genes) substantially,
while leaving its *predictive* embedding quality (the scIB suite) essentially
unchanged. Two grounding modes are shown side by side: literature-grounded
(citations enforced) and resource-unconstrained (noref, no provenance layer).

Layout: a 2x2 grid of grouped bars. Rows are the grounding mode, columns are
the epistemic vs. the predictive metric block. Each bar is the cohort mean over
seeds with an s.d. whisker, Stage A (grey) against Stage A+B (red).

Inputs:
    outs/embeddings/<dataset>/blueprint_<cohort>_<dataset>_seed*_summary.json
    outs/metrics/<dataset>/metrics.csv                 (scIB, label_col=celltype)
Output:
    outs/plots/blueprint_variability/<dataset>/refinement_ablation.{png,pdf}

Usage:
    python scripts/plots/plot_refinement_ablation.py --dataset human_pancreas
"""
from __future__ import annotations

import argparse
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

# (display label, Stage-A cohort, Stage-A+B cohort) per grounding mode. The
# cohort token is the model-slug fragment in the stem
# (blueprint_<cohort>_<dataset>_seed*); paired seed N of each cohort is the
# Stage-B descendant of Stage-A seed N (true paired before/after).
GROUNDING_PAIRS = {
    "human_pancreas": [
        ("Literature-grounded", "gpt_5_5_unrefined_high", "gpt_5_5_refine_high"),
        ("Resource-unconstrained", "gpt_5_5_noref_unrefined_high", "gpt_5_5_noref_refine_high"),
    ],
}

# (label, summary.json blueprint key); counts grow under refinement.
EPISTEMIC = [("resources", "n_resources"), ("DOIs", "n_dois_unique"),
             ("modules", "n_modules"), ("axes", "n_axes"), ("genes", "n_unique_genes")]
# (label, metrics.csv column); scIB sub-scores, all in [0, 1].
PREDICTIVE = [("Total", "composite_scib"), ("bio", "bio_score"),
              ("batch", "batch_score_scib"), ("ARI", "ari"), ("NMI", "nmi")]


def epistemic_frame(dataset: str, cohort: str) -> pd.DataFrame:
    """Per-seed epistemic counts from each run's blueprint summary."""
    emb = REPO_ROOT / "outs" / "embeddings" / dataset
    rows = []
    for f in sorted(emb.glob(f"blueprint_{cohort}_{dataset}_seed*_summary.json")):
        bp = json.loads(f.read_text())["blueprint"]
        rows.append({label: bp[key] for label, key in EPISTEMIC})
    if not rows:
        raise SystemExit(f"no summaries for cohort {cohort!r} under {emb}")
    return pd.DataFrame(rows)


def predictive_frame(metrics: pd.DataFrame, cohort: str) -> pd.DataFrame:
    """Per-seed scIB sub-scores for one cohort (celltype labels)."""
    m = metrics[(metrics["model"] == f"blueprint_{cohort}")
                & (metrics["label_col"] == "celltype")]
    if m.empty:
        raise SystemExit(f"no metrics rows for cohort {cohort!r}")
    return m[[key for _, key in PREDICTIVE]].rename(
        columns={key: label for label, key in PREDICTIVE})


def grouped_bars(ax, labels, a_mean, a_sd, b_mean, b_sd, *, fmt, top_pad):
    """Stage A (grey) vs Stage A+B (red) grouped bars with s.d. whiskers and
    value labels; returns the y-max used so the caller can place a title."""
    x = np.arange(len(labels))
    w = 0.38
    bars = [(-w / 2, a_mean, a_sd, A_GREY, "black", "normal"),
            (+w / 2, b_mean, b_sd, BP_RED, BP_RED, "bold")]
    for dx, mean, sd, fc, tc, weight in bars:
        ax.bar(x + dx, mean, w, color=fc, edgecolor="black",
               yerr=sd, capsize=2, error_kw=dict(lw=0.8, ecolor="#444444"))
        for xi, mu, er in zip(x, mean, sd):
            ax.text(xi + dx, mu + er + top_pad, fmt(mu), ha="center", va="bottom",
                    fontsize=7.5, color=tc, fontweight=weight)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=8.5)
    ymax = max((np.array(a_mean) + np.array(a_sd)).max(),
               (np.array(b_mean) + np.array(b_sd)).max())
    ax.set_ylim(0, ymax * 1.22)
    ax.margins(x=0.04)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", default="human_pancreas",
                   choices=sorted(GROUNDING_PAIRS))
    args = p.parse_args()

    pairs = GROUNDING_PAIRS[args.dataset]
    metrics = pd.read_csv(REPO_ROOT / "outs" / "metrics" / args.dataset / "metrics.csv")

    fig, axes = plt.subplots(len(pairs), 2, figsize=(11.0, 2.7 * len(pairs) + 0.6),
                             gridspec_kw={"width_ratios": [1.05, 1.0]})
    axes = np.atleast_2d(axes)

    for r, (mode, a_cohort, b_cohort) in enumerate(pairs):
        ea, eb = epistemic_frame(args.dataset, a_cohort), epistemic_frame(args.dataset, b_cohort)
        pa, pb = predictive_frame(metrics, a_cohort), predictive_frame(metrics, b_cohort)
        print(f"[{mode}] Stage A n={len(ea)} | Stage A+B n={len(eb)}")

        e_labels = [l for l, _ in EPISTEMIC]
        grouped_bars(axes[r, 0], e_labels,
                     ea[e_labels].mean().values, ea[e_labels].std().values,
                     eb[e_labels].mean().values, eb[e_labels].std().values,
                     fmt=lambda v: f"{v:.0f}" if v >= 10 else f"{v:.1f}",
                     top_pad=eb[e_labels].mean().max() * 0.012)
        axes[r, 0].set_ylabel("count")

        p_labels = [l for l, _ in PREDICTIVE]
        grouped_bars(axes[r, 1], p_labels,
                     pa[p_labels].mean().values, pa[p_labels].std().values,
                     pb[p_labels].mean().values, pb[p_labels].std().values,
                     fmt=lambda v: f"{v:.2f}", top_pad=0.01)
        axes[r, 1].set_ylim(0, 1.0)
        axes[r, 1].set_ylabel("scIB score")

        axes[r, 0].set_title(f"{mode}: specification coverage ↑", fontsize=11, loc="left")
        axes[r, 1].set_title(f"{mode}: embedding quality (scIB)", fontsize=11, loc="left")

    # Legend in the epistemic panel's empty upper-left (its bars grow left→right,
    # so the scIB panels' tall bars stay clear of it).
    handles = [plt.Rectangle((0, 0), 1, 1, fc=A_GREY, ec="black"),
               plt.Rectangle((0, 0), 1, 1, fc=BP_RED, ec="black")]
    axes[0, 0].legend(handles, ["Stage A (unrefined)", "Stage A+B (refined)"],
                      fontsize=9, frameon=False, loc="upper left")
    fig.tight_layout()

    out_dir = REPO_ROOT / "outs" / "plots" / "blueprint_variability" / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / "refinement_ablation"
    fig.savefig(f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"[ok] wrote {stem}.png / .pdf")



if __name__ == "__main__":
    main()
