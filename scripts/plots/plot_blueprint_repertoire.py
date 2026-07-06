"""
Blueprint repertoire & robustness figures — what the agent reuses across runs.

The three-way summary and the reproducibility figure give the *scalar* view of
blueprint variability (Jaccards, CKA, scIB dispersion). This figure gives the
*qualitative* view that those scalars summarise: which genes, modules, axes and
resources actually recur across replicate runs, and how stable the axis-kind /
resource-kind composition is. Two variants of a single selected blueprint (the
unrefined GPT-5.5/high cohort):

  repertoire_noref       — the resource-free (noref) variant. No citations, so:
      (a) axis kinds across runs  [per-run stacked bars]
      (b) top-20 genes   (c) top-20 modules   (d) top-20 axes   [recurrence bars]

  repertoire_resourced   — the literature-grounded variant. Adds an
      auditable citation layer, so:
      (a) axis kinds across runs       [per-run stacked bars]
      (b) resource kinds across runs   [per-run stacked bars]
      (c) top-20 genes   (d) top-20 modules
      (e) top-20 axes    (f) top-20 resources

Recurrence is the number of runs a concept appears in (out of n_runs), matched on the
expert-curated standardized names from configs/name_mappings/<dataset>/ (synonymous raw
labels mapped to a common concept; blank entries fall back to the canonical lowercased
form). Each recurrence panel is annotated with the cohort's off-diagonal mean pairwise
Jaccard ± 95% CI for that dimension, computed on the same standardized keys.

Freeform has no structured genes/modules/axes/resources to count, so this audit is
blueprint-only by construction — that asymmetry is the point (see the Results text).

Inputs (per cohort, from the blueprint-variability stage):
    outs/blueprint_variability/<dataset>/<cohort>/run_summary.csv
    outs/blueprint_variability/<dataset>/<cohort>/{gene,module,axis,resource}_recurrence.csv
    outs/blueprint_variability/<dataset>/cohort_metrics.csv   (off-diag Jaccards)

Outputs:
    outs/plots/blueprint_repertoire/<dataset>/repertoire_{noref,resourced}.{png,pdf}

Usage:
    python scripts/plots/plot_blueprint_repertoire.py --dataset human_pancreas
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
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")

# Match the three-way / reproducibility palette so the figures read as a set.
BLUEPRINT_BAR = "#1f77b4"
AXIS_KIND_COLORS = {"plain": "#9ecae1", "contrast": "#6a51a3", "gated": "#08519c"}
RES_KIND_COLORS = {"paper": "#2171b5", "atlas": "#41ab5d", "database": "#fd8d3c"}
TOP_N = 20


def cohort_dir(dataset: str, cohort: str) -> Path:
    return REPO_ROOT / "outs" / "blueprint_variability" / dataset / cohort


def jaccards_for(dataset: str, cohort: str) -> dict[str, tuple[float, float, float]]:
    """dimension -> (mean, ci_lo, ci_hi) off-diag Jaccard, from cohort_metrics.csv."""
    path = REPO_ROOT / "outs" / "blueprint_variability" / dataset / "cohort_metrics.csv"
    row = pd.read_csv(path).set_index("cohort").loc[cohort]
    out: dict[str, tuple[float, float, float]] = {}
    for dim, key in [("genes", "genes"), ("modules", "modnames"),
                     ("axes", "axisnames"), ("resources", "resources")]:
        out[dim] = (float(row[f"jaccard_{key}_off_diag"]),
                    float(row[f"jaccard_{key}_ci_lo"]),
                    float(row[f"jaccard_{key}_ci_hi"]))
    return out


def _jac_tag(jac: tuple[float, float, float]) -> str:
    m, lo, hi = jac
    return f"J̄ = {m:.2f} [{lo:.2f}, {hi:.2f}]"


# ---- per-run kind composition (stacked bars) -------------------------------
def kind_strip(ax, run_df: pd.DataFrame, count_cols: dict[str, str],
               colors: dict[str, str], title: str) -> None:
    runs = run_df.sort_values("run_id").reset_index(drop=True)
    x = np.arange(len(runs))
    bottom = np.zeros(len(runs))
    for kind, col in count_cols.items():
        vals = runs[col].to_numpy(dtype=float)
        ax.bar(x, vals, bottom=bottom, color=colors[kind], edgecolor="white",
               linewidth=0.5, label=kind)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels([str(i + 1) for i in range(len(runs))], fontsize=9)
    ax.set_xlabel("run", fontsize=10)
    ax.set_ylabel("count", fontsize=10)
    ax.set_ylim(0, bottom.max() * 1.30)               # headroom so legend clears bars
    ax.set_title(title, fontsize=12, loc="left", pad=8)
    ax.legend(fontsize=8, ncol=len(count_cols), frameon=False,
              loc="upper right")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


# ---- top-N recurrence (horizontal bars) ------------------------------------
def recurrence_panel(ax, df: pd.DataFrame, name_col: str, title: str,
                     jac: tuple[float, float, float], n_runs: int, *,
                     colors=None, kind_col: str | None = None,
                     trunc: int = 26) -> None:
    if df.empty:
        ax.text(0.5, 0.5, "none", ha="center", va="center", fontsize=10,
                color="#888888", transform=ax.transAxes)
        ax.set_title(f"{title}   ({_jac_tag(jac)})", fontsize=12, loc="left", pad=8)
        ax.axis("off")
        return
    top = (df.sort_values(["n_runs", name_col], ascending=[False, True])
             .head(TOP_N).iloc[::-1].reset_index(drop=True))
    y = np.arange(len(top))
    if kind_col is not None and colors is not None:
        bar_colors = [colors.get(str(k), "#bdbdbd") for k in top[kind_col]]
    else:
        bar_colors = BLUEPRINT_BAR if colors is None else colors
    ax.barh(y, top["n_runs"], color=bar_colors, edgecolor="white", linewidth=0.5)
    labels = [str(s)[:trunc] + ("…" if len(str(s)) > trunc else "")
              for s in top[name_col]]
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_ylim(-0.6, len(top) - 0.4)
    ax.set_xlim(0, n_runs)
    ax.set_xticks(range(0, n_runs + 1, max(1, n_runs // 5)))
    ax.set_xlabel("number of runs", fontsize=10)
    ax.set_title(f"{title}   ({_jac_tag(jac)})", fontsize=12, loc="left", pad=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    if kind_col is not None and colors is not None:
        present = [k for k in colors if k in set(top[kind_col].astype(str))]
        handles = [plt.Rectangle((0, 0), 1, 1, color=colors[k]) for k in present]
        ax.legend(handles, present, fontsize=7, frameon=False, loc="lower right")


# ---- gene recurrence spectrum (genes are too many / too universal for top-N) -
def gene_spectrum_panel(ax, df: pd.DataFrame, title: str,
                        jac: tuple[float, float, float], n_runs: int) -> None:
    counts = (df["n_runs"].value_counts()
              .reindex(range(1, n_runs + 1), fill_value=0).sort_index())
    x = counts.index.to_numpy()
    h = counts.to_numpy()
    ax.bar(x, h, color=BLUEPRINT_BAR, edgecolor="white", linewidth=0.5)
    for xi, hi in zip(x, h):
        if hi:
            ax.text(xi, hi, str(int(hi)), ha="center", va="bottom", fontsize=8)
    ax.set_xticks(range(1, n_runs + 1))
    ax.set_xlim(0.4, n_runs + 0.6)
    ax.set_ylim(0, h.max() * 1.20)
    ax.set_xlabel("number of runs a gene appears in", fontsize=10)
    ax.set_ylabel("number of genes", fontsize=10)
    ax.set_title(f"{title}   ({_jac_tag(jac)})", fontsize=12, loc="left", pad=8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def _load(cohort_path: Path):
    run_df = pd.read_csv(cohort_path / "run_summary.csv")
    rec = {k: pd.read_csv(cohort_path / f"{k}_recurrence.csv")
           for k in ("gene", "module", "axis", "resource")}
    return run_df, rec


def fig_noref(dataset: str, cohort: str, out_dir: Path) -> None:
    cpath = cohort_dir(dataset, cohort)
    run_df, rec = _load(cpath)
    jac = jaccards_for(dataset, cohort)

    n_runs = len(run_df)
    fig = plt.figure(figsize=(17.5, 9.0))
    # wide wspace: the structure-derived axis labels (panel d) are long and sit
    # to the left of their axis, so they need clearance from panel (c)'s bars.
    gs = GridSpec(2, 3, height_ratios=[1.0, 1.9], hspace=0.32, wspace=0.85,
                  figure=fig)
    kind_strip(fig.add_subplot(gs[0, :]), run_df,
               {"plain": "n_axes_plain", "contrast": "n_axes_contrast",
                "gated": "n_axes_gated"}, AXIS_KIND_COLORS,
               "(a) Axis kinds across runs")
    gene_spectrum_panel(fig.add_subplot(gs[1, 0]), rec["gene"],
                        "(b) Gene recurrence spectrum", jac["genes"], n_runs)
    recurrence_panel(fig.add_subplot(gs[1, 1]), rec["module"], "module_name",
                     "(c) Top-20 modules", jac["modules"], n_runs)
    recurrence_panel(fig.add_subplot(gs[1, 2]), rec["axis"], "axis_name",
                     "(d) Top-20 axes", jac["axes"], n_runs)
    _save(fig, out_dir / "repertoire_noref")


def fig_resourced(dataset: str, cohort: str, out_dir: Path) -> None:
    cpath = cohort_dir(dataset, cohort)
    run_df, rec = _load(cpath)
    jac = jaccards_for(dataset, cohort)

    n_runs = len(run_df)
    fig = plt.figure(figsize=(16.5, 14.0))
    # wide wspace so the long axis / resource labels (left of their axes) clear
    # the bars of the panel to their left.
    gs = GridSpec(4, 2, height_ratios=[0.8, 0.8, 1.9, 1.9], hspace=0.42,
                  wspace=0.5, figure=fig)
    kind_strip(fig.add_subplot(gs[0, :]), run_df,
               {"plain": "n_axes_plain", "contrast": "n_axes_contrast",
                "gated": "n_axes_gated"}, AXIS_KIND_COLORS,
               "(a) Axis kinds across runs")
    kind_strip(fig.add_subplot(gs[1, :]), run_df,
               {"paper": "n_resources_paper", "atlas": "n_resources_atlas",
                "database": "n_resources_database"}, RES_KIND_COLORS,
               "(b) Resource kinds across runs")
    gene_spectrum_panel(fig.add_subplot(gs[2, 0]), rec["gene"],
                        "(c) Gene recurrence spectrum", jac["genes"], n_runs)
    recurrence_panel(fig.add_subplot(gs[2, 1]), rec["module"], "module_name",
                     "(d) Top-20 modules", jac["modules"], n_runs)
    recurrence_panel(fig.add_subplot(gs[3, 0]), rec["axis"], "axis_name",
                     "(e) Top-20 axes", jac["axes"], n_runs)
    recurrence_panel(fig.add_subplot(gs[3, 1]), rec["resource"], "citation",
                     "(f) Top-20 resources", jac["resources"], n_runs,
                     colors=RES_KIND_COLORS, kind_col="kind", trunc=30)
    _save(fig, out_dir / "repertoire_resourced")


def _save(fig, stem: Path) -> None:
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(stem.with_suffix(f".{ext}"), dpi=200, bbox_inches="tight")
        print(f"[ok] wrote {stem.with_suffix('.' + ext)}", flush=True)
    plt.close(fig)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Blueprint repertoire/robustness figures.")
    p.add_argument("--dataset", required=True, choices=DATASET_CHOICES)
    p.add_argument("--noref-cohort", default="gpt_5_5_noref_unrefined_high",
                   help="cohort stem for the resource-free variant (sans _<dataset>)")
    p.add_argument("--resourced-cohort", default="gpt_5_5_unrefined_high",
                   help="cohort stem for the literature-grounded variant (sans _<dataset>)")
    p.add_argument("--out-dir", type=Path,
                   default=REPO_ROOT / "outs" / "plots" / "blueprint_repertoire")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir / args.dataset
    noref = f"{args.noref_cohort}_{args.dataset}"
    resourced = f"{args.resourced_cohort}_{args.dataset}"
    fig_noref(args.dataset, noref, out_dir)
    fig_resourced(args.dataset, resourced, out_dir)


if __name__ == "__main__":
    main()
