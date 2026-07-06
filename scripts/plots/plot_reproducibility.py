"""
Reproducibility companion figure to the blueprint-vs-freeform three-way summary.

The three-way summary (plot_three_way_summary.py) covers quality + structure and the
spec-level gene Jaccard (Level 1). This figure adds what it lacks: the *embedding*
and *conclusion* levels of reproducibility, and the conventional-method yardstick
(scVI, the only baseline with seed replicates). It reads the outputs of
scripts/embedding_reproducibility/metrics.py plus the spec gene-Jaccard from the
(corrected) blueprint-variability cohort metrics.

Panels:
  (a) "Specification churns, the embedding does not" — blueprint cohorts only.
      Per cohort: spec gene-Jaccard (Level 1) vs representation CKA (Level 2). The
      point: gene lists overlap only ~0.5, yet the
      embedding geometry is highly reproducible — spec variation does not propagate.
  (b) "Reproducibility yardstick" — all multi-seed cohorts, ordered
      scVI · scVI-nobatch · blueprint · freeform. Representation CKA (with run-level
      bootstrap CI) vs conclusion scIB ARI-vs-truth (mean ± SD across seeds, reused
      from the scIB metrics) — the SD whisker is the quality-reproducibility signal.
      The point: blueprint CKA matches/exceeds scVI's and far exceeds freeform's, and
      blueprint quality varies far less run-to-run than freeform. (The complementary
      pure cross-run agreement metric, cluster_ari, is in reproducibility.csv.)

Inputs:
    outs/embedding_reproducibility/<dataset>/reproducibility.csv
    outs/blueprint_variability/<dataset>/cohort_metrics.csv   (spec gene-Jaccard)

Outputs:
    outs/plots/reproducibility/<dataset>/reproducibility.{png,pdf}

Usage:
    python scripts/plots/plot_reproducibility.py --dataset human_pancreas
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

# Slightly larger text for the figure.
plt.rcParams.update({
    "font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12,
    "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 11,
})


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")

# Match the three-way summary palette so the two figures read as a set; scVI
# baselines get neutral greys.
BLUEPRINT_COLORS = ("#1f77b4", "#08519c", "#3182bd", "#6baed6")
FREEFORM_COLOR = {"constrained": "#ff7f0e", "naive": "#d62728", "noref": "#8c564b"}
BASELINE_COLOR = {"scvi": "#525252", "scvi_nobatch": "#969696"}


# Blueprint-arm labels aligned with the leaderboard convention
# ('BP <engine> <flags> [effort]'). The engine map is ordered so 'gpt_5_4_mini'
# is matched before 'gpt_5_4'. Flags: N (noref), U (unrefined, Stage A),
# R (refined, Stage A+B); non-high effort is appended verbatim.
_ENGINE_PRETTY = {"gpt_5_5": "GPT-5.5", "gpt_5_4_mini": "GPT-5.4-mini", "gpt_5_4": "GPT-5.4"}
_FLAG_PRETTY = {"noref": "N", "unrefined": "U", "refine": "R"}
_EFFORT_TOKENS = ("low", "medium", "high", "xhigh")


def arm_label(cohort: str) -> str:
    """Blueprint-arm label aligned with the leaderboard, e.g.
    blueprint_gpt_5_5_noref_unrefined_high -> 'BP GPT-5.5 N U',
    blueprint_gpt_5_4_mini_unrefined_medium -> 'BP GPT-5.4-mini U medium'."""
    if cohort.startswith("blueprint_"):
        s = cohort[len("blueprint_"):]
        eng = ""
        for k, v in _ENGINE_PRETTY.items():
            if k in s:
                eng, s = v, s.replace(k, "_").strip("_")
                break
        parts = [p for p in s.split("_") if p]
        flags = [_FLAG_PRETTY[p] for p in parts if p in _FLAG_PRETTY]
        effort = next((p for p in parts if p in _EFFORT_TOKENS and p != "high"), "")
        return " ".join(["BP"] + ([eng] if eng else []) + flags + ([effort] if effort else []))
    if cohort.startswith("freeform_"):
        v = cohort.split("_")[1]                     # data token; remap display only
        return "FF-" + ("literature-grounded" if v == "constrained" else v)
    return {"scvi": "scVI", "scvi_nobatch": "scVI-nobatch"}.get(cohort, cohort)


def arm_color(cohort: str, bp_index: int) -> str:
    if cohort.startswith("blueprint_"):
        return BLUEPRINT_COLORS[bp_index % len(BLUEPRINT_COLORS)]
    if cohort.startswith("freeform_"):
        return FREEFORM_COLOR.get(cohort.split("_")[1], "#7f7f7f")
    return BASELINE_COLOR.get(cohort, "#7f7f7f")


def spec_jaccard_by_cohort(dataset: str) -> dict[str, float]:
    """Map blueprint embedding-cohort stem -> spec gene-Jaccard, read from the
    corrected blueprint-variability cohort metrics. The variability cohort name is
    the embedding stem minus the 'blueprint_' prefix plus '_<dataset>'."""
    path = REPO_ROOT / "outs" / "blueprint_variability" / dataset / "cohort_metrics.csv"
    if not path.is_file():
        return {}
    vdf = pd.read_csv(path).set_index("cohort")["jaccard_genes_off_diag"]
    out: dict[str, float] = {}
    for stem in (c for c in vdf.index):
        out[f"blueprint_{stem.removesuffix('_' + dataset)}"] = float(vdf[stem])
    return out


def _style(ax, arms: list[str]) -> None:
    ax.set_xticks(range(len(arms)))
    ax.set_xticklabels([arm_label(a) for a in arms], fontsize=9, rotation=40, ha="right")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def panel_nonpropagation(ax, df: pd.DataFrame, dataset: str) -> None:
    bp = df[df["kind"] == "blueprint"].sort_values("cohort").reset_index(drop=True)
    spec = spec_jaccard_by_cohort(dataset)
    arms = bp["cohort"].tolist()
    x = np.arange(len(arms))
    w = 0.38
    spec_vals = [spec.get(c, np.nan) for c in arms]
    cka_vals = bp["cka_mean"].to_numpy()
    cka_err = np.vstack([cka_vals - bp["cka_ci_lo"], bp["cka_ci_hi"] - cka_vals])
    ax.bar(x - w / 2, spec_vals, w, color="#bdbdbd", edgecolor="white",
           label="spec gene-Jaccard (L1)")
    ax.bar(x + w / 2, cka_vals, w, yerr=cka_err, capsize=3,
           color=[arm_color(c, i) for i, c in enumerate(arms)], edgecolor="white",
           label="representation CKA (L2)")
    cka_hi = bp["cka_ci_hi"].to_numpy()
    for xi, sv, cv, hi in zip(x, spec_vals, cka_vals, cka_hi):
        if np.isfinite(sv):
            ax.text(xi - w / 2, sv + 0.02, f"{sv:.2f}", ha="center", va="bottom", fontsize=8)
        ax.text(xi + w / 2, hi + 0.02, f"{cv:.2f}", ha="center", va="bottom", fontsize=8)
    ax.set_ylim(0, 1.18)
    ax.set_ylabel("agreement across runs")
    ax.set_title("Specification churns, the embedding does not")
    ax.legend(loc="upper left", frameon=False, fontsize=10)
    _style(ax, arms)


def panel_yardstick(ax, df: pd.DataFrame) -> None:
    order = (["scvi", "scvi_nobatch"]
             + sorted(df[df["kind"] == "blueprint"]["cohort"])
             + sorted(df[df["kind"] == "freeform"]["cohort"]))
    arms = [c for c in order if c in set(df["cohort"])]
    sub = df.set_index("cohort").loc[arms]
    x = np.arange(len(arms))
    w = 0.38
    colors = [arm_color(c, i) for i, c in enumerate(arms)]
    cka = sub["cka_mean"].to_numpy()
    cka_err = np.vstack([cka - sub["cka_ci_lo"], sub["cka_ci_hi"] - cka])
    ax.bar(x - w / 2, cka, w, yerr=cka_err, capsize=3,
           color=colors, edgecolor="white", label="representation CKA (L2)")
    # L3: scIB ARI-vs-truth, mean ± SD across seeds (SD = quality reproducibility).
    ax.bar(x + w / 2, sub["scib_ari_mean"].to_numpy(), w,
           yerr=sub["scib_ari_sd"].to_numpy(), capsize=3, color="none",
           edgecolor=colors, hatch="////", linewidth=1.0,
           label="conclusion scIB ARI (L3, mean±SD)")
    # scVI reference band: its CKA CI, so blueprint/freeform read against it.
    sc = sub.loc["scvi"] if "scvi" in sub.index else None
    if sc is not None:
        ax.axhspan(sc["cka_ci_lo"], sc["cka_ci_hi"], color="#525252", alpha=0.10)
    ax.set_ylim(0, 1.18)
    ax.set_ylabel("agreement / score across runs")
    ax.set_title("Reproducibility yardstick — blueprint vs scVI vs free-form")
    ax.legend(loc="upper right", frameon=False, fontsize=10)
    _style(ax, arms)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Reproducibility companion figure.")
    p.add_argument("--dataset", required=True, choices=DATASET_CHOICES)
    p.add_argument("--out-dir", type=Path,
                   default=REPO_ROOT / "outs" / "plots" / "reproducibility")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    repro_csv = (REPO_ROOT / "outs" / "embedding_reproducibility" / args.dataset
                 / "reproducibility.csv")
    if not repro_csv.is_file():
        raise SystemExit(f"{repro_csv} not found — run embedding_reproducibility first.")
    df = pd.read_csv(repro_csv)

    fig, axes = plt.subplots(1, 2, figsize=(16.0, 6.2))
    panel_nonpropagation(axes[0], df, args.dataset)
    panel_yardstick(axes[1], df)
    fig.tight_layout()

    out_dir = args.out_dir / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = out_dir / "reproducibility"
    for ext in ("png", "pdf"):
        fig.savefig(stem.with_suffix(f".{ext}"), dpi=200, bbox_inches="tight")
        print(f"[ok] wrote {stem.with_suffix('.' + ext)}", flush=True)
    plt.close(fig)

    # Per-panel files, so the two panels can be lettered (a)/(b) separately. The
    # combined figure above is kept for quick inspection.
    figa, axa = plt.subplots(figsize=(8.0, 6.2))
    panel_nonpropagation(axa, df, args.dataset)
    figa.tight_layout()
    for ext in ("png", "pdf"):
        figa.savefig((out_dir / "reproducibility_a").with_suffix(f".{ext}"),
                     dpi=200, bbox_inches="tight")
    plt.close(figa)

    figb, axb = plt.subplots(figsize=(8.0, 6.2))
    panel_yardstick(axb, df)
    figb.tight_layout()
    for ext in ("png", "pdf"):
        figb.savefig((out_dir / "reproducibility_b").with_suffix(f".{ext}"),
                     dpi=200, bbox_inches="tight")
    plt.close(figb)
    print("[ok] wrote per-panel reproducibility_a / reproducibility_b", flush=True)


if __name__ == "__main__":
    main()
