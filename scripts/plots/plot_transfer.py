"""
Held-out reference->query transfer figure.

Reads outs/metrics/<dataset>/transfer_metrics.csv (written by
scripts/transfer/measure.py), aggregates mean +/- sd across seeds per
method-variant, and renders the frozen->surgery gradient explicitly:

  - blueprint                : single zero-shot anchor bar (no training).
  - scvi / expimap _frozen   : query mapped with NO fine-tune.
  - scvi / expimap _scarches : query mapped WITH scArches surgery.

Three panels: (1) kNN label-transfer accuracy + macro-F1, (2) scIB Total with
bio/batch sub-scores (same metrics as the main benchmark), (3) surgery
wall-clock (blueprint = 0, frozen ~ 0, scarches > 0).

Inputs:
    outs/metrics/<dataset>/transfer_metrics.csv

Outputs:
    outs/plots/transfer/<dataset>/transfer_eval.{svg,png}

Usage:
    python scripts/plots/plot_transfer.py --dataset mouse_cortex
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch

REPO_ROOT = Path(__file__).resolve().parents[2]

# canonical left->right order + display labels
VARIANT_ORDER = ["blueprint", "blueprint_noref", "scvi_frozen", "scvi_scarches",
                 "expimap_frozen", "expimap_scarches"]
VARIANT_LABEL = {
    "blueprint": "blueprint\n(literature-\ngrounded)", "blueprint_noref": "blueprint\n(noref)",
    "scvi_frozen": "scVI\nfrozen", "scvi_scarches": "scVI\nscArches",
    "expimap_frozen": "expiMap\nfrozen", "expimap_scarches": "expiMap\nscArches",
}
# color by family; surgery state distinguished by hatch
FAMILY_COLOR = {"blueprint": "#c1272d", "scvi": "#2a6f97", "expimap": "#3a8a3a"}
VARIANT_COLOR = {"blueprint": "#c1272d", "blueprint_noref": "#e07b39"}  # literature-grounded vs noref
HATCH = {"none": "", "frozen": "///", "scarches": ""}


def variant_color(variant, family):
    return VARIANT_COLOR.get(variant, FAMILY_COLOR[family])


def aggregate(df: pd.DataFrame) -> pd.DataFrame:
    """Mean +/- sd across seeds, one row per variant, in canonical order."""
    metrics = ["label_transfer_acc", "label_transfer_macro_f1",
               "composite_scib", "bio_score", "batch_score_scib", "surgery_secs"]
    g = df.groupby("variant")
    out = g[metrics].mean()
    sd = g[metrics].std().add_suffix("_sd")
    out = out.join(sd)
    out["family"] = g["method"].first()
    out["surgery"] = g["surgery"].first()
    out["n_seeds"] = g.size()
    present = [v for v in VARIANT_ORDER if v in out.index]
    return out.loc[present]


def _grouped_bars(ax, agg, cols, col_labels, title, ylabel):
    variants = list(agg.index)
    x = np.arange(len(variants))
    w = 0.8 / len(cols)
    top = 0.0
    # graded shade across the metric columns (e.g. Total / bio / batch) so adjacent
    # same-color columns are distinguishable, not just first-vs-rest
    alphas = np.linspace(1.0, 0.4, len(cols))
    for j, col in enumerate(cols):
        vals = agg[col].to_numpy(dtype=float)
        err = agg.get(f"{col}_sd")
        err = np.nan_to_num(err.to_numpy(dtype=float)) if err is not None else None
        colors = [variant_color(v, f) for v, f in zip(agg.index, agg["family"])]
        alpha = float(alphas[j])
        ax.bar(x + (j - (len(cols) - 1) / 2) * w, vals, w, yerr=err,
               capsize=2, color=colors, alpha=alpha,
               hatch=[HATCH[s] for s in agg["surgery"]],
               edgecolor="white", linewidth=0.6)
        for xi, v, e in zip(x + (j - (len(cols) - 1) / 2) * w, vals, (err if err is not None else [0]*len(vals))):
            if np.isfinite(v):
                ax.text(xi, v + e + 0.035 * max(vals.max(), 1e-9), f"{v:.2f}",
                        ha="center", va="bottom", fontsize=7.5, rotation=90)
                top = max(top, v + e)
    # headroom so the rotated value labels and the in-axes legend clear the bars
    ax.set_ylim(0, top * 1.36 if top > 0 else 1.0)
    ax.set_xticks(x)
    ax.set_xticklabels([VARIANT_LABEL.get(v, v) for v in variants], fontsize=9)
    ax.set_title(title, fontsize=12, pad=8)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.tick_params(axis="y", labelsize=9)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    if len(cols) > 1:
        # sit the metric legend in the clear band above the bars, not over them
        ax.legend(handles=[Patch(facecolor="0.4", alpha=float(alphas[j]), label=lab)
                           for j, lab in enumerate(col_labels)],
                  fontsize=8.5, frameon=False, loc="upper center", ncol=len(cols),
                  bbox_to_anchor=(0.5, 1.0), columnspacing=1.1, handlelength=1.3,
                  handletextpad=0.4)


def main() -> None:
    ap = argparse.ArgumentParser(description="Plot the held-out transfer Extended-Data figure.")
    ap.add_argument("--dataset", required=True)
    args = ap.parse_args()

    csv = REPO_ROOT / "outs" / "metrics" / args.dataset / "transfer_metrics.csv"
    if not csv.exists():
        raise FileNotFoundError(f"{csv} not found; run scripts/transfer/measure.py first.")
    df = pd.read_csv(csv)
    df = df[df["dataset"] == args.dataset]
    agg = aggregate(df)
    n_seeds = int(agg["n_seeds"].max())

    fig, axes = plt.subplots(1, 3, figsize=(15.5, 5.0))
    _grouped_bars(axes[0], agg, ["label_transfer_acc", "label_transfer_macro_f1"],
                  ["accuracy", "macro-F1"],
                  "Label transfer (reference -> query)", "score")
    _grouped_bars(axes[1], agg, ["composite_scib", "bio_score", "batch_score_scib"],
                  ["Total", "bio", "batch"],
                  "scIB on joint ref+query", "score")
    _grouped_bars(axes[2], agg, ["surgery_secs"], ["surgery wall-clock"],
                  "Query-side training cost", "seconds")
    axes[2].axhline(0, color="0.7", lw=0.6)

    # legend for family color + frozen hatch
    fam_handles = [Patch(facecolor=VARIANT_COLOR["blueprint"], label="blueprint (literature-grounded)"),
                   Patch(facecolor=VARIANT_COLOR["blueprint_noref"], label="blueprint (noref)"),
                   Patch(facecolor=FAMILY_COLOR["scvi"], label="scvi"),
                   Patch(facecolor=FAMILY_COLOR["expimap"], label="expimap"),
                   Patch(facecolor="0.5", hatch="///", label="frozen (no fine-tune)")]
    fig.legend(handles=fam_handles, loc="upper center", ncol=5, fontsize=9,
               frameon=False, bbox_to_anchor=(0.5, 1.0))
    ds = args.dataset.replace("_", " ")
    qk = df["query_values"].iloc[0] if "query_values" in df else "?"
    fig.suptitle(f"Held-out transfer — {ds} (query = {qk}; n={n_seeds} seeds)",
                 y=1.07, fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.93))

    out_dir = REPO_ROOT / "outs" / "plots" / "transfer" / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(out_dir / f"transfer_eval.{ext}", bbox_inches="tight", dpi=200)
    print(f"[ok] wrote {out_dir / 'transfer_eval.svg'} and .png")


if __name__ == "__main__":
    main()
