"""
Cross-dataset averaged bio-vs-batch frontier.

One marker per method = the unweighted mean over the three benchmark datasets of
its scIB bio-conservation (x) and batch-correction (y) sub-scores, each computed
at one cell-type resolution per dataset:

    mouse_cortex   -> celltype_harmonized   (fine subclass)
    human_pancreas -> celltype
    human_immune   -> final_annotation

Only methods present in ALL three datasets are plotted (so the average is over the
same three points for every method). The three underlying per-dataset points are
drawn faintly in the method's tier color (centroid + faint points), faint dashed
iso-composite lines (0.6*bio + 0.4*batch) let composite be read off position, and
the Pareto front over the centroids is highlighted.

The averaged value mixes datasets with different absolute score scales, so it is a
summary; per-dataset detail lives in the scIB tables (plot_scib_table.py).

Inputs:
    outs/metrics/<dataset>/metrics.csv   (written by embedding_metrics/metrics.py)

Outputs:
    outs/plots/metrics/frontier_avg.{png,pdf}

Usage:
    python scripts/plots/plot_frontier_avg.py
    python scripts/plots/plot_frontier_avg.py --overwrite
"""
from __future__ import annotations

import argparse
import importlib.util
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


def _load_plot_metrics():
    """Load scripts/plots/plot_metrics.py for its shared helpers (no side effects)."""
    path = Path(__file__).resolve().parent / "plot_metrics.py"
    spec = importlib.util.spec_from_file_location("plot_metrics", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


pm = _load_plot_metrics()

# One resolution per benchmark dataset (the fine cell-type label that is the
# honest cross-dataset analog; mouse_cortex's coarse class_harmonized is a supp).
DATASET_RES = {
    "mouse_cortex":   "celltype_harmonized",
    "human_pancreas": "celltype",
    "human_immune":   "final_annotation",
}

# Method -> tier, matching the prior-knowledge tiers. The
# data-driven tier (PCA, scVI with/without its batch key) fits on the data with no
# gene programs; the foundation tier adds large-corpus pretraining; the program-based
# tier are the marker-/program-based competitors (UCell and expiMap, both consuming
# explicit gene programs from a curated database); the agentic blueprint authors its
# programs from the literature. UCell and expiMap share the program-based class even
# though they differ in batch handling (UCell by construction, expiMap by correction),
# a distinction carried in the text rather than the color.
TIER_COLORS = {
    "datadriven":   "#4c78a8",   # PCA, scVI (de novo, with/without batch)
    "foundation":   "#9467bd",   # scVI-census, UCE
    "program":      "#2ca02c",   # UCell, expiMap (program-/marker-based competitors)
    "agentic":      "#d62728",   # blueprint
}
TIER_LABELS = {
    "datadriven":   "data-driven",
    "foundation":   "foundation",
    "program":      "program-based (UCell, expiMap)",
    "agentic":      "blueprint",
}
TIER_ORDER = ("datadriven", "foundation", "program", "agentic")


def tier_of(model: str) -> str:
    if model.startswith("blueprint_"):
        return "agentic"
    if model in ("scvi_census", "uce"):
        return "foundation"
    if model in ("ucell", "expimap"):
        return "program"
    # pca, scvi, scvi_nobatch
    return "datadriven"


def is_grounded_blueprint(model: str) -> bool:
    """The 'with references' blueprint: literature-grounded (every module cites a resource);
    the resource-unconstrained '_noref' cohort is the other grounding variant."""
    return model.startswith("blueprint_") and "noref" not in model


def marker_of(model: str) -> str:
    """Distinguish the two blueprint grounding variants: star = grounded (with
    references), diamond = resource-unconstrained (noref); circle = baseline."""
    if model.startswith("blueprint_"):
        return "*" if is_grounded_blueprint(model) else "D"
    return "o"


def marker_size(marker: str) -> float:
    """Per-marker area so star and diamond read at a comparable visual weight."""
    return {"*": 380, "D": 150}.get(marker, 90)


def short_label(model: str) -> str:
    """Compact point label for the frontier (shorter than prettify_model)."""
    if model.startswith("blueprint_"):
        return "BP" if is_grounded_blueprint(model) else "BP (noref)"
    return pm.prettify_model(model)


# Hand-tuned annotation offsets (points) pushing each label into nearby empty
# space so the dense cluster of near-coincident points stays legible. Two clusters
# need fanning: the marker/program cluster up top (scVI, expiMap, the two blueprint
# variants, UCell, all within bio~0.66-0.78 x batch~0.48-0.55) and the prior-free
# cloud around (0.69, 0.29). Top cluster: scVI -> NE, expiMap -> SE, UCell -> SW,
# blueprint noref -> N, blueprint literature-grounded -> S. Cloud: fanned N/E/S/W.
LABEL_OFFSETS = {
    "scvi":         (12, 6, "left"),     # NE (top-right of cluster)
    "expimap":      (11, -9, "left"),    # SE (into the open gap below-right)
    "ucell":        (-14, 4, "right"),   # W (leftmost cluster point), clears Blueprint
    "scvi_census":  (12, 7, "left"),     # upper-right
    "uce":          (-12, 7, "right"),   # upper-left
    "pca":          (12, -14, "left"),   # lower-right
    "scvi_nobatch": (-12, -14, "right"), # lower-left
}


def _offset_for(model: str) -> tuple[int, int, str]:
    if model in LABEL_OFFSETS:
        return LABEL_OFFSETS[model]
    if model.startswith("blueprint_"):
        # The two near-coincident blueprint markers split vertically: noref label
        # above, literature-grounded below; expiMap (SE) and UCell (SW) take the diagonals.
        return (0, 17, "center") if "noref" in model else (0, -20, "center")
    return (7, 6, "left")


def per_dataset_scores() -> pd.DataFrame:
    """Long frame: one row per (model, dataset) with bio_score / batch_score_scib.

    Aggregates seeds via plot_metrics.aggregate_seeds and keeps only the chosen
    resolution per dataset, and only rows with a finite batch score (a frontier
    point needs both axes).
    """
    rows = []
    for dataset, label_col in DATASET_RES.items():
        csv_path = REPO_ROOT / "outs" / "metrics" / dataset / "metrics.csv"
        if not csv_path.exists():
            raise FileNotFoundError(f"{csv_path} not found.")
        df = pd.read_csv(csv_path)
        agg = pm.aggregate_seeds(df)
        sub = agg[agg["label_col"] == label_col].copy()
        if sub.empty:
            raise ValueError(f"{dataset}: no rows for label_col={label_col!r}.")
        sub = sub[sub["batch_score_scib"].apply(np.isfinite)
                  & sub["bio_score"].apply(np.isfinite)]
        sub["dataset"] = dataset
        rows.append(sub[["model", "dataset", "bio_score", "batch_score_scib",
                         "composite_scib", "n_seeds"]])
    return pd.concat(rows, ignore_index=True)


# Benchmark cohort = the GPT-5.5-high blueprint (both grounding variants, unrefined)
# + baselines; the free-form, refine and gpt-5.4-mini cohorts are excluded from the
# frontier (mini lives only in human_pancreas, so leaving it in would add an
# off-cohort, duplicate-labelled point to the per-dataset panel).
R2_EXCLUDE = {"freeform", "refine", "mini"}


def is_r2_cohort(model: str) -> bool:
    return not (R2_EXCLUDE & set(model.split("_")))


def per_dataset_raw(dataset: str, label_col: str) -> pd.DataFrame:
    """Per-method mean/std of bio & batch across seeds for one dataset/resolution,
    restricted to the benchmark cohort (drops free-form/refine)."""
    df = pd.read_csv(REPO_ROOT / "outs" / "metrics" / dataset / "metrics.csv")
    df = df[(df["label_col"] == label_col) & df["batch_score_scib"].notna()
            & df["bio_score"].notna()]
    df = df[df["model"].map(is_r2_cohort)]
    g = df.groupby("model").agg(
        bio=("bio_score", "mean"), bio_sd=("bio_score", "std"),
        batch=("batch_score_scib", "mean"), batch_sd=("batch_score_scib", "std"),
        n=("seed", "size")).reset_index()
    return g


def draw_per_dataset(out_stem: Path) -> None:
    """Faceted bio-vs-batch frontier, one panel per benchmark dataset (un-averaged
    companion to the averaged frontier). Benchmark cohort only, seed error bars.
    Colored by the same method-class tiers as the averaged panel (TIER_COLORS)."""
    items = list(DATASET_RES.items())
    all_models = sorted({m for ds, lab in items
                         for m in per_dataset_raw(ds, lab)["model"]})
    # Extra height over the square panels reserves a clean band at the bottom for the
    # multi-row legend, so it never overlaps the per-panel x-axis labels.
    fig, axes = plt.subplots(1, len(items), figsize=(4.6 * len(items), 5.9),
                             squeeze=False)
    for ax, (ds, lab) in zip(axes[0], items):
        g = per_dataset_raw(ds, lab)
        for c in (0.4, 0.5, 0.6, 0.7, 0.8):
            bx = np.linspace(0, 1, 50); by = (c - 0.6 * bx) / 0.4
            m = (by >= 0) & (by <= 1)
            ax.plot(bx[m], by[m], ls="--", lw=0.6, color="#dddddd", zorder=0)
        for _, r in g.iterrows():
            mk = marker_of(r["model"]); col = TIER_COLORS[tier_of(r["model"])]
            ax.errorbar(r["bio"], r["batch"],
                        xerr=0 if np.isnan(r["bio_sd"]) else r["bio_sd"],
                        yerr=0 if np.isnan(r["batch_sd"]) else r["batch_sd"],
                        fmt=mk, ms=14 if mk == "*" else 9, color=col, ecolor=col,
                        mec="black", mew=0.5, capsize=2, lw=0.7, zorder=4)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal", adjustable="box")
        ax.set_title(f"{ds}\n({lab})", fontsize=10)
        ax.set_xlabel("Bio-conservation ↑")
        ax.grid(linestyle=":", linewidth=0.4, alpha=0.5)
    axes[0][0].set_ylabel("Batch-correction ↑")
    handles = [plt.Line2D([0], [0], marker=marker_of(m), linestyle="",
                          markerfacecolor=TIER_COLORS[tier_of(m)], markeredgecolor="black",
                          markeredgewidth=0.5, markersize=11 if marker_of(m) == "*" else 8,
                          label=short_label(m)) for m in all_models]
    fig.legend(handles=handles, loc="lower center", ncol=min(5, len(handles)),
               frameon=False, fontsize=8.5, bbox_to_anchor=(0.5, 0.01))
    fig.suptitle("Quality vs. batch robustness, per benchmark dataset", fontsize=12)
    # Reserve a generous bottom band for the 2-row legend (clear of the per-panel
    # x-axis labels) and a strip at top for the suptitle.
    fig.tight_layout(rect=(0, 0.22, 1, 0.95))
    pm._save(fig, out_stem)


def pareto_front(points: list[tuple[float, float]]) -> list[int]:
    """Indices of non-dominated points (both bio and batch higher = better)."""
    keep = []
    for i, (bi, ba) in enumerate(points):
        dominated = any(
            (bj >= bi and baj >= ba) and (bj > bi or baj > ba)
            for j, (bj, baj) in enumerate(points) if j != i
        )
        if not dominated:
            keep.append(i)
    return keep


def draw(long: pd.DataFrame, out_stem: Path) -> None:
    n_ds = long["dataset"].nunique()
    counts = long.groupby("model")["dataset"].nunique()
    methods = sorted(counts[counts == n_ds].index)  # present in ALL datasets
    dropped = sorted(counts[counts < n_ds].index)
    if dropped:
        print(f"[info] excluded (not in all {n_ds} datasets): {dropped}")

    # Centroid per method = unweighted mean over datasets.
    cent = (long[long["model"].isin(methods)]
            .groupby("model")
            .agg(bio=("bio_score", "mean"), batch=("batch_score_scib", "mean"),
                 n_seeds=("n_seeds", "sum"))
            .reset_index())

    fig, ax = plt.subplots(figsize=(6.6, 5.4))

    # Zoomed view into the data envelope, identical on both axes (square). The range
    # is set by the data: the lowest plotted point is batch ~0.235 (UCE) and the
    # highest is bio ~0.834 (scVI), so [0.20, 0.85] keeps every method in frame ---
    # including the collapsed prior-free baselines (batch 0.24-0.38), which are the
    # contrast --- while still zooming in vs the full [0,1] box for legibility.
    XLIM = YLIM = (0.20, 0.85)

    # Iso-composite reference lines: 0.6*bio + 0.4*batch = c. Label each where it
    # exits the visible window (rightmost in-view point), since the [0,1] corners are
    # now off-frame.
    bio_grid = np.linspace(0, 1, 200)
    for c in (0.4, 0.5, 0.6, 0.7, 0.8):
        batch_line = (c - 0.6 * bio_grid) / 0.4
        m = (batch_line >= 0) & (batch_line <= 1)
        ax.plot(bio_grid[m], batch_line[m], ls="--", lw=0.6, color="#cccccc",
                zorder=0)
        vis = m & (bio_grid >= XLIM[0]) & (bio_grid <= XLIM[1]) \
                & (batch_line >= YLIM[0]) & (batch_line <= YLIM[1])
        if vis.any():
            xi = bio_grid[vis][-1]
            yi = batch_line[vis][-1]
            ax.text(xi, yi, f" {c:.1f}", color="#999999", fontsize=8.5,
                    va="bottom", ha="right", zorder=0)

    # Pareto front over centroids (drawn only when >1 point is non-dominated; on
    # the averaged plot batch-keyed scVI typically dominates all, leaving a single
    # front point and no line — that single-point case is itself the honest result).
    pts = list(zip(cent["bio"], cent["batch"]))
    front = sorted(pareto_front(pts), key=lambda i: pts[i][0])
    drew_front = len(front) >= 2
    if drew_front:
        fx = [pts[i][0] for i in front]
        fy = [pts[i][1] for i in front]
        ax.plot(fx, fy, "-", lw=1.2, color="#444444", alpha=0.55, zorder=1,
                label="Pareto front")

    # Faint per-dataset points + centroid markers.
    for _, r in cent.iterrows():
        model = r["model"]
        tier = tier_of(model)
        color = TIER_COLORS[tier]
        mk = marker_of(model)
        big = marker_size(mk)
        # faint per-dataset points, lightly tethered to the centroid
        sub = long[long["model"] == model]
        for _, p in sub.iterrows():
            ax.plot([r["bio"], p["bio_score"]], [r["batch"], p["batch_score_scib"]],
                    color=color, lw=0.5, alpha=0.18, zorder=2)
            ax.scatter(p["bio_score"], p["batch_score_scib"], s=22, marker=mk,
                       color=color, alpha=0.30, edgecolors="none", zorder=2)
        ax.scatter(r["bio"], r["batch"], s=big, marker=mk, color=color,
                   edgecolors="black", linewidths=0.6, zorder=4)
        dx, dy, ha = _offset_for(model)
        ax.annotate(short_label(model), (r["bio"], r["batch"]),
                    textcoords="offset points", xytext=(dx, dy), fontsize=9.5,
                    color="#222", ha=ha, zorder=5)

    ax.set_xlabel("Bio-conservation score ↑")
    ax.set_ylabel("Batch-correction score ↑")
    ax.set_title("Quality vs. batch robustness (mean of 3 datasets)", fontsize=13)
    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.set_aspect("equal", adjustable="box")
    ax.grid(linestyle=":", linewidth=0.4, alpha=0.5)

    # Legend: color = method class; the two blueprint grounding variants get their
    # own marker (star = resource-unconstrained, diamond = literature-grounded).
    def _h(marker, facecolor, label, edge="black", ms=9):
        return plt.Line2D([0], [0], marker=marker, linestyle="",
                          markerfacecolor=facecolor, markeredgecolor=edge,
                          markeredgewidth=0.5, markersize=ms, label=label)

    handles = [
        _h("o", TIER_COLORS["datadriven"], "data-driven"),
        _h("o", TIER_COLORS["foundation"], "foundation"),
        _h("o", TIER_COLORS["program"], "program-based (UCell, expiMap)"),
        _h("*", TIER_COLORS["agentic"], "blueprint (literature-grounded)", ms=14),
        _h("D", TIER_COLORS["agentic"], "blueprint (noref)", ms=8),
    ]
    if drew_front:
        handles.append(plt.Line2D([0], [0], color="#444444", lw=1.2, alpha=0.55,
                                  label="Pareto front"))
    ax.legend(handles=handles, loc="upper left", frameon=False, fontsize=10,
              handletextpad=0.4, labelspacing=0.4)

    fig.tight_layout()
    pm._save(fig, out_stem)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Cross-dataset averaged bio-vs-batch frontier.")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--per-dataset", action="store_true",
                   help="Also emit the un-averaged per-dataset faceted frontier "
                        "(frontier_per_dataset.{png,pdf}).")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    matplotlib.rcParams.update(pm.PLOT_RC)
    # Larger fonts: this panel is shown at reduced size in a multi-panel layout,
    # so text needs to stay readable.
    matplotlib.rcParams.update({
        "font.size": 13, "axes.labelsize": 13, "axes.titlesize": 14,
        "xtick.labelsize": 12, "ytick.labelsize": 12,
    })

    out_dir = REPO_ROOT / "outs" / "plots" / "metrics"
    out_stem = out_dir / "frontier_avg"
    if out_stem.with_suffix(".png").exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite {out_stem}.png (pass --overwrite).")

    long = per_dataset_scores()
    draw(long, out_stem)
    print(f"[ok] wrote {out_stem}.{{png,pdf}}")

    pd_stem = None
    if args.per_dataset:
        pd_stem = out_dir / "frontier_per_dataset"
        draw_per_dataset(pd_stem)
        print(f"[ok] wrote {pd_stem}.{{png,pdf}}")



if __name__ == "__main__":
    main()
