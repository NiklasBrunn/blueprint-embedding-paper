"""
Figures from scripts/embedding_metrics/metrics.py outputs.

Reads a long-format metrics CSV (one row per (parquet, label_col)) and renders
three figures per (dataset, label_col):

  1. Heatmap   — scIB-style overview: rows = models (seeds aggregated to mean),
                 columns = metrics grouped by Bio / Batch / Composite.
  2. Bar grid  — one panel per metric; bar per model, error bars across seeds.
  3. Tradeoff  — bio score vs batch score scatter (one point per model);
                 only emitted when batch metrics are available.

Metrics are the scIB suite as computed by the scib-metrics Benchmarker
(bio: ARI, NMI, ASW-label, cLISI, isolated-label; batch: ASW-batch, iLISI, kBET,
graph connectivity, PCR). The composite is the package Total = 0.6*bio + 0.4*batch.

Figures land alongside the metrics tables to keep the project layout uniform.

Inputs:
    outs/metrics/<dataset>/metrics.csv

Outputs (per label_col):
    outs/plots/metrics/<dataset>/<label_col>/heatmap.{png,pdf}
    outs/plots/metrics/<dataset>/<label_col>/bars.{png,pdf}
    outs/plots/metrics/<dataset>/<label_col>/tradeoff.{png,pdf}   (when batch metrics present)

Usage:
    python scripts/plots/plot_metrics.py --dataset mouse_brain
    python scripts/plots/plot_metrics.py --dataset human_pancreas
    python scripts/plots/plot_metrics.py --dataset human_pancreas --label-cols celltype
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
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")

# Metric → (group, display label, direction). Direction is "↑" (higher=better)
# or "↓" (lower=better). All scib-metrics outputs are rescaled to [0,1] with
# higher=better, so every entry here is "↑" today; the column exists so future
# raw-distance metrics can declare "↓" without restructuring callers.
METRIC_SPEC: list[tuple[str, str, str, str]] = [
    ("ari",                "Bio",        "ARI",            "↑"),
    ("nmi",                "Bio",        "NMI",            "↑"),
    ("asw_label",          "Bio",        "ASW (label)",    "↑"),
    ("clisi",              "Bio",        "cLISI",          "↑"),
    ("isolated_label_f1",  "Bio",        "Isolated",       "↑"),
    # graph_connectivity is a BATCH-correction metric in scIB (Luecken 2022),
    # not bio — it belongs with the batch group here and in the composite.
    ("asw_batch",          "Batch",      "ASW (batch)",    "↑"),
    ("ilisi",              "Batch",      "iLISI",          "↑"),
    ("kbet_acceptance",    "Batch",      "kBET",           "↑"),
    ("graph_connectivity", "Batch",      "Graph conn",     "↑"),
    ("pcr",                "Batch",      "PCR",            "↑"),
    ("trajectory_corr",    "Trajectory", "stage corr",     "↑"),
]
BIO_METRICS   = [m for m, g, _, _ in METRIC_SPEC if g == "Bio"]
BATCH_METRICS = [m for m, g, _, _ in METRIC_SPEC if g == "Batch"]

# scIB scoring sets — MUST mirror SCIB_* in embedding_metrics/metrics.py (the
# canonical source). Used only as the fallback when metrics.csv lacks the
# precomputed bio_score/batch_score_scib/composite_scib columns (which the
# Benchmarker normally writes).
SCIB_BIO_METRICS:   list[str] = ["ari", "nmi", "asw_label", "clisi", "isolated_label_f1"]
SCIB_BATCH_METRICS: list[str] = ["asw_batch", "ilisi", "kbet_acceptance", "graph_connectivity", "pcr"]
# Batch metrics requiring a batch key — gate the batch score so metrics computable
# without a batch can't fabricate a batch score on batchless datasets.
BATCH_REQUIRED:     list[str] = ["asw_batch", "ilisi", "kbet_acceptance"]


def _annotate(label: str, direction: str) -> str:
    """Append the higher/lower-is-better arrow to a metric label."""
    return f"{label} {direction}"


# Display names for the `model` column (kept in sync with
# embedding_metrics/metrics.py:prettify_model). Used only for axis/legend labels;
# raw `model` stays the join/color key. noref->N, unrefined->U, refine->R.
_BASELINE_PRETTY = {
    "pca": "PCA", "scvi": "scVI", "scvi_census": "scVI-census",
    "scvi_nobatch": "scVI (no batch)", "uce": "UCE", "ucell": "UCell", "expimap": "expiMap",
}
_ENGINE_PRETTY = {"gpt_5_5": "GPT-5.5", "gpt_5": "GPT-5"}
_VARIANT_PRETTY = {"constrained": "literature-grounded", "naive": "naïve", "freeform": "freeform"}
_REFSTATE_PRETTY = {"unrefined": "U", "refine": "R"}
_EFFORT = ("low", "medium", "high", "xhigh")


def prettify_model(name: str) -> str:
    """Compact display label for a `model` value. Unknown names pass through."""
    if name in _BASELINE_PRETTY:
        return _BASELINE_PRETTY[name]
    if name.startswith("blueprint_"):
        family, s = "BP", name[len("blueprint_"):]
    elif name.startswith("freeform_"):
        family, s = "FF", name[len("freeform_"):]
    else:
        return name
    eng = ""
    for k, v in _ENGINE_PRETTY.items():
        if k in s:
            eng = v
            s = s.replace(k, "_").strip("_")
            break
    parts = [p for p in s.split("_") if p]
    effort = next((p for p in parts if p in _EFFORT), "")
    flags, variant = [], []
    for p in parts:
        if p in _EFFORT:
            continue
        if p in _REFSTATE_PRETTY:
            flags.append(_REFSTATE_PRETTY[p])
        elif p == "noref":
            flags.append("N")
        elif p == "userctx":
            flags.append("userctx")
        else:
            variant.append(_VARIANT_PRETTY.get(p, p))
    label = [family] + variant + ([eng] if eng else []) + flags
    if effort and effort != "high":
        label.append(effort)
    return " ".join(label)

# Stable colors for the long-established baselines so existing figures stay
# recognizable. Anything not in this map (e.g. blueprint_<model>_<reasoning>_<mode>
# variants) gets a deterministic slot from FALLBACK_PALETTE below, assigned by
# alphabetical position so colors are stable across reruns and across datasets.
MODEL_COLORS = {
    "pca":         "#888888",
    "scvi":        "#1f77b4",
    "scvi_census": "#75aadb",
    "expimap":     "#2ca02c",
    "ucell":       "#ff7f0e",
    "uce":         "#9467bd",
}

# 20-slot qualitative palette for unknown models. Hand-picked to be visually
# distinct from each other AND from the MODEL_COLORS values above (so adding
# baselines later won't collide). Designed to scale to ~20 blueprint variants.
FALLBACK_PALETTE = [
    "#d62728",  # red
    "#17becf",  # cyan
    "#bcbd22",  # olive
    "#e377c2",  # pink
    "#8c564b",  # brown
    "#393b79",  # navy
    "#637939",  # dark olive
    "#ad494a",  # brick
    "#a55194",  # mauve
    "#3182bd",  # steel blue
    "#31a354",  # leaf green
    "#e6550d",  # burnt orange
    "#756bb1",  # slate purple
    "#969696",  # mid grey
    "#fdae6b",  # sand
    "#9e9ac8",  # lavender
    "#c7e9c0",  # mint
    "#fdd0a2",  # peach
    "#dadaeb",  # lilac
    "#f7b6d2",  # rose
]

GROUP_COLORS = {"Bio": "#4c78a8", "Batch": "#e1812c", "Composite": "#54a24b", "Trajectory": "#9c27b0"}

DEFAULT_PARAMS = {
    "mouse_brain":         {"label_cols": ("primary_type", "subclass", "broad_type")},
    "human_pancreas":      {"label_cols": ("celltype",)},
    "mouse_gastrulation":  {"label_cols": ("celltype", "stage")},
    "mouse_cortex":        {"label_cols": ("celltype_harmonized", "class_harmonized")},
    "human_immune":        {"label_cols": ("final_annotation",)},
}

PLOT_RC = {
    "font.family": "DejaVu Sans",
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.facecolor": "white", "savefig.transparent": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
}

HEATMAP_CMAP = LinearSegmentedColormap.from_list(
    "scib_heat", ["#f7f7f7", "#cfe2f3", "#6ba6cd", "#2e5c8a", "#1a3552"], N=256,
)


def resolve_colors(models: list[str]) -> dict[str, str]:
    """Deterministic model → color map for a given set of models.

    Known baselines pull from MODEL_COLORS. Unknown models (e.g. blueprint_*
    variants) are sorted alphabetically and assigned consecutive slots from
    FALLBACK_PALETTE so the same model gets the same color across runs and
    across datasets, with no hash-based collisions.
    """
    out: dict[str, str] = {}
    unknowns = sorted(m for m in models if m not in MODEL_COLORS)
    for i, m in enumerate(unknowns):
        out[m] = FALLBACK_PALETTE[i % len(FALLBACK_PALETTE)]
    for m in models:
        if m in MODEL_COLORS:
            out[m] = MODEL_COLORS[m]
    return out


def aggregate_seeds(df: pd.DataFrame) -> pd.DataFrame:
    """Mean / std / n_seeds per (model, label_col) for every metric column.

    Metrics declared in METRIC_SPEC but absent from the CSV (e.g. older runs
    predating a column) are skipped to keep aggregation backwards-compatible.
    """
    metric_cols = [m for m, _, _, _ in METRIC_SPEC if m in df.columns]
    # Canonical per-row scores written by embedding_metrics/metrics.py (straight
    # from the Benchmarker aggregates). When present these are averaged over seeds
    # (they ARE the single source of truth); otherwise they are derived below.
    score_cols = [c for c in ("bio_score", "batch_score_scib", "composite_scib")
                  if c in df.columns]
    agg = df.groupby(["model", "label_col"], dropna=False).agg(
        **{f"{m}_mean": (m, "mean") for m in metric_cols},
        **{f"{m}_std":  (m, "std")  for m in metric_cols},
        **{c: (c, "mean") for c in score_cols},
        n_seeds=("seed", "size"),
    ).reset_index()

    canonical = {"bio_score", "batch_score_scib", "composite_scib"}
    if not canonical.issubset(agg.columns):
        # Fallback: derive from the scIB sets (composite is linear, so averaging
        # per-row equals deriving from per-metric means).
        bio_means = agg[[f"{m}_mean" for m in SCIB_BIO_METRICS if f"{m}_mean" in agg.columns]]
        agg["bio_score"] = bio_means.mean(axis=1, skipna=True)
        # Gate on batch-requiring metrics: some batch metrics are computable
        # without a batch, so on batchless datasets the batch score must be NaN.
        req_cols = [f"{m}_mean" for m in BATCH_REQUIRED if f"{m}_mean" in agg.columns]
        has_batch = (agg[req_cols].notna().any(axis=1) if req_cols
                     else pd.Series([False] * len(agg), index=agg.index))
        cols = [f"{m}_mean" for m in SCIB_BATCH_METRICS if f"{m}_mean" in agg.columns]
        batch = (agg[cols].mean(axis=1, skipna=True) if cols
                 else pd.Series([float("nan")] * len(agg), index=agg.index))
        batch = batch.where(has_batch, other=float("nan"))
        agg["batch_score_scib"] = batch
        has_b = batch.notna()
        agg["composite_scib"] = np.where(
            has_b, 0.6 * agg["bio_score"] + 0.4 * batch, agg["bio_score"],
        )
    return agg


def _save(fig, out_stem: Path, dpi: int = 600) -> None:
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_stem.with_suffix(f".{ext}"), dpi=dpi, bbox_inches="tight",
                    facecolor="white")
    plt.close(fig)


def draw_heatmap(agg_label: pd.DataFrame, dataset: str, label_col: str,
                 out_stem: Path) -> None:
    """scIB-style heatmap: models × metrics, sorted by composite score.

    On batchless datasets (no batch key → batch score NaN) every batch column is
    hidden, including graph_connectivity (computed without a batch but reported as
    a batch-correction result here), so no "Batch" section appears spuriously.
    """
    bs_col, comp_col = "batch_score_scib", "composite_scib"
    batchless = not agg_label[bs_col].notna().any()

    metric_cols = [m for m, _, _, _ in METRIC_SPEC]
    metric_lbls = [_annotate(lbl, d) for _, _, lbl, d in METRIC_SPEC]
    metric_grps = [g for _, g, _, _ in METRIC_SPEC]

    # Drop metric columns that are absent (older CSV), all-NaN, or any batch
    # metric on a batchless dataset.
    matrix_cols = [
        m for m in metric_cols
        if f"{m}_mean" in agg_label.columns
        and agg_label[f"{m}_mean"].notna().any()
        and not (batchless and m in BATCH_METRICS)
    ]
    matrix_lbls = [metric_lbls[metric_cols.index(m)] for m in matrix_cols]
    matrix_grps = [metric_grps[metric_cols.index(m)] for m in matrix_cols]

    extras: list[tuple[str, str, str]] = [("bio_score", "Bio", "Bio ↑")]
    if agg_label[bs_col].notna().any():
        extras.append((bs_col, "Batch", "Batch ↑"))
    extras.append((comp_col, "Composite", "Overall ↑"))
    matrix_cols += [c for c, _, _ in extras]
    matrix_lbls += [lbl for _, _, lbl in extras]
    matrix_grps += [g for _, g, _ in extras]

    ordered = agg_label.sort_values(comp_col, ascending=False).reset_index(drop=True)
    models = ordered["model"].tolist()
    n_seeds = ordered["n_seeds"].tolist()

    spec_metric_set = {m for m, _, _, _ in METRIC_SPEC}
    rows: list[list[float]] = []
    for _, row in ordered.iterrows():
        rows.append([
            row[f"{c}_mean"] if c in spec_metric_set else row[c]
            for c in matrix_cols
        ])
    M = np.array(rows, dtype=float)

    n_rows, n_cols = M.shape
    fig_w = max(5.0, 0.55 * n_cols + 2.6)
    fig_h = max(2.4, 0.45 * n_rows + 1.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    masked = np.ma.masked_invalid(M)
    im = ax.imshow(masked, cmap=HEATMAP_CMAP, vmin=0, vmax=1, aspect="auto")
    HEATMAP_CMAP.set_bad("#e8e8e8")

    for i in range(n_rows):
        for j in range(n_cols):
            v = M[i, j]
            txt = "—" if not np.isfinite(v) else f"{v:.2f}"
            colour = "white" if (np.isfinite(v) and v > 0.55) else "#222"
            ax.text(j, i, txt, ha="center", va="center", fontsize=7.2, color=colour)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(matrix_lbls, rotation=35, ha="right")
    ax.set_yticks(range(n_rows))
    ax.set_yticklabels([f"{prettify_model(m)} (n={ns})" for m, ns in zip(models, n_seeds)])
    ax.tick_params(axis="both", which="both", length=0)

    # Group bar above the column labels.
    bar_y = -0.6
    for j, grp in enumerate(matrix_grps):
        ax.add_patch(plt.Rectangle((j - 0.45, bar_y - 0.18), 0.9, 0.18,
                                   facecolor=GROUP_COLORS[grp], edgecolor="none",
                                   clip_on=False))
    seen: set[str] = set()
    legend_handles = []
    for grp in matrix_grps:
        if grp not in seen:
            legend_handles.append(Patch(color=GROUP_COLORS[grp], label=grp))
            seen.add(grp)
    ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
              frameon=False, fontsize=8, title="Metric group", title_fontsize=8.5)

    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.18)
    cbar.set_label("score (higher = better)", fontsize=8)
    cbar.ax.tick_params(labelsize=7.5)

    ax.set_title(f"{dataset} — {label_col}", pad=14)
    fig.tight_layout()
    _save(fig, out_stem)


def draw_bars(agg_label: pd.DataFrame, dataset: str, label_col: str,
              out_stem: Path) -> None:
    """Faceted bar grid — one panel per metric, models on x, error bars across seeds.

    Each panel ranks models by that panel's metric (best on the left, descending);
    color identity is shared across panels so a model is always the same color.
    """
    batchless = not agg_label["batch_score_scib"].notna().any()
    metric_cols = [
        (m, _annotate(lbl, d)) for m, _, lbl, d in METRIC_SPEC
        if f"{m}_mean" in agg_label.columns
        and agg_label[f"{m}_mean"].notna().any()
        and not (batchless and m in BATCH_METRICS)
    ]
    n = len(metric_cols)
    if n == 0:
        return

    ncols = 3
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.8 * ncols, 3.2 * nrows),
                             squeeze=False, sharex=False)

    # Per-panel ranking: each panel sorts models by that metric's mean (descending,
    # NaNs last) so bars descend monotonically left-to-right within the panel. Color
    # identity is preserved across panels via a single dataset-wide color map.
    color_map = resolve_colors(agg_label["model"].tolist())

    for k, (metric, label) in enumerate(metric_cols):
        ax = axes[k // ncols][k % ncols]
        panel = (agg_label
                 .sort_values(f"{metric}_mean", ascending=False, na_position="last")
                 .reset_index(drop=True))
        panel_models = panel["model"].tolist()
        means = panel[f"{metric}_mean"].to_numpy()
        # Single-seed methods have NaN std → fill 0.0 only for plotting math, but
        # render those bars hatched and without yerr so reviewers can see the
        # absence of uncertainty data rather than reading "perfectly reproducible".
        raw_stds = panel[f"{metric}_std"].to_numpy()
        stds = np.where(np.isnan(raw_stds), 0.0, raw_stds)
        n_seeds = panel["n_seeds"].to_numpy()
        is_single = n_seeds <= 1
        yerr = np.where(is_single, 0.0, stds)
        colors = [color_map[m] for m in panel_models]
        hatches = ["//" if s else "" for s in is_single]
        bars = ax.bar(range(len(panel_models)), means, yerr=yerr, capsize=2.5,
                      color=colors, edgecolor="black", linewidth=0.4)
        for bar, hatch in zip(bars, hatches):
            if hatch:
                bar.set_hatch(hatch)
        ax.set_title(label, fontsize=9.5)
        # iLISI values land in a narrow band well below 1.0 (typically 0.0–0.3
        # on real data), so a fixed [0, 1.0] axis crushes the bars. Use a
        # tight max here (still anchored at 0) so model differences are
        # readable; all other metrics keep the [0, 1.0] scale for cross-panel
        # comparability.
        if metric == "ilisi":
            top = float(np.nanmax(means + stds))
            ax.set_ylim(0.0, top * 1.10 if np.isfinite(top) and top > 0 else 1.0)
        else:
            ax.set_ylim(0, 1.0 if np.nanmax(means + stds) <= 1.0 else None)
        ax.set_xticks(range(len(panel_models)))
        ax.set_xticklabels([prettify_model(m) for m in panel_models],
                           rotation=35, ha="right", fontsize=7.5)
        ax.tick_params(axis="y", labelsize=7.5)
        ax.grid(axis="y", linestyle=":", linewidth=0.4, alpha=0.5)

    # Hide unused axes (when n_metrics is not a multiple of ncols).
    for k in range(n, nrows * ncols):
        axes[k // ncols][k % ncols].axis("off")

    fig.suptitle(f"{dataset} — {label_col}", y=1.0, fontsize=10.5)
    fig.tight_layout()
    _save(fig, out_stem)


def draw_tradeoff(agg_label: pd.DataFrame, dataset: str, label_col: str,
                  out_stem: Path) -> None:
    """Bio vs batch scatter (one point per model). No-op when batch is unavailable.

    Identification is via a side legend rather than inline annotations so labels
    don't overlap when several models cluster in the same region of the plot.
    """
    bs_col, comp_col = "batch_score_scib", "composite_scib"
    if not agg_label[bs_col].notna().any():
        return

    plotted = agg_label[agg_label[bs_col].apply(np.isfinite)].copy()
    if plotted.empty:
        return
    plotted = plotted.sort_values(comp_col, ascending=False).reset_index(drop=True)
    color_map = resolve_colors(plotted["model"].tolist())

    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    handles = []
    for _, row in plotted.iterrows():
        c = color_map[row["model"]]
        ax.errorbar(
            row["bio_score"], row[bs_col],
            xerr=np.nanstd([row[f"{m}_mean"] for m in BIO_METRICS if f"{m}_mean" in row]),
            yerr=np.nanstd([row[f"{m}_mean"] for m in BATCH_METRICS if f"{m}_mean" in row]),
            fmt="o", ms=8, color=c, ecolor=c, capsize=2.5, lw=0.8,
            mec="black", mew=0.6,
        )
        handles.append(plt.Line2D(
            [0], [0], marker="o", linestyle="", markerfacecolor=c,
            markeredgecolor="black", markeredgewidth=0.6, markersize=7,
            label=f"{prettify_model(row['model'])} (n={row['n_seeds']})",
        ))

    ax.set_xlim(0, 1.0)
    ax.set_ylim(0, 1.0)
    ax.set_xlabel("Bio-conservation score ↑ (mean of bio metrics)")
    ax.set_ylabel("Batch-correction score ↑ (mean of batch metrics)")
    ax.set_title(f"{dataset} — {label_col}\nbio vs batch tradeoff")
    ax.grid(linestyle=":", linewidth=0.4, alpha=0.5)
    ax.set_aspect("equal", adjustable="box")
    ax.legend(
        handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0),
        frameon=False, fontsize=8, title="Model (sorted by composite ↑)",
        title_fontsize=8.5, handletextpad=0.4, labelspacing=0.4,
        borderaxespad=0.0,
    )
    fig.tight_layout()
    _save(fig, out_stem)


def resolve_label_cols(agg: pd.DataFrame, args: argparse.Namespace) -> list[str]:
    if args.label_cols:
        requested = [c.strip() for c in args.label_cols.split(",") if c.strip()]
    else:
        requested = list(DEFAULT_PARAMS[args.dataset]["label_cols"])
    available = set(agg["label_col"].unique())
    missing = [c for c in requested if c not in available]
    if missing:
        raise ValueError(
            f"label_cols {missing} not present in metrics CSV "
            f"(available: {sorted(available)})."
        )
    return requested


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Plots from metrics.csv.")
    parser.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    parser.add_argument("--input", type=Path,
                        help="Override path to metrics.csv "
                             "(default: outs/metrics/<dataset>/metrics.csv).")
    parser.add_argument("--label-cols", type=str,
                        help="Comma-separated subset of label_col values to plot.")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    matplotlib.rcParams.update(PLOT_RC)

    csv_path = args.input or REPO_ROOT / "outs" / "metrics" / args.dataset / "metrics.csv"
    if not csv_path.exists():
        raise FileNotFoundError(
            f"{csv_path} not found. Run scripts/embedding_metrics/metrics.py first."
        )
    df = pd.read_csv(csv_path)
    if df.empty:
        raise ValueError(f"{csv_path} is empty.")

    agg = aggregate_seeds(df)
    label_cols = resolve_label_cols(agg, args)

    out_root = REPO_ROOT / "outs" / "plots" / "metrics" / args.dataset
    print(f"[ok] dataset={args.dataset} labels={label_cols} "
          f"models={sorted(agg['model'].unique())}")

    for label_col in label_cols:
        agg_label = agg[agg["label_col"] == label_col].copy()
        if agg_label.empty:
            print(f"[warn] no rows for label_col={label_col!r}; skipping")
            continue
        out_dir = out_root / label_col

        targets = [out_dir / "heatmap.png", out_dir / "bars.png", out_dir / "tradeoff.png"]
        existing = [p for p in targets if p.exists()]
        if existing and not args.overwrite:
            raise FileExistsError(
                "Refusing to overwrite (pass --overwrite):\n  "
                + "\n  ".join(str(p) for p in existing)
            )
        draw_heatmap(agg_label, args.dataset, label_col, out_dir / "heatmap")
        draw_bars(agg_label, args.dataset, label_col, out_dir / "bars")
        draw_tradeoff(agg_label, args.dataset, label_col, out_dir / "tradeoff")
        print(f"[ok] wrote {out_dir}/{{heatmap,bars,tradeoff}}.{{png,pdf}}")


if __name__ == "__main__":
    main()
