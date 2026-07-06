#!/usr/bin/env python3
"""Figures for the Cost / wall-clock / storage table.

Reads outs/resource_costs/<dataset>/costs.csv (produced by
scripts/resource_costs/collect.py) and joins outs/metrics/<dataset>/metrics.csv
for the quality axis. Writes to outs/plots/resource_costs/<dataset>/:

All the *time* axes EXCLUDE load_data_s (reading the input AnnData from disk) —
it's the same toll for every method, so it's noise in a comparison; the raw value
stays in costs.csv / costs_report.md. So "wall-clock" below means wall-clock − the
disk read (= obtain-model + build-embedding + other).

  1. tradeoff_total_vs_quality           — (left) bio-conservation score vs. total
     end-to-end wall-clock (disk read excluded); (right) vs. on-disk model
     footprint (log, humanised). All methods, identified by a side legend; marker
     shape distinguishes methods within a colour, with a [#] label keyed to the legend.
  2. tradeoff_obtain_model_vs_quality    — quality vs. the cost to *obtain a usable
     model* (log x; no disk read in this step anyway): de-novo VAEs / PCA = the
     fit/train; pretrained = checkpoint load + gene-space adaptation; agent methods
     = the LLM authoring (Stage A, + the optional Stage B refinement loop if
     --include-refine). Methods that acquire no model (UCell) and
     UCE (load+forward not separable) are auto-excluded (NaN x) and noted.
  3. tradeoff_build_embedding_vs_quality — quality vs. the cost to *construct the
     embedding for this dataset given the model*, disk read excluded (log x):
     trained models = the forward pass only (model already in memory); everyone
     else = the rest of the one-shot run.
  4. time_breakdown — one normalized-100% stacked bar per method (linear y; you
     can't stack on a log axis), segments = the three time buckets that sum to the
     run's wall-clock minus the disk read: obtain model · build embedding · other
     (the phase is shown by bar alpha; faint→solid). For the agent methods the bar
     spans Stage A (LLM authoring) [+ Stage B refinement if --include-refine] —
     these ARE the "obtain model" bucket — + the (disk-read-less) Embed run,
     which IS the "build embedding" bucket (a thin sliver; see agent_cost). The
     absolute total (= wall-clock − load_data, in seconds) is annotated above each
     bar; a T/L/—/A marker shows the model_source.
  5. memory_storage — three panels: (a) per-method RSS *trajectory* (4 faded→solid
     overlaid bars: after-load → after-model → after-embedding → peak; NOT
     stacked — peak RSS is a max, not a sum); (b) headline peak RSS (log y);
     (c) on-disk model footprint (log y, humanised — the "blueprint ≈ a few KB of
     script vs. multi-GB pretrained checkpoints" point; pretrained-checkpoint
     download is one-time and not shown as a per-run cost).
  6. agent_cost — agent cohorts only: token usage per run (Stage A by default;
     + Stage B refinement stacked when --include-refine) and wall-clock per cohort
     (Stage A | Embed build by default; Stage B refinement inserted between them
     when --include-refine).

Marker area ∝ embedding dimensionality; horizontal bars = per-seed IQR of the
cost; vertical bars = per-seed std of the quality metric. Style mirrors
scripts/plots/plot_metrics.py (DejaVu Sans 9 pt, no top/right spines, white
background, PNG + PDF at 600 dpi).

By default the optional Stage B refinement loop is *excluded* from every plot —
only Stage A (the populator) counts toward the blueprint methods' "obtain-model"
/ total costs, and Plot 6 (agent_cost) shows Stage A + Embed build only. Pass
--include-refine to add the Stage B refinement cost to every plot.

Usage:
    python scripts/plots/plot_resource_costs.py --dataset human_pancreas
    python scripts/plots/plot_resource_costs.py --dataset all --overwrite
    python scripts/plots/plot_resource_costs.py --dataset all --overwrite --include-refine
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib
matplotlib.use("Agg")

# Quiet the all-NaN/empty-slice noise from medians over not-yet-instrumented
# columns (the figures handle NaN explicitly).
warnings.filterwarnings("ignore", message="Mean of empty slice")
warnings.filterwarnings("ignore", message="All-NaN slice encountered")
warnings.filterwarnings("ignore", message="invalid value encountered in scalar divide")
import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from matplotlib.ticker import FuncFormatter

REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")
COSTS_ROOT = REPO_ROOT / "outs" / "resource_costs"
METRICS_ROOT = REPO_ROOT / "outs" / "metrics"
OUT_ROOT = REPO_ROOT / "outs" / "plots" / "resource_costs"

PRIMARY_LABEL = {"mouse_brain": "subclass", "human_pancreas": "celltype", "mouse_gastrulation": "celltype",
                 "mouse_cortex": "celltype_harmonized", "human_immune": "final_annotation"}
# Bio-conservation fallback set — mirrors SCIB_BIO_METRICS in
# embedding_metrics/metrics.py (graph_connectivity is a BATCH metric in scIB).
# Only used when the canonical per-row `bio_score` column is absent from
# metrics.csv (the Benchmarker normally writes it).
BIO_METRICS = ["ari", "nmi", "asw_label", "clisi", "isolated_label_f1"]

# Quality column on the cost-vs-quality y-axis. Default is the headline scIB
# composite (0.6*bio + 0.4*batch, the package Total; falls back to bio on batchless
# datasets). `--quality bio_score` reverts to biology-only. Set in main().
QUALITY_CHOICES = ("composite_scib", "bio_score")
_QUALITY = "composite_scib"

FAMILY_ORDER = ["conventional", "pretrained", "geneset", "blueprint", "freeform", "other"]
FAMILY_COLORS = {
    "conventional": "#1f77b4",  # blue   (data-driven baselines: PCA, scVI)
    "pretrained":   "#9467bd",  # purple (zero-shot foundation models: UCE, Census scVI)
    "geneset":      "#ff7f0e",  # orange (literature-informed gene-set methods: UCell, expiMap)
    "blueprint":    "#d62728",  # red    (agentic blueprint — the proposed method)
    "freeform":     "#8c564b",  # brown  (free-form agent baseline)
    "other":        "#888888",
}

# Distinct marker shapes so points of the same family colour are still
# distinguishable; cycled within each family. Paired with a numeric index drawn
# next to each point and echoed in the side legend.
MARKERS = ["o", "s", "^", "D", "v", "P", "X", "*", "h", "<", ">", "p"]

# All the cost plots exclude `load_data_s` (reading the input AnnData from disk):
# it's the same toll for every method, so it's noise in a comparison. The raw
# value stays in costs.csv / costs_report.md. ("Exclude" is exact for the model
# scripts + blueprint agents; for UCE the eval subprocess does its own internal
# read, which stays inside build_embedding.)
#
# The displayed time buckets (sum to walltime − load_data); rendered by bar alpha:
PHASE_KEYS    = ["obtain_model", "build_embedding", "other"]
PHASE_LABEL   = {"obtain_model": "obtain model  (de-novo: train · agent: Stage A — LLM authoring)",
                 "build_embedding": "build embedding  (forward/scoring · agent: the Embed run)",
                 "other": "other  (data prep · parquet write · glue)"}
PHASE_ALPHA   = {"obtain_model": 1.0, "build_embedding": 0.62, "other": 0.20}
# An agent bar's segments ARE the same buckets: Stage A [+ Stage B refinement when
# --include-refine] = "obtain model"; the Embed run (minus its read) = "build embedding".
AGENT_SEG_ALPHA = {"stage_a": PHASE_ALPHA["obtain_model"], "stage_b": PHASE_ALPHA["obtain_model"],
                   "embed": PHASE_ALPHA["build_embedding"]}
MODEL_SOURCE_MARK = {"trained": "T", "loaded": "L", "none": "—", "authored": "A", "?": "?", "": "?"}

PLOT_RC = {
    "font.family": "DejaVu Sans",
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.facecolor": "white", "savefig.transparent": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
}


# ----------------------------------------------------------------------- helpers
def _save(fig, out_stem: Path, dpi: int = 600) -> None:
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_stem.with_suffix(f".{ext}"), dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _human_bytes(x: float, _pos=None) -> str:
    if not np.isfinite(x) or x <= 0:
        return "0"
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if x < 1024 or unit == "TB":
            return f"{x:.0f}{unit}"
        x /= 1024.0
    return f"{x:.0f}TB"


# Method display labels aligned with the benchmark leaderboard. Baselines use their
# canonical casing; blueprints become 'BP <engine> <flags> [effort]' (engine map
# ordered so gpt_5_4_mini is matched before gpt_5_4; flags N=noref, U=unrefined,
# R=refined; non-high effort appended); free-form keeps its variant name.
_ENGINE_PRETTY = {"gpt_5_5": "GPT-5.5", "gpt_5_4_mini": "GPT-5.4-mini", "gpt_5_4": "GPT-5.4"}
_FLAG_PRETTY = {"noref": "N", "unrefined": "U", "refine": "R"}
_EFFORT_TOKENS = ("low", "medium", "high", "xhigh")
_BASELINE_PRETTY = {"scvi": "scVI", "scvi_nobatch": "scVI (no batch)",
                    "scvi_census": "scVI-census", "pca": "PCA", "uce": "UCE",
                    "ucell": "UCell", "expimap": "expiMap"}


def _short_method(m: str) -> str:
    if m in _BASELINE_PRETTY:
        return _BASELINE_PRETTY[m]
    if m.startswith("blueprint_"):
        s = m[len("blueprint_"):]
        eng = ""
        for k, v in _ENGINE_PRETTY.items():
            if k in s:
                eng, s = v, s.replace(k, "_").strip("_")
                break
        parts = [p for p in s.split("_") if p]
        flags = [_FLAG_PRETTY[p] for p in parts if p in _FLAG_PRETTY]
        effort = next((p for p in parts if p in _EFFORT_TOKENS and p != "high"), "")
        return " ".join(["BP"] + ([eng] if eng else []) + flags + ([effort] if effort else []))
    if m.startswith("freeform_"):
        v = m.split("_")[1]
        return "FF-" + ("literature-grounded" if v == "constrained" else v)
    return m


def _sb_tag(include_refine: bool) -> str:
    return "incl. Stage B refinement" if include_refine else "excl. Stage B (optional refinement loop)"


_NUM_COLS = ["walltime_s", "load_data_s", "obtain_model_s", "build_embedding_s", "other_s",
             "model_download_s", "agent_stage_a_s", "agent_stage_b_refine_s", "model_storage_bytes",
             "embedding_bytes", "embedding_dim", "peak_rss_mb", "rss_after_load_mb",
             "rss_after_obtain_mb", "rss_after_build_mb",
             "agent_tokens_stage_a", "agent_tokens_stage_b_refine", "agent_tokens_total"]


def _derive(df: pd.DataFrame, include_refine: bool = False) -> pd.DataFrame:
    """Add the analysis columns the figures use. All time axes EXCLUDE load_data_s
    (the shared disk-read toll — recorded raw in costs.csv): `walltime_xl_s` =
    walltime − load_data; `total_s` (end-to-end; for agents += Stage A [+ Stage B refine]);
    `obtain_s_plot` (cost to obtain a usable model — NaN where there is none — no
    disk read in this step anyway); `build_s_plot` (cost to construct the
    embedding given the model). Plus the four phase buckets coerced to numbers and
    the RSS trajectory (falling back to `peak` where not recorded)."""
    df = df.copy()
    for c in ("device", "host", "model_source"):
        if c not in df.columns:
            df[c] = ""
        df[c] = df[c].fillna("")
    for c in _NUM_COLS:
        if c not in df.columns:
            df[c] = np.nan
        df[c] = pd.to_numeric(df[c], errors="coerce")
    is_agent = df["family"].isin(["blueprint", "freeform"])
    sa = df["agent_stage_a_s"].fillna(0.0)
    sb = df["agent_stage_b_refine_s"].fillna(0.0) if include_refine else pd.Series(0.0, index=df.index)
    # the producing script's wall-clock minus the disk read (NaN load ⇒ unknown ⇒ 0)
    df["walltime_xl_s"] = df["walltime_s"].fillna(0.0) - df["load_data_s"].fillna(0.0)
    # 1) total end-to-end (sans disk read); for agents add the LLM authoring
    df["total_s"] = df["walltime_xl_s"]
    df.loc[is_agent, "total_s"] = df.loc[is_agent, "walltime_xl_s"] + sa[is_agent] + sb[is_agent]
    # 2) obtain-the-model cost (NaN ⇒ method has no model to obtain / not separable)
    obt = df["obtain_model_s"].copy()
    obt[df["model_source"] == "none"] = np.nan
    if include_refine:
        comb = df["agent_stage_a_s"].fillna(0.0) + df["agent_stage_b_refine_s"].fillna(0.0)
        comb[df["agent_stage_a_s"].isna() & df["agent_stage_b_refine_s"].isna()] = np.nan
        obt[is_agent] = comb[is_agent]
    else:
        obt[is_agent] = df.loc[is_agent, "agent_stage_a_s"]
    df["obtain_s_plot"] = pd.to_numeric(obt, errors="coerce")
    # 3) build-the-embedding cost given the model: trained ⇒ forward only (no disk
    #    read in it); everyone else ⇒ the whole one-shot run minus the disk read.
    bld = np.where(df["model_source"] == "trained", df["build_embedding_s"], df["walltime_xl_s"])
    df["build_s_plot"] = pd.to_numeric(pd.Series(bld, index=df.index), errors="coerce")
    df.loc[df["build_s_plot"].isna(), "build_s_plot"] = df.loc[df["build_s_plot"].isna(), "walltime_xl_s"]
    # RSS trajectory: fall back to peak where the snapshots weren't recorded.
    for c in ("rss_after_load_mb", "rss_after_obtain_mb", "rss_after_build_mb"):
        df[c] = df[c].where(df[c].notna(), df["peak_rss_mb"])
    return df


def _load_costs(dataset: str, include_refine: bool = False) -> pd.DataFrame:
    path = COSTS_ROOT / dataset / "costs.csv"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} not found — run scripts/resource_costs/collect.py --dataset {dataset} first.")
    return _derive(pd.read_csv(path), include_refine)


def _bio_score_by_method(dataset: str) -> pd.DataFrame:
    """Per-method bio-conservation score on the dataset's primary label, with std,
    from outs/metrics/<dataset>/metrics.csv. Prefers the canonical `bio_score`
    column (written by embedding_metrics/metrics.py); falls back to a mean over
    the scIB bio metrics for older CSVs."""
    path = METRICS_ROOT / dataset / "metrics.csv"
    if not path.exists():
        return pd.DataFrame(columns=["method", "bio_mean", "bio_std", "ari_mean", "n_metric_seeds"])
    md = pd.read_csv(path)
    lab = PRIMARY_LABEL.get(dataset)
    if "label_col" in md.columns and lab in set(md["label_col"].dropna()):
        md = md[md["label_col"] == lab]
    md = md.copy()
    # Quality column: the chosen `_QUALITY` (composite_scib by default), then
    # bio_score, then a mean over the scIB bio metrics (older CSVs).
    if _QUALITY in md.columns and md[_QUALITY].notna().any():
        md["bio"] = md[_QUALITY]
    elif "bio_score" in md.columns:
        md["bio"] = md["bio_score"]
    else:
        cols = [c for c in BIO_METRICS if c in md.columns]
        md["bio"] = md[cols].mean(axis=1, skipna=True) if cols else np.nan
    g = md.groupby("model").agg(
        bio_mean=("bio", "mean"), bio_std=("bio", "std"),
        ari_mean=("ari", "mean") if "ari" in md.columns else ("bio", "mean"),
        n_metric_seeds=("bio", "size"),
    ).reset_index().rename(columns={"model": "method"})
    return g


XCOLS = ("total_s", "obtain_s_plot", "build_s_plot")  # the three time-vs-quality x-axes


def _aggregate(df: pd.DataFrame, dataset: str) -> pd.DataFrame:
    """One row per method: medians (+ IQR for the three cost x-axes) + quality + device/source."""
    num = _NUM_COLS + list(XCOLS) + ["walltime_xl_s"]

    def _mode(s):
        s = s.dropna().astype(str)
        m = s[s != ""].mode()
        return str(m.iloc[0]) if len(m) else ""
    agg = df.groupby(["family", "method"], dropna=False).agg(
        n=("method", "size"),
        device=("device", _mode),
        model_source=("model_source", _mode),
        **{f"{c}_med": (c, "median") for c in num},
        **{f"{c}_lo": (c, lambda s: s.quantile(0.25)) for c in XCOLS},
        **{f"{c}_hi": (c, lambda s: s.quantile(0.75)) for c in XCOLS},
    ).reset_index()
    agg = agg.merge(_bio_score_by_method(dataset), on="method", how="left")
    agg["fam_order"] = agg["family"].map({f: i for i, f in enumerate(FAMILY_ORDER)}).fillna(99)
    return agg.sort_values(["fam_order", "method"]).reset_index(drop=True)


def _family_legend(families) -> list[Patch]:
    return [Patch(facecolor=FAMILY_COLORS.get(f, "#888888"), edgecolor="none", label=f)
            for f in FAMILY_ORDER if f in set(families)]


def _phase_legend_handles() -> list[Patch]:
    return [Patch(facecolor="#555555", alpha=PHASE_ALPHA[k], edgecolor="white", lw=0.4,
                  label=PHASE_LABEL[k]) for k in PHASE_KEYS]


def _msize(dim) -> float:
    if dim is None or not np.isfinite(dim) or dim <= 0:
        return 60.0
    return float(np.clip(30.0 + 55.0 * np.log10(dim), 35.0, 280.0))


# ------------------------------------------------------------------ scatter base
def _bio_ylabel(dataset: str) -> str:
    # Quality axis = the chosen `_QUALITY` column from metrics.csv.
    lab = PRIMARY_LABEL.get(dataset, "?")
    if _QUALITY == "composite_scib":
        return (f"scIB composite score ↑\n(0.6·bio + 0.4·batch, the package Total; "
                f"label = {lab})")
    return (f"Bio-conservation score ↑\n(mean of {len(BIO_METRICS)} scIB bio metrics: "
            f"ARI/NMI/ASW-label/cLISI/iso; label = {lab})")


def _ordered_methods(agg: pd.DataFrame) -> pd.DataFrame:
    """agg, sorted family-then-quality (NaN bio last) — the order points are
    numbered in and the side legend is listed in."""
    has_bio = agg["bio_mean"].notna().any()
    return agg.sort_values(["fam_order", ("bio_mean" if has_bio else "method")],
                           ascending=[True, not has_bio], na_position="last").reset_index(drop=True)


def _marker_map(agg: pd.DataFrame) -> dict[str, tuple[int, str, str]]:
    """{method: (index, marker_shape, colour)} — index runs over *all* methods in
    family-then-quality order (stable across panels & the legend); marker shape
    cycles within each family so same-coloured points stay distinguishable."""
    out: dict[str, tuple[int, str, str]] = {}
    per_fam: dict[str, int] = {}
    for i, (_, r) in enumerate(_ordered_methods(agg).iterrows(), start=1):
        fam = r["family"]
        k = per_fam.get(fam, 0)
        per_fam[fam] = k + 1
        out[r["method"]] = (i, MARKERS[k % len(MARKERS)], FAMILY_COLORS.get(fam, "#888888"))
    return out


def _scatter_panel(ax, agg: pd.DataFrame, dataset: str, *, xcol: str, xlabel: str, mmap: dict,
                   lo: str | None = None, hi: str | None = None, byte_x: bool = False) -> list[str]:
    has_bio = agg["bio_mean"].notna().any()
    drawn: list[str] = []
    for _, r in agg.iterrows():
        x = r.get(xcol)
        y = r["bio_mean"] if has_bio else r["n"]
        if x is None or not np.isfinite(x) or x <= 0 or not np.isfinite(y):
            continue
        idx, mk, c = mmap.get(r["method"], (0, "o", FAMILY_COLORS.get(r["family"], "#888888")))
        ax.scatter(x, y, s=_msize(r["embedding_dim_med"]), color=c, marker=mk, edgecolor="white",
                   linewidth=0.6, zorder=3, alpha=0.95)
        if has_bio and np.isfinite(r.get("bio_std", np.nan)) and r["bio_std"] > 0:
            ax.errorbar(x, y, yerr=r["bio_std"], color=c, alpha=0.45, capsize=2, lw=0.9, zorder=2)
        if lo and hi and np.isfinite(r.get(lo, np.nan)) and np.isfinite(r.get(hi, np.nan)) and r[hi] > r[lo] > 0:
            ax.plot([r[lo], r[hi]], [y, y], color=c, alpha=0.35, lw=0.9, zorder=2)
        t = ax.annotate(str(idx), (x, y), textcoords="offset points", xytext=(5, 4), fontsize=7.5,
                        fontweight="bold", color="#1a1a1a", zorder=5, annotation_clip=False)
        t.set_path_effects([pe.withStroke(linewidth=2.0, foreground="white")])
        drawn.append(r["method"])
    ax.set_xscale("log")
    if byte_x:
        ax.xaxis.set_major_formatter(FuncFormatter(_human_bytes))
    ax.set_xlabel(xlabel, fontsize=9)
    ax.set_ylabel(_bio_ylabel(dataset) if has_bio else "n runs", fontsize=10)
    ax.grid(True, which="both", axis="x", ls=":", lw=0.4, color="#cccccc", zorder=0)
    ax.grid(True, axis="y", ls=":", lw=0.4, color="#cccccc", zorder=0)
    return drawn


def _side_legend(ax, agg: pd.DataFrame, methods, mmap: dict, *,
                 title: str = "method  [#] = label on point") -> None:
    keep = set(methods)
    handles = []
    for _, r in _ordered_methods(agg).iterrows():
        if r["method"] not in keep:
            continue
        idx, mk, c = mmap.get(r["method"], (0, "o", FAMILY_COLORS.get(r["family"], "#888888")))
        n = int(r["n"]) if np.isfinite(r.get("n", np.nan)) else 0
        handles.append(plt.Line2D([0], [0], marker=mk, linestyle="", markerfacecolor=c,
                                  markeredgecolor="white", markeredgewidth=0.5, markersize=7.5,
                                  label=f"[{idx}] {_short_method(r['method'])} (n={n})"))
    if not handles:
        return
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.02, 1.0), frameon=False,
              fontsize=8, title=title, title_fontsize=8.5, labelspacing=0.4, handletextpad=0.5,
              borderaxespad=0.0)


# ------------------------------------------------------------------ figures 1-3
def fig_tradeoff_total(agg, dataset, out_stem, include_refine=False) -> None:
    mmap = _marker_map(agg)
    fig, axes = plt.subplots(1, 2, figsize=(14.0, 5.6), layout="constrained",
                             gridspec_kw={"width_ratios": [1, 1.25]})
    d1 = _scatter_panel(axes[0], agg, dataset, xcol="total_s_med", lo="total_s_lo", hi="total_s_hi", mmap=mmap,
                        xlabel="Total wall-clock to produce the embedding, disk read excluded (s, log)")
    d2 = _scatter_panel(axes[1], agg, dataset, xcol="model_storage_bytes_med", byte_x=True, mmap=mmap,
                        xlabel="Stored model footprint (log)")
    axes[0].set_title("(c) cost vs. quality", fontsize=11, loc="left")
    axes[1].set_title("(d) storage vs. quality", fontsize=11, loc="left")
    _side_legend(axes[1], agg, list(dict.fromkeys(d1 + d2)), mmap)
    _save(fig, out_stem)


def fig_tradeoff_obtain(agg, dataset, out_stem, include_refine=False) -> None:
    mmap = _marker_map(agg)
    fig, ax = plt.subplots(figsize=(8.8, 4.7), layout="constrained")
    agent_label = ("Stage A + Stage B refinement (LLM authors then iteratively refines the blueprint)" if include_refine
                   else "Stage A only (LLM authors the blueprint; the optional Stage B refinement loop is excluded)")
    drawn = _scatter_panel(ax, agg, dataset, xcol="obtain_s_plot_med", lo="obtain_s_plot_lo", hi="obtain_s_plot_hi", mmap=mmap,
                           xlabel="Cost to obtain a usable model (s, log) — no disk read in this step\n"
                                  f"de-novo VAEs / PCA: fit/train · pretrained: checkpoint load + gene-space adapt · agent: {agent_label}\n"
                                  "(not shown: UCE — load+forward not separable)")
    _side_legend(ax, agg, drawn, mmap)
    fig.suptitle(f"Model-acquisition cost vs. quality — {dataset}   (blueprint {_sb_tag(include_refine)}; "
                 "marker shape = method within a colour)", fontsize=10.0)
    _save(fig, out_stem)


def fig_tradeoff_build(agg, dataset, out_stem) -> None:
    mmap = _marker_map(agg)
    fig, ax = plt.subplots(figsize=(8.8, 4.7), layout="constrained")
    drawn = _scatter_panel(ax, agg, dataset, xcol="build_s_plot_med", lo="build_s_plot_lo", hi="build_s_plot_hi", mmap=mmap,
                           xlabel="Cost to construct the embedding for this dataset, given the model — disk read excluded (s, log)\n"
                                  "trained models (de-novo VAEs / PCA): forward pass only · "
                                  "pretrained / UCE: the rest of the one-shot run · agent: the Embed cell-scoring pass")
    _side_legend(ax, agg, drawn, mmap)
    fig.suptitle(f"Embedding-construction cost vs. quality — {dataset}   (disk read excluded · marker shape = method within a colour)", fontsize=10.0)
    _save(fig, out_stem)


# ------------------------------------------------------------------ figure 4
_PHASE_COL = {"obtain_model": "obtain_model_s_med", "build_embedding": "build_embedding_s_med",
              "other": "other_s_med"}


def _stack_segments(row, include_refine: bool):
    """Return [(label, alpha, seconds), ...] whose lengths sum to the method's
    wall-clock MINUS the disk read: the three phase buckets for local methods;
    Stage A [+ Stage B refine] + the (disk-read-less) Embed run for the agent methods."""
    def _v(col):
        try:
            x = float(row.get(col, np.nan))
        except (TypeError, ValueError):
            return 0.0
        return x if np.isfinite(x) else 0.0
    if row["family"] in ("blueprint", "freeform"):
        sa = _v("agent_stage_a_s_med")
        p5 = max(0.0, _v("walltime_s_med") - _v("load_data_s_med"))   # Embed run, sans its read
        segs = [("Stage A (LLM authoring)", AGENT_SEG_ALPHA["stage_a"], sa)]
        if include_refine:
            sb = _v("agent_stage_b_refine_s_med")
            if sb > 0:
                segs.append(("Stage B (LLM refinement)", AGENT_SEG_ALPHA["stage_b"], sb))
        segs.append(("Embed run", AGENT_SEG_ALPHA["embed"], p5))
        return segs
    # Local methods: bar = walltime − load_data; any bucket not recorded (stale
    # report) lands in `other` as a residual.
    tot = max(0.0, _v("walltime_s_med") - _v("load_data_s_med"))
    if tot <= 0:
        return []
    obt, bld = _v("obtain_model_s_med"), _v("build_embedding_s_med")
    vals = {"obtain_model": obt, "build_embedding": bld, "other": max(0.0, tot - obt - bld)}
    return [(PHASE_LABEL[k], PHASE_ALPHA[k], vals[k]) for k in PHASE_KEYS]


def fig_time_breakdown(agg: pd.DataFrame, dataset: str, out_stem: Path, include_refine: bool = False) -> None:
    methods = agg["method"].tolist()
    if not methods:
        return
    x = np.arange(len(methods))
    width = max(7.0, 0.62 * len(methods) + 2.0)
    fig, ax = plt.subplots(figsize=(width, 5.8), layout="constrained")
    for i, (_, r) in enumerate(agg.iterrows()):
        c = FAMILY_COLORS.get(r["family"], "#888888")
        segs = _stack_segments(r, include_refine)
        tot = sum(s for _, _, s in segs)
        if tot <= 0:
            continue
        bottom = 0.0
        for _, alpha, s in segs:
            if s <= 0:
                continue
            ax.bar(i, 100.0 * s / tot, 0.66, bottom=bottom, color=c, alpha=alpha, edgecolor="white", lw=0.4)
            bottom += 100.0 * s / tot
        # absolute total above the bar; model_source mark just under it
        lab = f"{tot:,.0f} s" if tot >= 1.0 else f"{tot:.2f} s"
        ax.text(i, 101.5, lab, ha="center", va="bottom", fontsize=6.6, rotation=0, color="#333333")
        ax.text(i, -3.0, MODEL_SOURCE_MARK.get(str(r.get("model_source") or ""), "?"),
                ha="center", va="top", fontsize=7.5, color=c, fontweight="bold")
    ax.set_ylim(0, 108)
    ax.set_ylabel("% of wall-clock  (disk read excluded)")
    ax.set_xticks(x)
    devs = [str(d) if (isinstance(d, str) and d) else "?" for d in agg["device"]]
    ax.set_xticklabels([f"{_short_method(m)}\n[{d}]" for m, d in zip(methods, devs)],
                       rotation=40, ha="right", fontsize=7)
    ax.set_axisbelow(True)
    ax.grid(True, axis="y", ls=":", lw=0.4, color="#cccccc")
    # Legend = the three buckets (alpha = phase). The agent segments ARE these
    # same buckets — Stage A [+ Stage B refinement] = "obtain model", the Embed run =
    # "build embedding" — so there's nothing extra to list. Placed to the right
    # (bars fill 0–100%, so there's no in-axes whitespace for it).
    ax.legend(handles=_phase_legend_handles(), fontsize=6.6,
              title="time bucket  (bar alpha = phase)", title_fontsize=7,
              frameon=False, loc="upper left", bbox_to_anchor=(1.005, 1.0), borderaxespad=0.0)
    fig.suptitle(f"Where the time goes, per method — {dataset}   (disk read excluded; blueprint {_sb_tag(include_refine)})", fontsize=10.0)
    fig.text(0.01, -0.02,
             "Each bar is normalized to 100% of that method's wall-clock MINUS the disk read (load_data_s — "
             "the same toll for everyone; recorded raw in costs.csv); the absolute total (seconds) printed above "
             "it is wall-clock − load_data.  Three buckets, same for every method: obtain model · build embedding · "
             "other.  For the agent methods, the **obtain-model** segment IS Stage A — the LLM authoring of the "
             "blueprint " + ("(+ Stage B, the LLM refinement loop — both shown as 'obtain model', split by the white "
             "edge) " if include_refine else "(the optional Stage B refinement loop is EXCLUDED — pass "
             "--include-refine; see agent_cost) ") + "— and the **build-embedding** segment is the "
             "Embed cell-scoring run (sans its read; a thin sliver here — its own breakdown is in agent_cost / "
             "the build-embedding tradeoff).  For de-novo VAEs / PCA, obtain-model = training/fit; "
             "for pretrained-load runs it's checkpoint load + gene-space adapt.  "
             "Mark under each bar: model_source — T=trained, L=loaded checkpoint, —=none, A=authored (LLM).  A "
             "faint full bar = the per-phase breakdown wasn't recorded for that run (it all falls in 'other'); "
             "re-run to populate.  Stacked bars use a linear axis (lengths must add up); absolute magnitudes vs. "
             "quality are in the three tradeoff scatters.",
             fontsize=6.4, color="#666666", va="top", ha="left", wrap=True)
    _save(fig, out_stem)


# ------------------------------------------------------------------ figure 5
def _bars(ax, x, heights, colors, *, log=True, width=0.8):
    h = np.nan_to_num(np.asarray(heights, dtype=float), nan=0.0)
    ax.bar(x, h, width, color=colors, alpha=0.95, edgecolor="white", lw=0.4)
    if log and (h > 0).any():
        ax.set_yscale("log")
    ax.grid(True, axis="y", which="both", ls=":", lw=0.4, color="#cccccc")


def fig_memory_storage(agg: pd.DataFrame, dataset: str, out_stem: Path, include_refine: bool = False) -> None:
    methods = agg["method"].tolist()
    if not methods:
        return
    x = np.arange(len(methods))
    colors = [FAMILY_COLORS.get(f, "#888888") for f in agg["family"]]
    width = max(7.0, 0.7 * len(methods) + 2.5)
    fig, axes = plt.subplots(3, 1, figsize=(width, 8.2), sharex=True, layout="constrained")

    # (a) RSS trajectory — 4 faded→solid overlaid bars per method (NOT stacked).
    ax = axes[0]
    traj = [("rss_after_load_mb_med", 0.30, "after load"), ("rss_after_obtain_mb_med", 0.50, "after model"),
            ("rss_after_build_mb_med", 0.72, "after embedding"), ("peak_rss_mb_med", 1.0, "peak")]
    for col, alpha, _ in traj:
        ax.bar(x, np.nan_to_num(agg[col].to_numpy(dtype=float), nan=0.0), 0.8,
               color=colors, alpha=alpha, edgecolor="white", lw=0.3)
    if (np.nan_to_num(agg["peak_rss_mb_med"].to_numpy(dtype=float)) > 0).any():
        ax.set_yscale("log")
    ax.grid(True, axis="y", which="both", ls=":", lw=0.4, color="#cccccc")
    ax.set_ylabel("RSS MB (log)")
    ax.set_title("(a) memory trajectory — overlaid bars, faint→solid: RSS after load · after model · "
                 "after embedding · peak  (a watermark, NOT a stack)", fontsize=7.6, pad=3, loc="left")
    ax.legend(handles=[Patch(facecolor="#555555", alpha=a, edgecolor="white", lw=0.3, label=lab)
                       for _, a, lab in traj], fontsize=6.4, frameon=False, loc="upper left", ncol=2,
              title="RSS milestone", title_fontsize=7)

    # (b) headline peak RSS.
    ax = axes[1]
    _bars(ax, x, agg["peak_rss_mb_med"].to_numpy(dtype=float), colors)
    ax.set_ylabel("peak RSS MB (log)")
    ax.set_title("(b) peak resident memory during the run  (\"fits on a 16 GB laptop?\")", fontsize=9, pad=3, loc="left")

    # (c) on-disk model footprint.
    ax = axes[2]
    _bars(ax, x, agg["model_storage_bytes_med"].to_numpy(dtype=float), colors)
    ax.yaxis.set_major_formatter(FuncFormatter(_human_bytes))
    ax.set_ylabel("stored model (log)")
    ax.set_title("(c) on-disk model footprint  (the populated .py script vs. trained / multi-GB pretrained checkpoints)",
                 fontsize=9, pad=3, loc="left")
    ax.set_xticks(x)
    devs = [str(d) if (isinstance(d, str) and d) else "?" for d in agg["device"]]
    ax.set_xticklabels([f"{_short_method(m)}\n[{d}]" for m, d in zip(methods, devs)],
                       rotation=40, ha="right", fontsize=7)

    fig.suptitle(f"Memory & storage footprint — {dataset}   (blueprint {_sb_tag(include_refine)}; "
                 "bar colour = method family)", fontsize=10.5)
    fig.text(0.01, -0.015,
             "(a) RSS is a high-water max, not a sum — these are watermarks at the phase boundaries, so they "
             "read as a trajectory (data read → model in memory → embedding built → peak); a single solid bar "
             "= the trajectory wasn't recorded for that run (re-run to populate).  (b) the peak number.  "
             "(c) for the agent methods this is the populated `.py` script (a few KB); for pretrained foundation "
             "models it's the `.pt` / checkpoint — additionally there is a one-time cold-cache *download* of "
             "that file (≈ the same size), not shown here as it isn't a per-run cost.",
             fontsize=6.4, color="#666666", va="top", ha="left", wrap=True)
    _save(fig, out_stem)


# ------------------------------------------------------------------ figure 6
def fig_agent_cost(df: pd.DataFrame, dataset: str, out_stem: Path, include_refine: bool = False) -> bool:
    a = df[df["family"].isin(["blueprint", "freeform"])].copy()
    if a.empty:
        return False
    cohorts = sorted(a["method"].unique())
    x = np.arange(len(cohorts))
    fig, axes = plt.subplots(1, 2, figsize=(max(10.0, 1.05 * len(cohorts) + 3.0), 5.6), layout="constrained")

    ax = axes[0]
    for i, m in enumerate(cohorts):
        sub = a[a["method"] == m]
        c = FAMILY_COLORS.get(sub["family"].iloc[0], "#888888")
        ma = sub["agent_tokens_stage_a"].dropna().median() if sub["agent_tokens_stage_a"].notna().any() else 0.0
        if include_refine:
            mb = sub["agent_tokens_stage_b_refine"].dropna().median() if sub["agent_tokens_stage_b_refine"].notna().any() else 0.0
            ax.bar(i, ma or 0.0, color=c, alpha=0.6, edgecolor="white", lw=0.4)
            ax.bar(i, mb or 0.0, bottom=ma or 0.0, color=c, alpha=0.95, edgecolor="white", lw=0.4)
            dots = sub["agent_tokens_total"].dropna()
        else:
            ax.bar(i, ma or 0.0, color=c, alpha=0.95, edgecolor="white", lw=0.4)
            dots = sub["agent_tokens_stage_a"].dropna()
        if len(dots):
            jit = (np.random.RandomState(i).rand(len(dots)) - 0.5) * 0.28
            ax.scatter(np.full(len(dots), i) + jit, dots, s=10, color="#222222", alpha=0.5, zorder=3)
    ax.set_ylabel("tokens / run")
    ax.set_title("(a) LLM tokens per pipeline", fontsize=11, loc="left")
    ax.set_xticks(x)
    ax.set_xticklabels([_short_method(m) for m in cohorts], rotation=40, ha="right", fontsize=9)
    ax.grid(True, axis="y", ls=":", lw=0.4, color="#cccccc")
    tokens_col = "agent_tokens_total" if include_refine else "agent_tokens_stage_a"
    if not a[tokens_col].notna().any():
        ax.text(0.5, 0.5, "no token data\n(claude text-mode emits no usage)", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="#999999")

    ax = axes[1]
    for i, m in enumerate(cohorts):
        sub = a[a["method"] == m]
        c = FAMILY_COLORS.get(sub["family"].iloc[0], "#888888")
        sa = sub["agent_stage_a_s"].median(skipna=True)
        sb = sub["agent_stage_b_refine_s"].median(skipna=True)
        p5 = sub["walltime_s"].median(skipna=True)
        sa, sb, p5 = (0.0 if not np.isfinite(v) else v for v in (sa, sb, p5))
        if include_refine:
            ax.bar(i, sa, color=c, alpha=0.55, edgecolor="white", lw=0.4)
            ax.bar(i, sb, bottom=sa, color=c, alpha=0.78, edgecolor="white", lw=0.4)
            ax.bar(i, p5, bottom=sa + sb, color=c, alpha=1.0, edgecolor="white", lw=0.4)
        else:
            ax.bar(i, sa, color=c, alpha=0.55, edgecolor="white", lw=0.4)
            ax.bar(i, p5, bottom=sa, color=c, alpha=1.0, edgecolor="white", lw=0.4)
        # Per-run dots (mirror panel a): total wall-clock per run, plotted for runs
        # with the authoring + build times recorded. Stage B is added when present
        # (NaN -> 0) so non-refine cohorts still get dots rather than being dropped
        # (they simply have no Stage B column).
        sb_run = sub["agent_stage_b_refine_s"].fillna(0.0) if include_refine else 0.0
        run_mask = sub["agent_stage_a_s"].notna() & sub["walltime_s"].notna()
        dots_b = (sub["agent_stage_a_s"] + sub["walltime_s"] + sb_run)[run_mask]
        if len(dots_b):
            jit = (np.random.RandomState(i).rand(len(dots_b)) - 0.5) * 0.28
            ax.scatter(np.full(len(dots_b), i) + jit, dots_b.to_numpy(),
                       s=10, color="#222222", alpha=0.5, zorder=3)
    ax.set_ylabel("wall-clock / run (s)")
    ax.set_title("(b) Wall-clock per pipeline", fontsize=11, loc="left")
    ax.set_xticks(x)
    ax.set_xticklabels([_short_method(m) for m in cohorts], rotation=40, ha="right", fontsize=9)
    ax.grid(True, axis="y", ls=":", lw=0.4, color="#cccccc")

    _save(fig, out_stem)
    return True


# --------------------------------------------------------------------------- main
_FIGNAMES = ("tradeoff_total_vs_quality", "tradeoff_obtain_model_vs_quality",
             "tradeoff_build_embedding_vs_quality", "time_breakdown", "memory_storage", "agent_cost")
_STALE = ("tradeoff_quality_vs_cost", "tradeoff_construct_embedding_vs_quality", "footprint_bars")


def run_dataset(dataset: str, overwrite: bool, include_refine: bool = False) -> None:
    df = _load_costs(dataset, include_refine)
    if df.empty:
        print(f"[warn] {dataset}: costs.csv is empty; nothing to plot.")
        return
    agg = _aggregate(df, dataset)
    out_dir = OUT_ROOT / dataset
    # When Stage B is included, figures get an `_incl_refine` suffix so they
    # sit alongside the default Stage-B-excluded set rather than replacing it.
    sfx = "_incl_refine" if include_refine else ""
    if out_dir.exists() and any(out_dir.glob(f"*{sfx}.png")) and not overwrite:
        raise FileExistsError(f"{out_dir} already has figures; pass --overwrite to replace them.")
    if not include_refine:
        for stem in _STALE:                              # superseded figure names
            (out_dir / f"{stem}.png").unlink(missing_ok=True)
            (out_dir / f"{stem}.pdf").unlink(missing_ok=True)
    fig_tradeoff_total(agg, dataset, out_dir / f"tradeoff_total_vs_quality{sfx}", include_refine)
    fig_tradeoff_obtain(agg, dataset, out_dir / f"tradeoff_obtain_model_vs_quality{sfx}", include_refine)
    fig_tradeoff_build(agg, dataset, out_dir / f"tradeoff_build_embedding_vs_quality{sfx}")
    fig_time_breakdown(agg, dataset, out_dir / f"time_breakdown{sfx}", include_refine)
    fig_memory_storage(agg, dataset, out_dir / f"memory_storage{sfx}", include_refine)
    drew = fig_agent_cost(df, dataset, out_dir / f"agent_cost{sfx}", include_refine)
    names = ", ".join(f"{n}{sfx}" for n in (_FIGNAMES if drew else _FIGNAMES[:-1])) + ("" if drew else " (no agent cohorts → skipped agent_cost)")
    try:
        loc = out_dir.relative_to(REPO_ROOT)
    except ValueError:
        loc = out_dir
    print(f"[ok] {dataset}: wrote {names} → {loc}/")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dataset", required=True, help=f"One of {DATASET_CHOICES}, or 'all'.")
    p.add_argument("--overwrite", action="store_true", help="Replace existing figures.")
    p.add_argument("--quality", choices=QUALITY_CHOICES, default="composite_scib",
                   help="Quality column on the cost-vs-quality y-axis "
                        "(default: composite_scib = the scIB Total; bio_score = biology only).")
    p.add_argument("--include-refine", action="store_true",
                   help="Include the optional Stage B refinement-loop cost in the blueprint methods' "
                        "obtain-model / total / time-breakdown (default: excluded; the agent_cost plot "
                        "always shows its full breakdown either way).")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    if args.dataset != "all" and args.dataset not in DATASET_CHOICES:
        print(f"--dataset must be one of {DATASET_CHOICES} or 'all' (got {args.dataset!r})", file=sys.stderr)
        sys.exit(1)
    matplotlib.rcParams.update(PLOT_RC)
    global _QUALITY
    _QUALITY = args.quality
    datasets = list(DATASET_CHOICES) if args.dataset == "all" else [args.dataset]
    for ds in datasets:
        try:
            run_dataset(ds, args.overwrite, args.include_refine)
        except FileNotFoundError as e:
            print(f"[skip] {ds}: {e}")


if __name__ == "__main__":
    main()
