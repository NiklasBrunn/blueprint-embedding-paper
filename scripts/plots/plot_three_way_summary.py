"""
Multi-arm summary: blueprint cohort(s) vs literature-grounded-freeform vs naive-freeform
on the five properties compared across the three arms.

Panels (2x3 grid):
  (a) Bio-conservation score          (scIB bio sub-score, boxplot, label_col-driven)
  (b) Batch-correction score          (scIB batch sub-score, boxplot)
  (c) Cross-seed gene Jaccard
  (d) Number of modules               (blueprint only, boxplot)
  (e) Number of axes / dims           (boxplot)
  (f) Citation density                (DOIs per run / script, boxplot)

Panels (a) and (b) are the two halves of the scIB composite (Luecken 2022):
each is the unweighted mean of the available metrics in its group, every metric
already scaled so higher = better. There is
no precomputed composite column in metrics.csv, so both are aggregated here.
Per scIB, graph connectivity is a *batch* metric and so belongs to panel (b).

Each panel shows every arm side-by-side for one freeform model tag and dataset.
The default is the classic three-way (one blueprint arm sharing the freeform
model tag); pass --blueprint-tags to compare several blueprint cohorts (e.g.
unrefined vs noref) against the same two freeform arms in one figure.

Inputs:
    work/<blueprint_tag>_<dataset>/<run>/summary.json
    work/freeform_<variant>_<model_tag>_<dataset>/<run>/.embed_status,
                                                       embed_freeform.py
    outs/embeddings/<dataset>/blueprint_<blueprint_tag>_seed*_summary.json
    outs/embeddings/<dataset>/freeform_<variant>_<model_tag>_<dataset>_seed*_status.json
    outs/metrics/<dataset>/metrics.csv

Outputs:
    outs/plots/freeform_baseline/three_way_summary_<model_tag>_<dataset>.{png,pdf}

Usage:
    # classic three-way (blueprint cohort tag == freeform model tag)
    python scripts/plots/plot_three_way_summary.py \
        --dataset mouse_brain --model-tag gpt_5_5_high --label-col subclass

    # several blueprint cohorts vs the two freeform arms
    python scripts/plots/plot_three_way_summary.py \
        --dataset human_pancreas --model-tag gpt_5_5_high \
        --blueprint-tags gpt_5_5_unrefined_high,gpt_5_5_noref_unrefined_high \
        --freeform-variants constrained,naive,noref
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import re
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
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation")
# Blueprint arms cycle through blues; the two freeform arms keep their fixed hues.
BLUEPRINT_COLORS = ("#1f77b4", "#08519c", "#3182bd", "#6baed6")
FREEFORM_COLOR = {"constrained": "#ff7f0e", "naive": "#d62728", "noref": "#8c564b"}
FREEFORM_LABEL = {"constrained": "FF-literature-grounded", "naive": "FF-naive", "noref": "FF-noref"}
DOI_PATTERN = re.compile(r"10\.\d{4,9}/[^\s\"',;)<>\[\]]+")

# scIB sub-score component metrics (all columns are pre-scaled higher = better).
# These are the metrics.csv columns belonging to each scIB group; graph
# connectivity is a batch metric per Luecken et al. (2022). The composite for a
# group is the unweighted row-wise mean over whichever of these columns the
# dataset actually populated (a dataset with no batch annotation contributes no
# batch columns, so panel (a) falls back to "no data").
# Mirrors SCIB_* in embedding_metrics/metrics.py: the row-wise mean over these
# columns equals the canonical bio_score / batch_score_scib the Benchmarker writes
# (kmeans NMI/ARI; cLISI in bio, PCR in batch).
BIO_METRICS = ("nmi", "ari", "asw_label", "clisi", "isolated_label_f1")
BATCH_METRICS = ("graph_connectivity", "asw_batch", "ilisi", "kbet_acceptance", "pcr")


# Engine / flag tokens for the blueprint label. Mirrors prettify_model in
# scripts/embedding_metrics/metrics.py so the legend matches the leaderboard,
# e.g. "BP GPT-5.5 U" / "BP GPT-5.5 N U".
_ENGINE_PRETTY = {"gpt_5_5": "GPT-5.5", "gpt_5_4_mini": "GPT-5.4-mini", "gpt_5_4": "GPT-5.4"}
_FLAG_PRETTY = {"noref": "N", "unrefined": "U", "refine": "R"}


def _arm_label(tag: str, model_tag: str) -> str:
    """Blueprint-arm label aligned with the leaderboard: 'BP <engine> <flags>',
    e.g. gpt_5_5_noref_unrefined_high -> 'BP GPT-5.5 N U'."""
    s = tag
    eng = ""
    for k, v in _ENGINE_PRETTY.items():
        if k in s:
            eng, s = v, s.replace(k, "_").strip("_")
            break
    flags = [_FLAG_PRETTY[p] for p in s.split("_") if p in _FLAG_PRETTY]
    return " ".join(["BP"] + ([eng] if eng else []) + flags)


def build_arms(model_tag: str, blueprint_tags: list[str],
               freeform_variants: list[str]) -> list[dict]:
    """Ordered arm specs: every blueprint cohort first, then one freeform arm
    per requested variant (all sharing the freeform model_tag). Each arm
    carries its key, kind, cohort tag, display label and colour."""
    arms: list[dict] = []
    for i, bt in enumerate(blueprint_tags):
        arms.append({
            "key": f"blueprint:{bt}",
            "kind": "blueprint",
            "tag": bt,
            "label": _arm_label(bt, model_tag),
            "color": BLUEPRINT_COLORS[i % len(BLUEPRINT_COLORS)],
        })
    for variant in freeform_variants:
        arms.append({
            "key": variant,
            "kind": "freeform",
            "variant": variant,
            "tag": model_tag,
            "label": FREEFORM_LABEL.get(variant, f"FF-{variant}"),
            "color": FREEFORM_COLOR.get(variant, "#7f7f7f"),
        })
    return arms


def cohort_path(arm: dict, dataset: str) -> Path:
    if arm["kind"] == "blueprint":
        return REPO_ROOT / "work" / f"{arm['tag']}_{dataset}"
    return REPO_ROOT / "work" / f"freeform_{arm['variant']}_{arm['tag']}_{dataset}"


def _ingest_summary_blueprint(data: dict, dois: list[int], dims: list[int],
                              mods: list[int], genes: list[set[str]]) -> None:
    """Append (n_dois_unique, n_axes, n_modules, gene set) from one
    summary.json's blueprint payload onto the running per-cohort accumulators."""
    bp = data.get("blueprint", {}) if isinstance(data, dict) else {}
    module_list = bp.get("modules", []) or []
    dois.append(int(bp.get("n_dois_unique", 0)))
    dims.append(int(bp.get("n_axes", 0)))
    mods.append(int(bp.get("n_modules", len(module_list))))
    gset: set[str] = set()
    for m in module_list:
        gset.update(g for g in m.get("genes", []) if isinstance(g, str))
    genes.append(gset)


def blueprint_stats(cohort_dir: Path, dataset: str | None = None) -> dict:
    """Blueprint per-run stats. Prefers the Embed summary.json mirror at
    outs/embeddings/<dataset>/blueprint_<cohort>_seed*_summary.json (written
    by agent_edit_blueprint{,_claude}.sh); falls back to walking the cohort's
    work/ runs for legacy cohorts predating the mirror convention.

    Same outs-first/work-fallback pattern as
    scripts/plots/umap_branch_interactive.py:find_blueprint_for_parquet."""
    cohort_name = cohort_dir.name
    mirrored: list[Path] = []
    if dataset:
        emb_dir = REPO_ROOT / "outs" / "embeddings" / dataset
        if emb_dir.is_dir():
            mirrored = sorted(emb_dir.glob(f"blueprint_{cohort_name}_seed*_summary.json"))

    dois: list[int] = []
    dims: list[int] = []
    mods: list[int] = []
    genes: list[set[str]] = []

    if mirrored:
        n_runs = n_fail = 0
        for sj in mirrored:
            n_runs += 1
            try:
                data = json.loads(sj.read_text())
            except Exception:
                n_fail += 1
                continue
            _ingest_summary_blueprint(data, dois, dims, mods, genes)
        return {"n_runs": n_runs, "n_failures": n_fail, "dois": dois,
                "dims": dims, "modules": mods, "genes": genes}

    if not cohort_dir.is_dir():
        return {"n_runs": 0, "n_failures": 0, "dois": [], "dims": [],
                "modules": [], "genes": []}
    n_runs = n_fail = 0
    for run in sorted(cohort_dir.glob("*_blueprint_run*")):
        sj = run / "summary.json"
        if not sj.is_file():
            continue
        n_runs += 1
        try:
            data = json.loads(sj.read_text())
        except Exception:
            n_fail += 1
            continue
        _ingest_summary_blueprint(data, dois, dims, mods, genes)
    return {"n_runs": n_runs, "n_failures": n_fail, "dois": dois,
            "dims": dims, "modules": mods, "genes": genes}


def _freeform_runs_from_outs(dataset: str, cohort_name: str) -> list[tuple[int, Path, Path | None]]:
    """List of (seed, status_json_path, script_path_or_None) for one
    freeform cohort, sourced from outs/. script_path is None when the
    .py wasn't mirrored (Stage A failure that left no script).

    Status mirror path: outs/embeddings/<dataset>/<cohort>_seed<N>_status.json
    Script mirror path: outs/freeform_scripts/<dataset>/<cohort>/embed_freeform_seed<N>.py
    Both written by Finalize of agent_edit_freeform.sh."""
    emb_dir = REPO_ROOT / "outs" / "embeddings" / dataset
    script_dir = REPO_ROOT / "outs" / "freeform_scripts" / dataset / cohort_name
    seed_re = re.compile(rf"^{re.escape(cohort_name)}_seed(\d+)_status\.json$")
    out: list[tuple[int, Path, Path | None]] = []
    if not emb_dir.is_dir():
        return out
    for sj in emb_dir.iterdir():
        m = seed_re.match(sj.name)
        if not m:
            continue
        seed = int(m.group(1))
        script = script_dir / f"embed_freeform_seed{seed}.py"
        out.append((seed, sj, script if script.is_file() else None))
    out.sort(key=lambda t: t[0])
    return out


def freeform_stats(cohort_dir: Path, dataset: str | None = None) -> dict:
    """Free-form per-run stats. Prefers outs/ mirrors written by Finalize of
    agent_edit_freeform.sh (status JSON + script archive); falls back to
    walking work/ for legacy cohorts predating the mirror convention.

    Same outs-first/work-fallback pattern as blueprint_stats() above."""
    cohort_name = cohort_dir.name
    mirrored = _freeform_runs_from_outs(dataset, cohort_name) if dataset else []

    if mirrored:
        n_runs = n_fail = 0
        dois: list[int] = []
        dims: list[int] = []
        genes: list[set[str]] = []
        for _seed, sj_path, script_path in mirrored:
            n_runs += 1
            try:
                status = json.loads(sj_path.read_text())
            except Exception:
                status = {}
            embed_status = status.get("embed_status", "")
            embed_detail = status.get("embed_detail", "")
            ok = embed_status == "OK"
            if not ok:
                n_fail += 1
            text = script_path.read_text(errors="replace") if script_path else ""
            dois.append(len(set(DOI_PATTERN.findall(text))))
            if ok and embed_detail:
                parquet = REPO_ROOT / embed_detail
                try:
                    import pyarrow.parquet as pq
                    names = set(pq.read_schema(parquet).names)
                    dims.append(sum(1 for n in names if n.startswith("dim_")))
                except Exception:
                    pass
            genes.append(set())  # populated downstream once vocabulary is known
        return {"n_runs": n_runs, "n_failures": n_fail, "dois": dois, "dims": dims,
                "modules": [], "genes": genes, "scripts": True}

    if not cohort_dir.is_dir():
        return {"n_runs": 0, "n_failures": 0, "dois": [], "dims": [],
                "modules": [], "genes": []}
    n_runs = n_fail = 0
    dois = []
    dims = []
    genes = []
    for run in sorted(cohort_dir.glob("*_freeform_run*")):
        n_runs += 1
        status = (run / ".embed_status").read_text().strip() if (run / ".embed_status").exists() else ""
        detail = (run / ".embed_detail").read_text().strip() if (run / ".embed_detail").exists() else ""
        ok = status == "OK"
        if not ok:
            n_fail += 1
        script = run / "embed_freeform.py"
        text = script.read_text(errors="replace") if script.is_file() else ""
        dois.append(len(set(DOI_PATTERN.findall(text))))
        if ok and detail:
            parquet = REPO_ROOT / detail
            try:
                import pyarrow.parquet as pq
                names = set(pq.read_schema(parquet).names)
                dims.append(sum(1 for n in names if n.startswith("dim_")))
            except Exception:
                pass
        genes.append(set())  # populated downstream once vocabulary is known
    return {"n_runs": n_runs, "n_failures": n_fail, "dois": dois, "dims": dims,
            "modules": [], "genes": genes, "scripts": True}


GENE_TOKEN_PATTERN = re.compile(r"\b[A-Z][A-Za-z0-9-]{1,15}\b")


def freeform_gene_sets(cohort_dir: Path, vocabulary: set[str],
                       dataset: str | None = None) -> list[set[str]]:
    """Gene tokens per run intersected with vocabulary. Prefers script mirrors
    under outs/freeform_scripts/<dataset>/<cohort>/; falls back to walking
    work/ for legacy cohorts. Order matches freeform_stats()."""
    cohort_name = cohort_dir.name
    mirrored = _freeform_runs_from_outs(dataset, cohort_name) if dataset else []
    out: list[set[str]] = []

    if mirrored:
        for _seed, _sj, script_path in mirrored:
            if script_path is None:
                out.append(set())
                continue
            try:
                text = script_path.read_text(errors="replace")
            except Exception:
                out.append(set())
                continue
            out.append(set(GENE_TOKEN_PATTERN.findall(text)) & vocabulary)
        return out

    for run in sorted(cohort_dir.glob("*_freeform_run*")):
        script = run / "embed_freeform.py"
        if not script.is_file():
            out.append(set())
            continue
        try:
            text = script.read_text(errors="replace")
        except Exception:
            out.append(set())
            continue
        out.append(set(GENE_TOKEN_PATTERN.findall(text)) & vocabulary)
    return out


def pairwise_jaccards(sets: list[set]) -> np.ndarray:
    """All pairwise gene Jaccard values across runs (the distribution behind
    the cross-seed reproducibility boxplot). Empty sets are dropped first."""
    sets = [s for s in sets if s]
    if len(sets) < 2:
        return np.array([])
    return np.array([len(a & b) / max(1, len(a | b))
                     for a, b in itertools.combinations(sets, 2)])


def mean_pairwise_jaccard(sets: list[set]) -> float:
    vals = pairwise_jaccards(sets)
    return float(vals.mean()) if vals.size else float("nan")


def collect(arms: list[dict], dataset: str) -> dict[str, dict]:
    """Per-arm stats keyed by arm['key']. The gene vocabulary used to score
    freeform gene reuse is the union of all blueprint arms' module genes."""
    out: dict[str, dict] = {}
    for arm in arms:
        if arm["kind"] == "blueprint":
            out[arm["key"]] = blueprint_stats(cohort_path(arm, dataset), dataset=dataset)
    vocab: set[str] = set()
    for arm in arms:
        if arm["kind"] == "blueprint":
            for g in out[arm["key"]]["genes"]:
                vocab |= g
    for arm in arms:
        if arm["kind"] != "freeform":
            continue
        cdir = cohort_path(arm, dataset)
        ff = freeform_stats(cdir, dataset=dataset)
        ff["genes"] = freeform_gene_sets(cdir, vocab, dataset=dataset) if vocab else []
        out[arm["key"]] = ff
    return out


def _arm_stem(arm: dict) -> str:
    return (f"blueprint_{arm['tag']}" if arm["kind"] == "blueprint"
            else f"freeform_{arm['variant']}_{arm['tag']}")


def composite_per_cohort(metrics_csv: Path, arms: list[dict], label_col: str,
                         metric_cols: tuple[str, ...]) -> dict[str, np.ndarray]:
    """Per-arm array of an scIB sub-score (one value per seed): the unweighted
    mean across the requested metric columns that the dataset actually
    populated. Columns absent from metrics.csv, or entirely NaN (e.g. batch
    metrics on a dataset with no batch annotation), are dropped; if none
    survive, the arm gets an empty array (panel falls back to "no data")."""
    if not metrics_csv.is_file():
        return {}
    df = pd.read_csv(metrics_csv)
    df = df[df["label_col"] == label_col]
    cols = [c for c in metric_cols if c in df.columns and df[c].notna().any()]
    out: dict[str, np.ndarray] = {}
    for arm in arms:
        sub = df[df["model"] == _arm_stem(arm)]
        if sub.empty or not cols:
            out[arm["key"]] = np.array([])
            continue
        # skipna=False: a seed is scored only if every surviving component is
        # present for it, so all arms average over the same metric panel.
        comp = sub[cols].astype(float).mean(axis=1, skipna=False).dropna()
        out[arm["key"]] = comp.to_numpy()
    return out


def _positions(arms: list[dict]) -> list[float]:
    return list(range(len(arms)))


def panel_jaccard(ax, arms: list[dict], stats: dict[str, dict]) -> None:
    """Boxplot of the per-pair gene Jaccard distribution across runs."""
    positions = _positions(arms)
    series = [pairwise_jaccards(stats[a["key"]]["genes"]) for a in arms]
    bp = ax.boxplot([s if s.size else np.array([np.nan]) for s in series],
                    positions=positions, widths=0.7, patch_artist=True, showmeans=True)
    for box, a in zip(bp["boxes"], arms):
        box.set_facecolor(a["color"])
        box.set_alpha(0.65)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Pairwise gene Jaccard")
    ax.set_title("(c) Cross-seed reproducibility")
    _annotate_axis(ax, arms)


def _box_panel(ax, arms: list[dict], stats: dict[str, dict], key: str,
               title: str, ylabel: str, log: bool = False) -> None:
    positions = _positions(arms)
    series = [np.asarray(stats[a["key"]][key], dtype=float) for a in arms]
    bp = ax.boxplot([s if s.size else np.array([np.nan]) for s in series],
                    positions=positions, widths=0.7, patch_artist=True, showmeans=True)
    for box, a in zip(bp["boxes"], arms):
        box.set_facecolor(a["color"])
        box.set_alpha(0.65)
    if log:
        ax.set_yscale("log")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    _annotate_axis(ax, arms)


def panel_score(ax, arms: list[dict], scores: dict[str, np.ndarray],
                title: str, ylabel: str) -> None:
    """Boxplot of an scIB sub-score (one point per seed), shared [0, 1] axis."""
    positions = _positions(arms)
    series = [scores.get(a["key"], np.array([])) for a in arms]
    if all(s.size == 0 for s in series):
        ax.text(0.5, 0.5, "no metrics rows", transform=ax.transAxes, ha="center")
        ax.set_title(f"{title} (no data)")
        _annotate_axis(ax, arms)
        return
    bp = ax.boxplot([s if s.size else np.array([np.nan]) for s in series],
                    positions=positions, widths=0.7, patch_artist=True, showmeans=True)
    for box, a in zip(bp["boxes"], arms):
        box.set_facecolor(a["color"])
        box.set_alpha(0.65)
    ax.set_ylim(0, 1)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    _annotate_axis(ax, arms)


def _annotate_axis(ax, arms: list[dict]) -> None:
    positions = _positions(arms)
    ax.set_xticks(positions)
    ax.set_xticklabels([a["label"] for a in arms], fontsize=10, rotation=20, ha="right")
    for tick, a in zip(ax.get_xticklabels(), arms):
        tick.set_color(a["color"])
        tick.set_fontweight("bold")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Blueprint vs freeform multi-arm summary.")
    p.add_argument("--dataset", required=True, choices=DATASET_CHOICES)
    p.add_argument("--model-tag", required=True,
                   help="Freeform model tag, e.g. gpt_5_5_high (model_slug + reasoning, joined by _).")
    p.add_argument("--blueprint-tags", default=None,
                   help="Comma-separated blueprint cohort tag(s), e.g. "
                        "gpt_5_5_unrefined_high,gpt_5_5_noref_unrefined_high. "
                        "Default: the --model-tag value (classic three-way).")
    p.add_argument("--freeform-variants", default="constrained,naive",
                   help="Comma-separated freeform variant(s) to include as arms "
                        "(constrained, naive, noref). Default: constrained,naive.")
    p.add_argument("--label-col", default=None,
                   help="Quality label column for panel (e). Default: subclass for "
                        "mouse_brain, celltype otherwise.")
    p.add_argument("--out-dir", type=Path,
                   default=REPO_ROOT / "outs" / "plots" / "freeform_baseline")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    label_col = args.label_col or ("subclass" if args.dataset == "mouse_brain" else "celltype")
    blueprint_tags = ([t.strip() for t in args.blueprint_tags.split(",") if t.strip()]
                      if args.blueprint_tags else [args.model_tag])
    freeform_variants = [v.strip() for v in args.freeform_variants.split(",") if v.strip()]
    arms = build_arms(args.model_tag, blueprint_tags, freeform_variants)

    stats = collect(arms, args.dataset)
    metrics_csv = REPO_ROOT / "outs" / "metrics" / args.dataset / "metrics.csv"
    bio = composite_per_cohort(metrics_csv, arms, label_col, BIO_METRICS)
    batch = composite_per_cohort(metrics_csv, arms, label_col, BATCH_METRICS)

    fig, axes = plt.subplots(2, 3, figsize=(1.8 * len(arms) + 7, 11.0))
    # Row 1 — the two scIB sub-scores + cross-seed reproducibility.
    panel_score(axes[0, 0], arms, bio,
                f"(a) Bio conservation — {label_col}", "Bio-conservation score")
    panel_score(axes[0, 1], arms, batch,
                "(b) Batch correction — scIB", "Batch-correction score")
    panel_jaccard(axes[0, 2], arms, stats)
    # Row 2 — blueprint structural size (modules / axes) + provenance.
    _box_panel(axes[1, 0], arms, stats, "modules",
               "(d) Number of modules", "n modules (blueprint only)")
    _box_panel(axes[1, 1], arms, stats, "dims",
               "(e) Number of axes", "n axes (BP) / dims (FF)")
    _box_panel(axes[1, 2], arms, stats, "dois",
               "(f) Citation density", "DOIs per run / script")

    handles = [plt.Rectangle((0, 0), 1, 1, color=a["color"], alpha=0.85) for a in arms]
    labels = [a["label"] for a in arms]
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.99),
               ncol=len(arms), frameon=False, fontsize=11)
    fig.tight_layout(rect=(0, 0.02, 1, 0.95))

    args.out_dir.mkdir(parents=True, exist_ok=True)
    out_stem = args.out_dir / f"three_way_summary_{args.model_tag}_{args.dataset}"
    for ext in ("png", "pdf"):
        path = out_stem.with_suffix(f".{ext}")
        fig.savefig(path, dpi=200, bbox_inches="tight")
        print(f"[ok] wrote {path}", flush=True)
    plt.close(fig)


if __name__ == "__main__":
    main()
