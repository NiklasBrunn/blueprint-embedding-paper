#!/usr/bin/env python3
"""Aggregate the Cost / wall-clock / storage table across every embedding method.

For a dataset, builds one long CSV (`outs/resource_costs/<dataset>/costs.csv`)
with one row per (method, seed-or-run), pulling from two sources:

  - Conventional / pretrained / gene-set models: the `## Resources` block written
    by each `scripts/models/<m>.py` into `outs/embeddings/<dataset>/<stem>_report.md`
    — walltime, load / obtain-model / build-embedding times, peak RSS, model
    storage, embedding bytes, embedding dim (legacy `## Stats` fit/transform
    reports are still read as a fallback). (`blueprint_*` / `freeform_*` parquets
    that also live under outs/embeddings/ are handled by the agent path instead.)
  - Agent runs (blueprint Stage A + the optional Stage B refinement loop, and
    the free-form baseline): `scripts/agentic_embedding/extract_cost.py`
    introspects each `work/<cohort>/<run>/` directory — token usage greped from
    the codex/claude transcripts, per-phase wall-clock (Stage A / Stage B
    refinement / Embed embedding-build), and on-disk footprint (the populated
    script = "the model", the output parquet).

Also writes a per-method summary (`costs_summary.csv`, n + median + IQR of every
numeric column) and a display-ready markdown table (`costs_report.md`). With
`--dataset all`, additionally writes `outs/resource_costs/costs_all.csv` and a
combined `costs_table.md`.

Usage:
    python scripts/resource_costs/collect.py --dataset human_pancreas
    python scripts/resource_costs/collect.py --dataset all --overwrite
"""
from __future__ import annotations

import argparse
import importlib.util
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
_THIS_SCRIPT_REPO = REPO_ROOT  # fixed; tests may rebind the data roots below.
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")
EMB_ROOT = REPO_ROOT / "outs" / "embeddings"
WORK_ROOT = REPO_ROOT / "work"
OUT_ROOT = REPO_ROOT / "outs" / "resource_costs"
METRICS_ROOT = REPO_ROOT / "outs" / "metrics"


def _rel(p: Path) -> Path:
    try:
        return p.relative_to(REPO_ROOT)
    except ValueError:
        return p

COLUMNS = [
    "method", "family", "dataset", "seed", "run_id", "cohort",
    "agent_model", "reasoning", "variant", "embedding_dim",
    # Where it ran:
    "device", "host",
    # Time taxonomy — the four buckets sum to walltime_s for every producing
    # script we own: load_data_s = read the input AnnData from disk; obtain_model_s
    # = get a usable model (train it / load a persisted-or-pretrained checkpoint +
    # adapt to the gene space — `model_source` tags which: trained / loaded / none
    # / authored; n/a/NaN when load+forward are not separable, e.g. UCE);
    # build_embedding_s = the forward pass / scoring; other_s = the residual (data
    # prep, parquet write, glue). For the agent methods this is the Embed build —
    # the *model-acquisition* cost is the LLM authoring (agent_stage_a_s, a
    # separate phase). model_download_s = one-time cold-cache fetch of a pretrained
    # checkpoint (0 when cached); never part of walltime in steady state.
    "walltime_s", "load_data_s", "obtain_model_s", "model_source", "build_embedding_s", "other_s",
    "model_download_s",
    # Memory: peak high-water RSS + a 4-point trajectory (after load / after model
    # / after embedding / peak). Memory is a max, not a sum — not decomposed.
    "peak_rss_mb", "rss_after_load_mb", "rss_after_obtain_mb", "rss_after_build_mb",
    # On-disk footprint:
    "model_storage_bytes", "embedding_bytes",
    # Agent-only: wall-clock + token cost of *authoring* the model (Stage A) and
    # the optional iterative refinement loop (Stage B, opt-in via --with-refine).
    "agent_stage_a_s", "agent_stage_b_refine_s",
    "agent_stage_b_iters", "agent_stage_b_final_rating",
    "agent_tokens_stage_a", "agent_tokens_stage_b_refine", "agent_tokens_total",
    # Provenance:
    "status", "source",
]
NUMERIC_COLS = [
    "embedding_dim", "walltime_s", "load_data_s", "obtain_model_s", "build_embedding_s",
    "other_s", "model_download_s", "peak_rss_mb", "rss_after_load_mb", "rss_after_obtain_mb",
    "rss_after_build_mb", "model_storage_bytes", "embedding_bytes",
    "agent_stage_a_s", "agent_stage_b_refine_s", "agent_stage_b_iters",
    "agent_tokens_stage_a", "agent_tokens_stage_b_refine", "agent_tokens_total",
]

# Pretrained foundation models: a report that predates the new `## Resources`
# block recorded their checkpoint-load-and-adapt time under "fit"; for those we
# infer model_source="loaded" (the new reports tag it explicitly).
_PRETRAINED = frozenset({"scvi_census", "uce"})
# De-novo parametric models we actually train on the data (PCA fits loadings too).
DE_NOVO = frozenset({"scvi", "scvi_nobatch", "expimap", "pca"})
# Non-parametric methods that acquire no model (UCell just reads a gene-set .gmt).
_NO_MODEL = frozenset({"ucell"})


def _infer_model_source(method: str) -> str:
    if method in _PRETRAINED:
        return "loaded"
    if method in DE_NOVO:
        return "trained"
    if method in _NO_MODEL:
        return "none"
    if method.startswith(("blueprint", "freeform")):
        return "authored"
    return "?"

# Method family for grouping / colour in plots; written as the `family` column
# of costs.csv and consumed by scripts/plots/plot_resource_costs.py.
_FAMILY = {
    "pca": "conventional",
    "scvi": "conventional", "scvi_nobatch": "conventional",
    "uce": "pretrained", "scvi_census": "pretrained",
    "ucell": "geneset", "expimap": "geneset",
}


def family_of(method: str) -> str:
    if method in _FAMILY:
        return _FAMILY[method]
    if method.startswith("blueprint"):
        return "blueprint"
    if method.startswith("freeform"):
        return "freeform"
    return "other"


# Cohort dir names.  Blueprint:  <model_slug>_<reasoning>_<dataset>
#                                (model_slug may carry `_userctx` for claude
#                                 no-isolate and/or `_refine` for --with-refine;
#                                 the non-greedy model group absorbs both, so
#                                 refined runs aggregate as a distinct cohort)
#                  Free-form:    freeform_<variant>_<model_slug>_<reasoning>_<dataset>
_R = r"(?:low|medium|high|xhigh)"
COHORT_RE = re.compile(rf"^(?P<model>.+?)_(?P<reasoning>{_R})_(?P<dataset>.+)$")
FREEFORM_COHORT_RE = re.compile(
    rf"^freeform_(?P<variant>naive|constrained|noref)_(?P<model>.+?)_(?P<reasoning>{_R})_(?P<dataset>.+)$"
)


def parse_cohort_name(name: str) -> dict | None:
    m = FREEFORM_COHORT_RE.match(name)
    if m:
        d = m.groupdict()
        return {
            "agent_model": d["model"], "reasoning": d["reasoning"],
            "dataset": d["dataset"], "variant": f"freeform_{d['variant']}",
            "method": f"freeform_{d['variant']}_{d['model']}_{d['reasoning']}",
        }
    m = COHORT_RE.match(name)
    if m:
        d = m.groupdict()
        return {
            "agent_model": d["model"], "reasoning": d["reasoning"],
            "dataset": d["dataset"], "variant": "blueprint",
            "method": f"blueprint_{d['model']}_{d['reasoning']}",
        }
    return None


# New `## Resources` block (scripts/models/_resource_meta.resource_lines). The
# parenthetical "what" labels can themselves contain parens (e.g. UCell's
# "...compute_ucell_scores (Mann-Whitney rank scoring)..."), so match up to the
# *value's* colon with `.*?:` rather than a non-nesting `\(...\)`.
RE_WALLTIME       = re.compile(r"^-\s*walltime:\s*([\d.]+)\s*s", re.M)
RE_LOAD_DATA      = re.compile(r"^-\s*load \(read input\b.*?:\s*(\S+)", re.M)
RE_MODEL_SOURCE   = re.compile(r"^-\s*obtain model — .*?\[(\w+)\]:", re.M)
RE_OBTAIN         = re.compile(r"^-\s*obtain model — .*?\[\w+\]:\s*(\S+)", re.M)
RE_BUILD          = re.compile(r"^-\s*build embedding\b.*?:\s*(\S+)", re.M)
RE_OTHER          = re.compile(r"^-\s*other \(.*?:\s*(\S+)", re.M)
RE_RSS_TRAJ       = re.compile(r"after-load ([\d.]+)\s*→\s*after-model ([\d.]+)\s*→\s*after-embedding ([\d.]+)\s*MB")
RE_MODEL_DOWNLOAD = re.compile(r"^-\s*model download:\s*([\d.]+)\s*s", re.M)
RE_PEAK_RSS       = re.compile(r"^-\s*peak RSS:\s*([\d.]+)\s*MB", re.M)
RE_MODEL_STORAGE  = re.compile(r"^-\s*model storage[^:]*:\s*(\d+)\s*bytes", re.M)
RE_EMB_BYTES      = re.compile(r"^-\s*embedding parquet:\s*(\d+)\s*bytes", re.M)
RE_EMB_DIM        = re.compile(r"embedding dim[^:]*:\s*(\d+)", re.M)  # also matches UCell's "(= embedding dim): N"
RE_DEVICE         = re.compile(r"^-\s*device:\s*(\S+)", re.M)
RE_HOST           = re.compile(r"^-\s*host:\s*(\S+)", re.M)
# Legacy `## Stats` block (reports that predate the taxonomy) — fallback only:
RE_FIT_OLD        = re.compile(r"^-\s*fit \([^)]*\):\s*(\S+)", re.M)
RE_TRANSFORM_OLD  = re.compile(r"^-\s*transform \([^)]*\):\s*(\S+)", re.M)


def _grp(rx: re.Pattern, text: str):
    m = rx.search(text)
    if not m:
        return None
    g = m.group(1)
    if g.lower() in ("n/a", "nan", "none"):
        return None
    try:
        return float(g) if "." in g else int(g)
    except ValueError:
        return None


def _grp_str(rx: re.Pattern, text: str):
    m = rx.search(text)
    g = m.group(1).strip("`'\"") if m else None
    return None if (g is None or g.lower() in ("n/a", "nan", "none", "?", "")) else g


def _parse_stem(stem: str, dataset: str) -> tuple[str, int | None] | None:
    """Recover (method, seed) from '<method>_<dataset>[_seed<N>]'."""
    if f"_{dataset}" not in stem:
        return None
    head, _, tail = stem.partition(f"_{dataset}")
    if not head:
        return None
    if not tail:
        return head, None
    if tail.startswith("_seed"):
        try:
            return head, int(tail[len("_seed"):])
        except ValueError:
            return head, None
    return head + tail, None


def conventional_rows(dataset: str) -> list[dict]:
    emb_dir = EMB_ROOT / dataset
    if not emb_dir.exists():
        return []
    rows: list[dict] = []
    for report in sorted(emb_dir.glob("*_report.md")):
        stem = report.name[: -len("_report.md")]
        if stem.startswith(("blueprint_", "freeform_")):
            continue  # agent embeddings — handled by agent_rows()
        parsed = _parse_stem(stem, dataset)
        if parsed is None:
            continue
        method, seed = parsed
        text = report.read_text(errors="ignore")
        if RE_MODEL_SOURCE.search(text):                 # new `## Resources` block
            model_source = _grp_str(RE_MODEL_SOURCE, text)
            obtain_s = _grp(RE_OBTAIN, text)             # None when "n/a"
            load_data_s, build_s, other_s = _grp(RE_LOAD_DATA, text), _grp(RE_BUILD, text), _grp(RE_OTHER, text)
            model_dl_s = _grp(RE_MODEL_DOWNLOAD, text)
            traj = RE_RSS_TRAJ.search(text)
            rss_load, rss_obt, rss_build = (tuple(float(g) for g in traj.groups()) if traj else (None, None, None))
        else:                                            # legacy report — only fit/transform/walltime
            model_source = _infer_model_source(method)
            obtain_s = _grp(RE_FIT_OLD, text)
            build_s = _grp(RE_TRANSFORM_OLD, text)
            load_data_s = other_s = model_dl_s = None
            rss_load = rss_obt = rss_build = None
        rows.append({
            "method": method, "family": family_of(method), "dataset": dataset,
            "seed": seed, "run_id": stem, "cohort": None,
            "agent_model": None, "reasoning": None, "variant": "conventional",
            "embedding_dim": _grp(RE_EMB_DIM, text),
            "device": _grp_str(RE_DEVICE, text), "host": _grp_str(RE_HOST, text),
            "walltime_s": _grp(RE_WALLTIME, text),
            "load_data_s": load_data_s, "obtain_model_s": obtain_s, "model_source": model_source,
            "build_embedding_s": build_s, "other_s": other_s, "model_download_s": model_dl_s,
            "peak_rss_mb": _grp(RE_PEAK_RSS, text),
            "rss_after_load_mb": rss_load, "rss_after_obtain_mb": rss_obt, "rss_after_build_mb": rss_build,
            "model_storage_bytes": _grp(RE_MODEL_STORAGE, text),
            "embedding_bytes": _grp(RE_EMB_BYTES, text),
            "agent_stage_a_s": None, "agent_stage_b_refine_s": None,
            "agent_stage_b_iters": None, "agent_stage_b_final_rating": None,
            "agent_tokens_stage_a": None, "agent_tokens_stage_b_refine": None,
            "agent_tokens_total": None,
            "status": "ok", "source": "report",
        })
    return rows


def _load_extract_cost():
    path = _THIS_SCRIPT_REPO / "scripts" / "agentic_embedding" / "extract_cost.py"
    spec = importlib.util.spec_from_file_location("extract_cost", path)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


def _refine_token_total(cost: dict) -> int | None:
    """Stage A + Stage B refinement token total. None when neither phase
    reported tokens (e.g. the claude text-mode track)."""
    ta = cost.get("agent_tokens_stage_a")
    tb = cost.get("agent_tokens_stage_b_refine")
    if ta is None and tb is None:
        return None
    return (ta or 0) + (tb or 0)


def agent_rows(dataset: str) -> list[dict]:
    if not WORK_ROOT.exists():
        return []
    ec = _load_extract_cost()
    rows: list[dict] = []
    for cohort_dir in sorted(WORK_ROOT.iterdir()):
        if not cohort_dir.is_dir() or not cohort_dir.name.endswith(f"_{dataset}"):
            continue
        # Exploratory local-model runs (OpenWebUI/Ollama/UFR) are excluded from
        # the aggregation for now; drop the guard to fold them back in.
        if "openwebui" in cohort_dir.name:
            continue
        meta = parse_cohort_name(cohort_dir.name)
        if meta is None:
            continue
        for run_dir in sorted(p for p in cohort_dir.iterdir()
                              if p.is_dir() and "_run" in p.name):
            c = ec._build_cost(run_dir)
            m = re.search(r"_run(\d+)$", run_dir.name)
            seed = int(m.group(1)) if m else None
            # Embed wall-clock: prefer the embedder's own [timing] line (whose
            # four buckets sum to it exactly); fall back to the orchestrator's.
            walltime = c.get("embed_walltime_s")
            if walltime is None:
                walltime = c.get("embed_s")
            rows.append({
                "method": meta["method"], "family": family_of(meta["method"]), "dataset": dataset,
                "seed": seed, "run_id": c["run_id"], "cohort": cohort_dir.name,
                "agent_model": meta["agent_model"], "reasoning": meta["reasoning"],
                "variant": meta["variant"],
                "embedding_dim": None,  # filled below from metrics.csv if available
                # Local Embed device (cpu); Stage A / Stage B run on the LLM provider's infra.
                "device": c.get("embed_device") or "cpu",
                "host": c.get("embed_host"),
                "walltime_s": walltime,
                # Embed buckets: load_data (read the input), obtain_model
                # (build_blueprint — ~ms; the *real* model-acquisition is Stage A,
                # see agent_stage_a_s), build_embedding (score every cell), other.
                "load_data_s": c.get("embed_load_data_s"),
                "obtain_model_s": c.get("embed_obtain_model_s"),
                "model_source": c.get("embed_model_source") or "authored",
                "build_embedding_s": c.get("embed_build_embedding_s"),
                "other_s": c.get("embed_other_s"),
                "model_download_s": c.get("embed_model_download_s", 0.0),
                "peak_rss_mb": c.get("embed_peak_rss_mb"),
                "rss_after_load_mb": c.get("embed_rss_after_load_mb"),
                "rss_after_obtain_mb": c.get("embed_rss_after_obtain_mb"),
                "rss_after_build_mb": c.get("embed_rss_after_build_mb"),
                "model_storage_bytes": c.get("model_storage_bytes"),  # the populated .py = "the model"
                "embedding_bytes": c.get("embedding_bytes"),
                "agent_stage_a_s": c.get("stage_a_s"),
                "agent_stage_b_refine_s": c.get("stage_b_s"),
                "agent_stage_b_iters": c.get("stage_b_iters"),
                "agent_stage_b_final_rating": c.get("stage_b_final_rating"),
                "agent_tokens_stage_a": c.get("agent_tokens_stage_a"),
                "agent_tokens_stage_b_refine": c.get("agent_tokens_stage_b_refine"),
                "agent_tokens_total": _refine_token_total(c),
                "status": (
                    f"refine={c.get('stage_b_status') or 'n/a'}; "
                    f"embed={c.get('embed_status') or 'n/a'}"
                ),
                "source": "cost",
            })
            # Stage A sibling row, when Embed also embedded the pre-refinement
            # snapshot (embedder_blueprint_unrefined.py). Shares Stage A authoring
            # cost with the refined row above; Stage B / refinement fields are
            # nulled because this cohort represents the blueprint *before* the
            # refine loop ran. The non-greedy cohort regex in every consumer
            # absorbs `_unrefined` into the model slot, so this row aggregates
            # as a distinct cohort under method = blueprint_<model>_unrefined.
            if c.get("embed_unrefined_walltime_s") is not None:
                u_method = meta["method"].replace("_refine", "_unrefined")
                u_agent_model = meta["agent_model"].replace("_refine", "_unrefined")
                u_walltime = (c.get("embed_unrefined_walltime_s")
                              if c.get("embed_unrefined_walltime_s") is not None
                              else c.get("embed_unrefined_s"))
                u_cost = {
                    "agent_tokens_stage_a": c.get("agent_tokens_stage_a"),
                    "agent_tokens_stage_b_refine": None,
                }
                rows.append({
                    "method": u_method, "family": family_of(u_method), "dataset": dataset,
                    "seed": seed, "run_id": c["run_id"], "cohort": cohort_dir.name,
                    "agent_model": u_agent_model, "reasoning": meta["reasoning"],
                    "variant": meta["variant"],
                    "embedding_dim": None,
                    "device": c.get("embed_unrefined_device") or "cpu",
                    "host": c.get("embed_unrefined_host"),
                    "walltime_s": u_walltime,
                    "load_data_s": c.get("embed_unrefined_load_data_s"),
                    "obtain_model_s": c.get("embed_unrefined_obtain_model_s"),
                    "model_source": c.get("embed_unrefined_model_source") or "authored",
                    "build_embedding_s": c.get("embed_unrefined_build_embedding_s"),
                    "other_s": c.get("embed_unrefined_other_s"),
                    "model_download_s": c.get("embed_unrefined_model_download_s", 0.0),
                    "peak_rss_mb": c.get("embed_unrefined_peak_rss_mb"),
                    "rss_after_load_mb": c.get("embed_unrefined_rss_after_load_mb"),
                    "rss_after_obtain_mb": c.get("embed_unrefined_rss_after_obtain_mb"),
                    "rss_after_build_mb": c.get("embed_unrefined_rss_after_build_mb"),
                    "model_storage_bytes": c.get("unrefined_model_storage_bytes"),
                    "embedding_bytes": c.get("embed_unrefined_bytes"),
                    # Stage A cost is shared (the authoring run produced both
                    # blueprints); refinement-loop fields are explicitly nulled.
                    "agent_stage_a_s": c.get("stage_a_s"),
                    "agent_stage_b_refine_s": None,
                    "agent_stage_b_iters": None,
                    "agent_stage_b_final_rating": None,
                    "agent_tokens_stage_a": c.get("agent_tokens_stage_a"),
                    "agent_tokens_stage_b_refine": None,
                    "agent_tokens_total": _refine_token_total(u_cost),
                    "status": (
                        "refine=n/a; "
                        f"embed={c.get('embed_unrefined_status') or 'n/a'}"
                    ),
                    "source": "cost",
                })
    return rows


def _fill_embedding_dim(df: pd.DataFrame, dataset: str) -> pd.DataFrame:
    """Backfill embedding_dim for rows that lack it (agent rows, older reports)
    from outs/metrics/<dataset>/metrics.csv when available."""
    mcsv = METRICS_ROOT / dataset / "metrics.csv"
    if not mcsv.exists():
        return df
    try:
        md = pd.read_csv(mcsv)
    except Exception:
        return df
    if "model" not in md.columns or "embedding_dim" not in md.columns:
        return df
    dim_by_model = (md.dropna(subset=["embedding_dim"])
                      .groupby("model")["embedding_dim"].median().to_dict())
    miss = df["embedding_dim"].isna()
    df.loc[miss, "embedding_dim"] = df.loc[miss, "method"].map(dim_by_model)
    return df


def build_dataset_df(dataset: str) -> pd.DataFrame:
    rows = conventional_rows(dataset) + agent_rows(dataset)
    df = pd.DataFrame(rows, columns=COLUMNS) if rows else pd.DataFrame(columns=COLUMNS)
    for c in NUMERIC_COLS:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = _fill_embedding_dim(df, dataset)
    if not df.empty:
        # Fill model_source where the (legacy) report didn't carry the tag.
        miss = df["model_source"].isna() | (df["model_source"].astype(str).str.strip() == "")
        df.loc[miss, "model_source"] = df.loc[miss, "method"].map(_infer_model_source)
    return df.sort_values(["family", "method", "seed"], na_position="last").reset_index(drop=True)


def _iqr(s: pd.Series) -> float:
    s = s.dropna()
    return float(s.quantile(0.75) - s.quantile(0.25)) if len(s) else float("nan")


def _mode_or_blank(s: pd.Series) -> str:
    s = s.dropna().astype(str)
    m = s.mode()
    return str(m.iloc[0]) if len(m) else ""


def summarize(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame()
    g = df.groupby(["family", "method", "dataset"], dropna=False)
    out = g.agg(
        n=("method", "size"),
        device=("device", _mode_or_blank),
        host=("host", _mode_or_blank),
        model_source=("model_source", _mode_or_blank),
        **{f"{c}_median": (c, "median") for c in NUMERIC_COLS},
        **{f"{c}_iqr": (c, _iqr) for c in NUMERIC_COLS},
    ).reset_index()
    return out.sort_values(["family", "method", "dataset"], na_position="last").reset_index(drop=True)


def _human_bytes(x) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "—"
    x = float(x)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if x < 1024 or unit == "TB":
            return f"{x:.0f} {unit}" if unit == "B" else f"{x:.1f} {unit}"
        x /= 1024.0
    return f"{x:.1f} TB"


def _fmt_s(x) -> str:
    return "—" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{float(x):.1f}"


def _fmt_n(x) -> str:
    return "—" if x is None or (isinstance(x, float) and np.isnan(x)) else f"{int(round(float(x)))}"


def write_markdown(summary: pd.DataFrame, path: Path, title: str) -> None:
    lines = [f"# {title}", "",
             "Per-method medians over seeds (IQR in `costs_summary.csv`). **Time taxonomy** — "
             "the four buckets sum to `walltime` for every script we own: **`load`** = read the "
             "input AnnData from disk; **`obtain-model`** = get a usable model — train it "
             "(de-novo VAEs, PCA loadings) or load a persisted/pretrained checkpoint + adapt it "
             "to the dataset's gene space (the `source` column tags which: *trained* / *loaded* / "
             "*none* — UCell acquires nothing — / *authored* — the LLM wrote the script); "
             "**`n/a`** when load-and-forward aren't separable upstream (UCE's eval subprocess); "
             "**`build-emb`** = the forward pass / scoring; **`other`** = the residual (data prep, "
             "parquet write, glue). For the **agent methods** `walltime` is just the Embed cell-"
             "scoring run — the *model-acquisition* cost is the LLM authoring, **`Stage A`** (+ the "
             "optional **`Stage B`** refinement loop), a separate phase, also reported in tokens. "
             "**`model dl`** = one-time cold-cache fetch of a pretrained checkpoint (0 when already "
             "on disk; never part of `walltime` in steady state). **Memory** (`peak RSS`) is a "
             "high-water max, not a sum — see `costs.csv` for the after-load → after-model → "
             "after-embedding trajectory. **`model storage`** = the persisted artifact on disk "
             "(scVI/expiMap checkpoint, UCE/Census `.pt`, gene-set `.gmt`, or — for the agent "
             "methods — the populated `.py` script). **`device`** = compute device (cpu/cuda/mps; "
             "cpu for the agents' Embed build). Blank cells = the producing script predates this "
             "instrumentation; re-run it to populate.", ""]
    if summary.empty:
        lines.append("_(no data found)_")
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("\n".join(lines) + "\n")
        return
    multi_ds = summary["dataset"].nunique() > 1
    hdr = ["method", "family"] + (["dataset"] if multi_ds else []) + [
        "device", "source", "n", "walltime s", "load s", "obtain-model s", "build-emb s", "other s",
        "model dl s", "peak RSS MB", "model storage", "emb size", "dim",
        "Stage A s", "Stage B refine s", "refine iters", "tokens A", "tokens B refine"]
    lines.append("| " + " | ".join(hdr) + " |")
    lines.append("|" + "|".join(["---"] * len(hdr)) + "|")
    for _, r in summary.iterrows():
        cells = [f"`{r['method']}`", str(r["family"])] + ([str(r["dataset"])] if multi_ds else []) + [
            (str(r["device"]) or "—"), (str(r.get("model_source") or "—")), _fmt_n(r["n"]),
            _fmt_s(r["walltime_s_median"]), _fmt_s(r["load_data_s_median"]),
            _fmt_s(r["obtain_model_s_median"]), _fmt_s(r["build_embedding_s_median"]),
            _fmt_s(r["other_s_median"]), _fmt_s(r["model_download_s_median"]),
            _fmt_s(r["peak_rss_mb_median"]), _human_bytes(r["model_storage_bytes_median"]),
            _human_bytes(r["embedding_bytes_median"]), _fmt_n(r["embedding_dim_median"]),
            _fmt_n(r["agent_stage_a_s_median"]), _fmt_n(r["agent_stage_b_refine_s_median"]),
            _fmt_n(r["agent_stage_b_iters_median"]),
            _fmt_n(r["agent_tokens_stage_a_median"]), _fmt_n(r["agent_tokens_stage_b_refine_median"])]
        lines.append("| " + " | ".join(cells) + " |")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _write_csv(df: pd.DataFrame, path: Path, overwrite: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(f"{path} exists; pass --overwrite to replace it.")
    df.to_csv(path, index=False)


def run_dataset(dataset: str, overwrite: bool) -> pd.DataFrame:
    df = build_dataset_df(dataset)
    out_dir = OUT_ROOT / dataset
    _write_csv(df, out_dir / "costs.csv", overwrite)
    summary = summarize(df)
    _write_csv(summary, out_dir / "costs_summary.csv", overwrite)
    write_markdown(summary, out_dir / "costs_report.md",
                   f"Cost / wall-clock / storage — {dataset}")
    print(f"[ok] {dataset}: {len(df)} rows over {df['method'].nunique()} methods "
          f"→ {_rel(out_dir)}/costs.csv")
    return df


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--dataset", required=True,
                   help=f"One of {DATASET_CHOICES}, or 'all'.")
    p.add_argument("--overwrite", action="store_true",
                   help="Replace existing CSV/markdown outputs.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    datasets = list(DATASET_CHOICES) if args.dataset == "all" else [args.dataset]
    if args.dataset != "all" and args.dataset not in DATASET_CHOICES:
        print(f"--dataset must be one of {DATASET_CHOICES} or 'all' (got {args.dataset!r})",
              file=sys.stderr)
        sys.exit(1)
    frames = []
    for ds in datasets:
        frames.append(run_dataset(ds, args.overwrite))
    if args.dataset == "all":
        combined = pd.concat([f for f in frames if not f.empty], ignore_index=True) if frames else pd.DataFrame(columns=COLUMNS)
        _write_csv(combined, OUT_ROOT / "costs_all.csv", args.overwrite)
        write_markdown(summarize(combined), OUT_ROOT / "costs_table.md",
                       "Cost / wall-clock / storage — all datasets")
        print(f"[ok] combined: {len(combined)} rows → {_rel(OUT_ROOT)}/costs_all.csv")


if __name__ == "__main__":
    main()
