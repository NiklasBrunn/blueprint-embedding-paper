"""
Blueprint variability metrics across replicate agent runs.

For every work/<cohort>/ parent dir matching --dataset, computes:

  Per-run scalars:
    - blueprint size (modules, axes, genes, resources, DOIs, redundancy)
    - axis-kind composition (plain / contrast / gated)
    - resource-kind composition (paper / atlas / database)
    - sanity flags (placeholders, orphaned modules, missing/unused resources,
      clip overrides)

  Per-item recurrence (modules / genes / axes / resources by DOI):
    - n_runs and fraction of runs containing the item
    - concentration: fraction in ≥80 % of runs, Gini of recurrence

  Pairwise Jaccard across runs (genes / DOIs / module names):
    - mean off-diagonal + bootstrap 95 % CI

  Within-name structural consistency:
    - axis kind agreement (scIB-style)
    - axis positive / gate / negatives stability
    - module gene-set Jaccard within a recurring module name
    - module role-set agreement

  Cohort-level summary (`cohort_summary.json` + cross-cohort `cohort_metrics.csv`):
    - means / stds / CVs of size metrics
    - mean kind / role / positive / gate / negatives agreements
    - bootstrap 95 % CIs on mean off-diag Jaccard
    - saturation slope at N (cumulative-union curves over random orderings)
    - structural stability composite ∈ [0,1] (equal-weighted, weights exposed;
      grounding-mode-invariant — excludes DOI-Jaccard, which is 0 for noref by
      construction; DOI overlap is reported separately as a grounding metric)
    - DOI-kind agreement, module-name near-duplicate count

Inputs:
    work/<cohort>/<runNN>/summary.json          (one per replicate)
    configs/name_mappings/<dataset>/{module,axis,resource}_names.csv
        (curated standardized-name maps; read for matching and
         merge-appended with any newly-seen raw names)

Outputs (per cohort):
    outs/blueprint_variability/<dataset>/<cohort>/
        run_summary.csv, modules_long.csv, axes_long.csv
        module_recurrence.csv  gene_recurrence.csv
        axis_recurrence.csv    resource_recurrence.csv
        module_geneset_stability.csv  module_role_consistency.csv
        axis_kind_consistency.csv     axis_structure_consistency.csv
        module_name_near_duplicates.csv
        unique_modules.csv  unique_axes.csv
        jaccard_genes.csv  jaccard_dois.csv  jaccard_resources.csv
        jaccard_module_names.csv  jaccard_axis_names.csv
        saturation.csv
        cohort_summary.json
        report.md

Cross-cohort (per dataset):
    outs/blueprint_variability/<dataset>/cohort_metrics.csv
    outs/blueprint_variability/<dataset>/name_mapping/extracted_{modules,axes,resources}.csv
    (side-effect) merge-appends newly-seen raw names into
        configs/name_mappings/<dataset>/{module,axis,resource}_names.csv

Usage:
    python scripts/blueprint_variability_metrics/metrics.py --dataset mouse_brain
    python scripts/blueprint_variability_metrics/metrics.py --dataset human_pancreas
    python scripts/blueprint_variability_metrics/metrics.py --dataset mouse_brain \
        --cohort gpt_5_5_high_mouse_brain

"""
from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")
WORK_ROOT = REPO_ROOT / "work"
# Human-curated name mappings (tracked, never clobbered). One dir per dataset
# holding module_names.csv / axis_names.csv / resource_names.csv; the human fills
# a `standardized_name` column that overrides lexical canonicalization when
# matching modules / axes / resources across runs (blank ⇒ canonical fallback).
CONFIG_ROOT = REPO_ROOT / "configs" / "name_mappings"

# Cohort dir naming: <model_slug>_<reasoning>_<dataset>.
# The <model> slot is non-greedy and stops at the first reasoning token,
# so optional model-slug markers like `_userctx` (claude no-isolate) and
# `_refine` (codex/claude --with-refine) are absorbed into the captured
# model name. That is by design: refined / no-isolate runs surface as
# distinct model labels in plots and aggregations, never commingled
# with their counterparts.
COHORT_RE = re.compile(
    r"^(?P<model>.+?)_(?P<reasoning>low|medium|high|xhigh)_(?P<dataset>.+)$"
)
DOI_PREFIX_RE = re.compile(r"^https?://(?:dx\.)?doi\.org/", re.IGNORECASE)

# Equal-weighted composite of *structural reproducibility* components. Exposed so
# a reviewer can reweight without forking.
#
# DOI-Jaccard is deliberately NOT a component: it is 0 for `noref` cohorts by
# construction (they carry no citations), so including it would penalise noref
# for being citation-free rather than for being unstable, making the composite
# non-comparable across grounding modes. DOI overlap is still reported on its own
# (`jaccard_dois_off_diag`, recurrence) as a *grounding* metric, not a stability one.
STABILITY_WEIGHTS: dict[str, float] = {
    "jaccard_genes_off_diag":              1.0,
    "jaccard_modnames_off_diag":           1.0,
    "mean_kind_agreement":                 1.0,
    "mean_role_agreement":                 1.0,
    "mean_within_module_geneset_jaccard":  1.0,
    "mean_axis_positive_agreement":        1.0,
}

DEFAULT_BOOTSTRAP = 1000
DEFAULT_SAT_ORDERINGS = 200
NEAR_DUP_MAX_EDIT = 2
RNG_SEED = 0

# Canonicalization for module / axis / gene names. Lowercase, then strip every
# non-alphanumeric character — collapses `l23_it` / `L23_IT` / `l2_3_it` / `L2/3 IT`
# all to `l23it`. Display tables retain the original spellings (mode of raw forms).
_CAN_RE = re.compile(r"[^a-z0-9]+")


# ============================================================
# Helpers
# ============================================================
def normalize_doi(s: str | None) -> str:
    if not s:
        return ""
    return DOI_PREFIX_RE.sub("", str(s).strip()).lower()


def canonicalize(s: str | None) -> str:
    if s is None:
        return ""
    return _CAN_RE.sub("", str(s).lower())


# ============================================================
# Human-curated name mapping (modules / axes / resources)
# ============================================================
def load_name_mapping(path: Path) -> dict[str, str]:
    """Read a curated mapping CSV → {key: standardized_name} for the rows the
    human has filled in (non-empty `standardized_name`). Key column is
    `resource_id` for resources, else `raw_name`. Missing file ⇒ empty mapping
    (matching then falls back to canonical / `resource_identifier`)."""
    if not path.exists():
        return {}
    df = pd.read_csv(path, dtype=str).fillna("")
    key_col = "resource_id" if "resource_id" in df.columns else "raw_name"
    out: dict[str, str] = {}
    for k, v in zip(df[key_col], df.get("standardized_name", "")):
        k, v = str(k), str(v).strip()
        if k and v:
            out[k] = v
    return out


def merge_append_mapping(path: Path, keys, key_col: str = "raw_name") -> None:
    """Create-or-merge the tracked mapping at `path`: add any newly-seen `keys`
    as rows with a blank `standardized_name`, **preserving every existing row
    verbatim** (human edits are never dropped or rewritten, even if a key goes
    stale). Sorted by the key column for stable git diffs."""
    cols = [key_col, "standardized_name", "notes"]
    if path.exists():
        existing = pd.read_csv(path, dtype=str).fillna("")
        for c in cols:
            if c not in existing.columns:
                existing[c] = ""
        existing = existing[cols]
    else:
        existing = pd.DataFrame(columns=cols)
    have = set(existing[key_col].astype(str))
    new_keys = sorted({str(k) for k in keys if str(k).strip()} - have)
    if new_keys:
        add = pd.DataFrame({key_col: new_keys, "standardized_name": "", "notes": ""})
        existing = pd.concat([existing, add], ignore_index=True)
    existing = existing.sort_values(key_col, kind="stable").reset_index(drop=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    existing.to_csv(path, index=False)


def make_resolver(mapping: dict[str, str]):
    """raw name → match key: `canonicalize(standardized)` when curated, else
    `canonicalize(raw)` (i.e. identical to the un-curated behavior). Canonicalizing
    the resolved name guarantees an empty mapping reproduces current outputs
    exactly and that two raws sharing a standardized name collide."""
    def resolve(raw: str | None) -> str:
        std = mapping.get(str(raw) if raw is not None else "", "")
        return canonicalize(std) if std else canonicalize(raw)
    return resolve


def make_resource_resolver(mapping: dict[str, str]):
    """resource dict → match key: `canonicalize(standardized)` when its
    identifier is curated, else `resource_identifier(r)` (DOI-else-key)."""
    def resolve(r: dict) -> str:
        rid = resource_identifier(r)
        std = mapping.get(rid, "")
        return canonicalize(std) if std else rid
    return resolve


def _std_display_map(mapping: dict[str, str]) -> dict[str, str]:
    """match key (`canonicalize(std)`) → standardized spelling, so curated names
    are shown verbatim in display columns / figures."""
    return {canonicalize(v): v for v in mapping.values() if v}


def _can_to_raw_map(raw_iter, can_iter) -> dict[str, str]:
    """canonical → most-common raw spelling. Ties broken by first-seen."""
    by_can: dict[str, list[str]] = {}
    for raw, can in zip(raw_iter, can_iter):
        if not can:
            continue
        by_can.setdefault(can, []).append(raw)
    return {can: Counter(forms).most_common(1)[0][0] for can, forms in by_can.items()}


def _add_display_column(rec_df: pd.DataFrame, canonical_col: str,
                        display_name: str, can_to_raw: dict[str, str]) -> pd.DataFrame:
    """Insert a `<display_name>` column with the most-common raw form right
    after `<canonical_col>`. Falls back to canonical when no raw seen."""
    if rec_df.empty:
        rec_df[display_name] = pd.Series(dtype=object)
        return rec_df
    rec_df = rec_df.copy()
    rec_df[display_name] = rec_df[canonical_col].map(can_to_raw).fillna(rec_df[canonical_col])
    cols = [c for c in rec_df.columns if c != display_name]
    idx = cols.index(canonical_col) + 1
    cols.insert(idx, display_name)
    return rec_df[cols]


def _collisions_count(raws_iter, cans_iter) -> int:
    """How many canonical keys map to ≥ 2 distinct raw spellings."""
    by_can: dict[str, set[str]] = {}
    for raw, can in zip(raws_iter, cans_iter):
        if not can:
            continue
        by_can.setdefault(can, set()).add(raw)
    return sum(1 for v in by_can.values() if len(v) > 1)


def _levenshtein(a: str, b: str, max_d: int = NEAR_DUP_MAX_EDIT) -> int:
    """Truncated Levenshtein; returns max_d+1 once exceeded."""
    if a == b:
        return 0
    if abs(len(a) - len(b)) > max_d:
        return max_d + 1
    n, m = len(a), len(b)
    if n == 0:
        return m
    if m == 0:
        return n
    prev = list(range(m + 1))
    for i, ca in enumerate(a, 1):
        curr = [i] + [0] * m
        row_min = i
        for j, cb in enumerate(b, 1):
            curr[j] = min(
                prev[j] + 1,
                curr[j - 1] + 1,
                prev[j - 1] + (0 if ca == cb else 1),
            )
            if curr[j] < row_min:
                row_min = curr[j]
        if row_min > max_d:
            return max_d + 1
        prev = curr
    return prev[m]


def gini(values) -> float:
    v = np.asarray(values, dtype=float)
    if v.size == 0:
        return float("nan")
    if v.sum() == 0:
        return 0.0
    v = np.sort(v)
    n = v.size
    cum = np.cumsum(v)
    return float((n + 1 - 2 * cum.sum() / cum[-1]) / n)


def cv(s: pd.Series) -> float:
    if s.empty:
        return float("nan")
    m = float(s.mean())
    if not np.isfinite(m) or m == 0.0:
        return float("nan")
    sd = float(s.std(ddof=1)) if len(s) > 1 else 0.0
    return sd / m


# ============================================================
# Discovery + loading
# ============================================================
def discover_cohorts(dataset: str, cohort: str | None) -> list[Path]:
    if not WORK_ROOT.exists():
        raise FileNotFoundError(f"work root not found: {WORK_ROOT}")
    candidates = [
        p for p in sorted(WORK_ROOT.iterdir())
        if p.is_dir() and p.name.endswith(f"_{dataset}")
        # Exploratory local-model runs (OpenWebUI/Ollama/UFR) are excluded from
        # the aggregation for now; drop the guard to fold them back in.
        and "openwebui" not in p.name
    ]
    if cohort is not None:
        candidates = [p for p in candidates if p.name == cohort]
        if not candidates:
            raise FileNotFoundError(
                f"cohort {cohort!r} not found under {WORK_ROOT}")
    return candidates


def discover_runs(cohort_dir: Path) -> list[Path]:
    return sorted(p for p in cohort_dir.iterdir()
                  if p.is_dir() and (p / "summary.json").is_file())


def load_summary(run_dir: Path) -> dict | None:
    try:
        with (run_dir / "summary.json").open() as fh:
            return json.load(fh)
    except (OSError, json.JSONDecodeError) as exc:
        print(f"[warn] could not read {run_dir.name}/summary.json: {exc}",
              file=sys.stderr)
        return None


def parse_cohort_name(name: str) -> dict[str, str]:
    m = COHORT_RE.match(name)
    if m is None:
        return {"model": "", "reasoning": "", "dataset": ""}
    return m.groupdict()


def _as_seq(v) -> list:
    """Scenario fields that may be a str, a list, or None -> a list (empty if None).

    `species` is array-typed in current scenario configs but older
    `summary.json` files on disk still carry it as a bare string; this
    accepts both. Used wherever such a field needs to be a list (or a
    hashable tuple, or a "|"-joined label).
    """
    if v is None:
        return []
    return list(v) if isinstance(v, list) else [v]


def scenario_key(scenario: dict) -> tuple:
    # species / tissue are array-typed; celltype / condition are scalar-or-null
    # (and `species` may also arrive as a bare string in older summaries).
    # Normalise each to a tuple so the returned key stays hashable — it is used
    # as a dict key when grouping summaries by scenario.
    return (
        tuple(_as_seq(scenario.get("species"))),
        tuple(_as_seq(scenario.get("tissue"))),
        tuple(_as_seq(scenario.get("celltype"))),
        scenario.get("age"),
        scenario.get("sex"),
        tuple(_as_seq(scenario.get("condition"))),
    )


# ============================================================
# DataFrame builders
# ============================================================
def build_run_df(summaries: list[dict]) -> pd.DataFrame:
    rows = []
    for s in summaries:
        bp = s.get("blueprint") or {}
        sc = s.get("scenario") or {}
        sa = s.get("sanity") or {}
        akb = bp.get("n_axes_by_kind") or {}
        rkb = bp.get("n_resources_by_kind") or {}
        n_axes = int(bp.get("n_axes", 0) or 0)
        n_total = int(bp.get("n_total_genes", 0) or 0)
        n_unique = int(bp.get("n_unique_genes", 0) or 0)
        n_resources = int(bp.get("n_resources", 0) or 0)
        modules = bp.get("modules") or []
        axes = bp.get("axes") or []
        gpm = [int(m.get("n_genes", 0) or 0) for m in modules]
        rkpm = [int(m.get("n_resource_keys", 0) or 0) for m in modules]
        apm = [int(m.get("n_axes_using", 0) or 0) for m in modules]
        npa = [len(a.get("negatives") or []) for a in axes]
        n_paper = int(rkb.get("paper", 0) or 0)
        n_atlas = int(rkb.get("atlas", 0) or 0)
        n_database = int(rkb.get("database", 0) or 0)
        rows.append({
            "run_id": s.get("run_id"),
            "species": "|".join(_as_seq(sc.get("species"))),
            "tissue": "|".join(sc.get("tissue") or []),
            "celltype": sc.get("celltype"),
            "scenario_key": str(scenario_key(sc)),
            "n_modules": int(bp.get("n_modules", 0) or 0),
            "n_axes": n_axes,
            "n_axes_plain":    int(akb.get("plain", 0) or 0),
            "n_axes_contrast": int(akb.get("contrast", 0) or 0),
            "n_axes_gated":    int(akb.get("gated", 0) or 0),
            "frac_plain":    (akb.get("plain", 0) / n_axes) if n_axes else 0.0,
            "frac_contrast": (akb.get("contrast", 0) / n_axes) if n_axes else 0.0,
            "frac_gated":    (akb.get("gated", 0) / n_axes) if n_axes else 0.0,
            "n_total_genes":  n_total,
            "n_unique_genes": n_unique,
            "gene_redundancy": (n_total / n_unique) if n_unique else float("nan"),
            "n_resources":          n_resources,
            "n_resources_paper":    n_paper,
            "n_resources_atlas":    n_atlas,
            "n_resources_database": n_database,
            "frac_resources_paper":    (n_paper / n_resources) if n_resources else 0.0,
            "frac_resources_atlas":    (n_atlas / n_resources) if n_resources else 0.0,
            "frac_resources_database": (n_database / n_resources) if n_resources else 0.0,
            "n_dois_unique": int(bp.get("n_dois_unique", 0) or 0),
            "n_orphaned_modules":          len(sa.get("orphaned_modules") or []),
            "n_modules_missing_resources": len(sa.get("modules_missing_resources") or []),
            "n_unused_resources":          len(sa.get("unused_resources") or []),
            "placeholder_strings_present": bool(sa.get("placeholder_strings_present")),
            "set_clip_overridden":         bool(sa.get("set_clip_overridden")),
            "mean_genes_per_module":   float(np.mean(gpm)) if gpm else 0.0,
            "median_genes_per_module": float(np.median(gpm)) if gpm else 0.0,
            "max_genes_per_module":    int(np.max(gpm)) if gpm else 0,
            "mean_resources_per_module": float(np.mean(rkpm)) if rkpm else 0.0,
            "mean_axes_per_module":      float(np.mean(apm)) if apm else 0.0,
            "mean_negatives_per_axis":   float(np.mean(npa)) if npa else 0.0,
        })
    return pd.DataFrame(rows).sort_values("run_id").reset_index(drop=True)


def build_module_df(summaries: list[dict],
                    mod_resolver=canonicalize) -> pd.DataFrame:
    """Per-run-per-module rows. Both raw and canonical names/genes are kept:

      - `module_name`     — original spelling (display)
      - `module_name_can` — match key: `mod_resolver(name)` (curated standardized
                            name when available, else lexical `canonicalize`)
      - `genes`           — original gene list (display)
      - `genes_can`       — canonical gene list (matching key; genes untouched)
    """
    rows = []
    for s in summaries:
        run_id = s.get("run_id")
        sa = s.get("sanity") or {}
        orphans = set(sa.get("orphaned_modules") or [])
        bp = s.get("blueprint") or {}
        key_to_doi = {
            r.get("key"): normalize_doi(r.get("doi"))
            for r in (bp.get("resources") or [])
            if r.get("key")
        }
        for m in bp.get("modules") or []:
            name_raw = m.get("name") or ""
            name_can = mod_resolver(name_raw)
            roles = sorted({u.get("role") for u in (m.get("used_in_axes") or [])
                            if u.get("role")})
            rkeys = list(m.get("resource_keys") or [])
            rdois = [d for d in (key_to_doi.get(k, "") for k in rkeys) if d]
            genes_raw = list(m.get("genes") or [])
            genes_can = [canonicalize(g) for g in genes_raw]
            rows.append({
                "run_id":          run_id,
                "module_name":     name_raw,
                "module_name_can": name_can,
                "n_genes":         int(m.get("n_genes", 0) or 0),
                "n_resource_keys": int(m.get("n_resource_keys", 0) or 0),
                "n_axes_using":    int(m.get("n_axes_using", 0) or 0),
                "is_orphan":       name_raw in orphans,
                "genes":           genes_raw,
                "genes_can":       genes_can,
                "resource_keys":   rkeys,
                "resource_dois":   rdois,
                "roles":           "|".join(roles),
            })
    return pd.DataFrame(rows)


def build_axis_df(summaries: list[dict],
                  axis_resolver=canonicalize,
                  mod_resolver=canonicalize) -> pd.DataFrame:
    """Per-run-per-axis rows with raw + match-key companion columns. Axis names
    use `axis_resolver`; the module-referencing fields (positive / gate /
    negatives) use `mod_resolver` so they stay consistent with module keys."""
    rows = []
    for s in summaries:
        run_id = s.get("run_id")
        for a in (s.get("blueprint") or {}).get("axes") or []:
            ax_raw = a.get("name") or ""
            pos_raw = a.get("positive") or ""
            gate_raw = a.get("gate") or ""
            negs_raw = list(a.get("negatives") or [])
            rows.append({
                "run_id":           run_id,
                "axis_name":        ax_raw,
                "axis_name_can":    axis_resolver(ax_raw),
                "kind":             a.get("kind"),
                "positive":         pos_raw,
                "positive_can":     mod_resolver(pos_raw),
                "negatives":        negs_raw,
                "negatives_can":    [mod_resolver(x) for x in negs_raw],
                "n_negatives":      len(negs_raw),
                "negative_weight":  a.get("negative_weight"),
                "gate":             gate_raw,
                "gate_can":         mod_resolver(gate_raw) if gate_raw else "",
            })
    return pd.DataFrame(rows)


# ============================================================
# Pairwise Jaccard
# ============================================================
def _union_per_run(module_df: pd.DataFrame, col: str, *,
                   flatten: bool) -> dict[str, set]:
    """Per run, union the values in `col`. `col` should be a canonical column
    (e.g. `module_name_can`, `genes_can`, `resource_dois`) so set operations
    measure agreement on concepts, not on punctuation/case variants."""
    out: dict[str, set] = {}
    for rid, sub in module_df.groupby("run_id"):
        s: set = set()
        for v in sub[col]:
            if flatten:
                s.update(v)
            else:
                s.add(v)
        out[rid] = s
    return out


def jaccard_matrix(sets_by_run: dict[str, set]) -> pd.DataFrame:
    runs = sorted(sets_by_run)
    n = len(runs)
    mat = np.zeros((n, n))
    for i, ra in enumerate(runs):
        for j, rb in enumerate(runs):
            a, b = sets_by_run[ra], sets_by_run[rb]
            u = a | b
            mat[i, j] = (len(a & b) / len(u)) if u else 0.0
    return pd.DataFrame(mat, index=runs, columns=runs)


def _off_diag(jacc: pd.DataFrame) -> np.ndarray:
    if jacc.shape[0] < 2:
        return np.array([])
    arr = jacc.values
    return arr[~np.eye(arr.shape[0], dtype=bool)]


def bootstrap_offdiag_ci(jacc: pd.DataFrame, rng: np.random.Generator,
                         n_boot: int = DEFAULT_BOOTSTRAP) -> dict[str, float]:
    """Mean off-diagonal + 95 % CI from `n_boot` resamples of the pair distribution."""
    off = _off_diag(jacc)
    if off.size == 0:
        return {"mean": float("nan"), "ci_lo": float("nan"), "ci_hi": float("nan"),
                "min": float("nan"), "max": float("nan")}
    means = np.empty(n_boot, dtype=float)
    n = off.size
    for k in range(n_boot):
        means[k] = float(rng.choice(off, size=n, replace=True).mean())
    return {
        "mean":  float(off.mean()),
        "ci_lo": float(np.percentile(means, 2.5)),
        "ci_hi": float(np.percentile(means, 97.5)),
        "min":   float(off.min()),
        "max":   float(off.max()),
    }


# ============================================================
# Recurrence
# ============================================================
def _recurrence(records: list[tuple[str, str]], n_runs: int,
                label_col: str) -> pd.DataFrame:
    by_item: dict[str, set[str]] = {}
    for run_id, item in records:
        by_item.setdefault(item, set()).add(run_id)
    rows = [{
        label_col: item,
        "n_runs": len(runs),
        "fraction_runs": (len(runs) / n_runs) if n_runs else 0.0,
        "run_ids": "|".join(sorted(runs)),
    } for item, runs in by_item.items()]
    if not rows:
        return pd.DataFrame(columns=[label_col, "n_runs", "fraction_runs", "run_ids"])
    return (pd.DataFrame(rows)
            .sort_values(["n_runs", label_col], ascending=[False, True])
            .reset_index(drop=True))


def module_recurrence(module_df: pd.DataFrame, n_runs: int,
                      can_to_raw: dict[str, str]) -> pd.DataFrame:
    rec = _recurrence(list(zip(module_df["run_id"], module_df["module_name_can"])),
                      n_runs, "module_name_can")
    return _add_display_column(rec, "module_name_can", "module_name", can_to_raw)


def gene_recurrence(module_df: pd.DataFrame, n_runs: int,
                    can_to_raw: dict[str, str]) -> pd.DataFrame:
    seen: set[tuple[str, str]] = set()
    for run_id, genes_can in zip(module_df["run_id"], module_df["genes_can"]):
        seen.update((run_id, g) for g in genes_can if g)
    rec = _recurrence(list(seen), n_runs, "gene_can")
    return _add_display_column(rec, "gene_can", "gene", can_to_raw)


def axis_recurrence(axis_df: pd.DataFrame, n_runs: int,
                    can_to_raw: dict[str, str]) -> pd.DataFrame:
    rec = _recurrence(list(zip(axis_df["run_id"], axis_df["axis_name_can"])),
                      n_runs, "axis_name_can")
    return _add_display_column(rec, "axis_name_can", "axis_name", can_to_raw)


def resource_identifier(r: dict) -> str:
    """Stable cross-run identifier for a resource: its normalized DOI if it has
    one (papers, atlases), else its canonicalized key (databases such as KEGG /
    MSigDB / Reactome carry no DOI). Lets database resources participate in
    recurrence / Jaccard rather than being silently dropped; key-based matching
    is run-local and therefore lexical (a conservative lower bound), as for
    module / axis names."""
    doi = normalize_doi(r.get("doi"))
    if doi:
        return doi
    key = canonicalize(r.get("key") or r.get("citation") or "")
    return f"key:{key}" if key else ""


def _resources_by_run(summaries: list[dict],
                      res_resolver=resource_identifier) -> dict[str, set]:
    """Per run, the set of declared-catalog resource match keys. Default key is
    DOI-else-canonical-key (`resource_identifier`); a curated `res_resolver`
    collapses standardized variants (mostly DOI-less databases)."""
    out: dict[str, set] = {}
    for s in summaries:
        ids = {res_resolver(r)
               for r in (s.get("blueprint") or {}).get("resources") or []}
        out[s.get("run_id")] = {i for i in ids if i}
    return out


def resource_recurrence(n_runs: int, summaries: list[dict],
                        res_resolver=resource_identifier) -> pd.DataFrame:
    """Recurrence over the declared resource catalog, keyed on the resource match
    key (DOI-else-key, or a curated standardized name) so every resource kind
    (paper / atlas / database) is represented. Same universe as the per-run
    resource-kind counts (`n_resources_by_kind`)."""
    seen: set[tuple[str, str]] = set()
    citation_by_id: dict[str, str] = {}
    kind_by_id: dict[str, str] = {}
    doi_by_id: dict[str, str] = {}
    for s in summaries:
        run_id = s.get("run_id")
        for r in (s.get("blueprint") or {}).get("resources") or []:
            rid = res_resolver(r)
            if not rid:
                continue
            seen.add((run_id, rid))
            citation_by_id.setdefault(rid, r.get("citation") or "")
            kind_by_id.setdefault(rid, r.get("kind") or "")
            doi_by_id.setdefault(rid, normalize_doi(r.get("doi")))
    rec = _recurrence(list(seen), n_runs, "resource_id")
    if not rec.empty:
        rec["citation"] = rec["resource_id"].map(lambda i: citation_by_id.get(i, ""))
        rec["kind"]     = rec["resource_id"].map(lambda i: kind_by_id.get(i, ""))
        rec["doi"]      = rec["resource_id"].map(lambda i: doi_by_id.get(i, ""))
    return rec


def recurrence_concentration(recur_df: pd.DataFrame, n_runs: int) -> dict[str, float]:
    if recur_df.empty or n_runs < 1:
        return {"frac_in_ge80pct": 0.0, "gini": float("nan"),
                "n_universal": 0, "n_majority": 0, "n_singleton": 0,
                "n_unique_total": 0}
    threshold = max(1, int(np.ceil(0.8 * n_runs)))
    return {
        "frac_in_ge80pct": float((recur_df["n_runs"] >= threshold).mean()),
        "gini":            gini(recur_df["n_runs"].to_numpy()),
        "n_universal":     int((recur_df["n_runs"] == n_runs).sum()),
        "n_majority":      int((recur_df["n_runs"] >= (n_runs + 1) // 2).sum()),
        "n_singleton":     int((recur_df["n_runs"] == 1).sum()),
        "n_unique_total":  int(len(recur_df)),
    }


# ============================================================
# Within-name structural consistency
# ============================================================
def module_geneset_stability(module_df: pd.DataFrame,
                              mod_can_to_raw: dict[str, str]) -> pd.DataFrame:
    """Within-module-name pairwise canonical-gene Jaccard, for modules that
    recur in ≥ 2 runs. Grouping is on canonical module name; gene sets use
    canonical genes so case/punctuation variants do not deflate Jaccard."""
    rows = []
    for can, sub in module_df.groupby("module_name_can"):
        if len(sub) < 2:
            continue
        gene_sets = [set(g) for g in sub["genes_can"]]
        jaccards = []
        for a, b in combinations(gene_sets, 2):
            u = a | b
            jaccards.append((len(a & b) / len(u)) if u else 1.0)
        union = set().union(*gene_sets)
        sizes = [len(g) for g in gene_sets]
        rows.append({
            "module_name_can":  can,
            "module_name":      mod_can_to_raw.get(can, can),
            "n_runs_seen":      len(gene_sets),
            "mean_jaccard":     float(np.mean(jaccards)),
            "min_jaccard":      float(np.min(jaccards)),
            "max_jaccard":      float(np.max(jaccards)),
            "n_unique_genes":   len(union),
            "mean_n_genes":     float(np.mean(sizes)),
        })
    if not rows:
        return pd.DataFrame(columns=[
            "module_name_can", "module_name", "n_runs_seen", "mean_jaccard",
            "min_jaccard", "max_jaccard", "n_unique_genes", "mean_n_genes",
        ])
    return (pd.DataFrame(rows)
            .sort_values(["mean_jaccard", "n_runs_seen"], ascending=[True, False])
            .reset_index(drop=True))


def module_role_consistency(module_df: pd.DataFrame,
                            mod_can_to_raw: dict[str, str]) -> pd.DataFrame:
    rows = []
    for can, sub in module_df.groupby("module_name_can"):
        roles = [r if r else "(none)" for r in sub["roles"].tolist()]
        c = Counter(roles)
        modal, modal_n = c.most_common(1)[0]
        rows.append({
            "module_name_can":  can,
            "module_name":      mod_can_to_raw.get(can, can),
            "n_runs_seen":      len(roles),
            "modal_roles":      modal,
            "role_agreement":   modal_n / len(roles),
            "n_distinct_role_sets": len(c),
            "role_distribution": "|".join(f"{r}:{n}" for r, n in c.most_common()),
        })
    if not rows:
        return pd.DataFrame(columns=[
            "module_name_can", "module_name", "n_runs_seen", "modal_roles",
            "role_agreement", "n_distinct_role_sets", "role_distribution",
        ])
    return (pd.DataFrame(rows)
            .sort_values(["role_agreement", "n_runs_seen", "module_name_can"],
                         ascending=[True, False, True])
            .reset_index(drop=True))


def axis_kind_consistency(axis_df: pd.DataFrame,
                          axis_can_to_raw: dict[str, str]) -> pd.DataFrame:
    rows = []
    for can, sub in axis_df.groupby("axis_name_can"):
        kinds = sub["kind"].tolist()
        c = Counter(kinds)
        modal_kind, modal_n = c.most_common(1)[0]
        rows.append({
            "axis_name_can":    can,
            "axis_name":        axis_can_to_raw.get(can, can),
            "n_runs_seen":      len(kinds),
            "modal_kind":       modal_kind,
            "kind_agreement":   modal_n / len(kinds),
            "n_distinct_kinds": len(c),
            "kinds":            "|".join(f"{k}:{v}" for k, v in c.most_common()),
        })
    if not rows:
        return pd.DataFrame(columns=[
            "axis_name_can", "axis_name", "n_runs_seen", "modal_kind",
            "kind_agreement", "n_distinct_kinds", "kinds",
        ])
    return (pd.DataFrame(rows)
            .sort_values(["kind_agreement", "n_runs_seen", "axis_name_can"],
                         ascending=[True, False, True])
            .reset_index(drop=True))


def axis_structure_consistency(axis_df: pd.DataFrame,
                               axis_can_to_raw: dict[str, str],
                               mod_can_to_raw: dict[str, str]) -> pd.DataFrame:
    """Per axis_name (canonical): positive / gate stability + within-name
    Jaccard of the canonical negative-module set across runs."""
    rows = []
    for can, sub in axis_df.groupby("axis_name_can"):
        n = len(sub)
        positives = sub["positive_can"].fillna("").tolist()
        c_pos = Counter(positives)
        modal_pos, modal_pos_n = c_pos.most_common(1)[0]
        gates = sub["gate_can"].fillna("").tolist()
        c_gate = Counter(gates)
        modal_gate, modal_gate_n = c_gate.most_common(1)[0]
        neg_sets = [set(x for x in xs if x) for xs in sub["negatives_can"]]
        if len(neg_sets) >= 2:
            jaccards = []
            for a, b in combinations(neg_sets, 2):
                u = a | b
                jaccards.append((len(a & b) / len(u)) if u else 1.0)
            mean_neg = float(np.mean(jaccards))
        else:
            mean_neg = float("nan")
        rows.append({
            "axis_name_can":            can,
            "axis_name":                axis_can_to_raw.get(can, can),
            "n_runs_seen":              n,
            "modal_positive_can":       modal_pos or "",
            "modal_positive":           mod_can_to_raw.get(modal_pos, modal_pos) if modal_pos else "",
            "positive_agreement":       modal_pos_n / n,
            "n_distinct_positives":     len(c_pos),
            "modal_gate_can":           modal_gate or "",
            "modal_gate":               mod_can_to_raw.get(modal_gate, modal_gate) if modal_gate else "",
            "gate_agreement":           modal_gate_n / n,
            "mean_negatives_jaccard":   mean_neg,
        })
    if not rows:
        return pd.DataFrame(columns=[
            "axis_name_can", "axis_name", "n_runs_seen",
            "modal_positive_can", "modal_positive",
            "positive_agreement", "n_distinct_positives",
            "modal_gate_can", "modal_gate", "gate_agreement",
            "mean_negatives_jaccard",
        ])
    return (pd.DataFrame(rows)
            .sort_values(["positive_agreement", "n_runs_seen", "axis_name_can"],
                         ascending=[True, False, True])
            .reset_index(drop=True))


# ============================================================
# Sanity
# ============================================================
def module_name_near_duplicates(module_df: pd.DataFrame,
                                mod_can_to_raw: dict[str, str],
                                max_d: int = NEAR_DUP_MAX_EDIT) -> pd.DataFrame:
    """Pairs of distinct *canonical* module names within edit distance ≤ max_d.

    Trivial spelling drift (e.g. `l23_it` vs `L23_IT`, `non_neuronal_core` vs
    `nonneuronal_core`) is already merged at canonicalization. What survives
    here is real edit-distance proximity on canonical forms — typos, suffixes
    like `panneuron` vs `panneurons`, or distinct concepts that happen to
    share strings (e.g. `l4it` vs `l5it`). Pairs with `co_occurring_runs == 0`
    are higher-signal drift candidates than those with `co_occurring_runs > 0`.
    """
    by_can: dict[str, set[str]] = {}
    for run_id, can in zip(module_df["run_id"], module_df["module_name_can"]):
        if not can:
            continue
        by_can.setdefault(can, set()).add(run_id)
    cans = sorted(by_can)
    rows = []
    for a, b in combinations(cans, 2):
        d = _levenshtein(a, b, max_d=max_d)
        if d <= max_d:
            rows.append({
                "name_a_can":        a,
                "name_b_can":        b,
                "name_a":            mod_can_to_raw.get(a, a),
                "name_b":            mod_can_to_raw.get(b, b),
                "edit_distance":     d,
                "n_runs_a":          len(by_can[a]),
                "n_runs_b":          len(by_can[b]),
                "co_occurring_runs": len(by_can[a] & by_can[b]),
            })
    if not rows:
        return pd.DataFrame(columns=[
            "name_a_can", "name_b_can", "name_a", "name_b",
            "edit_distance", "n_runs_a", "n_runs_b", "co_occurring_runs",
        ])
    return (pd.DataFrame(rows)
            .sort_values(["edit_distance", "name_a_can"], ascending=[True, True])
            .reset_index(drop=True))


def doi_kind_agreement(summaries: list[dict]) -> float:
    """Per DOI, fraction of runs that agree on its `kind`. Mean across DOIs."""
    by_doi: dict[str, list[str]] = {}
    for s in summaries:
        for r in (s.get("blueprint") or {}).get("resources") or []:
            doi = normalize_doi(r.get("doi"))
            if doi:
                by_doi.setdefault(doi, []).append(r.get("kind") or "")
    if not by_doi:
        return float("nan")
    agreements = []
    for kinds in by_doi.values():
        c = Counter(kinds)
        modal_n = c.most_common(1)[0][1]
        agreements.append(modal_n / len(kinds))
    return float(np.mean(agreements))


def unique_module_names_table(module_df: pd.DataFrame, n_runs: int,
                              top_genes_k: int = 8) -> pd.DataFrame:
    """One row per canonical module name, expanded for human inspection.

    Columns:
      - `module_name_canonical`: matching key (lowercased + alnum-only)
      - `n_raw_spellings`:       distinct raw spellings observed for this concept
      - `raw_spellings`:         "name (count) | name (count) | …" sorted by frequency
      - `n_runs`, `fraction_runs`: support
      - `modal_role_set`:        most-common role string (positive/gate/negative)
      - `n_distinct_genes`:      union of canonical gene set across runs
      - `top_genes`:             representative genes (mode of raw spellings) for
                                 the top-k canonical genes by run frequency

    Curation of the matching key now lives in the tracked per-dataset mapping
    (`configs/name_mappings/<dataset>/module_names.csv`), not here — this table is
    a per-cohort spelling inventory only.
    """
    rows: list[dict] = []
    for can, sub in module_df.groupby("module_name_can"):
        if not can:
            continue
        raw_counter = Counter(sub["module_name"])
        raw_spellings = " | ".join(f"{r} ({n})" for r, n in raw_counter.most_common())
        roles_counter = Counter(r if r else "(none)" for r in sub["roles"])
        modal_role, _ = roles_counter.most_common(1)[0]
        gene_counter: Counter = Counter()
        gene_to_raw: dict[str, Counter] = {}
        for genes, genes_can in zip(sub["genes"], sub["genes_can"]):
            gene_counter.update(genes_can)
            for g, gc in zip(genes, genes_can):
                gene_to_raw.setdefault(gc, Counter()).update([g])
        top_genes = []
        for gc, _ in gene_counter.most_common(top_genes_k):
            if gc in gene_to_raw and gene_to_raw[gc]:
                top_genes.append(gene_to_raw[gc].most_common(1)[0][0])
        n_seen = len(sub)
        rows.append({
            "module_name_canonical": can,
            "n_raw_spellings":       len(raw_counter),
            "raw_spellings":         raw_spellings,
            "n_runs":                n_seen,
            "fraction_runs":         n_seen / n_runs if n_runs else 0.0,
            "modal_role_set":        modal_role,
            "n_distinct_genes":      len(gene_counter),
            "top_genes":             ", ".join(top_genes),
        })
    cols = ["module_name_canonical", "n_raw_spellings", "raw_spellings", "n_runs",
            "fraction_runs", "modal_role_set", "n_distinct_genes", "top_genes"]
    if not rows:
        return pd.DataFrame(columns=cols)
    return (pd.DataFrame(rows, columns=cols)
            .sort_values(["n_runs", "module_name_canonical"], ascending=[False, True])
            .reset_index(drop=True))


def unique_axis_names_table(axis_df: pd.DataFrame, n_runs: int,
                            mod_can_to_raw: dict[str, str]) -> pd.DataFrame:
    """One row per canonical axis name; columns mirror unique_module_names_table."""
    rows: list[dict] = []
    for can, sub in axis_df.groupby("axis_name_can"):
        if not can:
            continue
        raw_counter = Counter(sub["axis_name"])
        raw_spellings = " | ".join(f"{r} ({n})" for r, n in raw_counter.most_common())
        kind_counter = Counter(sub["kind"])
        modal_kind, _ = kind_counter.most_common(1)[0]
        pos_counter  = Counter(p for p in sub["positive_can"] if p)
        modal_pos_can  = pos_counter.most_common(1)[0][0]  if pos_counter  else ""
        gate_counter = Counter(g for g in sub["gate_can"] if g)
        modal_gate_can = gate_counter.most_common(1)[0][0] if gate_counter else ""
        rows.append({
            "axis_name_canonical":  can,
            "n_raw_spellings":      len(raw_counter),
            "raw_spellings":        raw_spellings,
            "n_runs":               len(sub),
            "fraction_runs":        len(sub) / n_runs if n_runs else 0.0,
            "modal_kind":           modal_kind,
            "modal_positive":       mod_can_to_raw.get(modal_pos_can, modal_pos_can),
            "modal_gate":           mod_can_to_raw.get(modal_gate_can, modal_gate_can),
        })
    cols = ["axis_name_canonical", "n_raw_spellings", "raw_spellings", "n_runs",
            "fraction_runs", "modal_kind", "modal_positive", "modal_gate"]
    if not rows:
        return pd.DataFrame(columns=cols)
    return (pd.DataFrame(rows, columns=cols)
            .sort_values(["n_runs", "axis_name_canonical"], ascending=[False, True])
            .reset_index(drop=True))


def top_hubs(axis_df: pd.DataFrame, mod_can_to_raw: dict[str, str],
             k: int = 5) -> dict[str, list[dict]]:
    """Most-reused modules across the cohort: as positive, gate, and negative.

    Counts are over (run × axis), so a module that sits as `positive` in 8
    axes across 30 runs scores up to 240. Useful for "which gene module is
    most central to the embedding" — human-readable hub identification.
    """
    pos_counts  = Counter(p for p in axis_df["positive_can"] if p)
    gate_counts = Counter(g for g in axis_df["gate_can"] if g)
    neg_counts  = Counter(n for negs in axis_df["negatives_can"] for n in negs if n)

    def _fmt(counter: Counter) -> list[dict]:
        return [{
            "module":        mod_can_to_raw.get(can, can),
            "module_can":    can,
            "n_axes_using":  int(n),
        } for can, n in counter.most_common(k)]

    return {
        "top_positives": _fmt(pos_counts),
        "top_gates":     _fmt(gate_counts),
        "top_negatives": _fmt(neg_counts),
    }


def top_resources(res_recur: pd.DataFrame, k: int = 5) -> list[dict]:
    """Top-k cited DOIs across the cohort (citation trail)."""
    if res_recur.empty:
        return []
    out = []
    for _, r in res_recur.head(k).iterrows():
        out.append({
            "doi":           str(r["doi"]),
            "citation":      str(r.get("citation", "") or ""),
            "kind":          str(r.get("kind", "") or ""),
            "n_runs":        int(r["n_runs"]),
            "fraction_runs": float(r["fraction_runs"]),
        })
    return out


# ============================================================
# Saturation
# ============================================================
def saturation_curve(items_by_run: dict[str, set], rng: np.random.Generator,
                     n_orderings: int = DEFAULT_SAT_ORDERINGS) -> dict[str, np.ndarray]:
    """Mean / SD of |union| after k runs over `n_orderings` random orderings."""
    runs = sorted(items_by_run)
    n = len(runs)
    if n == 0:
        return {"x": np.array([]), "mean": np.array([]), "std": np.array([])}
    curves = np.zeros((n_orderings, n), dtype=float)
    for k in range(n_orderings):
        order = rng.permutation(n)
        union: set = set()
        for i, idx in enumerate(order):
            union |= items_by_run[runs[idx]]
            curves[k, i] = len(union)
    return {
        "x":    np.arange(1, n + 1, dtype=float),
        "mean": curves.mean(axis=0),
        "std":  curves.std(axis=0),
    }


def saturation_slope_at_n(curve: dict[str, np.ndarray]) -> float:
    """Absolute new-item gain per added run at the boundary k = N."""
    m = curve.get("mean")
    if m is None or m.size < 2:
        return float("nan")
    return float(m[-1] - m[-2])


# ============================================================
# Composite
# ============================================================
def stability_composite(components: dict[str, float]) -> float:
    val = 0.0
    w_sum = 0.0
    for k, w in STABILITY_WEIGHTS.items():
        v = components.get(k, float("nan"))
        if isinstance(v, (int, float)) and np.isfinite(v):
            val += w * float(v)
            w_sum += w
    return val / w_sum if w_sum > 0 else float("nan")


# ============================================================
# Dataset-wide name universe (for the shared curated mappings)
# ============================================================
def gather_dataset_summaries(dataset: str) -> list[tuple[str, dict]]:
    """Every (cohort_name, summary) across ALL cohorts of the dataset, so the
    shared per-dataset name mapping stays complete regardless of any --cohort
    filter."""
    tagged: list[tuple[str, dict]] = []
    for cd in discover_cohorts(dataset, None):
        for rd in discover_runs(cd):
            s = load_summary(rd)
            if s is not None:
                tagged.append((cd.name, s))
    return tagged


def collect_name_context(tagged: list[tuple[str, dict]]):
    """Dataset-wide module / axis / resource universe with curation context.
    Returns (modules_df, axes_df, resources_df); one row per distinct raw name
    (resource_id for resources), with run/cohort support and a few descriptive
    columns the human reads when filling `standardized_name`."""
    mod: dict[str, dict] = {}
    ax: dict[str, dict] = {}
    res: dict[str, dict] = {}

    def _bump(store: dict, key: str, cohort: str, run_id) -> dict:
        e = store.setdefault(key, {"runs": set(), "cohorts": set()})
        e["runs"].add((cohort, run_id))
        e["cohorts"].add(cohort)
        return e

    for cohort, s in tagged:
        run_id = s.get("run_id")
        bp = s.get("blueprint") or {}
        for m in bp.get("modules") or []:
            name = m.get("name") or ""
            if name:
                e = _bump(mod, name, cohort, run_id)
                e.setdefault("genes", Counter()).update(m.get("genes") or [])
        for a in bp.get("axes") or []:
            name = a.get("name") or ""
            if name:
                e = _bump(ax, name, cohort, run_id)
                e.setdefault("kinds", Counter())[a.get("kind") or ""] += 1
                e.setdefault("positives", Counter())[a.get("positive") or ""] += 1
        for r in bp.get("resources") or []:
            rid = resource_identifier(r)
            if rid:
                e = _bump(res, rid, cohort, run_id)
                e.setdefault("citation", r.get("citation") or "")
                e.setdefault("kind", r.get("kind") or "")
                e.setdefault("doi", normalize_doi(r.get("doi")))

    def _common(key_col: str, key: str, e: dict) -> dict:
        return {key_col: key, "n_runs": len(e["runs"]),
                "n_cohorts": len(e["cohorts"]),
                "cohorts": ", ".join(sorted(e["cohorts"]))}

    mod_rows = [{**_common("raw_name", name, e), "canonical": canonicalize(name),
                 "top_genes": ", ".join(g for g, _ in e["genes"].most_common(8))}
                for name, e in mod.items()]
    ax_rows = [{**_common("raw_name", name, e), "canonical": canonicalize(name),
                "modal_kind": (e["kinds"].most_common(1)[0][0] if e["kinds"] else ""),
                "modal_positive": (e["positives"].most_common(1)[0][0] if e["positives"] else "")}
               for name, e in ax.items()]
    res_rows = [{**_common("resource_id", rid, e), "citation": e.get("citation", ""),
                 "kind": e.get("kind", ""), "doi": e.get("doi", "")}
                for rid, e in res.items()]

    mod_df = pd.DataFrame(mod_rows, columns=["raw_name", "canonical", "n_runs",
                                             "n_cohorts", "cohorts", "top_genes"])
    ax_df = pd.DataFrame(ax_rows, columns=["raw_name", "canonical", "n_runs",
                                           "n_cohorts", "cohorts", "modal_kind",
                                           "modal_positive"])
    res_df = pd.DataFrame(res_rows, columns=["resource_id", "n_runs", "n_cohorts",
                                             "cohorts", "citation", "kind", "doi"])
    for d, key in [(mod_df, "raw_name"), (ax_df, "raw_name"), (res_df, "resource_id")]:
        if not d.empty:
            d.sort_values(["n_runs", key], ascending=[False, True], inplace=True)
    return (mod_df.reset_index(drop=True), ax_df.reset_index(drop=True),
            res_df.reset_index(drop=True))


# ============================================================
# Cohort driver
# ============================================================
def compute_cohort(cohort_dir: Path, dataset: str, out_root: Path,
                   overwrite: bool, rng: np.random.Generator,
                   mod_map: dict | None = None, axis_map: dict | None = None,
                   res_map: dict | None = None) -> dict | None:
    cohort_name = cohort_dir.name
    out_dir = out_root / cohort_name
    cohort_summary_path = out_dir / "cohort_summary.json"
    if cohort_summary_path.exists() and not overwrite:
        raise FileExistsError(
            f"{cohort_summary_path} exists. Pass --overwrite to replace it.")

    run_dirs = discover_runs(cohort_dir)
    summaries: list[dict] = []
    for d in run_dirs:
        s = load_summary(d)
        if s is not None:
            summaries.append(s)
    if not summaries:
        print(f"[warn] {cohort_name}: no summary.json files found; skipping",
              file=sys.stderr)
        return None

    # If multiple scenario groups are present, keep the largest and warn.
    groups: dict[tuple, list[dict]] = {}
    for s in summaries:
        groups.setdefault(scenario_key(s.get("scenario") or {}), []).append(s)
    if len(groups) > 1:
        chosen = max(groups, key=lambda k: len(groups[k]))
        kept = len(groups[chosen])
        dropped = len(summaries) - kept
        print(f"[warn] {cohort_name}: {len(groups)} scenario groups present; "
              f"keeping the largest (n={kept}, dropped={dropped})",
              file=sys.stderr)
        summaries = groups[chosen]

    out_dir.mkdir(parents=True, exist_ok=True)

    mod_map  = mod_map or {}
    axis_map = axis_map or {}
    res_map  = res_map or {}
    mod_resolver  = make_resolver(mod_map)
    axis_resolver = make_resolver(axis_map)
    res_resolver  = make_resource_resolver(res_map)

    run_df    = build_run_df(summaries)
    module_df = build_module_df(summaries, mod_resolver)
    axis_df   = build_axis_df(summaries, axis_resolver, mod_resolver)
    n_runs = len(run_df)

    # A cohort whose runs declare no gene modules (e.g. non-blueprint baseline
    # runs that still drop a summary.json) carries no blueprint-variability
    # signal and would crash the recurrence/Jaccard machinery below on the
    # missing `module_name` column. Skip it, mirroring the no-summaries skip.
    if module_df.empty:
        print(f"[warn] {cohort_name}: no gene modules in any run; skipping",
              file=sys.stderr)
        return None

    # Match key → display spelling. Curated standardized names are shown verbatim
    # (`_std_display_map`); otherwise the most-common raw form is used, so output
    # tables / plots stay readable even though matching is on the resolved key.
    mod_can_to_raw  = {**_can_to_raw_map(module_df["module_name"], module_df["module_name_can"]),
                       **_std_display_map(mod_map)}
    axis_can_to_raw = {**_can_to_raw_map(axis_df["axis_name"], axis_df["axis_name_can"]),
                       **_std_display_map(axis_map)}
    gene_pairs_raw, gene_pairs_can = [], []
    for raws, cans in zip(module_df["genes"], module_df["genes_can"]):
        gene_pairs_raw.extend(raws)
        gene_pairs_can.extend(cans)
    gene_can_to_raw = _can_to_raw_map(gene_pairs_raw, gene_pairs_can)

    mod_recur = module_recurrence(module_df, n_runs, mod_can_to_raw)
    gen_recur = gene_recurrence(module_df,   n_runs, gene_can_to_raw)
    res_recur = resource_recurrence(n_runs, summaries, res_resolver)
    ax_recur  = axis_recurrence(axis_df,    n_runs, axis_can_to_raw)

    geneset_stab = module_geneset_stability(module_df, mod_can_to_raw)
    role_cons    = module_role_consistency(module_df, mod_can_to_raw)
    kind_cons    = axis_kind_consistency(axis_df, axis_can_to_raw)
    struct_cons  = axis_structure_consistency(axis_df, axis_can_to_raw, mod_can_to_raw)
    near_dup     = module_name_near_duplicates(module_df, mod_can_to_raw)
    unique_modules_tbl = unique_module_names_table(module_df, n_runs)
    unique_axes_tbl    = unique_axis_names_table(axis_df, n_runs, mod_can_to_raw)

    genes_by_run    = _union_per_run(module_df, "genes_can",       flatten=True)
    dois_by_run     = _union_per_run(module_df, "resource_dois",   flatten=True)
    modnames_by_run = _union_per_run(module_df, "module_name_can", flatten=False)
    axisnames_by_run = _union_per_run(axis_df,  "axis_name_can",   flatten=False)
    resids_by_run    = _resources_by_run(summaries, res_resolver)  # all kinds incl. databases

    jacc_genes     = jaccard_matrix(genes_by_run)
    jacc_dois      = jaccard_matrix(dois_by_run)       # papers + atlases (DOI only)
    jacc_resources = jaccard_matrix(resids_by_run)     # all kinds incl. databases
    jacc_modnames  = jaccard_matrix(modnames_by_run)
    jacc_axisnames = jaccard_matrix(axisnames_by_run)

    boot_genes    = bootstrap_offdiag_ci(jacc_genes,     rng)
    boot_dois     = bootstrap_offdiag_ci(jacc_dois,      rng)
    boot_resources = bootstrap_offdiag_ci(jacc_resources, rng)
    boot_modnames = bootstrap_offdiag_ci(jacc_modnames,  rng)
    # Axis-name set Jaccard: reported for completeness (genes/modules/axes/resources
    # all get a pairwise overlap), but intentionally NOT in STABILITY_WEIGHTS — the
    # composite already carries axis structure via mean_axis_positive_agreement, and
    # adding this term would shift the established composite values.
    boot_axisnames = bootstrap_offdiag_ci(jacc_axisnames, rng)

    sat_genes   = saturation_curve(genes_by_run,    rng)
    sat_modules = saturation_curve(modnames_by_run, rng)
    sat_dois    = saturation_curve(dois_by_run,     rng)

    # ---- per-cohort tables
    run_df.to_csv(out_dir / "run_summary.csv", index=False)
    md = module_df.copy()
    md["genes"]         = md["genes"].apply(lambda xs: "|".join(xs))
    md["genes_can"]     = md["genes_can"].apply(lambda xs: "|".join(xs))
    md["resource_keys"] = md["resource_keys"].apply(lambda xs: "|".join(xs))
    md["resource_dois"] = md["resource_dois"].apply(lambda xs: "|".join(xs))
    md.to_csv(out_dir / "modules_long.csv", index=False)
    ad = axis_df.copy()
    ad["negatives"]     = ad["negatives"].apply(lambda xs: "|".join(xs))
    ad["negatives_can"] = ad["negatives_can"].apply(lambda xs: "|".join(xs))
    ad.to_csv(out_dir / "axes_long.csv", index=False)
    mod_recur.to_csv(out_dir / "module_recurrence.csv",         index=False)
    gen_recur.to_csv(out_dir / "gene_recurrence.csv",           index=False)
    ax_recur .to_csv(out_dir / "axis_recurrence.csv",           index=False)
    res_recur.to_csv(out_dir / "resource_recurrence.csv",       index=False)
    geneset_stab.to_csv(out_dir / "module_geneset_stability.csv", index=False)
    role_cons   .to_csv(out_dir / "module_role_consistency.csv",  index=False)
    kind_cons   .to_csv(out_dir / "axis_kind_consistency.csv",    index=False)
    struct_cons .to_csv(out_dir / "axis_structure_consistency.csv", index=False)
    near_dup    .to_csv(out_dir / "module_name_near_duplicates.csv", index=False)
    unique_modules_tbl.to_csv(out_dir / "unique_modules.csv", index=False)
    unique_axes_tbl   .to_csv(out_dir / "unique_axes.csv",    index=False)
    jacc_genes    .to_csv(out_dir / "jaccard_genes.csv")
    jacc_dois     .to_csv(out_dir / "jaccard_dois.csv")
    jacc_resources.to_csv(out_dir / "jaccard_resources.csv")
    jacc_modnames .to_csv(out_dir / "jaccard_module_names.csv")
    jacc_axisnames.to_csv(out_dir / "jaccard_axis_names.csv")
    sat_rows: list[dict] = []
    for kind, curve in [("genes", sat_genes), ("modules", sat_modules), ("dois", sat_dois)]:
        for x, mu, sd in zip(curve["x"], curve["mean"], curve["std"]):
            sat_rows.append({"kind": kind, "k": int(x),
                             "mean": float(mu), "std": float(sd)})
    pd.DataFrame(sat_rows, columns=["kind", "k", "mean", "std"]).to_csv(
        out_dir / "saturation.csv", index=False)

    # ---- cohort-level scalars
    parts = parse_cohort_name(cohort_name)
    size_cols = ["n_modules", "n_axes", "n_unique_genes", "n_resources",
                 "n_dois_unique", "gene_redundancy", "mean_genes_per_module"]
    size_stats: dict[str, dict[str, float]] = {}
    for col in size_cols:
        s = run_df[col]
        size_stats[col] = {
            "mean": float(s.mean()),
            "std":  float(s.std(ddof=1)) if len(s) > 1 else 0.0,
            "min":  float(s.min()),
            "max":  float(s.max()),
            "cv":   cv(s),
        }

    rec_modules   = recurrence_concentration(mod_recur, n_runs)
    rec_axes      = recurrence_concentration(ax_recur,  n_runs)
    rec_genes     = recurrence_concentration(gen_recur, n_runs)
    rec_resources = recurrence_concentration(res_recur, n_runs)

    mean_kind_agreement = (float(kind_cons["kind_agreement"].mean())
                           if not kind_cons.empty else float("nan"))
    mean_role_agreement = (float(role_cons["role_agreement"].mean())
                           if not role_cons.empty else float("nan"))
    mean_geneset_jac    = (float(geneset_stab["mean_jaccard"].mean())
                           if not geneset_stab.empty else float("nan"))
    if not struct_cons.empty:
        mean_axis_pos  = float(struct_cons["positive_agreement"].mean())
        mean_axis_gate = float(struct_cons["gate_agreement"].mean())
        mean_axis_neg  = float(struct_cons["mean_negatives_jaccard"].mean(skipna=True))
    else:
        mean_axis_pos = mean_axis_gate = mean_axis_neg = float("nan")

    # DOI-Jaccard (boot_dois) is intentionally excluded — see STABILITY_WEIGHTS.
    composite_inputs = {
        "jaccard_genes_off_diag":              boot_genes["mean"],
        "jaccard_modnames_off_diag":           boot_modnames["mean"],
        "mean_kind_agreement":                 mean_kind_agreement,
        "mean_role_agreement":                 mean_role_agreement,
        "mean_within_module_geneset_jaccard":  mean_geneset_jac,
        "mean_axis_positive_agreement":        mean_axis_pos,
    }
    composite = stability_composite(composite_inputs)

    sat_slope = {
        "genes":   saturation_slope_at_n(sat_genes),
        "modules": saturation_slope_at_n(sat_modules),
        "dois":    saturation_slope_at_n(sat_dois),
    }

    # Coverage: total unique items observed across the entire cohort.
    coverage = {
        "n_unique_modules":   rec_modules["n_unique_total"],
        "n_unique_axes":      rec_axes["n_unique_total"],
        "n_unique_genes":     rec_genes["n_unique_total"],
        "n_unique_resources": rec_resources["n_unique_total"],
    }

    # Resource-kind composition: mean / SD of paper / atlas / database fractions
    # across runs. Captures the "what does the agent cite" mix.
    def _series_mean_std(s: pd.Series) -> dict[str, float]:
        return {
            "mean": float(s.mean()),
            "std":  float(s.std(ddof=1)) if len(s) > 1 else 0.0,
        }
    resource_kind_composition = {
        "frac_paper":    _series_mean_std(run_df["frac_resources_paper"]),
        "frac_atlas":    _series_mean_std(run_df["frac_resources_atlas"]),
        "frac_database": _series_mean_std(run_df["frac_resources_database"]),
    }

    hubs = top_hubs(axis_df, mod_can_to_raw, k=5)
    top_res = top_resources(res_recur, k=5)

    # Canonicalization merges: how many concepts had ≥ 2 distinct raw spellings
    # across runs, collapsed by `canonicalize()` before the rest of the pipeline.
    # High values indicate the agent isn't picking a single canonical name.
    n_merged_modules = _collisions_count(module_df["module_name"], module_df["module_name_can"])
    n_merged_axes    = _collisions_count(axis_df["axis_name"],     axis_df["axis_name_can"])
    n_merged_genes   = _collisions_count(gene_pairs_raw,           gene_pairs_can)

    # n_drift: near-duplicate canonical pairs that never co-occur (higher-signal
    # drift; see module_name_near_duplicates docstring).
    n_drift = (int((near_dup["co_occurring_runs"] == 0).sum())
               if not near_dup.empty else 0)
    sanity = {
        "frac_with_placeholders":         float(run_df["placeholder_strings_present"].mean()),
        "mean_orphaned_modules":          float(run_df["n_orphaned_modules"].mean()),
        "mean_modules_missing_resources": float(run_df["n_modules_missing_resources"].mean()),
        "mean_unused_resources":          float(run_df["n_unused_resources"].mean()),
        "frac_clip_overridden":           float(run_df["set_clip_overridden"].mean()),
        "n_modules_merged_by_canonicalization": n_merged_modules,
        "n_axes_merged_by_canonicalization":    n_merged_axes,
        "n_genes_merged_by_canonicalization":   n_merged_genes,
        "n_module_name_drift_candidates": n_drift,
        "n_module_name_near_duplicates":  int(len(near_dup)),
        "doi_kind_agreement":             doi_kind_agreement(summaries),
    }

    cohort_summary = {
        "cohort": cohort_name,
        **parts,
        "n_runs": n_runs,
        "scenario": run_df.iloc[0][["species", "tissue", "celltype"]].to_dict(),
        "run_ids": run_df["run_id"].tolist(),
        "size": size_stats,
        "jaccard": {
            "genes":     boot_genes,
            "dois":      boot_dois,         # papers + atlases (DOI-keyed grounding)
            "resources": boot_resources,    # all kinds incl. databases (key-else-DOI)
            "modnames":  boot_modnames,
            "axisnames": boot_axisnames,
        },
        "structural_consistency": {
            "mean_kind_agreement":                mean_kind_agreement,
            "mean_role_agreement":                mean_role_agreement,
            "mean_within_module_geneset_jaccard": mean_geneset_jac,
            "mean_axis_positive_agreement":       mean_axis_pos,
            "mean_axis_gate_agreement":           mean_axis_gate,
            "mean_axis_negatives_jaccard":        mean_axis_neg,
        },
        "recurrence": {
            "modules":   rec_modules,
            "axes":      rec_axes,
            "genes":     rec_genes,
            "resources": rec_resources,
        },
        "saturation_slope_at_N":      sat_slope,
        "coverage":                   coverage,
        "resource_kind_composition":  resource_kind_composition,
        "top_hubs":                   hubs,
        "top_resources":              top_res,
        "stability_composite":        composite,
        "stability_weights":          STABILITY_WEIGHTS,
        "sanity":                     sanity,
    }
    with cohort_summary_path.open("w") as fh:
        json.dump(cohort_summary, fh, indent=2, default=str)
        fh.write("\n")
    write_report(cohort_summary, out_dir / "report.md",
                 geneset_stab, struct_cons, near_dup)

    # ---- flat row for cross-cohort cohort_metrics.csv
    flat: dict[str, Any] = {"cohort": cohort_name, **parts, "n_runs": n_runs}
    for col, st in size_stats.items():
        flat[f"{col}_mean"] = st["mean"]
        flat[f"{col}_std"]  = st["std"]
        flat[f"{col}_cv"]   = st["cv"]
    for k, b in [("genes", boot_genes), ("dois", boot_dois),
                 ("resources", boot_resources), ("modnames", boot_modnames),
                 ("axisnames", boot_axisnames)]:
        flat[f"jaccard_{k}_off_diag"] = b["mean"]
        flat[f"jaccard_{k}_ci_lo"]    = b["ci_lo"]
        flat[f"jaccard_{k}_ci_hi"]    = b["ci_hi"]
    flat.update({
        "mean_kind_agreement":                mean_kind_agreement,
        "mean_role_agreement":                mean_role_agreement,
        "mean_within_module_geneset_jaccard": mean_geneset_jac,
        "mean_axis_positive_agreement":       mean_axis_pos,
        "mean_axis_gate_agreement":           mean_axis_gate,
        "mean_axis_negatives_jaccard":        mean_axis_neg,
    })
    for kind, rec in [("modules", rec_modules), ("axes", rec_axes),
                      ("genes", rec_genes), ("resources", rec_resources)]:
        flat[f"frac_{kind}_ge80pct"] = rec["frac_in_ge80pct"]
        flat[f"gini_{kind}"]         = rec["gini"]
        flat[f"n_universal_{kind}"]  = rec["n_universal"]
    flat["saturation_slope_genes"]   = sat_slope["genes"]
    flat["saturation_slope_modules"] = sat_slope["modules"]
    flat["saturation_slope_dois"]    = sat_slope["dois"]
    # Coverage + resource-kind + top hub + top resource (cross-cohort comparable)
    flat["coverage_n_unique_modules"]   = coverage["n_unique_modules"]
    flat["coverage_n_unique_axes"]      = coverage["n_unique_axes"]
    flat["coverage_n_unique_genes"]     = coverage["n_unique_genes"]
    flat["coverage_n_unique_resources"] = coverage["n_unique_resources"]
    flat["frac_resources_paper_mean"]    = resource_kind_composition["frac_paper"]["mean"]
    flat["frac_resources_atlas_mean"]    = resource_kind_composition["frac_atlas"]["mean"]
    flat["frac_resources_database_mean"] = resource_kind_composition["frac_database"]["mean"]
    flat["top_positive_module"] = hubs["top_positives"][0]["module"] if hubs["top_positives"] else ""
    flat["top_gate_module"]     = hubs["top_gates"][0]["module"]     if hubs["top_gates"]     else ""
    flat["top_resource_doi"]      = top_res[0]["doi"]      if top_res else ""
    flat["top_resource_citation"] = top_res[0]["citation"] if top_res else ""
    flat["top_resource_n_runs"]   = top_res[0]["n_runs"]   if top_res else 0
    flat.update(sanity)
    flat["stability_composite"] = composite
    return flat


# ============================================================
# Reporting
# ============================================================
def _fmt_size(d: dict[str, float]) -> str:
    return (f"{d['mean']:.2f} ± {d['std']:.2f}  "
            f"[{d['min']:.2f}, {d['max']:.2f}], CV={d['cv']:.3f}")


def _fmt_jacc(d: dict[str, float]) -> str:
    if not np.isfinite(d.get("mean", float("nan"))):
        return "n/a"
    return f"{d['mean']:.3f}  (95% CI [{d['ci_lo']:.3f}, {d['ci_hi']:.3f}])"


def write_report(cs: dict, path: Path,
                 geneset_stab, struct_cons, near_dup) -> None:
    L: list[str] = []
    L.append(f"# Blueprint variability — {cs['cohort']}  ({cs['n_runs']} runs)\n")
    sc = cs["scenario"]
    L.append(f"**Scenario:** species={sc.get('species')}, tissue={sc.get('tissue')}, "
             f"celltype={sc.get('celltype')}\n")
    L.append(f"**Stability composite:** **{cs['stability_composite']:.3f}**  "
             f"_(equal-weighted; components in `cohort_summary.json`)_\n")

    L.append("## Blueprint size (mean ± SD [min, max], CV)\n")
    for k in ["n_modules", "n_axes", "n_unique_genes", "n_resources",
              "n_dois_unique", "gene_redundancy", "mean_genes_per_module"]:
        L.append(f"- **{k}**: {_fmt_size(cs['size'][k])}")

    L.append("\n## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)\n")
    L.append(f"- **genes**:        {_fmt_jacc(cs['jaccard']['genes'])}")
    L.append(f"- **DOIs**:         {_fmt_jacc(cs['jaccard']['dois'])}")
    L.append(f"- **module names**: {_fmt_jacc(cs['jaccard']['modnames'])}")

    L.append("\n## Structural consistency (per-name mean agreement, ∈ [0,1])\n")
    sc2 = cs["structural_consistency"]
    L.append(f"- mean axis-kind agreement:        {sc2['mean_kind_agreement']:.3f}")
    L.append(f"- mean axis positive agreement:    {sc2['mean_axis_positive_agreement']:.3f}")
    L.append(f"- mean axis gate agreement:        {sc2['mean_axis_gate_agreement']:.3f}")
    L.append(f"- mean axis negatives Jaccard:     {sc2['mean_axis_negatives_jaccard']:.3f}")
    L.append(f"- mean module role-set agreement:  {sc2['mean_role_agreement']:.3f}")
    L.append(f"- mean within-module gene Jaccard: {sc2['mean_within_module_geneset_jaccard']:.3f}")

    L.append("\n## Recurrence (universal / majority / singleton / unique-total; "
             "frac in ≥80 %; Gini)\n")
    for kind in ["modules", "axes", "genes", "resources"]:
        r = cs["recurrence"][kind]
        L.append(f"- **{kind}**: "
                 f"{r['n_universal']} / {r['n_majority']} / {r['n_singleton']} / "
                 f"{r['n_unique_total']}; "
                 f"≥80 %: {r['frac_in_ge80pct']*100:.1f} %; "
                 f"Gini: {r['gini']:.3f}")

    L.append("\n## Saturation slope at N "
             "(absolute new-item gain when adding the N-th run)\n")
    s = cs["saturation_slope_at_N"]
    L.append(f"- genes:    {s['genes']:.2f}")
    L.append(f"- modules:  {s['modules']:.2f}")
    L.append(f"- DOIs:     {s['dois']:.2f}")

    cov = cs["coverage"]
    L.append("\n## Cohort coverage (totals across runs)\n")
    L.append(f"- unique modules:   **{cov['n_unique_modules']}**")
    L.append(f"- unique axes:      **{cov['n_unique_axes']}**")
    L.append(f"- unique genes:     **{cov['n_unique_genes']}**")
    L.append(f"- unique resources (DOI): **{cov['n_unique_resources']}**")

    rkc = cs["resource_kind_composition"]
    L.append("\n## Resource-kind composition (mean fraction per run ± SD)\n")
    L.append(f"- paper:    {rkc['frac_paper']['mean']:.3f} ± {rkc['frac_paper']['std']:.3f}")
    L.append(f"- atlas:    {rkc['frac_atlas']['mean']:.3f} ± {rkc['frac_atlas']['std']:.3f}")
    L.append(f"- database: {rkc['frac_database']['mean']:.3f} ± {rkc['frac_database']['std']:.3f}")

    th = cs["top_hubs"]
    L.append("\n## Top hubs (most-reused modules across all runs × axes)\n")

    def _hub_block(label: str, items: list) -> None:
        L.append(f"### As {label}")
        if not items:
            L.append("- _(none)_")
            return
        for r in items:
            L.append(f"- `{r['module']}` — used in {r['n_axes_using']} (run × axis) slots")

    _hub_block("positive (axis identity carrier)", th["top_positives"])
    _hub_block("gate (axis activation switch)",    th["top_gates"])
    _hub_block("negative (axis contrast partner)", th["top_negatives"])

    tr = cs["top_resources"]
    if tr:
        L.append("\n## Top resources (most-cited DOIs)\n")
        for r in tr:
            cite = r["citation"] or r["doi"]
            L.append(f"- {cite}  (`{r['doi']}`, {r['kind']}) — "
                     f"{r['n_runs']}/{cs['n_runs']} runs "
                     f"({r['fraction_runs']*100:.0f} %)")

    sn = cs["sanity"]
    L.append("\n## Sanity\n")
    L.append(f"- runs with placeholders:                {sn['frac_with_placeholders']*100:.1f} %")
    L.append(f"- mean orphaned modules / run:           {sn['mean_orphaned_modules']:.2f}")
    L.append(f"- mean modules missing resources / run:  {sn['mean_modules_missing_resources']:.2f}")
    L.append(f"- mean unused resources / run:           {sn['mean_unused_resources']:.2f}")
    L.append(f"- runs with clip override:               {sn['frac_clip_overridden']*100:.1f} %")
    L.append("- canonicalization merges (concepts with ≥ 2 distinct raw spellings): "
             f"modules={sn['n_modules_merged_by_canonicalization']}, "
             f"axes={sn['n_axes_merged_by_canonicalization']}, "
             f"genes={sn['n_genes_merged_by_canonicalization']}")
    L.append(f"- post-canonical near-duplicate pairs (edit ≤ 2): "
             f"{sn['n_module_name_near_duplicates']}; "
             f"of which {sn['n_module_name_drift_candidates']} never co-occur "
             "(higher-signal drift candidates).")
    L.append(f"- DOI-kind agreement:                    {sn['doi_kind_agreement']:.3f}")

    drift = (near_dup[near_dup["co_occurring_runs"] == 0]
             if not near_dup.empty else near_dup)
    if not drift.empty:
        L.append("\n### Post-canonical near-duplicate names (review — possible typos)\n")
        for _, r in drift.head(15).iterrows():
            L.append(f"- `{r['name_a']}` ↔ `{r['name_b']}` "
                     f"(canonical `{r['name_a_can']}` / `{r['name_b_can']}`, "
                     f"edit={r['edit_distance']}, "
                     f"runs_a={r['n_runs_a']}, runs_b={r['n_runs_b']})")

    if not struct_cons.empty:
        L.append("\n## Top axes by structural drift (low = unstable)\n")
        L.append("| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |")
        L.append("|---|---|---|---|---|")
        for _, r in struct_cons.head(15).iterrows():
            mn = r["mean_negatives_jaccard"]
            mn_s = f"{mn:.2f}" if np.isfinite(mn) else "—"
            L.append(f"| `{r['axis_name']}` | {int(r['n_runs_seen'])} | "
                     f"{r['positive_agreement']:.2f} | "
                     f"{r['gate_agreement']:.2f} | {mn_s} |")

    if not geneset_stab.empty:
        L.append("\n## Least-stable modules (within-name gene Jaccard)\n")
        for _, r in geneset_stab.head(10).iterrows():
            L.append(f"- `{r['module_name']}`: mean Jaccard {r['mean_jaccard']:.3f} "
                     f"(seen in {int(r['n_runs_seen'])} runs)")

    path.write_text("\n".join(L) + "\n")


# ============================================================
# Main
# ============================================================
def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description="Blueprint variability metrics across agent replicates per cohort.")
    ap.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    ap.add_argument("--cohort", type=str, default=None,
                    help="Restrict to a single work/<cohort>/ dir name.")
    ap.add_argument("--overwrite", action="store_true")
    return ap.parse_args()


def main() -> None:
    args = parse_args()
    cohort_dirs = discover_cohorts(args.dataset, args.cohort)
    if not cohort_dirs:
        raise SystemExit(f"No cohorts found for --dataset {args.dataset}.")

    out_root = REPO_ROOT / "outs" / "blueprint_variability" / args.dataset
    out_root.mkdir(parents=True, exist_ok=True)

    # ---- shared per-dataset name mappings (modules / axes / resources).
    # Scan ALL cohorts (independent of --cohort) so the curated mapping stays
    # complete; merge-append new names into the tracked configs/ files (human
    # edits preserved), then write a context snapshot under outs/ to curate
    # against. Blank standardized_name ⇒ canonical / resource_identifier fallback.
    tagged = gather_dataset_summaries(args.dataset)
    mod_ctx, axis_ctx, res_ctx = collect_name_context(tagged)
    cfg_dir = CONFIG_ROOT / args.dataset
    mod_map_path  = cfg_dir / "module_names.csv"
    axis_map_path = cfg_dir / "axis_names.csv"
    res_map_path  = cfg_dir / "resource_names.csv"
    merge_append_mapping(mod_map_path,  mod_ctx["raw_name"],     key_col="raw_name")
    merge_append_mapping(axis_map_path, axis_ctx["raw_name"],    key_col="raw_name")
    merge_append_mapping(res_map_path,  res_ctx["resource_id"],  key_col="resource_id")
    mod_map  = load_name_mapping(mod_map_path)
    axis_map = load_name_mapping(axis_map_path)
    res_map  = load_name_mapping(res_map_path)

    name_out = out_root / "name_mapping"
    name_out.mkdir(parents=True, exist_ok=True)
    mod_ctx["standardized_resolved"]  = mod_ctx["raw_name"].map(make_resolver(mod_map))
    axis_ctx["standardized_resolved"] = axis_ctx["raw_name"].map(make_resolver(axis_map))
    res_ctx["standardized_resolved"]  = res_ctx["resource_id"].map(
        lambda i: canonicalize(res_map[i]) if res_map.get(i) else i)
    mod_ctx.to_csv(name_out / "extracted_modules.csv",   index=False)
    axis_ctx.to_csv(name_out / "extracted_axes.csv",     index=False)
    res_ctx.to_csv(name_out / "extracted_resources.csv", index=False)
    n_curated = (sum(bool(v) for v in mod_map.values()),
                 sum(bool(v) for v in axis_map.values()),
                 sum(bool(v) for v in res_map.values()))
    print(f"[ok] name mappings ({cfg_dir}): "
          f"modules {len(mod_ctx)} names ({n_curated[0]} curated), "
          f"axes {len(axis_ctx)} ({n_curated[1]} curated), "
          f"resources {len(res_ctx)} ({n_curated[2]} curated)")

    flats: list[dict] = []
    for cd in cohort_dirs:
        # Per-cohort RNG: rerunning on the same cohort gives deterministic outputs.
        rng = np.random.default_rng(RNG_SEED)
        flat = compute_cohort(cd, args.dataset, out_root, args.overwrite, rng,
                              mod_map=mod_map, axis_map=axis_map, res_map=res_map)
        if flat is None:
            continue
        flats.append(flat)
        print(f"[ok] {cd.name}: n_runs={flat['n_runs']} "
              f"composite={flat['stability_composite']:.3f} "
              f"jaccard_genes={flat['jaccard_genes_off_diag']:.3f} "
              f"[{flat['jaccard_genes_ci_lo']:.3f}, {flat['jaccard_genes_ci_hi']:.3f}] "
              f"kind_agreement={flat['mean_kind_agreement']:.3f}")

    if flats:
        cohort_csv = out_root / "cohort_metrics.csv"
        if cohort_csv.exists() and not args.overwrite:
            raise FileExistsError(
                f"{cohort_csv} exists. Pass --overwrite to replace it.")
        df = pd.DataFrame(flats).sort_values("cohort").reset_index(drop=True)
        df.to_csv(cohort_csv, index=False)
        print(f"[ok] wrote {cohort_csv} ({len(flats)} cohorts)")


if __name__ == "__main__":
    main()
