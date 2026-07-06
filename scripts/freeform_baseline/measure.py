"""
Measure free-form agent baseline runs and compare to the constrained
blueprint cohort. For each free-form workdir, computes:

  - hard-failure: was the script written, did it run, did it produce a
    parquet with cell_id + dim_*?
  - dimensionality variance across seeds
  - provenance density: count of DOI-shaped strings in embed_freeform.py
  - gene-list extraction: regex heuristic on uppercase tokens, filtered
    against gene-symbol vocabularies discovered in the cohort itself
    (genes mentioned in any blueprint-cohort summary.json for the same
    dataset, used as a vocabulary anchor)
  - cross-seed gene Jaccard within the free-form cohort and vs the
    matching blueprint cohort

Embedding quality itself is scored separately by scripts/embedding_metrics/metrics.py
(from the parquets under outs/embeddings/<dataset>/, stem freeform_<...>_seed<N>);
this script does not compute it.

Inputs:
    work/freeform_<cohort>/<run>/{embed_freeform.py, .embed_status, ...}
    work/<blueprint_cohort>/<run>/summary.json   (for vocabulary + Jaccard)
    outs/metrics/<dataset>/metrics.csv           (for quality comparison)

Outputs:
    outs/freeform_baseline/<freeform_cohort>/measurement.csv
    outs/freeform_baseline/<freeform_cohort>/measurement_summary.md

Usage:
    python scripts/freeform_baseline/measure.py \
        --freeform-cohort freeform_gpt_5_4_mini_medium_human_pancreas \
        --blueprint-cohort gpt_5_4_mini_medium_human_pancreas \
        --dataset human_pancreas
"""
from __future__ import annotations

import argparse
import csv
import itertools
import json
import re
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]
DOI_PATTERN = re.compile(r"10\.\d{4,9}/[^\s\"',;)<>\[\]]+")
# Match both human ALL-CAPS (INS, NKX6-1) and mouse CamelCase (Pvalb, Mybpc1)
# gene symbols. The intersection with the cohort vocabulary filters non-gene
# uppercase tokens (Python constants, class names) downstream.
GENE_TOKEN_PATTERN = re.compile(r"\b[A-Z][A-Za-z0-9-]{1,15}\b")


def load_runs(cohort_dir: Path, run_glob: str) -> list[Path]:
    return sorted(d for d in cohort_dir.glob(run_glob) if d.is_dir())


def cohort_genes(cohort_dir: Path) -> set[str]:
    """Union of all genes appearing in any summary.json under this cohort."""
    out: set[str] = set()
    for run in load_runs(cohort_dir, "*"):
        sj = run / "summary.json"
        if not sj.is_file():
            continue
        try:
            data = json.loads(sj.read_text())
        except Exception:
            continue
        for m in data.get("blueprint", {}).get("modules", []):
            out.update(g for g in m.get("genes", []) if isinstance(g, str))
    return out


def freeform_run_genes(script_path: Path, vocabulary: set[str]) -> set[str]:
    """Extract gene tokens from a free-form script by intersecting uppercase
    tokens with the supplied vocabulary."""
    if not script_path.is_file():
        return set()
    try:
        text = script_path.read_text(errors="replace")
    except Exception:
        return set()
    candidates = set(GENE_TOKEN_PATTERN.findall(text))
    return candidates & vocabulary


def parquet_columns(parquet_path: Path) -> tuple[bool, int]:
    """Return (schema_ok, n_dim_columns). schema_ok = has cell_id + at least
    one dim_*. Robust to parquet read errors (returns (False, 0))."""
    if not parquet_path.is_file():
        return False, 0
    try:
        import pyarrow.parquet as pq
        schema = pq.read_schema(parquet_path)
    except Exception:
        return False, 0
    names = set(schema.names)
    dim_cols = sorted(n for n in names if n.startswith("dim_"))
    return ("cell_id" in names and bool(dim_cols)), len(dim_cols)


def jaccard(a: set, b: set) -> float:
    if not a and not b:
        return float("nan")
    return len(a & b) / max(1, len(a | b))


def mean_pairwise_jaccard(sets: list[set]) -> float:
    sets = [s for s in sets if s]
    if len(sets) < 2:
        return float("nan")
    pairs = list(itertools.combinations(sets, 2))
    return sum(jaccard(a, b) for a, b in pairs) / len(pairs)


def cross_pair_mean_jaccard(a_sets: list[set], b_sets: list[set]) -> float:
    a_sets = [s for s in a_sets if s]
    b_sets = [s for s in b_sets if s]
    if not a_sets or not b_sets:
        return float("nan")
    n = len(a_sets) * len(b_sets)
    return sum(jaccard(a, b) for a in a_sets for b in b_sets) / n


def per_run_record(run_dir: Path, vocabulary: set[str]) -> dict:
    name = run_dir.name
    script = run_dir / "embed_freeform.py"
    stage_a_rc = (run_dir / ".stage_a_rc").read_text().strip() if (run_dir / ".stage_a_rc").exists() else ""
    embed_status = (run_dir / ".embed_status").read_text().strip() if (run_dir / ".embed_status").exists() else ""
    embed_detail = (run_dir / ".embed_detail").read_text().strip() if (run_dir / ".embed_detail").exists() else ""

    script_present = script.is_file() and script.stat().st_size > 0
    script_size = script.stat().st_size if script_present else 0

    text = script.read_text(errors="replace") if script_present else ""
    n_dois = len(set(DOI_PATTERN.findall(text)))
    genes = freeform_run_genes(script, vocabulary)

    parquet_path = None
    schema_ok, n_dims = False, 0
    if embed_status == "OK" and embed_detail:
        candidate = REPO_ROOT / embed_detail
        if candidate.is_file():
            parquet_path = candidate
            schema_ok, n_dims = parquet_columns(parquet_path)

    return {
        "run":              name,
        "stage_a_rc":       stage_a_rc,
        "script_present":   int(script_present),
        "script_size":      script_size,
        "n_dois":           n_dois,
        "n_genes_extracted": len(genes),
        "embed_status":    embed_status,
        "embed_detail":    embed_detail,
        "parquet_path":     str(parquet_path.relative_to(REPO_ROOT)) if parquet_path else "",
        "parquet_schema_ok": int(schema_ok),
        "n_dim_columns":    n_dims,
        "_genes":           genes,
    }


def blueprint_run_genes(run_dir: Path) -> set[str]:
    sj = run_dir / "summary.json"
    if not sj.is_file():
        return set()
    try:
        data = json.loads(sj.read_text())
    except Exception:
        return set()
    out: set[str] = set()
    for m in data.get("blueprint", {}).get("modules", []):
        out.update(g for g in m.get("genes", []) if isinstance(g, str))
    return out


def quality_table(metrics_csv: Path, free_model: str, blueprint_model: str) -> pd.DataFrame:
    if not metrics_csv.is_file():
        return pd.DataFrame()
    df = pd.read_csv(metrics_csv)
    sub = df[df["model"].isin([free_model, blueprint_model])].copy()
    if sub.empty:
        return sub
    return sub


def render_summary(
    free_records: list[dict], free_jaccard_self: float, free_vs_bp_jaccard: float,
    bp_jaccard_self: float, free_dims: list[int], bp_dims: list[int],
    free_dois: list[int], bp_dois_per_run: list[int],
    quality_df: pd.DataFrame, args: argparse.Namespace,
) -> str:
    n_total = len(free_records)
    n_script = sum(r["script_present"] for r in free_records)
    n_ok = sum(1 for r in free_records if r["embed_status"] == "OK" and r["parquet_schema_ok"])
    fail_rate = 1.0 - (n_ok / n_total) if n_total else float("nan")

    def stats(xs):
        if not xs:
            return "—"
        s = pd.Series(xs)
        return f"{s.median():.2f} (min {s.min()}, max {s.max()}, sd {s.std():.2f})"

    lines = [
        f"# Free-form baseline measurement — `{args.freeform_cohort}`",
        "",
        f"- **Free-form cohort**: `{args.freeform_cohort}`",
        f"- **Blueprint reference cohort**: `{args.blueprint_cohort}`",
        f"- **Dataset**: `{args.dataset}`",
        f"- **Runs scanned**: {n_total}",
        "",
        "## Stage A and Finalize outcomes",
        "",
        f"- Scripts written:        {n_script}/{n_total}",
        f"- Successful parquets:    {n_ok}/{n_total}",
        f"- Hard-failure rate:      {fail_rate:.2f}",
        "",
        "## Cross-seed reproducibility (gene set)",
        "",
        f"- Free-form pairwise Jaccard (mean):   {free_jaccard_self:.3f}",
        f"- Blueprint pairwise Jaccard (mean):   {bp_jaccard_self:.3f}",
        f"- Free-form × blueprint mean Jaccard:  {free_vs_bp_jaccard:.3f}",
        "",
        "## Dimensionality (parquet `dim_*` columns)",
        "",
        f"- Free-form dim:  {stats(free_dims)}",
        f"- Blueprint dim:  {stats(bp_dims)}",
        "",
        "## Provenance density (DOI-shaped tokens per script)",
        "",
        f"- Free-form (in `embed_freeform.py`):     {stats(free_dois)}",
        f"- Blueprint (per run, sum across resources): {stats(bp_dois_per_run)}",
        "",
        "## Embedding quality (from outs/metrics/<dataset>/metrics.csv)",
        "",
    ]
    if quality_df.empty:
        lines.append("_metrics.csv not found or no rows for these models. "
                     "Run `scripts/embedding_metrics/metrics.py --dataset <X>` "
                     "after Finalize completes._")
    else:
        cols = ["model", "label_col", "ari", "nmi", "clisi",
                "isolated_label_f1", "ilisi", "kbet_acceptance"]
        present = [c for c in cols if c in quality_df.columns]
        agg = (
            quality_df.groupby(["model", "label_col"])[
                [c for c in present if c not in ("model", "label_col")]
            ].median().reset_index()
        )
        lines.append("Median across seeds:")
        lines.append("")
        lines.append("| " + " | ".join(present) + " |")
        lines.append("|" + "|".join(["---"] * len(present)) + "|")
        for _, r in agg.iterrows():
            cells = []
            for c in present:
                v = r[c]
                cells.append(f"{v:.3f}" if isinstance(v, float) and pd.notna(v) else str(v))
            lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Measure free-form agent baseline.")
    p.add_argument("--freeform-cohort", required=True,
                   help="Cohort dir name under work/, e.g. "
                        "freeform_gpt_5_4_mini_medium_human_pancreas")
    p.add_argument("--blueprint-cohort", required=True,
                   help="Matching blueprint cohort under work/ for vocabulary "
                        "and Jaccard reference, e.g. "
                        "gpt_5_4_mini_medium_human_pancreas")
    p.add_argument("--dataset", required=True,
                   help="mouse_brain | human_pancreas; selects metrics.csv path.")
    p.add_argument("--out-dir", type=Path,
                   default=None,
                   help="Override output directory (default outs/freeform_baseline/<cohort>)")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    free_dir = REPO_ROOT / "work" / args.freeform_cohort
    bp_dir = REPO_ROOT / "work" / args.blueprint_cohort
    if not free_dir.is_dir():
        raise FileNotFoundError(f"Free-form cohort not found: {free_dir}")
    if not bp_dir.is_dir():
        raise FileNotFoundError(f"Blueprint cohort not found: {bp_dir}")

    out_dir = args.out_dir or REPO_ROOT / "outs" / "freeform_baseline" / args.freeform_cohort
    out_dir.mkdir(parents=True, exist_ok=True)

    vocabulary = cohort_genes(bp_dir)
    print(f"[ok] vocabulary size from blueprint cohort: {len(vocabulary)}", flush=True)

    free_records: list[dict] = []
    for run in load_runs(free_dir, "*_freeform_run*"):
        free_records.append(per_run_record(run, vocabulary))
    print(f"[ok] free-form runs scanned: {len(free_records)}", flush=True)

    bp_genes = [blueprint_run_genes(r) for r in load_runs(bp_dir, "*_blueprint_run*")]
    bp_genes = [s for s in bp_genes if s]
    free_genes = [r["_genes"] for r in free_records if r["_genes"]]

    free_jaccard_self = mean_pairwise_jaccard(free_genes)
    bp_jaccard_self = mean_pairwise_jaccard(bp_genes)
    free_vs_bp_jaccard = cross_pair_mean_jaccard(free_genes, bp_genes)

    free_dims = [r["n_dim_columns"] for r in free_records if r["n_dim_columns"] > 0]
    bp_dims: list[int] = []
    for run in load_runs(bp_dir, "*_blueprint_run*"):
        sj = run / "summary.json"
        if sj.is_file():
            try:
                data = json.loads(sj.read_text())
                bp_dims.append(int(data.get("blueprint", {}).get("n_axes", 0)))
            except Exception:
                pass
    free_dois = [r["n_dois"] for r in free_records]
    bp_dois_per_run: list[int] = []
    for run in load_runs(bp_dir, "*_blueprint_run*"):
        sj = run / "summary.json"
        if sj.is_file():
            try:
                data = json.loads(sj.read_text())
                bp_dois_per_run.append(int(data.get("blueprint", {}).get("n_dois_unique", 0)))
            except Exception:
                pass

    metrics_csv = REPO_ROOT / "outs" / "metrics" / args.dataset / "metrics.csv"
    suffix = f"_{args.dataset}"
    free_model = args.freeform_cohort
    if free_model.endswith(suffix):
        free_model = free_model[: -len(suffix)]
    bp_model = args.blueprint_cohort
    if bp_model.endswith(suffix):
        bp_model = bp_model[: -len(suffix)]
    bp_model = "blueprint_" + bp_model
    quality_df = quality_table(metrics_csv, free_model=free_model, blueprint_model=bp_model)

    csv_path = out_dir / "measurement.csv"
    with csv_path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow([
            "run", "stage_a_rc", "script_present", "script_size",
            "n_dois", "n_genes_extracted",
            "embed_status", "parquet_schema_ok", "n_dim_columns",
            "parquet_path",
        ])
        for r in free_records:
            writer.writerow([
                r["run"], r["stage_a_rc"], r["script_present"], r["script_size"],
                r["n_dois"], r["n_genes_extracted"],
                r["embed_status"], r["parquet_schema_ok"], r["n_dim_columns"],
                r["parquet_path"],
            ])

    summary = render_summary(
        free_records, free_jaccard_self, free_vs_bp_jaccard, bp_jaccard_self,
        free_dims, bp_dims, free_dois, bp_dois_per_run, quality_df, args,
    )
    (out_dir / "measurement_summary.md").write_text(summary)

    print(f"[ok] wrote {csv_path}", flush=True)
    print(f"[ok] wrote {out_dir / 'measurement_summary.md'}", flush=True)


if __name__ == "__main__":
    main()
