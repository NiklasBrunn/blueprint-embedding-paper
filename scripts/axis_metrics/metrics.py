"""
Per-axis biological-quality metrics for blueprint embeddings.

For each Embed blueprint parquet, this script computes two views:

  1. AUROC of each axis-direction's score vs the binary label
     `adata.obs[label_col] == resolved_target`, where `resolved_target`
     is the cell-type label the axis claims to mark (resolved via the
     bidirectional token-subset auto-matcher + alias-config sidecar).
     Plain axes → one row; gated axes → one row (positive only;
     the gate is scope, not target); contrast axes → one row for
     the positive module + one row for each negative module (scored
     against `-axis_value`).
  2. Per-cohort label-coverage summary: for each cell-type label
     value, which axis (if any) scores it best (`best_auroc`), how
     many axes target it (`n_axes_targeting`), and the mean AUROC
     across those targeting axes.

Multi-axis-on-same-label is not penalized in `axis_metrics.csv`;
`label_coverage.csv` is the label-side view that doesn't double-count.

Inputs:
    outs/embeddings/<dataset>/blueprint_*.parquet
    outs/embeddings/<dataset>/blueprint_*_summary.json
    data/preprocessed/<dataset>_pp.h5ad
    configs/axis_label_aliases.jsonc

Outputs (under outs/axis_metrics/<dataset>/):
    axis_metrics.csv         one row per (parquet, axis, direction, label_col)
    label_coverage.csv       one row per (parquet, label_col, label_value)
    <stem>_axis_metrics.md   per-parquet markdown companion
    unmapped_axes.log        axes/directions that resolved to no label

Usage:
    python scripts/axis_metrics/metrics.py --dataset human_pancreas --overwrite
    python scripts/axis_metrics/metrics.py --dataset mouse_brain    --skip-existing
    python scripts/axis_metrics/metrics.py --dataset human_pancreas \\
        --input outs/embeddings/human_pancreas/blueprint_X_seed1.parquet
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[2]

DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")

# Label columns per dataset. Kept in sync with scripts/embedding_metrics/metrics.py.
LABEL_COLS = {
    "mouse_brain":         ("primary_type", "subclass", "broad_type"),
    "human_pancreas":      ("celltype",),
    "mouse_gastrulation":  ("celltype", "stage"),
    "mouse_cortex":        ("celltype_harmonized", "class_harmonized"),
    "human_immune":        ("final_annotation",),
}

DEFAULT_ALIAS_CONFIG = REPO_ROOT / "configs" / "axis_label_aliases.jsonc"

AXIS_METRICS_COLS = [
    "model", "dataset", "seed",
    "axis_name", "axis_kind", "direction", "module_name",
    "label_col", "label_value", "match_source",
    "auroc", "n_pos", "n_neg",
]

LABEL_COVERAGE_COLS = [
    "model", "dataset", "seed", "label_col", "label_value",
    "n_axes_targeting", "best_auroc", "best_axis_name",
    "best_module_name", "best_direction", "mean_auroc", "n_cells",
]


# --------------------------------------------------------------------
# Pure-Python helpers (no heavy deps)
# --------------------------------------------------------------------


def parse_stem(stem: str, dataset: str) -> tuple[str, int | None]:
    """Recover (model, seed) from '<model>_<dataset>[_seed<N>]'."""
    if f"_{dataset}" not in stem:
        raise ValueError(f"Stem {stem!r} does not contain dataset {dataset!r}.")
    head, _, tail = stem.partition(f"_{dataset}")
    if not head:
        raise ValueError(f"Empty model name parsed from stem {stem!r}.")
    if not tail:
        return head, None
    if tail.startswith("_seed"):
        return head, int(tail[len("_seed"):])
    return head + tail, None


def discover_parquets(emb_dir: Path, override: Path | None) -> list[Path]:
    """Return blueprint_*.parquet files (skipping UMAP coord siblings)."""
    if override is not None:
        if not override.exists():
            raise FileNotFoundError(f"--input not found: {override}")
        return [override.resolve()]
    if not emb_dir.exists():
        raise FileNotFoundError(f"Embeddings directory not found: {emb_dir}")
    parquets = [
        p for p in sorted(emb_dir.glob("blueprint_*.parquet"))
        if not p.name.endswith("_umap_coords.parquet")
    ]
    if not parquets:
        raise FileNotFoundError(f"No blueprint parquets in {emb_dir}")
    return parquets


def load_alias_config(path: Path) -> dict:
    """Strip // line comments and parse the jsonc alias file. Returns {} on missing."""
    if not path.exists():
        return {}
    text = path.read_text()
    # Naive regex is safe: the alias file never embeds // inside a string.
    stripped = re.sub(r"//.*", "", text)
    return json.loads(stripped)


# Auto-matching: bidirectional token-subset rule.
_SEP_CHARS = " -/+()"


def normalize(s: str) -> set[str]:
    """Lowercase; treat separators as `_`; split on `_`; drop empties."""
    s = s.lower()
    for c in _SEP_CHARS:
        s = s.replace(c, "_")
    return {t for t in s.split("_") if t}


def auto_match(name: str, label_values: list[str]) -> str | None:
    """Match `name` against any of `label_values` via token-subset.

    Direction priority: (a) label tokens ⊆ name tokens (label fits inside
    module/axis name); (b) name tokens ⊆ label tokens (module fits inside
    a multi-word label). Within each direction, ties are broken by
    most-specific (largest token set for direction (a); smallest for (b));
    unresolved ties → None.
    """
    name_toks = normalize(name)
    if not name_toks:
        return None
    # Direction (a).
    a_cands: list[tuple[str, int]] = []
    for lab in label_values:
        lab_toks = normalize(lab)
        if lab_toks and lab_toks <= name_toks:
            a_cands.append((lab, len(lab_toks)))
    if len(a_cands) == 1:
        return a_cands[0][0]
    if len(a_cands) >= 2:
        a_cands.sort(key=lambda t: -t[1])
        if a_cands[0][1] > a_cands[1][1]:
            return a_cands[0][0]
        return None
    # Direction (b).
    b_cands: list[tuple[str, int]] = []
    for lab in label_values:
        lab_toks = normalize(lab)
        if lab_toks and name_toks <= lab_toks:
            b_cands.append((lab, len(lab_toks)))
    if len(b_cands) == 1:
        return b_cands[0][0]
    if len(b_cands) >= 2:
        b_cands.sort(key=lambda t: t[1])
        if b_cands[0][1] < b_cands[1][1]:
            return b_cands[0][0]
    return None


def resolve_label(
    module_name: str,
    axis_name: str,
    label_col: str,
    label_values: list[str],
    aliases: dict,
    dataset: str,
) -> tuple[str | None, str]:
    """Return (label_value, match_source). match_source ∈
    {module_auto, module_alias, axis_auto, axis_alias, none}."""
    dataset_aliases = (aliases.get(dataset, {}) or {}).get(label_col, {}) or {}
    # 1. module name auto.
    m = auto_match(module_name, label_values)
    if m is not None:
        return m, "module_auto"
    # 2. module name alias.
    if module_name.lower() in dataset_aliases:
        v = dataset_aliases[module_name.lower()]
        if v in label_values:
            return v, "module_alias"
    # 3. axis name auto (fallback when module name doesn't carry the identity).
    if axis_name and axis_name != module_name:
        a = auto_match(axis_name, label_values)
        if a is not None:
            return a, "axis_auto"
        if axis_name.lower() in dataset_aliases:
            v = dataset_aliases[axis_name.lower()]
            if v in label_values:
                return v, "axis_alias"
    return None, "none"


# --------------------------------------------------------------------
# I/O
# --------------------------------------------------------------------


def load_summary(summary_path: Path) -> dict:
    if not summary_path.exists():
        raise FileNotFoundError(f"summary.json missing: {summary_path}")
    return json.loads(summary_path.read_text())


def load_embedding(parquet_path: Path) -> pd.DataFrame:
    """Read cell_id + dim_* columns; return a DataFrame indexed by cell_id."""
    df = pd.read_parquet(parquet_path)
    if "cell_id" not in df.columns:
        raise ValueError(f"{parquet_path} lacks cell_id column.")
    df = df.set_index("cell_id")
    dim_cols = [c for c in df.columns if c.startswith("dim_")]
    if not dim_cols:
        raise ValueError(f"{parquet_path} has no dim_* columns.")
    return df[dim_cols]


def load_anndata(dataset: str):
    """Load the preprocessed AnnData."""
    import anndata as ad
    pp_path = REPO_ROOT / "data" / "preprocessed" / f"{dataset}_pp.h5ad"
    if not pp_path.exists():
        raise FileNotFoundError(f"preprocessed h5ad missing: {pp_path}")
    return ad.read_h5ad(pp_path)


def labels_for(adata, label_col: str) -> pd.Series:
    """Return a Series of `label_col` indexed by adata.obs_names (str)."""
    if label_col not in adata.obs.columns:
        return None
    return adata.obs[label_col].astype(str).copy()


# --------------------------------------------------------------------
# Computation
# --------------------------------------------------------------------


def compute_auroc(score: np.ndarray, y_true: np.ndarray) -> float:
    """sklearn AUROC."""
    from sklearn.metrics import roc_auc_score
    if y_true.sum() == 0 or y_true.sum() == len(y_true):
        return float("nan")
    return float(roc_auc_score(y_true, score))


# --------------------------------------------------------------------
# Per-parquet driver
# --------------------------------------------------------------------


def process_parquet(
    parquet: Path,
    dataset: str,
    label_cols: tuple[str, ...],
    adata,
    aliases: dict,
    unmapped_log: list[str],
) -> tuple[list[dict], list[dict]]:
    """Compute axis_metrics and label_coverage rows for one parquet.

    Returns (axis_rows, coverage_rows).
    """
    stem = parquet.stem
    model, seed = parse_stem(stem, dataset)
    summary_path = parquet.with_name(f"{stem}_summary.json")
    summary = load_summary(summary_path)

    emb = load_embedding(parquet)
    # Align by cell_id.
    obs_names = adata.obs_names.astype(str)
    common = obs_names.intersection(emb.index)
    if len(common) == 0:
        raise ValueError(f"No cell_id overlap between {parquet} and adata.")
    emb = emb.loc[common]
    cell_idx = obs_names.get_indexer(common)
    # Subviews of adata that match cell-by-cell with emb.
    adata_aligned = adata[cell_idx, :]

    bp = summary["blueprint"]
    modules = {m["name"]: m for m in bp["modules"]}
    axes = bp["axes"]
    dim_cols = [f"dim_{i}" for i in range(len(axes))]
    if dim_cols and dim_cols[-1] not in emb.columns:
        # Older parquets may have fewer dims than axes when build_blueprint
        # changed since the parquet was written. Truncate to what's on disk.
        present_dims = [c for c in emb.columns if c.startswith("dim_")]
        n_present = len(present_dims)
        if n_present < len(axes):
            axes = axes[:n_present]
            dim_cols = present_dims

    axis_rows: list[dict] = []

    # Precompute per-label-col label arrays + uniques.
    label_data: dict[str, tuple[pd.Series, list[str]]] = {}
    for lc in label_cols:
        s = labels_for(adata_aligned, lc)
        if s is None:
            continue
        label_data[lc] = (s.reset_index(drop=True), sorted(s.unique().tolist()))

    # ---- per-axis AUROC ----
    for dim_idx, axis in enumerate(axes):
        axis_name = axis["name"]
        axis_kind = axis.get("kind", "plain")
        col = f"dim_{dim_idx}"
        if col not in emb.columns:
            continue
        score = emb[col].to_numpy(dtype=float)

        # List of (module_name, direction, sign) tuples.
        directions: list[tuple[str, str, int]] = []
        if axis.get("positive"):
            directions.append((axis["positive"], "positive", +1))
        for neg in axis.get("negatives") or []:
            directions.append((neg, "negative", -1))
        if axis_kind == "gated":
            # gated: the gate is scope, not a target; only score positive.
            directions = [d for d in directions if d[1] == "positive"]

        for module_name, direction, sign in directions:
            signed_score = score if sign > 0 else -score
            for lc, (lab_series, lab_uniques) in label_data.items():
                target, src = resolve_label(
                    module_name, axis_name, lc, lab_uniques, aliases, dataset
                )
                if target is None:
                    unmapped_log.append(
                        f"{stem}\taxis={axis_name}\tdirection={direction}"
                        f"\tmodule={module_name}\tlabel_col={lc}"
                    )
                    continue
                y_true = (lab_series.values == target)
                au = compute_auroc(signed_score, y_true.astype(bool))
                axis_rows.append({
                    "model": model, "dataset": dataset, "seed": seed,
                    "axis_name": axis_name, "axis_kind": axis_kind,
                    "direction": direction, "module_name": module_name,
                    "label_col": lc, "label_value": target,
                    "match_source": src,
                    "auroc": au,
                    "n_pos": int(y_true.sum()), "n_neg": int((~y_true).sum()),
                })

    # ---- label coverage rollup from axis_rows ----
    coverage_rows: list[dict] = []
    if axis_rows:
        df = pd.DataFrame(axis_rows)
        for (lc, lv), g in df.groupby(["label_col", "label_value"]):
            best = g.loc[g["auroc"].idxmax()] if g["auroc"].notna().any() else None
            coverage_rows.append({
                "model": model, "dataset": dataset, "seed": seed,
                "label_col": lc, "label_value": lv,
                "n_axes_targeting": int(len(g)),
                "best_auroc": float(best["auroc"]) if best is not None else float("nan"),
                "best_axis_name": str(best["axis_name"]) if best is not None else "",
                "best_module_name": str(best["module_name"]) if best is not None else "",
                "best_direction": str(best["direction"]) if best is not None else "",
                "mean_auroc": float(g["auroc"].mean()) if g["auroc"].notna().any() else float("nan"),
                "n_cells": int(g["n_pos"].iloc[0]),
            })

    return axis_rows, coverage_rows


# --------------------------------------------------------------------
# CSV merge / write
# --------------------------------------------------------------------


def merge_csv(csv_path: Path, new_rows: list[dict], key_cols: tuple[str, ...]) -> pd.DataFrame:
    """Replace rows whose key matches a new row; append the rest."""
    new_df = pd.DataFrame(new_rows)
    if not csv_path.exists() or new_df.empty:
        return new_df
    existing = pd.read_csv(csv_path)
    if existing.empty:
        return new_df
    mask = pd.Series(False, index=existing.index)
    for _, row in new_df.iterrows():
        m = pd.Series(True, index=existing.index)
        for k in key_cols:
            val = row[k]
            if pd.isna(val):
                m &= existing[k].isna()
            else:
                m &= (existing[k] == val)
        mask |= m
    return pd.concat([existing[~mask], new_df], ignore_index=True)


def write_csv(path: Path, rows: list[dict], cols: list[str], full_sweep: bool,
              overwrite: bool, key_cols: tuple[str, ...]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if full_sweep:
        if path.exists() and not overwrite:
            raise FileExistsError(f"{path} exists. Pass --overwrite to replace it.")
        df = pd.DataFrame(rows, columns=cols)
    else:
        df = merge_csv(path, rows, key_cols)
        df = df.reindex(columns=cols)
    df = df.sort_values([c for c in ["model", "seed", "label_col", "axis_name", "module_name", "label_value"]
                         if c in df.columns], na_position="last").reset_index(drop=True)
    df.to_csv(path, index=False)


def _md_table(df: pd.DataFrame, float_cols: tuple[str, ...] = ()) -> str:
    """Render a small DataFrame as a pipe-formatted markdown table.

    Self-contained (no `tabulate` dependency, which pandas' to_markdown needs).
    """
    cols = list(df.columns)
    def fmt(c, v):
        if pd.isna(v):
            return ""
        if c in float_cols:
            return f"{float(v):.3f}"
        return str(v)
    header = "| " + " | ".join(cols) + " |"
    sep = "|" + "|".join("---" for _ in cols) + "|"
    body = "\n".join(
        "| " + " | ".join(fmt(c, row[c]) for c in cols) + " |"
        for _, row in df.iterrows()
    )
    return "\n".join([header, sep, body])


def write_companion_md(stem: str, dataset: str, out_dir: Path, axis_rows: list[dict],
                       coverage_rows: list[dict]) -> None:
    md_path = out_dir / f"{stem}_axis_metrics.md"
    lines: list[str] = []
    lines.append(f"# Axis metrics — {stem}")
    lines.append("")
    lines.append(f"- Timestamp: `{datetime.now().isoformat(timespec='seconds')}`")
    lines.append(f"- Dataset:   `{dataset}`")
    lines.append(f"- n_axis_rows: {len(axis_rows)}; "
                 f"n_coverage_rows: {len(coverage_rows)}")
    lines.append("")
    if axis_rows:
        lines.append("## Per-axis AUROC (top 20 by AUROC)")
        df = pd.DataFrame(axis_rows).sort_values("auroc", ascending=False).head(20)
        cols = ["axis_name", "direction", "module_name", "label_col",
                "label_value", "auroc", "match_source"]
        lines.append(_md_table(df[cols], float_cols=("auroc",)))
        lines.append("")
    if coverage_rows:
        lines.append("## Per-label best AUROC")
        df = pd.DataFrame(coverage_rows).sort_values(
            ["label_col", "best_auroc"], ascending=[True, False])
        cols = ["label_col", "label_value", "n_axes_targeting",
                "best_auroc", "best_axis_name", "best_direction", "mean_auroc"]
        lines.append(_md_table(df[cols], float_cols=("best_auroc", "mean_auroc")))
        lines.append("")
    md_path.write_text("\n".join(lines) + "\n")


# --------------------------------------------------------------------
# Main
# --------------------------------------------------------------------


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dataset", required=True, choices=DATASET_CHOICES)
    p.add_argument("--input", type=Path, default=None,
                   help="Single parquet (overrides discovery).")
    p.add_argument("--alias-config", type=Path, default=DEFAULT_ALIAS_CONFIG)
    p.add_argument("--label-cols", default=None,
                   help="Comma-separated override (defaults to LABEL_COLS[dataset]).")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--skip-existing", action="store_true",
                   help="Skip parquets whose <stem>_axis_metrics.md exists.")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    label_cols = (tuple(args.label_cols.split(",")) if args.label_cols
                  else LABEL_COLS[args.dataset])
    out_dir = REPO_ROOT / "outs" / "axis_metrics" / args.dataset
    out_dir.mkdir(parents=True, exist_ok=True)
    emb_dir = REPO_ROOT / "outs" / "embeddings" / args.dataset

    aliases = load_alias_config(args.alias_config)
    print(f"[ok] alias config: {args.alias_config}  "
          f"(datasets: {sorted(aliases.keys())})")

    parquets = discover_parquets(emb_dir, args.input)
    if args.skip_existing:
        parquets = [p for p in parquets
                    if not (out_dir / f"{p.stem}_axis_metrics.md").exists()]
    if not parquets:
        print("[ok] no parquets to process (all have _axis_metrics.md).")
        return
    print(f"[ok] dataset={args.dataset} files={len(parquets)} "
          f"labels={list(label_cols)}")

    adata = load_anndata(args.dataset)
    print(f"[ok] adata: {adata.n_obs} cells × {adata.n_vars} genes")

    unmapped: list[str] = []
    axis_rows_all: list[dict] = []
    coverage_rows_all: list[dict] = []

    for parquet in parquets:
        try:
            a_rows, c_rows = process_parquet(
                parquet, args.dataset, label_cols, adata, aliases, unmapped
            )
        except Exception as exc:
            print(f"[fail] {parquet.name}: {exc}", file=sys.stderr)
            continue
        axis_rows_all.extend(a_rows)
        coverage_rows_all.extend(c_rows)
        write_companion_md(parquet.stem, args.dataset, out_dir,
                           a_rows, c_rows)
        print(f"[ok] {parquet.name}: axis={len(a_rows)} "
              f"coverage={len(c_rows)}")

    full_sweep = (args.input is None) and (not args.skip_existing)
    write_csv(out_dir / "axis_metrics.csv", axis_rows_all, AXIS_METRICS_COLS,
              full_sweep, args.overwrite,
              key_cols=("model", "seed", "axis_name", "direction", "module_name",
                        "label_col"))
    write_csv(out_dir / "label_coverage.csv", coverage_rows_all, LABEL_COVERAGE_COLS,
              full_sweep, args.overwrite,
              key_cols=("model", "seed", "label_col", "label_value"))
    log_path = out_dir / "unmapped_axes.log"
    if unmapped:
        log_path.write_text("\n".join(unmapped) + "\n")
        print(f"[ok] {len(unmapped)} unmapped axis-direction-label_col rows → {log_path.name}")
    elif log_path.exists():
        log_path.unlink()

    print(f"[ok] wrote {out_dir / 'axis_metrics.csv'} "
          f"({len(axis_rows_all)} rows)")
    print(f"[ok] wrote {out_dir / 'label_coverage.csv'} "
          f"({len(coverage_rows_all)} rows)")


if __name__ == "__main__":
    main()
