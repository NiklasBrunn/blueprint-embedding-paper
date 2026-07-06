"""
Full scIB numbers table: per-(dataset, method) mean ± s.d. across seeds.

Reads the three benchmark metrics.csv (the verified Benchmarker output) and emits,
for the benchmark cohort (blueprint grounding variants + baselines; free-form and
refine cohorts excluded):

  * outs/metrics/r2_summary_table.csv   — long, ALL scIB sub-metrics + aggregates
                                          + participation ratio, mean & s.d.
  * outs/metrics/r2_metrics_table.tex   — a \\footnotesize LaTeX table* (selected
                                          columns).

Each cell is mean ± s.d. across the cohort's seeds (deterministic baselines show
a single value, no ±). Rows are grouped by dataset/resolution and sorted by the
scIB Total within each group.

Usage:
    python scripts/plots/build_metrics_table.py
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]


def _load_metrics_module():
    path = REPO_ROOT / "scripts" / "embedding_metrics" / "metrics.py"
    spec = importlib.util.spec_from_file_location("embedding_metrics", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


M = _load_metrics_module()

# (dataset, label_col, display name) at the reported resolutions.
GROUPS = [
    ("mouse_cortex",   "celltype_harmonized", "mouse\\_cortex (subclass)"),
    ("mouse_cortex",   "class_harmonized",    "mouse\\_cortex (class)"),
    ("human_pancreas", "celltype",            "human\\_pancreas"),
    ("human_immune",   "final_annotation",    "human\\_immune"),
]
# "mini" drops the gpt-5.4-mini engine-ablation cohorts (human_pancreas only)
# from the model-comparison table.
R2_EXCLUDE = {"freeform", "refine", "mini"}

ALL_METRICS = ["bio_score", "batch_score_scib", "composite_scib",
               "ari", "nmi", "asw_label", "clisi", "isolated_label_f1",
               "asw_batch", "ilisi", "kbet_acceptance", "graph_connectivity",
               "pcr", "participation_ratio"]
# LaTeX table layout: every scIB sub-metric, grouped Bio | Batch (canonical scib
# ordering), each group closed by its aggregate sub-score, then Total and PR. The
# group spans drive the \cmidrule headers in write_tex.
BIO_COLS = [("ari", "ARI"), ("nmi", "NMI"), ("asw_label", r"ASW$_\mathrm{L}$"),
            ("clisi", "cLISI"), ("isolated_label_f1", "iso"),
            ("bio_score", r"\textbf{Bio}")]
BATCH_COLS = [("asw_batch", r"ASW$_\mathrm{B}$"), ("ilisi", "iLISI"),
              ("kbet_acceptance", "kBET"), ("graph_connectivity", "GC"),
              ("pcr", "PCR"), ("batch_score_scib", r"\textbf{Batch}")]
TAIL_COLS = [("composite_scib", r"\textbf{Total}"), ("participation_ratio", "PR")]
TEX_METRICS = BIO_COLS + BATCH_COLS + TAIL_COLS


def is_r2(model: str) -> bool:
    return not (R2_EXCLUDE & set(model.split("_")))


def aggregate() -> pd.DataFrame:
    rows = []
    for ds, lab, disp in GROUPS:
        df = pd.read_csv(REPO_ROOT / "outs" / "metrics" / ds / "metrics.csv")
        df = df[(df["label_col"] == lab) & df["model"].map(is_r2)]
        for model, g in df.groupby("model"):
            rec = {"dataset": ds, "label_col": lab, "display": disp,
                   "model": model, "pretty": M.prettify_model(model),
                   "n_seeds": len(g)}
            for m in ALL_METRICS:
                vals = g[m].dropna()
                rec[f"{m}_mean"] = vals.mean() if len(vals) else np.nan
                rec[f"{m}_sd"] = vals.std(ddof=1) if len(vals) > 1 else np.nan
            rows.append(rec)
    return pd.DataFrame(rows)


def _cell(mean: float, sd: float, prec: int = 2, bold: bool = False) -> str:
    if np.isnan(mean):
        return "--"
    m = f"\\textbf{{{mean:.{prec}f}}}" if bold else f"{mean:.{prec}f}"
    if np.isnan(sd):
        return m
    return f"{m}\\,$\\pm$\\,{sd:.{prec}f}"


def write_csv(agg: pd.DataFrame) -> Path:
    cols = ["dataset", "label_col", "model", "pretty", "n_seeds"]
    cols += [f"{m}_{s}" for m in ALL_METRICS for s in ("mean", "sd")]
    out = REPO_ROOT / "outs" / "metrics" / "r2_summary_table.csv"
    agg[cols].to_csv(out, index=False)
    return out


def write_tex(agg: pd.DataFrame) -> Path:
    n_metric = len(TEX_METRICS)
    n_bio, n_batch = len(BIO_COLS), len(BATCH_COLS)
    col_spec = "l" + "r" * n_metric           # Method + every metric (no Dataset col)
    n_cols = n_metric + 1
    # Grouped header: spanned "Bio conservation" / "Batch correction" rules over
    # their sub-metric columns; Total and PR stand alone.
    grp = (r" & \multicolumn{%d}{c}{Bio conservation} & "
           r"\multicolumn{%d}{c}{Batch correction} & & \\" % (n_bio, n_batch))
    cmid = (r"\cmidrule(lr){2-%d}\cmidrule(lr){%d-%d}"
            % (1 + n_bio, 2 + n_bio, 1 + n_bio + n_batch))
    head = "Method & " + " & ".join(h for _, h in TEX_METRICS) + r" \\"
    lines = [
        r"\begin{sidewaystable}[p]", r"\centering", r"\scriptsize",
        r"\setlength{\tabcolsep}{3.2pt}",
        r"\caption{\textbf{Blueprints reach competitive quality and are "
        r"batch-invariant by construction.} Full scIB benchmark numbers (mean $\pm$ "
        r"s.d.\ across seeds). Every scIB sub-metric, grouped into the bio-conservation "
        r"block (ARI, NMI, ASW$_\mathrm{L}$ $=$ ASW-label, cLISI, iso $=$ "
        r"isolated-label $F_1$) and the batch-correction block (ASW$_\mathrm{B}$ $=$ "
        r"ASW-batch, iLISI, kBET, GC $=$ graph connectivity, PCR), each closed by its "
        r"aggregate sub-score (\textbf{Bio}, \textbf{Batch}). \textbf{Total} "
        r"$=0.6\,$Bio$+0.4\,$Batch (the scIB package Total) and PR $=$ participation "
        r"ratio (effective dimensionality). Rows are the two blueprint grounding "
        r"variants (\texttt{BP GPT-5.5 U} $=$ literature-grounded, "
        r"\texttt{BP GPT-5.5 N U} $=$ resource-unconstrained, where U $=$ unrefined "
        r"and N $=$ \texttt{noref}), the program- and marker-based comparators "
        r"(UCell, expiMap), and the unsupervised and foundation baselines, grouped by "
        r"benchmark dataset and cell-type resolution and sorted by Total. The best "
        r"mean per column within each dataset block is in bold (PR excluded, as it is "
        r"not a quality score). Deterministic methods (PCA, scVI-census, UCE, UCell) "
        r"are scored once and carry no s.d.}",
        r"\label{supptab:r2numbers}",
        r"\begin{tabular}{" + col_spec + "}", r"\toprule", grp, cmid, head,
    ]
    # Participation ratio is effective dimensionality, not a quality score, so it
    # has no "winner" to bold; every other column is higher-is-better.
    bold_metrics = [m for m, _ in TEX_METRICS if m != "participation_ratio"]
    prec = {"participation_ratio": 1}
    for gi, (ds, lab, disp) in enumerate(GROUPS):
        sub = agg[(agg.dataset == ds) & (agg.label_col == lab)].copy()
        sub = sub.sort_values("composite_scib_mean", ascending=False)
        # Per metric column, the model with the highest mean wins (bold).
        winners = {m: sub[f"{m}_mean"].idxmax() for m in bold_metrics
                   if sub[f"{m}_mean"].notna().any()}
        lines.append(r"\midrule")
        lines.append(r"\multicolumn{%d}{l}{\textbf{%s}}\\[1pt]" % (n_cols, disp))
        for ridx, r in sub.iterrows():
            cells = [_cell(r[f"{m}_mean"], r[f"{m}_sd"], prec.get(m, 2),
                           bold=(winners.get(m) == ridx))
                     for m, _ in TEX_METRICS]
            lines.append(f"{r['pretty']} & " + " & ".join(cells) + r" \\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{sidewaystable}"]
    out = REPO_ROOT / "outs" / "metrics" / "r2_metrics_table.tex"
    out.write_text("\n".join(lines) + "\n")
    return out


def main() -> None:
    agg = aggregate()
    csv = write_csv(agg)
    tex = write_tex(agg)
    print(f"[ok] wrote {csv}")
    print(f"[ok] wrote {tex}  ({len(agg)} rows)")


if __name__ == "__main__":
    main()
