"""Bake webapp/eta_table.json from outs/resource_costs/*/costs.csv.

Key: "<model_key>|<effort>|<refine 0/1>" -> {stage_a_s, stage_b_s, n}.
model_key normalizes model names (gpt-5.5 -> gpt_5_5). Pools over datasets and
noref/ref-resource variants. Run from the repo root with the webapp venv:

    webapp/.venv/bin/python webapp/build_eta_table.py
"""
from __future__ import annotations
import json
from pathlib import Path
import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
COSTS = sorted((REPO_ROOT / "outs" / "resource_costs").glob("*/costs.csv"))
OUT = Path(__file__).resolve().parent / "eta_table.json"

EFFORTS = {"low", "medium", "high", "xhigh"}
FLAGS = {"noref", "unrefined", "refine"}


def parse_cohort(cohort: str):
    """Parse a blueprint cohort name into (model_key, effort, refine)|None.

    Cohort format is "<model>_<flags>_<effort>_<dataset>", e.g.
    "gpt_5_5_noref_refine_high_human_pancreas". The effort token (low/medium/
    high/xhigh) splits the model+flag prefix from the trailing dataset suffix;
    no dataset name contains an effort token, so the first match is the effort.
    """
    toks = cohort.split("_")
    eff_idx = next((i for i, t in enumerate(toks) if t in EFFORTS), None)
    if eff_idx is None:
        return None
    effort = toks[eff_idx]
    prefix = toks[:eff_idx]                  # model tokens + flags
    refine = "refine" in prefix              # "unrefined" != "refine"
    model_key = "_".join(t for t in prefix if t not in FLAGS)
    if not model_key:
        return None
    return model_key, effort, int(refine)


def main():
    rows = []
    for csv in COSTS:
        df = pd.read_csv(csv)
        if "cohort" not in df or "agent_stage_a_s" not in df:
            continue
        if "variant" in df:
            df = df[df["variant"] == "blueprint"]
        for _, r in df.iterrows():
            parsed = parse_cohort(str(r["cohort"]))
            if not parsed:
                continue
            model_key, effort, refine = parsed
            rows.append({
                "model_key": model_key, "effort": effort, "refine": refine,
                "stage_a_s": r.get("agent_stage_a_s"),
                "stage_b_s": r.get("agent_stage_b_refine_s"),
            })
    table = {}
    if rows:
        agg = pd.DataFrame(rows)
        for (mk, eff, rf), g in agg.groupby(["model_key", "effort", "refine"]):
            key = f"{mk}|{eff}|{rf}"
            table[key] = {
                "stage_a_s": float(g["stage_a_s"].median(skipna=True)),
                "stage_b_s": (float(g["stage_b_s"].median(skipna=True))
                              if g["stage_b_s"].notna().any() else None),
                "n": int(g["stage_a_s"].notna().sum()),
            }
    OUT.write_text(json.dumps(table, indent=2, sort_keys=True))
    print(f"wrote {OUT} with {len(table)} keys")


if __name__ == "__main__":
    main()
