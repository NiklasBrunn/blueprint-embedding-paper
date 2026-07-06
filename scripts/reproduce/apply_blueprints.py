#!/usr/bin/env python3
"""Apply the shipped blueprints to preprocessed data, rebuilding outs/embeddings/.

Consumes the authored blueprints committed under
``outs/blueprints/<dataset>/<cohort>/`` (so no agent re-authoring is needed) and
writes the canonical embedding artifacts that the metric and plot scripts read
from ``outs/embeddings/<dataset>/``:

    blueprint_<cohort>_<dataset>_seed<N>.parquet          cell x axis table
    blueprint_<cohort>_<dataset>_seed<N>_coverage.json    written by the embedder
    blueprint_<cohort>_<dataset>_seed<N>_summary.json     copied from the blueprint

The output stem mirrors the Embed step of ``analysis/agent_edit_blueprint.sh``:
the cohort directory name is the model slug, so
``outs/blueprints/mouse_cortex/gpt_5_5_unrefined_high/embedder_blueprint_seed2.py``
becomes ``blueprint_gpt_5_5_unrefined_high_mouse_cortex_seed2.parquet`` (model
``blueprint_gpt_5_5_unrefined_high`` in outs/metrics/<dataset>/metrics.csv).

Two modes:

  (default)        Run each blueprint's own CLI against
                   data/preprocessed/<dataset>_pp.h5ad to build the parquet +
                   coverage sidecar, then stage the shipped summary alongside it.
  --metadata-only  Copy only the shipped _summary.json / _coverage.json into the
                   outs/embeddings layout (no data, no parquet). Enough for the
                   metadata-only figures (blueprint_card, three_way_summary) that
                   the light reproduce_figures.sh renders without embeddings.

With --freeform the free-form baseline scripts under
``outs/freeform_scripts/<dataset>/<cohort>/`` are applied the same way (their
cohort name already carries the dataset, so the stem is ``<cohort>_seed<N>``).

Usage:
    python scripts/reproduce/apply_blueprints.py                    # all datasets, full apply
    python scripts/reproduce/apply_blueprints.py --dataset human_pancreas
    python scripts/reproduce/apply_blueprints.py --metadata-only    # stage sidecars only
    python scripts/reproduce/apply_blueprints.py --freeform         # also apply free-form scripts
"""
from __future__ import annotations

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
BLUEPRINTS_ROOT = REPO_ROOT / "outs" / "blueprints"
FREEFORM_ROOT = REPO_ROOT / "outs" / "freeform_scripts"
EMBEDDINGS_ROOT = REPO_ROOT / "outs" / "embeddings"

_SEED_RE = re.compile(r"_seed(\d+)\.py$")


def _pp_input(dataset: str) -> Path:
    return REPO_ROOT / "data" / "preprocessed" / f"{dataset}_pp.h5ad"


def _seed_of(script: Path) -> int | None:
    m = _SEED_RE.search(script.name)
    return int(m.group(1)) if m else None


def _discover(root: Path, datasets: list[str] | None) -> list[str]:
    """Dataset subdirs present under a shipped-artifact root, optionally filtered."""
    if not root.is_dir():
        return []
    found = sorted(p.name for p in root.iterdir() if p.is_dir())
    if datasets:
        found = [d for d in found if d in datasets]
    return found


def _stage_sidecars(src_stem: Path, out_stem: Path, overwrite: bool) -> list[str]:
    """Copy the <src_stem>_{summary,coverage}.json sidecars to <out_stem>_*.json."""
    staged = []
    for kind in ("summary", "coverage"):
        src = src_stem.with_name(f"{src_stem.name}_{kind}.json")
        if not src.exists():
            continue
        dst = out_stem.with_name(f"{out_stem.name}_{kind}.json")
        if dst.exists() and not overwrite:
            continue
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(src, dst)
        staged.append(kind)
    return staged


def _run_embedder(script: Path, dataset: str, parquet: Path, overwrite: bool) -> bool:
    """Run a blueprint / free-form script's CLI to build one parquet. Returns ok."""
    if parquet.exists() and not overwrite:
        print(f"  [skip] {parquet.relative_to(REPO_ROOT)} (exists)")
        return True
    inp = _pp_input(dataset)
    if not inp.exists():
        print(f"  [error] missing input {inp.relative_to(REPO_ROOT)} — run "
              f"`bash analysis/run_data_pipeline.sh {dataset}` first.", file=sys.stderr)
        return False
    parquet.parent.mkdir(parents=True, exist_ok=True)
    cmd = [sys.executable, str(script), "--input", str(inp), "--output", str(parquet)]
    if "embed_freeform" in script.name:                     # free-form contract wants --seed
        seed = _seed_of(script)
        cmd += ["--seed", str(seed if seed is not None else 0)]
    print(f"  [embed] {script.relative_to(REPO_ROOT)} -> {parquet.relative_to(REPO_ROOT)}")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(f"  [error] {script.name} failed (rc={proc.returncode})\n{proc.stderr}",
              file=sys.stderr)
        return False
    return True


def apply_blueprint_cohort(cohort_dir: Path, dataset: str, *, metadata_only: bool,
                           overwrite: bool) -> tuple[int, int]:
    """Apply (or stage) every seed of one blueprint cohort. Returns (ok, total)."""
    cohort = cohort_dir.name
    ok = total = 0
    for script in sorted(cohort_dir.glob("embedder_blueprint_seed*.py")):
        seed = _seed_of(script)
        if seed is None:
            continue
        total += 1
        stem = f"blueprint_{cohort}_{dataset}_seed{seed}"
        out_stem = EMBEDDINGS_ROOT / dataset / stem
        src_stem = script.with_name(script.stem)           # .../embedder_blueprint_seed<N>
        if metadata_only:
            _stage_sidecars(src_stem, out_stem, overwrite)
            ok += 1
            continue
        if _run_embedder(script, dataset, out_stem.with_suffix(".parquet"), overwrite):
            # The embedder writes <stem>_coverage.json itself; stage the (structural,
            # data-independent) summary.json from the shipped blueprint next to it.
            _stage_sidecars(src_stem, out_stem, overwrite=True)
            ok += 1
    return ok, total


def apply_freeform_cohort(cohort_dir: Path, dataset: str, overwrite: bool) -> tuple[int, int]:
    """Apply every seed of one free-form cohort. Returns (ok, total)."""
    cohort = cohort_dir.name                                # already ends with <dataset>
    ok = total = 0
    for script in sorted(cohort_dir.glob("embed_freeform_seed*.py")):
        seed = _seed_of(script)
        if seed is None:
            continue
        total += 1
        parquet = EMBEDDINGS_ROOT / dataset / f"{cohort}_seed{seed}.parquet"
        if _run_embedder(script, dataset, parquet, overwrite):
            # Status mirror consumed by plot_three_way_summary's freeform_stats
            # (embed_status / embed_detail); the Finalize step of the freeform
            # orchestrator writes the same file during real runs.
            status = parquet.with_name(f"{cohort}_seed{seed}_status.json")
            status.write_text(json.dumps(
                {"embed_status": "OK",
                 "embed_detail": str(parquet.relative_to(REPO_ROOT))}, indent=2) + "\n")
            ok += 1
    return ok, total


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--dataset", action="append", dest="datasets",
                   help="Restrict to this dataset (repeatable). Default: all shipped.")
    p.add_argument("--metadata-only", action="store_true",
                   help="Only stage shipped summary/coverage sidecars (no data, no parquet).")
    p.add_argument("--freeform", action="store_true",
                   help="Also apply the free-form baseline scripts under outs/freeform_scripts/.")
    p.add_argument("--overwrite", action="store_true",
                   help="Rebuild/restage even when the target already exists.")
    args = p.parse_args()

    datasets = args.datasets
    total_ok = total_all = 0

    for dataset in _discover(BLUEPRINTS_ROOT, datasets):
        cohorts = sorted(p for p in (BLUEPRINTS_ROOT / dataset).iterdir() if p.is_dir())
        if not cohorts:
            continue
        mode = "stage" if args.metadata_only else "embed"
        print(f"[{dataset}] {mode} {len(cohorts)} blueprint cohort(s)")
        for cohort_dir in cohorts:
            ok, tot = apply_blueprint_cohort(cohort_dir, dataset,
                                             metadata_only=args.metadata_only,
                                             overwrite=args.overwrite)
            total_ok += ok
            total_all += tot

    if args.freeform and not args.metadata_only:
        for dataset in _discover(FREEFORM_ROOT, datasets):
            cohorts = sorted(p for p in (FREEFORM_ROOT / dataset).iterdir() if p.is_dir())
            if not cohorts:
                continue
            print(f"[{dataset}] embed {len(cohorts)} free-form cohort(s)")
            for cohort_dir in cohorts:
                ok, tot = apply_freeform_cohort(cohort_dir, dataset, args.overwrite)
                total_ok += ok
                total_all += tot

    verb = "staged" if args.metadata_only else "applied"
    print(f"\n[done] {verb} {total_ok}/{total_all} run(s).")
    if total_ok != total_all:
        sys.exit(1)


if __name__ == "__main__":
    main()
