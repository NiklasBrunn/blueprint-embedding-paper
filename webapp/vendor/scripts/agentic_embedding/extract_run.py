#!/usr/bin/env python3
"""Extract a structured summary from one blueprint-run directory.

Usage:
    python analysis/extract_run.py <run_dir>

Reads <run_dir>/embedder_blueprint.py + <run_dir>/scenario_config.jsonc,
introspects the agent's edited build_blueprint(), and writes
<run_dir>/summary.json with module/axis/reference structure for
cross-run analysis. No data embedding is performed.
"""
from __future__ import annotations

import argparse
import datetime as _dt
import importlib.util
import json
import re
import sys
from pathlib import Path
from typing import Any


PLACEHOLDER_TOKEN = "<<EXAMPLE_"
DEFAULT_CLIP = 4.0


def _load_blueprint_module(run_dir: Path, run_id: str, source: str = "embedder_blueprint.py"):
    script_path = run_dir / source
    # Stem differentiates module cache entries when the same run is introspected
    # against both the refined and unrefined blueprints in one process.
    module_name = f"_blueprint_{run_id}_{Path(source).stem}"
    spec = importlib.util.spec_from_file_location(module_name, script_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not create import spec for {script_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def _load_scenario(run_dir: Path) -> dict[str, Any]:
    text = (run_dir / "scenario_config.jsonc").read_text()
    # Strip line comments (// ...). The scenario file is small and
    # never embeds // inside a string literal in practice.
    stripped = re.sub(r"//.*", "", text)
    return json.loads(stripped)


def _axis_kind(axis) -> str:
    if axis.gate is not None:
        return "gated"
    if axis.negatives:
        return "contrast"
    return "plain"


def _module_role_in_axis(axis, module_name: str) -> list[str]:
    roles: list[str] = []
    if axis.positive == module_name:
        roles.append("positive")
    if module_name in axis.negatives:
        roles.append("negative")
    if axis.gate == module_name:
        roles.append("gate")
    return roles


def _build_summary(
    run_dir: Path,
    source: str = "embedder_blueprint.py",
    coverage_path: Path | None = None,
) -> dict[str, Any]:
    run_id = run_dir.name
    module = _load_blueprint_module(run_dir, run_id, source)
    blueprint = module.build_blueprint("extract")

    modules = blueprint.modules        # dict[str, ModuleSpec]
    axes = blueprint.axes              # tuple[AxisSpec, ...]
    resources = blueprint.resources    # dict[str, Resource]

    total_genes = sum(len(m.genes) for m in modules.values())
    unique_genes: set[str] = set()
    for m in modules.values():
        unique_genes.update(m.genes)

    # resources cited per module (for the cross-reference back into resource records)
    resource_to_modules: dict[str, list[str]] = {k: [] for k in resources}
    for m in modules.values():
        for k in m.resource_keys:
            if k in resource_to_modules:
                resource_to_modules[k].append(m.name)

    n_axes_by_kind = {"plain": 0, "contrast": 0, "gated": 0}
    for a in axes:
        n_axes_by_kind[_axis_kind(a)] += 1

    n_resources_by_kind = {"paper": 0, "atlas": 0, "database": 0}
    for r in resources.values():
        n_resources_by_kind[r.kind] += 1

    module_records: list[dict[str, Any]] = []
    orphaned_modules: list[str] = []
    modules_missing_resources: list[str] = []
    for m in modules.values():
        used_in: list[dict[str, str]] = []
        for a in axes:
            for role in _module_role_in_axis(a, m.name):
                used_in.append({"axis": a.name, "role": role})
        if not used_in:
            orphaned_modules.append(m.name)
        rkeys = m.resource_keys
        if not rkeys:
            modules_missing_resources.append(m.name)
        module_records.append(
            {
                "name": m.name,
                "genes": list(m.genes),
                "n_genes": len(m.genes),
                "resource_keys": list(rkeys),
                "n_resource_keys": len(rkeys),
                "used_in_axes": used_in,
                "n_axes_using": len(used_in),
            }
        )

    axis_records: list[dict[str, Any]] = []
    for a in axes:
        axis_records.append(
            {
                "name": a.name,
                "kind": _axis_kind(a),
                "positive": a.positive,
                "negatives": list(a.negatives),
                "negative_weight": a.negative_weight,
                "gate": a.gate,
            }
        )

    resource_records: list[dict[str, Any]] = []
    unused_resources: list[str] = []
    for r in resources.values():
        users = sorted(resource_to_modules.get(r.key, []))
        if not users:
            unused_resources.append(r.key)
        resource_records.append(
            {
                "key": r.key,
                "kind": r.kind,
                "citation": r.citation,
                "doi": r.doi,
                "used_in_modules": users,
            }
        )

    placeholder_present = any(
        PLACEHOLDER_TOKEN in n
        for n in list(modules.keys()) + [a.name for a in axes]
    )

    clip = blueprint._clip
    set_clip_overridden = clip != DEFAULT_CLIP

    summary: dict[str, Any] = {
        "run_id": run_id,
        "extracted_at": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        "scenario": _load_scenario(run_dir),
        "blueprint": {
            "name": blueprint.name,
            "clip": clip,
            "n_modules": len(modules),
            "n_axes": len(axes),
            "n_axes_by_kind": n_axes_by_kind,
            "n_resources": len(resources),
            "n_resources_by_kind": n_resources_by_kind,
            "n_dois_unique": len({r.doi for r in resources.values() if r.doi}),
            "n_total_genes": total_genes,
            "n_unique_genes": len(unique_genes),
            "modules": module_records,
            "axes": axis_records,
            "resources": resource_records,
        },
        "sanity": {
            "placeholder_strings_present": placeholder_present,
            "orphaned_modules": orphaned_modules,
            "modules_missing_resources": modules_missing_resources,
            "unused_resources": sorted(unused_resources),
            "set_clip_overridden": set_clip_overridden,
        },
    }

    # Optional dataset-coverage merge. The embedder writes a sidecar
    # next to the output parquet during build_embedding; when the
    # caller threads its path in, fold per-module coverage fields into
    # blueprint.modules and add a top-level dataset_coverage header.
    # Absent sidecar → summary is unchanged.
    if coverage_path is not None and coverage_path.exists():
        coverage = json.loads(coverage_path.read_text())
        cov_by_name = {m["name"]: m for m in coverage.get("modules", [])}
        for rec in summary["blueprint"]["modules"]:
            cov = cov_by_name.get(rec["name"])
            if cov is None:
                continue
            rec["n_genes_present"] = cov["n_genes_present"]
            rec["n_genes_missing"] = cov["n_genes_missing"]
            rec["missing_genes"] = cov["missing_genes"]
            rec["fraction_nonzero_cells"] = cov["fraction_nonzero_cells"]
        summary["dataset_coverage"] = {
            "h5ad_path": coverage.get("h5ad_path"),
            "layer": coverage.get("layer"),
            "n_cells": coverage.get("n_cells"),
            "n_vars": coverage.get("n_vars"),
        }

        # Truthful clip: the sidecar records the clip actually applied to
        # this embedding, which may be a CLI --clip override that
        # build_blueprint() alone does not reflect (extract_run re-instantiates
        # the blueprint without the flag).
        if "clip" in coverage:
            summary["blueprint"]["clip"] = coverage["clip"]
            summary["sanity"]["set_clip_overridden"] = (
                coverage["clip"] != DEFAULT_CLIP
            )
        summary["blueprint"]["n_min_genes"] = coverage.get("n_min_genes", 0)

        # Prune axes to those actually emitted, in dim order, so that
        # dim_i <-> axes[i] holds for every positional downstream consumer
        # (axis_metrics, plot_axis_*, temporal_metrics, umap_branch_interactive).
        # Only rewrites when something was excluded; the default n_min_genes=0
        # keeps every axis, leaving the summary untouched.
        cov_axes = coverage.get("axes")
        if cov_axes is not None:
            excluded = [a["name"] for a in cov_axes if not a.get("kept")]
            if excluded:
                kept = sorted(
                    (a for a in cov_axes if a.get("kept")),
                    key=lambda a: a["dim"],
                )
                by_name = {r["name"]: r for r in summary["blueprint"]["axes"]}
                pruned = [by_name[a["name"]] for a in kept if a["name"] in by_name]
                kinds = {"plain": 0, "contrast": 0, "gated": 0}
                for r in pruned:
                    kinds[r["kind"]] += 1
                summary["blueprint"]["n_axes_declared"] = len(
                    summary["blueprint"]["axes"]
                )
                summary["blueprint"]["excluded_axes"] = excluded
                summary["blueprint"]["axes"] = pruned
                summary["blueprint"]["n_axes"] = len(pruned)
                summary["blueprint"]["n_axes_by_kind"] = kinds
    return summary


def _first_line(p: Path) -> str | None:
    """First non-empty line of a file, or None."""
    try:
        txt = p.read_text(errors="ignore")
    except OSError:
        return None
    for line in txt.splitlines():
        if line.strip():
            return line.strip()
    return None


def _render_axis_table(summary: dict[str, Any]) -> str:
    """The `## Embedding dimensions` markdown table — one row per axis, in dim
    order. Column 3 lists the axis's modules: positive first, then negatives,
    then `gate=<module>` for a gated axis. This is the exact format
    blueprint_tree.py and umap_dimension_plots.py parse."""
    lines = ["## Embedding dimensions", "", "| Axis | Kind | Modules |", "|---|---|---|"]
    for axis in summary["blueprint"]["axes"]:
        mods: list[str] = []
        if axis["positive"]:
            mods.append(axis["positive"])
        mods.extend(axis["negatives"])
        if axis["gate"]:
            mods.append(f"gate={axis['gate']}")
        lines.append(f"| {axis['name']} | {axis['kind']} | {', '.join(mods)} |")
    return "\n".join(lines)


def _render_report(summary: dict[str, Any], run_dir: Path, meta: dict[str, str],
                   source: str = "embedder_blueprint.py") -> str:
    """Render the <stem>_report.md companion. Carries the run metadata header +
    the `## Embedding dimensions` table, both derived from the introspected
    blueprint — no fragile shell heredoc, no plain-text --describe dump."""
    bp = summary["blueprint"]
    # When introspecting the pre-refinement snapshot (embedder_blueprint_unrefined.py),
    # the report describes Stage A regardless of whether the surrounding workdir
    # also went on to run Stage B.
    unrefined_source = source.endswith("_unrefined.py")
    with_refine = (run_dir / ".stage_b_iters").exists() and not unrefined_source
    stage = "Stage A + B refined" if with_refine else "Stage A only"
    by_kind = bp["n_axes_by_kind"]
    lines = [f"# Blueprint embedding ({stage})", ""]
    lines.append(f"- run: `{run_dir.name}`")
    if meta.get("dataset"):
        lines.append(f"- dataset: {meta['dataset']}")
    if meta.get("model"):
        eff = f" ({meta['reasoning']})" if meta.get("reasoning") else ""
        lines.append(f"- model: {meta['model']}{eff}")
    if with_refine:
        iters = _first_line(run_dir / ".stage_b_iters") or "?"
        rating = _first_line(run_dir / ".stage_b_final_rating") or "?"
        status = _first_line(run_dir / ".stage_b_status") or "?"
        lines.append(f"- refinement: {iters} iter(s), final rating=`{rating}` (status: `{status}`)")
    if meta.get("input"):
        lines.append(f"- input: `{meta['input']}`")
    if meta.get("parquet"):
        lines.append(f"- parquet: `{meta['parquet']}`")
    lines.append(
        f"- modules: {bp['n_modules']} · axes: {bp['n_axes']} "
        f"(plain {by_kind['plain']}, contrast {by_kind['contrast']}, gated {by_kind['gated']}) "
        f"· resources: {bp['n_resources']}"
    )
    lines.append("")
    lines.append(_render_axis_table(summary))
    lines.append("")
    return "\n".join(lines)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "run_dir",
        type=Path,
        help="Path to a work/<date>_blueprint_runNN/ directory.",
    )
    # Optional report-companion emission. When --report-out is given, the
    # `## Embedding dimensions` report markdown is written there in addition to
    # summary.json. The remaining flags only populate the report's metadata
    # header (the orchestrator passes them; they are not needed for summary.json).
    parser.add_argument("--report-out", type=Path, default=None,
                        help="Also write the <stem>_report.md companion to this path.")
    parser.add_argument("--dataset", default="", help="Report header: dataset name.")
    parser.add_argument("--model", default="", help="Report header: agent model.")
    parser.add_argument("--reasoning", default="", help="Report header: reasoning effort.")
    parser.add_argument("--input", default="", help="Report header: embedding input path.")
    parser.add_argument("--parquet", default="", help="Report header: output parquet path.")
    # Which blueprint .py to introspect inside <run_dir>. Default is the canonical
    # (post-refinement) blueprint; pass `embedder_blueprint_unrefined.py` to derive
    # summary + report from the Stage A snapshot instead.
    parser.add_argument("--blueprint-source", default="embedder_blueprint.py",
                        help="Blueprint .py file inside <run_dir> to introspect.")
    parser.add_argument("--summary-out", default=None,
                        help="Path to write the summary JSON. Defaults to "
                             "<run_dir>/summary.json for the canonical source and "
                             "<run_dir>/summary_unrefined.json for the unrefined source.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    if not run_dir.is_dir():
        print(f"not a directory: {run_dir}", file=sys.stderr)
        sys.exit(1)
    for required in (args.blueprint_source, "scenario_config.jsonc"):
        if not (run_dir / required).is_file():
            print(f"missing required file: {run_dir / required}", file=sys.stderr)
            sys.exit(1)

    coverage_path: Path | None = None
    if args.parquet:
        parquet_path = Path(args.parquet)
        if not parquet_path.is_absolute():
            # The orchestrator passes --parquet as a repo-relative path
            # (REPO_ROOT-stripped); reconstruct REPO_ROOT from run_dir
            # (work/<cohort>/<run>/). Absolute paths pass through.
            parquet_path = run_dir.parents[2] / parquet_path
        coverage_path = parquet_path.parent / f"{parquet_path.stem}_coverage.json"

    summary = _build_summary(run_dir, args.blueprint_source, coverage_path)
    if args.summary_out is not None:
        out_path = Path(args.summary_out)
        if not out_path.is_absolute():
            out_path = run_dir / out_path
    elif args.blueprint_source.endswith("_unrefined.py"):
        out_path = run_dir / "summary_unrefined.json"
    else:
        out_path = run_dir / "summary.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=False)
        handle.write("\n")
    print(f"Wrote {out_path}")

    if args.report_out is not None:
        meta = {
            "dataset": args.dataset, "model": args.model, "reasoning": args.reasoning,
            "input": args.input, "parquet": args.parquet,
        }
        args.report_out.parent.mkdir(parents=True, exist_ok=True)
        args.report_out.write_text(_render_report(summary, run_dir, meta, args.blueprint_source))
        print(f"Wrote {args.report_out}")


if __name__ == "__main__":
    main()
