#!/usr/bin/env python3
"""Detect unauthorized removals between two blueprint versions.

Usage:
    python regression_guard.py <workdir> <before.py> <after.py> [--report PATH]

Loads both blueprint files via importlib, calls ``build_blueprint("guard")``
on each, and diffs the (resource, module, axis) key sets. If ``after`` is
missing any keys not listed in the sanctioned-removal set from
``reports/report.json``, prints the unauthorized removals and exits 1.

Sanctioned removal format in report.json: each issue may carry a
``remove: ["<kind>:<name>", ...]`` array. ``<kind>`` ∈ {resource, module, axis}.
Whitespace around entries is tolerated; case is preserved (matched exactly).

Exits 0 if no unauthorized removals (including the case of zero removals).
"""
from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path


def _load_blueprint(source: Path, tag: str):
    module_name = f"_blueprint_guard_{tag}"
    spec = importlib.util.spec_from_file_location(module_name, source)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load spec for {source}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module.build_blueprint("guard")


def _keys(blueprint) -> set[str]:
    keys: set[str] = set()
    for k in blueprint.resources:
        keys.add(f"resource:{k}")
    for m in blueprint.modules:
        keys.add(f"module:{m}")
    for axis in blueprint.axes:
        keys.add(f"axis:{axis.name}")
    return keys


def _sanctioned(report_path: Path) -> set[str]:
    if not report_path.exists():
        return set()
    try:
        data = json.loads(report_path.read_text())
    except Exception:
        return set()
    sanctioned: set[str] = set()
    for issue in data.get("issues", []) or []:
        for entry in issue.get("remove", []) or []:
            if isinstance(entry, str) and entry.strip():
                sanctioned.add(entry.strip())
    return sanctioned


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workdir", type=Path)
    parser.add_argument("before", type=Path, help="pre-Generator blueprint snapshot")
    parser.add_argument("after", type=Path, help="post-Generator blueprint")
    parser.add_argument("--report", type=Path, default=None,
                        help="defaults to <workdir>/reports/report.json")
    args = parser.parse_args()

    report_path = args.report or (args.workdir / "reports" / "report.json")
    before = _load_blueprint(args.before, "before")
    after = _load_blueprint(args.after, "after")
    removed = sorted(_keys(before) - _keys(after))
    sanctioned = _sanctioned(report_path)
    unauthorized = [k for k in removed if k not in sanctioned]
    if unauthorized:
        print("REGRESSION-DETECTED")
        for u in unauthorized:
            print(f"  unauthorized removal: {u}")
        return 1
    print("REGRESSION-CLEAN")
    return 0


if __name__ == "__main__":
    sys.exit(main())
