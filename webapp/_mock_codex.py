#!/usr/bin/env python3
"""Mock codex for token-free e2e: populate ./embedder_blueprint.py by copying a
real archived blueprint. Activated by the driver when WEBAPP_MOCK_CODEX=1."""
import shutil
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    src = next(
        (REPO_ROOT / "outs" / "blueprints").rglob("embedder_blueprint_seed*.py"),
        None,
    )
    if src is None:
        print("mock_codex: no archived blueprint to copy", file=sys.stderr)
        return 1
    shutil.copy(src, Path.cwd() / "embedder_blueprint.py")
    print(f"mock_codex: populated embedder_blueprint.py from {src}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
