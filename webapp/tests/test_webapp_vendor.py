"""Standalone vendoring: webapp/vendor/ holds byte-identical copies of the
canonical pipeline sources so the app runs without the parent repo. These tests
guard against (a) a missing vendored file and (b) drift from the repo original
(when the repo is present). In a true standalone deployment the repo originals
are absent and the drift check skips."""
import filecmp
from pathlib import Path

import pytest

WEBAPP = Path(__file__).resolve().parents[1]
REPO = WEBAPP.parent
VENDOR = WEBAPP / "vendor"

VENDORED = [
    "scripts/agentic_embedding/embedder_blueprint.py",
    "scripts/agentic_embedding/embedder_blueprint_noref.py",
    "scripts/agentic_embedding/extract_run.py",
    "scripts/agentic_embedding/regression_guard.py",
    "scripts/agentic_embedding/extract_cost.py",
    "prompt_templates/agent_prompt_blueprint_edit.md",
    "prompt_templates/agent_prompt_blueprint_edit_noref.md",
    "prompt_templates/agent_prompt_blueprint_discriminate.md",
    "prompt_templates/agent_prompt_blueprint_discriminate_noref.md",
    "prompt_templates/agent_prompt_blueprint_generate.md",
    "prompt_templates/agent_prompt_blueprint_generate_noref.md",
]


def test_all_sources_vendored():
    missing = [rel for rel in VENDORED if not (VENDOR / rel).exists()]
    assert not missing, f"vendor/ missing {missing}; run webapp/sync_vendor.sh"


@pytest.mark.parametrize("rel", VENDORED)
def test_vendor_matches_repo(rel):
    repo_file = REPO / rel
    if not repo_file.exists():
        pytest.skip(f"repo original absent (standalone deployment): {rel}")
    assert filecmp.cmp(repo_file, VENDOR / rel, shallow=False), (
        f"vendor/{rel} drifted from the repo — run webapp/sync_vendor.sh")


def test_src_root_prefers_vendor():
    from webapp.pipeline import SRC_ROOT, VENDOR as PVENDOR
    assert SRC_ROOT == PVENDOR  # vendor/ exists -> it is preferred over the repo
