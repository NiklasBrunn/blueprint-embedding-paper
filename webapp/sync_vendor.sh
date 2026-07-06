#!/usr/bin/env bash
# Copy the canonical pipeline sources from the parent repo into webapp/vendor/,
# mirroring the repo's directory layout, so the standalone app can resolve them
# WITHOUT the repo present. Run this in-repo whenever the originals change; the
# drift-guard test (test_webapp_vendor.py) fails if vendor/ falls out of sync.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${HERE}/.." && pwd)"
VENDOR="${HERE}/vendor"

FILES=(
  "scripts/agentic_embedding/embedder_blueprint.py"
  "scripts/agentic_embedding/embedder_blueprint_noref.py"
  "scripts/agentic_embedding/extract_run.py"
  "scripts/agentic_embedding/regression_guard.py"
  "scripts/agentic_embedding/extract_cost.py"
  "prompt_templates/agent_prompt_blueprint_edit.md"
  "prompt_templates/agent_prompt_blueprint_edit_noref.md"
  "prompt_templates/agent_prompt_blueprint_discriminate.md"
  "prompt_templates/agent_prompt_blueprint_discriminate_noref.md"
  "prompt_templates/agent_prompt_blueprint_generate.md"
  "prompt_templates/agent_prompt_blueprint_generate_noref.md"
)

for rel in "${FILES[@]}"; do
  src="${REPO_ROOT}/${rel}"
  dst="${VENDOR}/${rel}"
  if [ ! -f "${src}" ]; then
    echo "MISSING in repo: ${rel}" >&2
    exit 1
  fi
  mkdir -p "$(dirname "${dst}")"
  cp "${src}" "${dst}"
  echo "vendored ${rel}"
done
echo "Done. Vendored ${#FILES[@]} files into ${VENDOR}"
