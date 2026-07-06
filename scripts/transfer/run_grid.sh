#!/usr/bin/env bash
# Sequential driver for the held-out reference-to-query transfer grid, one run at
# a time so the query-side surgery wall-clock is not contaminated by concurrent
# CPU load. Covers mouse_cortex and human_pancreas, seeds 1-5.
#
# Prerequisites: `conda activate agentic-embeddings` (see environment.yml).
set -u
REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

MC_BP="outs/blueprints/mouse_cortex/gpt_5_5_unrefined_high/embedder_blueprint_seed"
HP_BP="outs/blueprints/human_pancreas/gpt_5_5_unrefined_high/embedder_blueprint_seed"

run () {  # ds seed blueprint_prefix
  local ds="$1" s="$2" bp="$3"
  echo "=== RUN ds=${ds} seed=${s} $(date +%H:%M:%S) ==="
  python scripts/transfer/measure.py --dataset "${ds}" \
    --methods blueprint,scvi,expimap --blueprint-run "${bp}${s}.py" --seed "${s}"
  echo "=== DONE ds=${ds} seed=${s} rc=$? $(date +%H:%M:%S) ==="
}

for s in 1 2 3 4 5; do run mouse_cortex "${s}" "${MC_BP}"; done
for s in 1 2 3 4 5; do run human_pancreas "${s}" "${HP_BP}"; done
echo "=== GRID COMPLETE $(date +%H:%M:%S) ==="
