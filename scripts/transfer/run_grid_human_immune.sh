#!/usr/bin/env bash
# Sequential driver for the human_immune held-out transfer (one run at a time so
# surgery wall-clock stays clean). Query = held-out batches Oetjen_U (bone marrow)
# + Sun_sample4_TC (PBMC). Both blueprint arms (resourced + noref) are produced per run.
# Prerequisites: `conda activate agentic-embeddings` (see environment.yml).
set -u
REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$REPO_ROOT"

C="outs/blueprints/human_immune/gpt_5_5_unrefined_high/embedder_blueprint_seed"
N="outs/blueprints/human_immune/gpt_5_5_noref_unrefined_high/embedder_blueprint_seed"

for s in 1 2 3 4 5; do
  echo "=== RUN ds=human_immune seed=${s} $(date +%H:%M:%S) ==="
  python scripts/transfer/measure.py --dataset human_immune \
    --methods blueprint,scvi,expimap \
    --blueprint-run "${C}${s}.py" --blueprint-noref-run "${N}${s}.py" --seed "${s}"
  echo "=== DONE ds=human_immune seed=${s} rc=$? $(date +%H:%M:%S) ==="
done
echo "=== HUMAN_IMMUNE COMPLETE $(date +%H:%M:%S) ==="
