#!/usr/bin/env bash
# Reproduce the manuscript's quantitative and structural figures from the results
# already committed under outs/ (metric CSVs + authored blueprints). This path does
# NOT download data, run baselines, or re-author blueprints — it re-plots from the
# shipped artifacts, so it runs in minutes with only the conda env.
#
# Covered here: the scIB frontier + leaderboard tables, the blueprint-vs-freeform
# three-way summary, cross-run reproducibility, resource cost, the refinement
# ablations, the developmental-tree temporal panels + axis tree, the held-out
# transfer bars, and the authored-blueprint cards.
#
# NOT covered (they need the raw data + regenerated embeddings): the UMAP grids,
# marker-gene dot plots, axis-faithfulness matrices, per-axis UMAPs, and branch
# filmstrips. For those, and for a full end-to-end rebuild, run:
#     bash analysis/reproduce_all.sh
#
# Prerequisites: `conda activate agentic-embeddings` (see environment.yml).
# Figures are written under outs/plots/.
#
#   bash analysis/reproduce_figures.sh          # regenerate all light-path figures
#   bash analysis/reproduce_figures.sh --help
set -uo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

case "${1:-}" in
  -h|--help)
    sed -n '2,25p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    exit 0 ;;
esac

PY="python"
PLOTS="scripts/plots"
ok=0; fail=0; failed_cmds=()

run() {  # run a plot step; keep going on failure and tally the result
  printf '\n>>> %s\n' "$*"
  if "$@"; then ok=$((ok + 1)); else
    fail=$((fail + 1)); failed_cmds+=("$*"); printf '[FAILED] %s\n' "$*" >&2
  fi
}

# --- guard: the shipped metric CSVs must be present -------------------------
if [[ ! -f outs/metrics/mouse_cortex/metrics.csv ]]; then
  echo "error: shipped results not found under outs/metrics/." >&2
  echo "This script re-plots the committed CSVs; for a full rebuild from data run" >&2
  echo "  bash analysis/reproduce_all.sh" >&2
  exit 1
fi

# --- stage blueprint metadata into the outs/embeddings/ layout --------------
# Copies the shipped summary/coverage sidecars so the metadata figures
# (blueprint_card, three_way_summary, axis tree) find them without embeddings.
run "$PY" scripts/reproduce/apply_blueprints.py --metadata-only

# --- scIB benchmark: averaged frontier + per-dataset + leaderboard tables ---
run "$PY" "$PLOTS/plot_frontier_avg.py" --per-dataset --overwrite
run "$PY" "$PLOTS/build_metrics_table.py"
for ds in mouse_cortex human_pancreas human_immune; do
  # plot_scib_table needs scib-metrics; if absent it fails here but build_metrics_table
  # above still carries the same numbers.
  run "$PY" "$PLOTS/plot_scib_table.py" --dataset "$ds" --overwrite
done

# --- structured-framework analyses (human pancreas) -------------------------
run "$PY" "$PLOTS/plot_reproducibility.py" --dataset human_pancreas
run "$PY" "$PLOTS/plot_three_way_summary.py" --dataset human_pancreas --model-tag gpt_5_5_high
run "$PY" "$PLOTS/plot_blueprint_repertoire.py" --dataset human_pancreas
run "$PY" "$PLOTS/plot_resource_costs.py" --dataset human_pancreas --overwrite
run "$PY" "$PLOTS/plot_refinement_ablation.py" --dataset human_pancreas

# --- developmental-tree showcase (mouse gastrulation) -----------------------
run "$PY" "$PLOTS/plot_temporal_metrics.py" --dataset mouse_gastrulation
run "$PY" "$PLOTS/plot_refinement_effect.py" --dataset mouse_gastrulation
run "$PY" "$PLOTS/blueprint_tree.py" --dataset mouse_gastrulation --overwrite

# --- held-out transfer + authored-blueprint cards --------------------------
for ds in mouse_cortex human_pancreas human_immune; do
  run "$PY" "$PLOTS/plot_transfer.py" --dataset "$ds"
  run "$PY" "$PLOTS/blueprint_card.py" --dataset "$ds" --overwrite
done

# --- summary ----------------------------------------------------------------
printf '\n=== reproduce_figures: %d ok, %d failed ===\n' "$ok" "$fail"
if (( fail > 0 )); then
  printf 'failed steps:\n'; printf '  %s\n' "${failed_cmds[@]}"
  echo "(plot_scib_table.py needs scib-metrics installed; see environment.yml)"
fi
echo "Figures written under outs/plots/."
echo "Embedding-visualization figures (UMAPs, dot plots, axis faithfulness) require"
echo "  bash analysis/reproduce_all.sh"
(( fail == 0 ))
