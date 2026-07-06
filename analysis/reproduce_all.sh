#!/usr/bin/env bash
# Full end-to-end reproduction of the manuscript results, WITHOUT re-running the
# agentic authoring pipeline. It consumes the authored blueprints already committed
# under outs/blueprints/, applies them to the datasets, runs every baseline,
# recomputes the metric CSVs, and renders every figure — including the
# embedding-visualization figures that the light analysis/reproduce_figures.sh
# cannot make (UMAP grids, marker dot plots, axis-faithfulness matrices, branch
# filmstrips).
#
# This is the HEAVY path. Prerequisites beyond the conda env (environment.yml):
#   - Several GB of dataset downloads (run_data_pipeline.sh).
#   - The foundation-model baselines: the UCE checkpoint (~5.7 GB) and the
#     CZ CELLxGENE Census cache (scvi_census). scvi-tools + pytorch.
#   - Time: hours, mostly baseline training + the held-out transfer grid.
#
# Not recomputed (kept as shipped): outs/resource_costs/*.csv — the agent token/
# wall-clock telemetry is captured during authoring (needs work/ run transcripts),
# which this script deliberately does not reproduce.
#
# Prerequisites: `conda activate agentic-embeddings`.
#
# Stages run in order; pass any subset of flags to run only those (default: all):
#   --data        download + QC + preprocess the datasets
#   --embed       apply the shipped blueprints (+ free-form scripts) to the data
#   --baselines   PCA / scVI / scVI-census / UCE / UCell / expiMap
#   --metrics     recompute scIB / axis / temporal / variability / reproducibility
#   --plots       render every figure (light set + embedding-visualization figures)
# Other options:
#   --datasets "a b c"   restrict to these datasets (default: all four)
#   --overwrite          rebuild artifacts that already exist
#   -h | --help
#
#   bash analysis/reproduce_all.sh                       # everything
#   bash analysis/reproduce_all.sh --embed --metrics --plots   # data already present
#   bash analysis/reproduce_all.sh --plots --datasets "mouse_cortex"
set -uo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

ALL_DATASETS=(mouse_cortex human_pancreas human_immune mouse_gastrulation)
BENCH=(mouse_cortex human_pancreas human_immune)      # the scIB benchmark trio
PY="python"; PLOTS="scripts/plots"; MODELS="scripts/models"
DO_DATA=0 DO_EMBED=0 DO_BASELINES=0 DO_METRICS=0 DO_PLOTS=0 ANY=0
OVERWRITE=""; DATASETS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --data)      DO_DATA=1; ANY=1 ;;
    --embed)     DO_EMBED=1; ANY=1 ;;
    --baselines) DO_BASELINES=1; ANY=1 ;;
    --metrics)   DO_METRICS=1; ANY=1 ;;
    --plots)     DO_PLOTS=1; ANY=1 ;;
    --overwrite) OVERWRITE="--overwrite" ;;
    --datasets)  shift; read -r -a DATASETS <<< "${1:-}" ;;
    -h|--help)   sed -n '2,38p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'; exit 0 ;;
    *) echo "unknown arg: $1" >&2; exit 2 ;;
  esac
  shift
done
if (( ANY == 0 )); then DO_DATA=1; DO_EMBED=1; DO_BASELINES=1; DO_METRICS=1; DO_PLOTS=1; fi
(( ${#DATASETS[@]} )) || DATASETS=("${ALL_DATASETS[@]}")

# Intersection helper: is $1 in the selected DATASETS?
selected() { local d; for d in "${DATASETS[@]}"; do [[ "$d" == "$1" ]] && return 0; done; return 1; }
in_bench() { local d; for d in "${BENCH[@]}"; do [[ "$d" == "$1" ]] && return 0; done; return 1; }

ok=0; fail=0; failed_cmds=()
run() {
  printf '\n>>> %s\n' "$*"
  if "$@"; then ok=$((ok + 1)); else
    fail=$((fail + 1)); failed_cmds+=("$*"); printf '[FAILED] %s\n' "$*" >&2
  fi
}
stage() { printf '\n========== %s ==========\n' "$1"; }

# --- Stage 0: data ----------------------------------------------------------
if (( DO_DATA )); then
  stage "data pipeline"
  run bash analysis/run_data_pipeline.sh "${DATASETS[@]}"
fi

# --- Stage 1: apply shipped blueprints (+ free-form) ------------------------
if (( DO_EMBED )); then
  stage "apply blueprints -> outs/embeddings/"
  ds_flags=(); for d in "${DATASETS[@]}"; do ds_flags+=(--dataset "$d"); done
  run "$PY" scripts/reproduce/apply_blueprints.py --freeform "${ds_flags[@]}" ${OVERWRITE}
fi

# --- Stage 2: baselines -----------------------------------------------------
if (( DO_BASELINES )); then
  stage "baseline embeddings"
  for d in "${DATASETS[@]}"; do
    if in_bench "$d"; then                                   # full comparator set
      run "$PY" "$MODELS/pca.py"         --dataset "$d"
      run "$PY" "$MODELS/scvi.py"        --dataset "$d"
      run "$PY" "$MODELS/scvi.py"        --dataset "$d" --batch-key none   # scvi_nobatch
      run "$PY" "$MODELS/scvi_census.py" --dataset "$d"
      run "$PY" "$MODELS/uce.py"         --dataset "$d"
      run "$PY" "$MODELS/ucell.py"       --dataset "$d"
      run "$PY" "$MODELS/expimap.py"     --dataset "$d"
    elif [[ "$d" == "mouse_gastrulation" ]]; then            # DPT-substrate baselines
      run "$PY" "$MODELS/pca.py"         --dataset "$d"
      run "$PY" "$MODELS/scvi.py"        --dataset "$d" --batch-key none   # scvi_nobatch
      run "$PY" "$MODELS/scvi_census.py" --dataset "$d"
    fi
  done
fi

# --- Stage 3: recompute metric CSVs -----------------------------------------
if (( DO_METRICS )); then
  stage "metrics"
  for d in "${DATASETS[@]}"; do
    run "$PY" scripts/embedding_metrics/metrics.py --dataset "$d" ${OVERWRITE}
    run "$PY" scripts/axis_metrics/metrics.py      --dataset "$d" ${OVERWRITE}
  done
  for d in "${DATASETS[@]}"; do
    in_bench "$d" || continue
    run "$PY" scripts/blueprint_variability_metrics/metrics.py --dataset "$d" ${OVERWRITE}
    run "$PY" scripts/embedding_reproducibility/metrics.py     --dataset "$d" ${OVERWRITE}
  done
  if selected mouse_gastrulation; then
    run "$PY" scripts/temporal_metrics/metrics.py             --dataset mouse_gastrulation
    run "$PY" scripts/temporal_metrics/pseudotime_substrate.py --dataset mouse_gastrulation
  fi
  # Held-out reference-to-query transfer (trains scVI/expiMap on the reference).
  if selected mouse_cortex || selected human_pancreas; then
    run bash scripts/transfer/run_grid.sh
  fi
  if selected human_immune; then
    run bash scripts/transfer/run_grid_human_immune.sh
  fi
  echo "note: outs/resource_costs/*.csv is authoring-time telemetry and is kept as shipped."
fi

# --- Stage 4: figures -------------------------------------------------------
if (( DO_PLOTS )); then
  stage "figures — light (CSV-driven) set"
  run bash analysis/reproduce_figures.sh

  stage "figures — embedding-visualization set"
  # Interpretability / faithfulness (mouse cortex showcase).
  run "$PY" "$PLOTS/plot_axis_faithfulness.py" --dataset mouse_cortex ${OVERWRITE}
  run "$PY" "$PLOTS/plot_gene_dotplot.py"      --dataset mouse_cortex ${OVERWRITE}
  run "$PY" "$PLOTS/plot_axis_assignment.py"   --dataset mouse_cortex ${OVERWRITE}
  run "$PY" "$PLOTS/umap_dimension_plots.py"   --dataset mouse_cortex ${OVERWRITE}
  # UMAP coords (cached) + per-method grids for every dataset.
  for d in "${DATASETS[@]}"; do
    run "$PY" "$PLOTS/umap_plots.py" --dataset "$d"
    in_bench "$d" && run "$PY" "$PLOTS/umap_grid.py"        --dataset "$d" ${OVERWRITE}
    in_bench "$d" && run "$PY" "$PLOTS/plot_transfer_umap.py" --dataset "$d"
  done
  # Gastrulation branch filmstrips.
  if selected mouse_gastrulation; then
    run "$PY" "$PLOTS/umap_branch_filmstrip.py" --dataset mouse_gastrulation \
        --input outs/embeddings/mouse_gastrulation/blueprint_gpt_5_5_refine_high_mouse_gastrulation_seed3.parquet
  fi
fi

# --- summary ----------------------------------------------------------------
printf '\n=== reproduce_all: %d ok, %d failed ===\n' "$ok" "$fail"
if (( fail > 0 )); then
  printf 'failed steps:\n'; printf '  %s\n' "${failed_cmds[@]}"
fi
echo "Figures under outs/plots/, recomputed CSVs under outs/metrics|axis_metrics|..."
(( fail == 0 ))
