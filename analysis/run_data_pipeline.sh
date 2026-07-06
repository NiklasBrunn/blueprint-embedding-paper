#!/usr/bin/env bash
# Run the data pipeline end-to-end: load raw -> QC -> preprocess.
# Idempotent: source downloads are cached; raw h5ad and downstream artifacts
# are re-written each run.
#
# Prerequisites: `conda activate agentic-embeddings` (see environment.yml).
#
# Available datasets:
#   mouse_brain         Tasic et al. 2016, V1 cortex SMART-Seq
#   human_pancreas      scIB integration benchmark (Luecken et al. 2022, 6 studies / 9 batches)
#   mouse_gastrulation  Pijuan-Sala et al. 2019, gastrulation atlas E6.5-E8.5 (cellxgene census)
#   mouse_cortex        Integrated three-source adult-mouse isocortex (Tasic 2016 SMART-Seq +
#                       Tasic 2018 SMART-Seq v4 + Zeisel 2018 10x cortex slice), harmonized to
#                       the Tasic-2018 subclass vocabulary for cross-protocol integratability.
#   human_immune        scIB "Immune (human)" integration benchmark (Luecken et al. 2022) —
#                       PBMC + bone marrow across many studies / chemistries (10x, Smart-seq2),
#                       10 batches in obs["batch"], labels in obs["final_annotation"].
#
# Dataset selection (positional args; choose any subset, default = all):
#   bash analysis/run_data_pipeline.sh                                                                      # all datasets (default)
#   bash analysis/run_data_pipeline.sh mouse_brain                                                          # one only
#   bash analysis/run_data_pipeline.sh mouse_brain human_pancreas mouse_gastrulation mouse_cortex human_immune  # explicit all
#
# Force re-download of source files (otherwise cached files in data/raw are reused):
#   bash analysis/run_data_pipeline.sh --force-download                       # both
#   bash analysis/run_data_pipeline.sh --force-download human_pancreas        # one
#
#   bash analysis/run_data_pipeline.sh --help                     # this message
set -euo pipefail

REPO_ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

FORCE_DOWNLOAD=0
DATASETS=()
for arg in "$@"; do
    case "$arg" in
        --force-download)                                       FORCE_DOWNLOAD=1 ;;
        mouse_brain|human_pancreas|mouse_gastrulation|mouse_cortex|human_immune)    DATASETS+=("$arg") ;;
        -h|--help)
            cat <<EOF
Usage: $(basename "$0") [--force-download] [mouse_brain] [human_pancreas] [mouse_gastrulation] [mouse_cortex] [human_immune]

Runs load -> QC -> preprocess for the named datasets
(default: all five). With --force-download, source files are re-fetched even
if cached locally (note: mouse_gastrulation uses a cellxgene_census cache and
--force-download is ignored for it).
mouse_cortex downloads ~2.5 GB (Zeisel l5_all.loom) on first use and caches a
filtered isocortex subset for subsequent reruns.
EOF
            exit 0 ;;
        *) echo "unknown arg: $arg" >&2; exit 2 ;;
    esac
done
[ "${#DATASETS[@]}" -eq 0 ] && DATASETS=(mouse_brain human_pancreas mouse_gastrulation mouse_cortex human_immune)

LOAD_FLAGS=(--overwrite)
[ "$FORCE_DOWNLOAD" -eq 1 ] && LOAD_FLAGS+=(--force-download)

stage() { printf "\n=== %s ===\n" "$1"; }

for ds in "${DATASETS[@]}"; do
    stage "[$ds] 1/3 load raw"
    python "scripts/IO/load_${ds}.py" "${LOAD_FLAGS[@]}"

    stage "[$ds] 2/3 qc"
    python scripts/data_qc/qc.py --dataset "$ds"

    stage "[$ds] 3/3 preprocess"
    python scripts/data_preprocessing/preprocess.py --dataset "$ds"
done

printf "\n[ok] pipeline complete for: %s\n" "${DATASETS[*]}"
