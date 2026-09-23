#!/usr/bin/env bash
# ============================================================
# WEBAPP-VENDORED COPY of analysis/agent_edit_blueprint.sh.
# This file lives under webapp/ and is the only edited copy; the
# original analysis/ script is never touched. SCRIPT_DIR/.. still
# resolves to the repo root, so the canonical prompt templates and
# embedder_blueprint.py are found read-only.
#
# It adds caller overrides so the web app can drive the pipeline
# headless against an arbitrary scenario + input:
#   --scenario-file PATH   use this scenario_config.jsonc (instead of
#                          configs/scenario_config_<dataset>.jsonc).
#                          When given, --dataset is optional (defaults
#                          to the `webapp` cohort tag).
#   --input PATH           embed against this h5ad (instead of
#                          data/preprocessed/<dataset>_pp.h5ad).
#   --work-root PATH       place workdirs under this root (instead of
#                          <repo>/work).
#   --non-interactive      do not prompt for the parallel-run count;
#                          take it from --n-parallel (default 1).
#   --n-parallel N         parallel-run count for --non-interactive.
#   --print-config         print the resolved {scenario_file,
#                          blueprint_input, work_root, n_parallel,
#                          dataset, model} JSON and exit 0, before any
#                          workdir creation or codex execution.
# Both `--flag value` and `--flag=value` forms are accepted.
# ============================================================
#
# Two-stage blueprint pipeline:
#   Stage A populates,
#   Stage B iteratively refines via a Discriminator/Generator loop
#     (opt-in via --with-refine).
# A final embedding phase runs each accepted blueprint against the
# dataset's data to produce a comparable parquet under outs/embeddings/.
# All codex session transcripts and the embedder's stdout log live
# under <workdir>/logs/ — see project_workdir_logs_convention.
#
# WORKFLOW STEPS (per session, run sequentially; Stage A, Stage B and Embed
# are parallel across workdirs internally):
#
#   Stage A    — populator. Each workdir gets a fresh copy of
#                embedder_blueprint.py + scenario_config + edit prompt.
#                Codex edits build_blueprint() to declare resources, gene
#                modules, and axes. Workdirs whose codex exec exits != 0
#                are flagged A-FAIL and skipped from later steps.
#   artifact staging  — Per-flag artifact copy into A-OK workdirs only:
#                  · --with-refine → discriminate + generate prompts;
#                                    reports/ subdirectory created.
#   Extract (Stage A)  — Early extraction. extract_run.py runs on every workdir
#                and writes summary.json from the populator-only state.
#                Extract (Stage B) overwrites this later if Stage B runs.
#   Stage B    — iterative refinement, isolated per-iter
#                context). Skipped unless --with-refine. Snapshots
#                embedder_blueprint.py to embedder_blueprint_unrefined.py
#                once at loop start. Each iter runs in a workdir
#                stripped of prior reports / prior refinement logs /
#                the Stage A transcript / the Stage A snapshot
#                (all moved into .refine_hidden/), so the freshly
#                spawned codex session sees a clean slate. Per iter
#                (up to --max-iterations, default 3):
#                  · Discriminator (fresh codex session) judges the
#                    current blueprint and writes
#                    reports/report.{md,json} (no round tag) with
#                    rating ∈ {fine, minor, substantial}.
#                  · If rating == "fine" → loop exits (no Generator).
#                  · Else if iter < max → embedder_blueprint_unrefined.py
#                    is restored (audit trail), the pre-Generator
#                    blueprint is snapshotted, the Generator rewrites
#                    the blueprint, and regression_guard.py compares
#                    pre vs post key sets against the report's
#                    `remove` directives. Unauthorized removal →
#                    revert and re-spawn once with a corrective
#                    addendum; second regression marks the iter
#                    REGRESSION-BLOCKED and the loop continues.
#                  · Else → MAX-ITER exit.
#                Post-loop, .refine_hidden/iter*/ contents are
#                restored to reports/report_iter<NN>.{md,json} and
#                logs/codex_<role>_<NN>.log; reports/INDEX.md is
#                reconstructed from the per-iter JSONs. Records final
#                iter + rating + wall-clock to dotfiles.
#   Extract (Stage B)    — Late extraction. Re-runs extract_run.py on every workdir,
#                overwriting Extract (Stage A)'s summary.json with one that
#                reflects Stage B refinement edits. Skipped when
#                --with-refine is off (Extract (Stage A)'s summary is already
#                final in that case).
#   Finalize    — Session-log finalization. Appends finished_at and a
#                per-workdir status table.
#   Embed    — Embedding (default on; opt out with --no-embed). Runs
#                each accepted workdir's populated embedder_blueprint.py
#                against data/preprocessed/<dataset>_pp.h5ad and writes
#                the parquet plus a report companion directly under
#                outs/embeddings/<dataset>/, named
#                blueprint_<model_slug>_<reasoning>_<dataset>_seed<N>.parquet
#                so it sits alongside other models' embeddings and
#                aggregates as replicates in metrics/plot scripts.
#                The <model_slug> always carries a refinement-status tag:
#                `_unrefined` for Stage-A-only runs, `_refine` when
#                --with-refine is set. Eligibility: Stage A must be OK,
#                and Stage B (if enabled) must not have hard-failed.
#                MAX-ITER refinement plateaus are still eligible.
#                When --with-refine is on, the Stage A snapshot
#                (embedder_blueprint_unrefined.py) is *also* embedded
#                into a sibling parquet under the `_unrefined` cohort
#                slug. Per-workdir Stage A artifacts: logs/embed_unrefined.log,
#                .embed_unrefined_{secs,status,detail,rc}, plus a
#                Stage-A-only _report.md + _summary.json in outs/embeddings/
#                and an archived embedder_blueprint_seed<N>.py under
#                outs/blueprints/<dataset>/<model>_unrefined_<reasoning>/.
#   Cost    — Cost extraction. extract_cost.py writes cost.json
#                per workdir (tokens, per-phase wall-clock, footprint).
#
# HOW TO RUN (from the repo root). --dataset, --model and --reasoning-effort
# are required; they feed the work/<...> parent dir name and Embed output
# layout. Other flags are optional.
#
#   bash analysis/agent_edit_blueprint.sh \
#       --dataset mouse_brain \
#       --model gpt-5.5 --reasoning-effort high
#
#   # Cheaper sweep:
#   bash analysis/agent_edit_blueprint.sh \
#       --dataset human_pancreas \
#       --model gpt-5.6-luna --reasoning-effort medium
#
#   # Add Stage B refinement (Discriminator/Generator loop, default cap 3 iters):
#   bash analysis/agent_edit_blueprint.sh \
#       --dataset mouse_brain \
#       --model gpt-5.5 --reasoning-effort high \
#       --with-refine
#
#   # Tighter refinement cap:
#   bash analysis/agent_edit_blueprint.sh \
#       --dataset mouse_brain \
#       --model gpt-5.5 --reasoning-effort high \
#       --with-refine --max-iterations 3
#
#   # Skip the dataset-level embedding (no parquet written):
#   bash analysis/agent_edit_blueprint.sh \
#       --dataset mouse_brain \
#       --model gpt-5.5 --reasoning-effort high --no-embed
#
#   # Resource-unconstrained variant (no resources in the prompt
#   # or the API). Cohort tag _noref_unrefined; Stage A only.
#   bash analysis/agent_edit_blueprint.sh \
#       --dataset mouse_brain \
#       --model gpt-5.5 --reasoning-effort high --noref
#
#   # See all flags:
#   bash analysis/agent_edit_blueprint.sh --help
#
# FLAGS:
#   --dataset NAME             REQUIRED. mouse_brain | human_pancreas | mouse_gastrulation | mouse_cortex | human_immune.
#   --model SLUG               REQUIRED. e.g. gpt-5.5, gpt-5.6-luna.
#   --reasoning-effort LEVEL   REQUIRED. low | medium | high | xhigh.
#   --with-refine              opt INTO Stage B (Discriminator/Generator refinement).
#                              When set, MODEL_SLUG carries a `_refine` suffix; when
#                              unset, it carries `_unrefined`. The two refinement-status
#                              tags keep refined and Stage-A-only cohorts from commingling
#                              in metrics.py / extract_cost.py / plot scripts.
#   --max-iterations N         maximum Discriminator rounds (default: 3).
#                              Only matters when --with-refine is set.
#   --no-embed                 skip Embed. Default: Embed runs.
#   --noref                    resource-unconstrained variant. Swaps the populator
#                              prompt + blueprint template (and, with --with-refine,
#                              the discriminate/generate prompts) so the agent MAY
#                              but need not cite literature. MODEL_SLUG carries
#                              `_noref_unrefined` (Stage A only) or `_noref_refine`
#                              (with --with-refine; sibling snapshot stays in the
#                              `_noref_unrefined` cohort, mirroring the constrained
#                              `_unrefined` philosophy).
#   -h, --help                 print usage and exit.
#
# Notes: Web search and shell network access are always on for both
# stages (not configurable). Stage B's Discriminator uses web access
# for general biology lookups only.
# Resolved knobs are recorded in <PARENT_DIR>/sessions/<session_id>.log.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Canonical pipeline sources: prefer the vendored copies under webapp/vendor/
# (standalone), falling back to the parent repo when running in-tree.
VENDOR_ROOT="${SCRIPT_DIR}/vendor"
if [ -d "${VENDOR_ROOT}" ]; then SRC_ROOT="${VENDOR_ROOT}"; else SRC_ROOT="${REPO_ROOT}"; fi

SOURCE_BLUEPRINT="${SRC_ROOT}/scripts/agentic_embedding/embedder_blueprint.py"
SOURCE_BLUEPRINT_NOREF="${SRC_ROOT}/scripts/agentic_embedding/embedder_blueprint_noref.py"
SOURCE_PROMPT_EDIT="${SRC_ROOT}/prompt_templates/agent_prompt_blueprint_edit.md"
SOURCE_PROMPT_EDIT_NOREF="${SRC_ROOT}/prompt_templates/agent_prompt_blueprint_edit_noref.md"
SOURCE_PROMPT_DISCRIMINATE="${SRC_ROOT}/prompt_templates/agent_prompt_blueprint_discriminate.md"
SOURCE_PROMPT_DISCRIMINATE_NOREF="${SRC_ROOT}/prompt_templates/agent_prompt_blueprint_discriminate_noref.md"
SOURCE_PROMPT_GENERATE="${SRC_ROOT}/prompt_templates/agent_prompt_blueprint_generate.md"
SOURCE_PROMPT_GENERATE_NOREF="${SRC_ROOT}/prompt_templates/agent_prompt_blueprint_generate_noref.md"
SOURCE_EXTRACT_SCRIPT="${SRC_ROOT}/scripts/agentic_embedding/extract_run.py"
WORK_ROOT="${REPO_ROOT}/work"

DEFAULT_MAX_REFINE_ITERATIONS=3
DATASET_CHOICES=(mouse_brain human_pancreas mouse_gastrulation mouse_cortex human_immune webapp)
REASONING_CHOICES=(low medium high xhigh)

usage() {
  cat <<'EOF'
Usage: agent_edit_blueprint.sh --dataset NAME --model SLUG --reasoning-effort LEVEL
                               [--with-refine] [--max-iterations N]
                               [--noref]
                               [--no-embed] [-h|--help]

Default behaviour: Stage A only, resource-constrained. Opt into the
Stage B refinement loop with --with-refine. Opt into the
resource-unconstrained variant with --noref. The two flags compose.

Required:
  --dataset NAME             mouse_brain | human_pancreas | mouse_gastrulation | mouse_cortex | human_immune.
  --model SLUG               Codex model slug, e.g. gpt-5.5, gpt-5.6-luna.
  --reasoning-effort LEVEL   One of: low, medium, high, xhigh.

Optional:
  --with-refine              Opt INTO Stage B (Discriminator/Generator refinement loop).
                             When set, MODEL_SLUG carries `_refine`; when unset,
                             `_unrefined`. Combines with --noref → `_noref_refine`
                             (sibling Stage A snapshot lands in the `_noref_unrefined` cohort).
  --max-iterations N         Maximum Discriminator rounds when --with-refine is set
                             (default: 3). The loop exits early on rating "fine".
  --noref                    Resource-unconstrained variant: swaps to
                             embedder_blueprint_noref.py + agent_prompt_blueprint_edit_noref.md
                             (and the *_noref discriminate/generate prompts when
                             --with-refine is also set). MODEL_SLUG carries
                             `_noref_unrefined` (or `_noref_refine` with --with-refine).
  --no-embed                 Skip Embed. Default: Embed runs.
  -h, --help                 Print this message and exit.

webapp-vendored overrides (this copy only):
  --scenario-file PATH       Use this scenario_config.jsonc instead of the
                             canonical configs/scenario_config_<dataset>.jsonc.
                             When set, --dataset is optional (cohort tag `webapp`).
  --input PATH               Embed against this h5ad instead of
                             data/preprocessed/<dataset>_pp.h5ad.
  --work-root PATH           Place workdirs under this root instead of <repo>/work.
  --non-interactive          Don't prompt for the parallel-run count; use --n-parallel.
  --n-parallel N             Parallel-run count when --non-interactive (default 1).
  --print-config             Print resolved config JSON and exit 0 (no execution).

Output layout:
  Workdirs land under  work/<model_slug>_{unrefined|refine|noref_unrefined|noref_refine}_<reasoning>_<dataset>/<YYMMDD>_blueprint_runNN/
  All session transcripts + embed.log live under <workdir>/logs/.
  Session log lives at work/<...>/sessions/<session_id>.log
  Embed parquets at  outs/embeddings/<dataset>/blueprint_<model_slug>_{unrefined|refine|noref_unrefined|noref_refine}_<reasoning>_<dataset>_seed<N>.parquet

Web search is forced on for both stages (not configurable).
The number of parallel runs is asked interactively after parsing.
EOF
}

DATASET=""
MODEL=""
REASONING_EFFORT=""
WITH_REFINE="off"
MAX_REFINE_ITERATIONS="${DEFAULT_MAX_REFINE_ITERATIONS}"
EMBED="on"
NOREF="off"

# webapp-vendored overrides (see header). All empty/off by default so the
# script behaves identically to the original 6-dataset path unless a caller
# supplies them.
SCENARIO_FILE=""
BLUEPRINT_INPUT_OVERRIDE=""
WORK_ROOT_OVERRIDE=""
NON_INTERACTIVE="off"
N_PARALLEL_FLAG=""
PRINT_CONFIG="off"

while (( $# > 0 )); do
  case "$1" in
    --dataset)
      [[ $# -ge 2 ]] || { printf -- '--dataset requires a value\n' >&2; exit 1; }
      DATASET="$2"
      shift 2
      ;;
    --dataset=*)
      DATASET="${1#--dataset=}"
      shift
      ;;
    --model)
      [[ $# -ge 2 ]] || { printf -- '--model requires a value\n' >&2; exit 1; }
      MODEL="$2"
      shift 2
      ;;
    --model=*)
      MODEL="${1#--model=}"
      shift
      ;;
    --reasoning-effort)
      [[ $# -ge 2 ]] || { printf -- '--reasoning-effort requires a value\n' >&2; exit 1; }
      REASONING_EFFORT="$2"
      shift 2
      ;;
    --reasoning-effort=*)
      REASONING_EFFORT="${1#--reasoning-effort=}"
      shift
      ;;
    --with-refine)
      WITH_REFINE="on"
      shift
      ;;
    --no-refine)
      WITH_REFINE="off"
      shift
      ;;
    --max-iterations)
      [[ $# -ge 2 ]] || { printf -- '--max-iterations requires a value\n' >&2; exit 1; }
      MAX_REFINE_ITERATIONS="$2"
      shift 2
      ;;
    --max-iterations=*)
      MAX_REFINE_ITERATIONS="${1#--max-iterations=}"
      shift
      ;;
    --no-embed)
      EMBED="off"
      shift
      ;;
    --embed)
      EMBED="on"
      shift
      ;;
    --noref)
      NOREF="on"
      shift
      ;;
    --scenario-file)
      [[ $# -ge 2 ]] || { printf -- '--scenario-file requires a value\n' >&2; exit 1; }
      SCENARIO_FILE="$2"
      shift 2
      ;;
    --scenario-file=*)
      SCENARIO_FILE="${1#--scenario-file=}"
      shift
      ;;
    --input)
      [[ $# -ge 2 ]] || { printf -- '--input requires a value\n' >&2; exit 1; }
      BLUEPRINT_INPUT_OVERRIDE="$2"
      shift 2
      ;;
    --input=*)
      BLUEPRINT_INPUT_OVERRIDE="${1#--input=}"
      shift
      ;;
    --work-root)
      [[ $# -ge 2 ]] || { printf -- '--work-root requires a value\n' >&2; exit 1; }
      WORK_ROOT_OVERRIDE="$2"
      shift 2
      ;;
    --work-root=*)
      WORK_ROOT_OVERRIDE="${1#--work-root=}"
      shift
      ;;
    --non-interactive)
      NON_INTERACTIVE="on"
      shift
      ;;
    --n-parallel)
      [[ $# -ge 2 ]] || { printf -- '--n-parallel requires a value\n' >&2; exit 1; }
      N_PARALLEL_FLAG="$2"
      shift 2
      ;;
    --n-parallel=*)
      N_PARALLEL_FLAG="${1#--n-parallel=}"
      shift
      ;;
    --print-config)
      PRINT_CONFIG="on"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Unknown argument: %s\n' "$1" >&2
      usage >&2
      exit 1
      ;;
  esac
done

# webapp variant: when a caller supplies its own scenario file, the fixed
# dataset name is not required — default the cohort tag to `webapp` (a valid
# DATASET_CHOICES entry) so the rest of the pipeline still resolves.
[[ -z "${DATASET}" && -n "${SCENARIO_FILE}" ]] && DATASET="webapp"

missing=()
[[ -z "${DATASET}" ]]          && missing+=(--dataset)
[[ -z "${MODEL}" ]]            && missing+=(--model)
[[ -z "${REASONING_EFFORT}" ]] && missing+=(--reasoning-effort)
if (( ${#missing[@]} > 0 )); then
  printf 'Missing required flag(s): %s\n' "${missing[*]}" >&2
  usage >&2
  exit 1
fi

valid_choice() {
  local needle="$1"; shift
  local item
  for item in "$@"; do [[ "${item}" == "${needle}" ]] && return 0; done
  return 1
}

if ! valid_choice "${DATASET}" "${DATASET_CHOICES[@]}"; then
  printf -- '--dataset must be one of: %s (got: %s)\n' \
    "${DATASET_CHOICES[*]}" "${DATASET}" >&2
  exit 1
fi
if ! valid_choice "${REASONING_EFFORT}" "${REASONING_CHOICES[@]}"; then
  printf -- '--reasoning-effort must be one of: %s (got: %s)\n' \
    "${REASONING_CHOICES[*]}" "${REASONING_EFFORT}" >&2
  exit 1
fi

if ! [[ "${MAX_REFINE_ITERATIONS}" =~ ^[0-9]+$ ]] || (( MAX_REFINE_ITERATIONS < 1 )); then
  printf -- '--max-iterations must be a positive integer (got: %s)\n' \
    "${MAX_REFINE_ITERATIONS}" >&2
  exit 1
fi

# Swap the populator prompt + blueprint template under --noref so the
# agent edits the resource-free variant. When --with-refine is
# additionally set, swap the Discriminator/Generator prompts too so the
# whole loop runs under the no-ref contract.
if [[ "${NOREF}" == "on" ]]; then
  SOURCE_BLUEPRINT="${SOURCE_BLUEPRINT_NOREF}"
  SOURCE_PROMPT_EDIT="${SOURCE_PROMPT_EDIT_NOREF}"
  if [[ "${WITH_REFINE}" == "on" ]]; then
    SOURCE_PROMPT_DISCRIMINATE="${SOURCE_PROMPT_DISCRIMINATE_NOREF}"
    SOURCE_PROMPT_GENERATE="${SOURCE_PROMPT_GENERATE_NOREF}"
  fi
fi

BLUEPRINT_INPUT="${REPO_ROOT}/data/preprocessed/${DATASET}_pp.h5ad"

SOURCE_SCENARIO="${REPO_ROOT}/configs/scenario_config_${DATASET}.jsonc"

# webapp-vendored overrides: let a caller point at its own scenario file and
# input h5ad instead of the canonical per-dataset locations. Applied right
# after the defaults are assigned so every downstream use (required-file
# checks, workdir copy, Embed embed) sees the override.
[[ -n "${SCENARIO_FILE}" ]]           && SOURCE_SCENARIO="${SCENARIO_FILE}"
[[ -n "${BLUEPRINT_INPUT_OVERRIDE}" ]] && BLUEPRINT_INPUT="${BLUEPRINT_INPUT_OVERRIDE}"

# model_slug: gpt-5.5 -> gpt_5_5, gpt-5.6-luna -> gpt_5_6_luna.
MODEL_SLUG="${MODEL//./_}"
MODEL_SLUG="${MODEL_SLUG//-/_}"
# Always append a cohort-identity tag so every cohort surfaces with a
# stable identity in metrics.py, costs.csv, and the plot scripts (mirrors
# the `_userctx` marker on the claude variant). The tag is
# `[_noref]_{unrefined|refine}` — the grounding-mode word (present only
# under --noref) followed by the refinement-state word, so every cohort
# spells out both. Four combinations:
#   default               → `_unrefined`        (Stage A only, constrained)
#   --with-refine         → `_refine`           (sibling `_unrefined` snapshot)
#   --noref               → `_noref_unrefined`  (Stage A only, no resources)
#   --noref --with-refine → `_noref_refine`     (sibling `_noref_unrefined` snapshot)
# Pure-Stage-A runs land at the exact same artifact layout as the Stage A
# snapshot half of the matching --with-refine run (constrained `_unrefined`
# == refine snapshot; noref `_noref_unrefined` == noref-refine snapshot).
# This keeps the refine-from-existing orchestrator
# (analysis/agent_refine_blueprint.sh) a no-rename additive operation, and
# lets aggregators pair Stage A vs Stage A+B per grounding mode.
if [[ "${NOREF}" == "on" ]]; then
  if [[ "${WITH_REFINE}" == "on" ]]; then
    MODEL_SLUG="${MODEL_SLUG}_noref_refine"
  else
    MODEL_SLUG="${MODEL_SLUG}_noref_unrefined"
  fi
elif [[ "${WITH_REFINE}" == "on" ]]; then
  MODEL_SLUG="${MODEL_SLUG}_refine"
else
  MODEL_SLUG="${MODEL_SLUG}_unrefined"
fi
# Sibling slug for the Stage A snapshot of a --with-refine run: swap the
# trailing `_refine` for `_unrefined`, preserving any `_noref` grounding
# word (`_noref_refine` → `_noref_unrefined`, `_refine` → `_unrefined`).
# For Stage-A-only runs MODEL_SLUG already ends in `_unrefined`, so this is
# a no-op and UNREFINED_SLUG == MODEL_SLUG (no separate snapshot is emitted;
# the Embed emit_unrefined guard stays 0). Kept defined unconditionally to
# satisfy set -u where unrefined_stem is evaluated.
UNREFINED_SLUG="${MODEL_SLUG%_refine}"
[[ "${UNREFINED_SLUG}" == "${MODEL_SLUG}" ]] || UNREFINED_SLUG="${UNREFINED_SLUG}_unrefined"

# webapp-vendored override: redirect the work-root before PARENT_DIR /
# SESSION_LOG_DIR are derived from it.
[[ -n "${WORK_ROOT_OVERRIDE}" ]] && WORK_ROOT="${WORK_ROOT_OVERRIDE}"

PARENT_NAME="${MODEL_SLUG}_${REASONING_EFFORT}_${DATASET}"
PARENT_DIR="${WORK_ROOT}/${PARENT_NAME}"
SESSION_LOG_DIR="${PARENT_DIR}/sessions"

CODEX_ARGS=(
  --skip-git-repo-check
  -s workspace-write
  -c "tools.web_search=true"
  -c "sandbox_workspace_write.network_access=true"
  -m "${MODEL}"
  -c "model_reasoning_effort=\"${REASONING_EFFORT}\""
)
CODEX_ARGS+=(--disable fast_mode)
# Isolation: hide the user's Codex plugins, apps, memories and skills from the agent.
CODEX_ARGS+=(--disable plugins --disable apps --disable memories -c "skills.include_instructions=false")

required_files=(
  "${SOURCE_BLUEPRINT}"
  "${SOURCE_PROMPT_EDIT}"
  "${SOURCE_SCENARIO}"
  "${SOURCE_EXTRACT_SCRIPT}"
)
if [[ "${WITH_REFINE}" == "on" ]]; then
  required_files+=("${SOURCE_PROMPT_DISCRIMINATE}" "${SOURCE_PROMPT_GENERATE}")
fi
for required_file in "${required_files[@]}"; do
  if [[ ! -f "${required_file}" ]]; then
    printf 'Required file not found: %s\n' "${required_file}" >&2
    exit 1
  fi
done

# webapp-vendored: resolve the parallel-run count from a flag when running
# non-interactively (the web app drives the script headless); otherwise fall
# back to the original interactive prompt.
if [[ "${NON_INTERACTIVE}" == "on" ]]; then
  N_PARALLEL="${N_PARALLEL_FLAG:-1}"
else
  printf 'How many parallel runs? [default: 1]: '
  read -r N_PARALLEL_INPUT
  N_PARALLEL="${N_PARALLEL_INPUT:-1}"
fi
if ! [[ "${N_PARALLEL}" =~ ^[0-9]+$ ]] || (( N_PARALLEL < 1 )); then
  printf 'Invalid input: must be a positive integer (got: %s)\n' "${N_PARALLEL}" >&2
  exit 1
fi

# webapp-vendored: dump the resolved config and exit before any workdir
# creation or codex execution. All overrides + N_PARALLEL are resolved above.
if [[ "${PRINT_CONFIG}" == "on" ]]; then
  printf '{"scenario_file":"%s","blueprint_input":"%s","work_root":"%s","n_parallel":%s,"dataset":"%s","model":"%s"}\n' \
    "${SOURCE_SCENARIO}" "${BLUEPRINT_INPUT}" "${WORK_ROOT}" "${N_PARALLEL}" "${DATASET}" "${MODEL}"
  exit 0
fi

mkdir -p "${PARENT_DIR}"
mkdir -p "${SESSION_LOG_DIR}"

RUN_DATE="$(date '+%y%m%d')"
SESSION_ID="${RUN_DATE}_$(date '+%H%M%S')_$$"
SESSION_LOG="${SESSION_LOG_DIR}/${SESSION_ID}.log"
NEXT_RUN=1

shopt -s nullglob
# Seeds (= runNN) stay unique per cohort across dates.
for existing_dir in "${PARENT_DIR}/"*_blueprint_run*; do
  existing_name="$(basename "${existing_dir}")"
  if [[ "${existing_name}" =~ ^[0-9]{6}_blueprint_run([0-9]{2})$ ]]; then
    run_number=$((10#${BASH_REMATCH[1]}))
    if (( run_number >= NEXT_RUN )); then
      NEXT_RUN=$((run_number + 1))
    fi
  fi
done
shopt -u nullglob

LAST_RUN=$((NEXT_RUN + N_PARALLEL - 1))
if (( LAST_RUN > 99 )); then
  printf 'Too many runs for %s in %s; would allocate run %d which exceeds 99.\n' \
    "${RUN_DATE}" "${PARENT_NAME}" "${LAST_RUN}" >&2
  exit 1
fi

EDIT_PROMPT="$(<"${SOURCE_PROMPT_EDIT}")"
if [[ "${WITH_REFINE}" == "on" ]]; then
  DISCRIMINATE_PROMPT="$(<"${SOURCE_PROMPT_DISCRIMINATE}")"
  GENERATE_PROMPT="$(<"${SOURCE_PROMPT_GENERATE}")"
fi

WORKDIRS=()
for ((i = 0; i < N_PARALLEL; i++)); do
  run_num=$((NEXT_RUN + i))
  workdir="${PARENT_DIR}/${RUN_DATE}_blueprint_run$(printf '%02d' "${run_num}")"
  mkdir -p "${workdir}"
  mkdir -p "${workdir}/logs"
  cp "${SOURCE_BLUEPRINT}"   "${workdir}/embedder_blueprint.py"
  cp "${SOURCE_SCENARIO}"    "${workdir}/scenario_config.jsonc"
  cp "${SOURCE_PROMPT_EDIT}" "${workdir}/agent_prompt_blueprint_edit.md"
  WORKDIRS+=("${workdir}")
  printf 'Prepared working directory: %s\n' "${workdir}"
done

{
  printf 'session_id: %s\n' "${SESSION_ID}"
  printf 'started_at: %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')"
  printf 'n_parallel: %d\n' "${N_PARALLEL}"
  printf 'dataset: %s\n'                "${DATASET}"
  printf 'codex_model: %s\n'            "${MODEL}"
  printf 'codex_reasoning_effort: %s\n' "${REASONING_EFFORT}"
  printf 'codex_fast_mode: off\n'
  printf 'with_refine: %s\n'            "${WITH_REFINE}"
  printf 'max_refine_iterations: %d\n'  "${MAX_REFINE_ITERATIONS}"
  printf 'embed: %s\n'                  "${EMBED}"
  printf 'noref: %s\n'                  "${NOREF}"
  printf 'codex_web_search: always-on\n'
  printf 'codex_network_access: always-on\n'
  printf 'codex_sandbox: workspace-write (isolated: no plugins/apps/memories/skills)\n'
  printf 'parent_dir: %s\n' "${PARENT_NAME}"
  printf 'workdirs:\n'
  for w in "${WORKDIRS[@]}"; do
    printf '  - %s\n' "$(basename "${w}")"
  done
} > "${SESSION_LOG}"

for workdir in "${WORKDIRS[@]}"; do
  {
    printf 'session_id: %s\n' "${SESSION_ID}"
    printf 'session_log: %s\n' "${SESSION_LOG}"
    printf 'siblings:\n'
    for w in "${WORKDIRS[@]}"; do
      [[ "${w}" == "${workdir}" ]] && continue
      printf '  - %s\n' "$(basename "${w}")"
    done
  } > "${workdir}/session.txt"
done

printf 'Session %s; shared log: %s\n' "${SESSION_ID}" "${SESSION_LOG}"

# ============================================================
# Stage A — populator, parallel across workdirs
# ============================================================

printf 'Stage A (populator) for %d workdir(s).\n' "${N_PARALLEL}"
pids=()
for workdir in "${WORKDIRS[@]}"; do
  (
    cd "${workdir}"
    log="${workdir}/logs/codex_edit.log"
    printf '[%s] [%s] stage_a starting\n' "$(date '+%H:%M:%S')" "$(basename "${workdir}")" >> "${SESSION_LOG}"
    set +e
    _t0=$(date +%s)
    # webapp-vendored: when WEBAPP_MOCK_CODEX=1, skip the real codex CLI and
    # populate the blueprint from a real archived one (token-free e2e).
    if [[ "${WEBAPP_MOCK_CODEX:-}" == "1" ]]; then
      python "${SCRIPT_DIR}/_mock_codex.py" >"${log}" 2>&1
      rc=$?
    else
      codex exec ${CODEX_ARGS[@]+"${CODEX_ARGS[@]}"} "${EDIT_PROMPT}" </dev/null >"${log}" 2>&1
      rc=$?
    fi
    _t1=$(date +%s)
    set -e
    # A clean exit that left the template untouched is a failure.
    (( rc == 0 )) && cmp -s "${SOURCE_BLUEPRINT}" embedder_blueprint.py && rc=3
    echo "$((_t1 - _t0))" > "${workdir}/.stage_a_secs"
    printf '[%s] [%s] stage_a finished (exit=%d)\n' "$(date '+%H:%M:%S')" "$(basename "${workdir}")" "${rc}" >> "${SESSION_LOG}"
    echo "${rc}" > "${workdir}/.stage_a_rc"
    exit "${rc}"
  ) &
  pids+=("$!")
done

wait

# Determine surviving workdirs (Stage A exited 0)
SURVIVORS=()
A_FAIL=()
for workdir in "${WORKDIRS[@]}"; do
  rc=$(cat "${workdir}/.stage_a_rc")
  if [[ "${rc}" == "0" ]]; then
    SURVIVORS+=("${workdir}")
  else
    A_FAIL+=("${workdir}")
    printf '[%s] [%s] A-FAIL (rc=%s); skipping later steps\n' \
        "$(date '+%H:%M:%S')" "$(basename "${workdir}")" "${rc}" >> "${SESSION_LOG}"
  fi
done

if (( ${#SURVIVORS[@]} == 0 )); then
  printf 'All Stage A runs failed; no later steps to run.\n' >&2
fi

# ============================================================
# artifact staging: artifact copy into surviving workdirs only.
# - --with-refine → Discriminator + Generator prompts, reports/ dir.
# ============================================================

if (( ${#SURVIVORS[@]} > 0 )); then
  if [[ "${WITH_REFINE}" == "on" ]]; then
    for workdir in "${SURVIVORS[@]}"; do
      mkdir -p "${workdir}/reports"
      cp "${SOURCE_PROMPT_DISCRIMINATE}" "${workdir}/agent_prompt_blueprint_discriminate.md"
      cp "${SOURCE_PROMPT_GENERATE}"     "${workdir}/agent_prompt_blueprint_generate.md"
    done
  fi
fi

# ============================================================
# Extract (Stage A): early extraction — populator-only summary, overwritten
# by Extract (Stage B) once Stage B refinement edits the blueprint.
# Runs on every workdir (including A-FAIL ones; failures tolerated).
# ============================================================

printf 'Extract (Stage A): early extraction (after Stage A).\n'
for workdir in "${WORKDIRS[@]}"; do
  cp "${SOURCE_EXTRACT_SCRIPT}" "${workdir}/"
  python "${workdir}/extract_run.py" "${workdir}" || true
done

# ============================================================
# Stage B (Discriminator/Generator refinement loop).
# Skipped unless --with-refine.
#
# Each iteration runs in an isolated context: prior reports, prior
# refinement logs, the Stage A transcript (logs/codex_edit.log), and
# the embedder_blueprint_unrefined.py snapshot are moved into
# .refine_hidden/ so the freshly spawned agent sees a clean slate.
# After the loop completes, everything is restored with round-tagged
# filenames under reports/ and logs/ so downstream consumers
# (extract_cost.py, plot scripts) read the conventional layout.
# regression_guard.py protects the Generator step: unauthorized
# removals are reverted + the Generator re-spawned once with a
# corrective addendum.
# ============================================================

if [[ "${WITH_REFINE}" == "on" ]] && (( ${#SURVIVORS[@]} > 0 )); then
  printf 'Stage B (refinement, isolated per-iter context) for %d workdir(s); max %d iter(s).\n' \
      "${#SURVIVORS[@]}" "${MAX_REFINE_ITERATIONS}"
  REGRESSION_GUARD="${SRC_ROOT}/scripts/agentic_embedding/regression_guard.py"
  pids=()
  for workdir in "${SURVIVORS[@]}"; do
    (
      cd "${workdir}"
      cp embedder_blueprint.py embedder_blueprint_unrefined.py
      [[ -f summary.json ]] && cp summary.json summary_unrefined.json
      mkdir -p .refine_hidden/snapshot .refine_hidden/stage_a reports logs
      # Park the Stage A codex transcript so refinement agents never read it.
      [[ -f logs/codex_edit.log ]] && mv logs/codex_edit.log .refine_hidden/stage_a/codex_edit.log

      iter=1
      iters_completed=0
      final_rating="?"
      loop_secs=0
      status="FAIL"
      while (( iter <= MAX_REFINE_ITERATIONS )); do
        nn=$(printf '%02d' "${iter}")
        iter_dir=".refine_hidden/iter${nn}"
        mkdir -p "${iter_dir}"

        # ---- Discriminator: hide every prior-round artifact ----
        [[ -f embedder_blueprint_unrefined.py ]] && mv embedder_blueprint_unrefined.py .refine_hidden/snapshot/
        [[ -f summary_unrefined.json          ]] && mv summary_unrefined.json          .refine_hidden/snapshot/
        rm -f reports/report.md reports/report.json

        disc_log="${iter_dir}/discriminate.log"
        printf '[%s] [%s] stage_b iter=%s discriminator starting (isolated)\n' \
            "$(date '+%H:%M:%S')" "$(basename "${workdir}")" "${nn}" >> "${SESSION_LOG}"
        set +e
        _t0=$(date +%s)
        codex exec ${CODEX_ARGS[@]+"${CODEX_ARGS[@]}"} "${DISCRIMINATE_PROMPT}" </dev/null \
          >"${disc_log}" 2>&1
        disc_rc=$?
        _t1=$(date +%s)
        set -e
        loop_secs=$((loop_secs + _t1 - _t0))

        rating=$(python - <<'PYEOF' 2>/dev/null
import json
try:
    r = json.load(open("reports/report.json"))["rating"]
    if r not in ("fine", "minor", "substantial"):
        r = "substantial"
except Exception:
    r = "substantial"
print(r)
PYEOF
)
        rating="${rating:-substantial}"
        final_rating="${rating}"
        iters_completed=$((iters_completed + 1))
        printf '[%s] [%s] stage_b iter=%s discriminator finished (rc=%d) rating=%s\n' \
            "$(date '+%H:%M:%S')" "$(basename "${workdir}")" \
            "${nn}" "${disc_rc}" "${rating}" >> "${SESSION_LOG}"

        # Stash the Discriminator's output for post-loop restoration.
        [[ -f reports/report.md   ]] && cp reports/report.md   "${iter_dir}/report.md"
        [[ -f reports/report.json ]] && cp reports/report.json "${iter_dir}/report.json"

        if [[ "${rating}" == "fine" ]]; then
          status="OK"
          break
        fi
        if (( iter == MAX_REFINE_ITERATIONS )); then
          status="MAX-ITER"
          break
        fi

        # ---- Generator: restore the Stage A audit-trail snapshot ----
        [[ -f .refine_hidden/snapshot/embedder_blueprint_unrefined.py ]] && \
          mv .refine_hidden/snapshot/embedder_blueprint_unrefined.py .
        [[ -f .refine_hidden/snapshot/summary_unrefined.json ]] && \
          mv .refine_hidden/snapshot/summary_unrefined.json .
        cp embedder_blueprint.py "${iter_dir}/blueprint_before_gen.py"

        gen_attempt=1
        gen_max_attempts=2
        gen_prompt_used="${GENERATE_PROMPT}"
        while (( gen_attempt <= gen_max_attempts )); do
          gen_log="${iter_dir}/generate_attempt$(printf '%02d' "${gen_attempt}").log"
          printf '[%s] [%s] stage_b iter=%s generator attempt=%d starting\n' \
              "$(date '+%H:%M:%S')" "$(basename "${workdir}")" \
              "${nn}" "${gen_attempt}" >> "${SESSION_LOG}"
          set +e
          _t2=$(date +%s)
          codex exec ${CODEX_ARGS[@]+"${CODEX_ARGS[@]}"} "${gen_prompt_used}" </dev/null \
            >"${gen_log}" 2>&1
          gen_rc=$?
          _t3=$(date +%s)
          set -e
          loop_secs=$((loop_secs + _t3 - _t2))

          set +e
          guard_out=$(python "${REGRESSION_GUARD}" \
              "${workdir}" \
              "${iter_dir}/blueprint_before_gen.py" \
              embedder_blueprint.py 2>&1)
          guard_rc=$?
          set -e
          printf '=== attempt %02d ===\n%s\n' "${gen_attempt}" "${guard_out}" \
            >> "${iter_dir}/regression_guard.log"

          if (( guard_rc == 0 )); then
            printf '[%s] [%s] stage_b iter=%s generator attempt=%d finished (rc=%d) regression=clean\n' \
                "$(date '+%H:%M:%S')" "$(basename "${workdir}")" \
                "${nn}" "${gen_attempt}" "${gen_rc}" >> "${SESSION_LOG}"
            break
          fi

          unauthorized=$(printf '%s\n' "${guard_out}" | sed -n 's/^  unauthorized removal: //p')
          printf '[%s] [%s] stage_b iter=%s generator attempt=%d REGRESSION-DETECTED; reverting\n' \
              "$(date '+%H:%M:%S')" "$(basename "${workdir}")" \
              "${nn}" "${gen_attempt}" >> "${SESSION_LOG}"
          cp "${iter_dir}/blueprint_before_gen.py" embedder_blueprint.py
          if (( gen_attempt == gen_max_attempts )); then
            printf '[%s] [%s] stage_b iter=%s REGRESSION-BLOCKED after %d attempts; continuing loop\n' \
                "$(date '+%H:%M:%S')" "$(basename "${workdir}")" \
                "${nn}" "${gen_attempt}" >> "${SESSION_LOG}"
            break
          fi
          gen_prompt_used="${GENERATE_PROMPT}

# CORRECTIVE ADDENDUM (your prior attempt was reverted)

The harness reverted your previous attempt because it removed these
keys without authorization from \`reports/report.json\`:

${unauthorized}

Re-edit \`embedder_blueprint.py\` to address the report's findings
WITHOUT removing those keys. If you need to alter one of them,
edit the call in place (swap markers, change DOIs, restructure
arguments) rather than deleting it. Removals are sanctioned ONLY
when the report's \`remove\` array lists the key explicitly."
          gen_attempt=$((gen_attempt + 1))
        done

        # Hide the snapshot again before the next iter's Discriminator.
        [[ -f embedder_blueprint_unrefined.py ]] && mv embedder_blueprint_unrefined.py .refine_hidden/snapshot/
        [[ -f summary_unrefined.json          ]] && mv summary_unrefined.json          .refine_hidden/snapshot/
        rm -f reports/report.md reports/report.json
        iter=$((iter + 1))
      done

      echo "${loop_secs}"       > "${workdir}/.stage_b_secs"
      echo "${iters_completed}" > "${workdir}/.stage_b_iters"
      echo "${final_rating}"    > "${workdir}/.stage_b_final_rating"
      case "${status}" in
        OK)        echo "OK ${iters_completed} ${final_rating}"        > "${workdir}/.stage_b_status" ;;
        MAX-ITER)  echo "MAX-ITER ${iters_completed} ${final_rating}"  > "${workdir}/.stage_b_status" ;;
        *)         echo "FAIL ${iters_completed}"                       > "${workdir}/.stage_b_status" ;;
      esac

      # ---- Post-loop restoration: bring everything back into the
      #      conventional tagged layout for downstream consumers. ----
      mkdir -p reports logs
      [[ -f .refine_hidden/stage_a/codex_edit.log ]] && \
        mv .refine_hidden/stage_a/codex_edit.log logs/codex_edit.log
      [[ -f .refine_hidden/snapshot/embedder_blueprint_unrefined.py ]] && \
        mv .refine_hidden/snapshot/embedder_blueprint_unrefined.py .
      [[ -f .refine_hidden/snapshot/summary_unrefined.json ]] && \
        mv .refine_hidden/snapshot/summary_unrefined.json .
      shopt -s nullglob
      for iter_dir in .refine_hidden/iter*; do
        [[ -d "${iter_dir}" ]] || continue
        nn=$(basename "${iter_dir}" | sed 's/^iter//')
        [[ -f "${iter_dir}/report.md"            ]] && cp "${iter_dir}/report.md"            "reports/report_iter${nn}.md"
        [[ -f "${iter_dir}/report.json"          ]] && cp "${iter_dir}/report.json"          "reports/report_iter${nn}.json"
        [[ -f "${iter_dir}/discriminate.log"     ]] && cp "${iter_dir}/discriminate.log"     "logs/codex_discriminate_${nn}.log"
        [[ -f "${iter_dir}/regression_guard.log" ]] && cp "${iter_dir}/regression_guard.log" "logs/regression_guard_${nn}.log"
        for gen_attempt_log in "${iter_dir}"/generate_attempt*.log; do
          [[ -f "${gen_attempt_log}" ]] || continue
          attempt=$(basename "${gen_attempt_log}" .log | sed 's/^generate_attempt//')
          if [[ "${attempt}" == "01" ]]; then
            cp "${gen_attempt_log}" "logs/codex_generate_${nn}.log"
          else
            cp "${gen_attempt_log}" "logs/codex_generate_${nn}_attempt${attempt}.log"
          fi
        done
      done
      shopt -u nullglob
      # Reconstruct INDEX.md from the per-iter JSONs.
      python - <<'PYEOF' >/dev/null 2>&1 || true
import json
from pathlib import Path
reports = sorted(Path("reports").glob("report_iter*.json"))
lines = []
for r in reports:
    nn = r.stem.replace("report_iter", "")
    try:
        d = json.loads(r.read_text())
    except Exception:
        d = {}
    rating = d.get("rating", "?")
    summary = (d.get("summary") or "").splitlines()[0][:200]
    lines.append(f"- iter {nn} — {rating} — {summary}")
if lines:
    Path("reports/INDEX.md").write_text("\n".join(lines) + "\n")
PYEOF

      printf '[%s] [%s] stage_b done (%s; iters=%d; final=%s)\n' \
          "$(date '+%H:%M:%S')" "$(basename "${workdir}")" \
          "${status}" "${iters_completed}" "${final_rating}" >> "${SESSION_LOG}"
    ) &
    pids+=("$!")
  done
  for pid in "${pids[@]}"; do
    wait "${pid}" || true
  done
fi

# ============================================================
# Extract (Stage B): post-run extraction — overwrites the Extract (Stage A) summary
# with one reflecting the final blueprint on disk (Stage B refinement
# edits). Skipped when --with-refine is off (the Extract (Stage A) summary is
# already final for a Stage-A-only run).
# ============================================================

if [[ "${WITH_REFINE}" == "on" ]]; then
  printf 'Extract (Stage B): extracting blueprint summaries.\n'
  for workdir in "${WORKDIRS[@]}"; do
    python "${workdir}/extract_run.py" "${workdir}" || true
  done
fi

# ============================================================
# Finalize: session log finalization
# ============================================================

{
  printf 'finished_at: %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')"
  printf 'results:\n'
  for workdir in "${WORKDIRS[@]}"; do
    name=$(basename "${workdir}")
    a_rc=$(cat "${workdir}/.stage_a_rc" 2>/dev/null || echo "?")
    if [[ "${a_rc}" != "0" ]]; then
      printf '  %-40s A-FAIL  attempts=0\n' "${name}"
      continue
    fi
    parts=()
    if [[ "${WITH_REFINE}" == "on" ]]; then
      b_status=$(cat "${workdir}/.stage_b_status" 2>/dev/null || echo "?")
      parts+=("refine=${b_status}")
    fi
    if (( ${#parts[@]} == 0 )); then
      printf '  %-40s STAGE-A-ONLY  attempts=0\n' "${name}"
    else
      printf '  %-40s %s\n' "${name}" "${parts[*]}"
    fi
  done
} >> "${SESSION_LOG}"

# ============================================================
# Embed: embedding — run the populated blueprint against the
# dataset's data input and write a parquet (+ report companion)
# directly under outs/embeddings/<dataset>/ as
# blueprint_<model_slug>[_refine]_<reasoning>_<dataset>_seed<N>.{parquet,_report.md}.
# Skipped entirely under --no-embed.
# ============================================================

embed_log_line() {
  printf '[%s] [%s] %s\n' "$(date '+%H:%M:%S')" "$1" "$2" >> "${SESSION_LOG}"
}

if [[ "${EMBED}" == "off" ]]; then
  printf 'Embed: disabled (--no-embed).\n'
  printf 'embed: disabled (--no-embed)\n' >> "${SESSION_LOG}"
else
  EMBED_OUT_DIR="${REPO_ROOT}/outs/embeddings/${DATASET}"

  if [[ ! -f "${BLUEPRINT_INPUT}" ]]; then
    printf '[warn] Embed skipped: input not found: %s\n' "${BLUEPRINT_INPUT}" >&2
    printf '[warn] Run analysis/run_data_pipeline.sh %s first to build it.\n' "${DATASET}" >&2
    printf 'embed: skipped (input not found: %s)\n' "${BLUEPRINT_INPUT#${REPO_ROOT}/}" >> "${SESSION_LOG}"
  else
    mkdir -p "${EMBED_OUT_DIR}"
    printf 'Embed: embedding accepted blueprints onto %s\n' "${BLUEPRINT_INPUT#${REPO_ROOT}/}"
    printf 'embed_started_at: %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" >> "${SESSION_LOG}"
    printf 'embed_results:\n' >> "${SESSION_LOG}"

    # Per-workdir status (SKIP-FILTERED|SKIP-EXISTS|OK|FAIL) and detail are
    # written to <workdir>/.embed_status and .embed_detail. The .embed_*
    # dotfile names are kept for backward compatibility with extract_cost.py
    # and other readers; only the human-facing phase number changed (5→6).
    pids=()

    for workdir in "${WORKDIRS[@]}"; do
      name=$(basename "${workdir}")
      a_rc=$(cat "${workdir}/.stage_a_rc" 2>/dev/null || echo "?")
      eligible=0
      if [[ "${a_rc}" == "0" ]]; then
        eligible=1
        if [[ "${WITH_REFINE}" == "on" ]]; then
          b_status_line=$(cat "${workdir}/.stage_b_status" 2>/dev/null || echo "")
          case "${b_status_line}" in
            OK*|MAX-ITER*) ;;  # both produce a usable refined blueprint
            *) eligible=0 ;;
          esac
        fi
      fi

      if (( eligible == 0 )); then
        echo "SKIP-FILTERED" > "${workdir}/.embed_status"
        : > "${workdir}/.embed_detail"
        embed_log_line "${name}" "embed SKIP-FILTERED"
        continue
      fi

      run_nn=""
      if [[ "${name}" =~ blueprint_run([0-9]{2})$ ]]; then
        run_nn="${BASH_REMATCH[1]}"
      else
        echo "FAIL" > "${workdir}/.embed_status"
        echo "cannot parse run number from ${name}" > "${workdir}/.embed_detail"
        embed_log_line "${name}" "embed FAIL (cannot parse run number)"
        continue
      fi

      seed_n=$((10#${run_nn}))
      stem="blueprint_${MODEL_SLUG}_${REASONING_EFFORT}_${DATASET}_seed${seed_n}"
      parquet="${EMBED_OUT_DIR}/${stem}.parquet"
      report="${EMBED_OUT_DIR}/${stem}_report.md"

      # Stage A sibling stem (only meaningful when --with-refine snapshotted the
      # pre-refinement blueprint as embedder_blueprint_unrefined.py). The Stage A
      # parquet sits next to the refined one and is recognised as a distinct
      # cohort by every downstream consumer (non-greedy regex absorbs _unrefined).
      unrefined_src="${workdir}/embedder_blueprint_unrefined.py"
      unrefined_stem="blueprint_${UNREFINED_SLUG}_${REASONING_EFFORT}_${DATASET}_seed${seed_n}"
      unrefined_parquet="${EMBED_OUT_DIR}/${unrefined_stem}.parquet"
      unrefined_report="${EMBED_OUT_DIR}/${unrefined_stem}_report.md"
      emit_unrefined=0
      if [[ -f "${unrefined_src}" && ! -e "${unrefined_parquet}" ]]; then
        emit_unrefined=1
      fi

      # Skip the workdir only when both the refined and the Stage A parquets are
      # already on disk (or no Stage A is expected). Lets a re-invocation backfill
      # only the missing one.
      if [[ -e "${parquet}" ]] && (( emit_unrefined == 0 )); then
        echo "SKIP-EXISTS" > "${workdir}/.embed_status"
        echo "${parquet#${REPO_ROOT}/}" > "${workdir}/.embed_detail"
        embed_log_line "${name}" "embed SKIP-EXISTS -> ${parquet#${REPO_ROOT}/}"
        continue
      fi

      (
        rc=0
        if [[ -e "${parquet}" ]]; then
          echo "SKIP-EXISTS" > "${workdir}/.embed_status"
          echo "${parquet#${REPO_ROOT}/}" > "${workdir}/.embed_detail"
        else
          log="${workdir}/logs/embed.log"
          set +e
          _t0=$(date +%s)
          python "${workdir}/embedder_blueprint.py" \
              --input  "${BLUEPRINT_INPUT}" \
              --output "${parquet}" \
              >"${log}" 2>&1
          rc=$?
          _t1=$(date +%s)
          set -e
          echo "$((_t1 - _t0))" > "${workdir}/.embed_secs"
          if (( rc == 0 )); then
            # Archive the populated blueprint into a nested layout so downstream
            # tools (umap_branch_interactive, blueprint_tree, audits) can locate
            # it by (dataset, model, seed) without globbing work/...
            #   outs/blueprints/<dataset>/<model_slug>_<reasoning>/embedder_blueprint_seed<N>.py
            blueprint_model_dir="${REPO_ROOT}/outs/blueprints/${DATASET}/${MODEL_SLUG}_${REASONING_EFFORT}"
            blueprint_dst="${blueprint_model_dir}/embedder_blueprint_seed${seed_n}.py"
            mkdir -p "${blueprint_model_dir}"
            cp "${workdir}/embedder_blueprint.py" "${blueprint_dst}"
            # Render the companion (metadata header + `## Embedding dimensions`
            # table) from the introspected blueprint via extract_run.py. This
            # also re-extracts summary.json from the final (refined) blueprint.
            python "${workdir}/extract_run.py" "${workdir}" \
                --report-out "${report}" \
                --dataset "${DATASET}" \
                --model "${MODEL}" --reasoning "${REASONING_EFFORT}" \
                --input "${BLUEPRINT_INPUT#${REPO_ROOT}/}" \
                --parquet "${parquet#${REPO_ROOT}/}" \
                >/dev/null 2>&1 \
              || printf '# Blueprint embedding\n\n(report generation failed)\n' > "${report}"
            # Cache the structured blueprint summary alongside the parquet so
            # downstream plot/metric scripts have a parse-free source for the
            # axis / module / resource structure.
            if [[ -s "${workdir}/summary.json" ]]; then
              cp "${workdir}/summary.json" "${EMBED_OUT_DIR}/${stem}_summary.json"
            fi
            echo "OK" > "${workdir}/.embed_status"
            echo "${parquet#${REPO_ROOT}/}" > "${workdir}/.embed_detail"
          else
            echo "FAIL" > "${workdir}/.embed_status"
            echo "rc=${rc} log=${workdir#${REPO_ROOT}/}/logs/embed.log" > "${workdir}/.embed_detail"
          fi
          echo "${rc}" > "${workdir}/.embed_rc"
        fi

        # Stage A sibling embed. Independent of the refined branch — embeds
        # embedder_blueprint_unrefined.py (snapshotted at refine-loop entry) into
        # a parquet under the `_unrefined` cohort name, archives the Stage A
        # blueprint, mirrors summary_unrefined.json, and renders a Stage-A-only
        # report. Failure here does not affect the subshell rc (refined and
        # Stage A are independent products of the same workdir).
        if (( emit_unrefined == 1 )); then
          unrefined_log="${workdir}/logs/embed_unrefined.log"
          set +e
          _u0=$(date +%s)
          python "${unrefined_src}" \
              --input  "${BLUEPRINT_INPUT}" \
              --output "${unrefined_parquet}" \
              >"${unrefined_log}" 2>&1
          urc=$?
          _u1=$(date +%s)
          set -e
          echo "$((_u1 - _u0))" > "${workdir}/.embed_unrefined_secs"
          if (( urc == 0 )); then
            unrefined_model_dir="${REPO_ROOT}/outs/blueprints/${DATASET}/${UNREFINED_SLUG}_${REASONING_EFFORT}"
            mkdir -p "${unrefined_model_dir}"
            cp "${unrefined_src}" "${unrefined_model_dir}/embedder_blueprint_seed${seed_n}.py"
            python "${workdir}/extract_run.py" "${workdir}" \
                --blueprint-source embedder_blueprint_unrefined.py \
                --report-out "${unrefined_report}" \
                --dataset "${DATASET}" \
                --model "${MODEL}" --reasoning "${REASONING_EFFORT}" \
                --input "${BLUEPRINT_INPUT#${REPO_ROOT}/}" \
                --parquet "${unrefined_parquet#${REPO_ROOT}/}" \
                >/dev/null 2>&1 \
              || printf '# Blueprint embedding (Stage A only)\n\n(report generation failed)\n' > "${unrefined_report}"
            if [[ -s "${workdir}/summary_unrefined.json" ]]; then
              cp "${workdir}/summary_unrefined.json" "${EMBED_OUT_DIR}/${unrefined_stem}_summary.json"
            fi
            echo "OK" > "${workdir}/.embed_unrefined_status"
            echo "${unrefined_parquet#${REPO_ROOT}/}" > "${workdir}/.embed_unrefined_detail"
          else
            echo "FAIL" > "${workdir}/.embed_unrefined_status"
            echo "rc=${urc} log=${workdir#${REPO_ROOT}/}/logs/embed_unrefined.log" > "${workdir}/.embed_unrefined_detail"
          fi
          echo "${urc}" > "${workdir}/.embed_unrefined_rc"
        fi

        exit "${rc}"
      ) &
      pids+=("$!")
    done

    for pid in "${pids[@]}"; do
      wait "${pid}" || true
    done

    # Per-workdir log lines now that all subshells finished. Reads the
    # status/detail files the subshells just wrote.
    for workdir in "${WORKDIRS[@]}"; do
      name=$(basename "${workdir}")
      status=$(cat "${workdir}/.embed_status" 2>/dev/null || echo "SKIP-FILTERED")
      detail=$(cat "${workdir}/.embed_detail" 2>/dev/null || true)
      case "${status}" in
        OK)            embed_log_line "${name}" "embed OK -> ${detail}" ;;
        FAIL)          embed_log_line "${name}" "embed FAIL ${detail}" ;;
        SKIP-EXISTS|SKIP-FILTERED) ;;  # already logged inline above
      esac
      if [[ -f "${workdir}/.embed_unrefined_status" ]]; then
        ustatus=$(cat "${workdir}/.embed_unrefined_status" 2>/dev/null || true)
        udetail=$(cat "${workdir}/.embed_unrefined_detail" 2>/dev/null || true)
        case "${ustatus}" in
          OK)   embed_log_line "${name}" "embed stage_a OK -> ${udetail}" ;;
          FAIL) embed_log_line "${name}" "embed stage_a FAIL ${udetail}" ;;
        esac
      fi
    done

    {
      for workdir in "${WORKDIRS[@]}"; do
        name=$(basename "${workdir}")
        status=$(cat "${workdir}/.embed_status" 2>/dev/null || echo "SKIP-FILTERED")
        detail=$(cat "${workdir}/.embed_detail" 2>/dev/null || true)
        combined="${detail}"
        if [[ -f "${workdir}/.embed_unrefined_status" ]]; then
          ustatus=$(cat "${workdir}/.embed_unrefined_status" 2>/dev/null || true)
          udetail=$(cat "${workdir}/.embed_unrefined_detail" 2>/dev/null || true)
          if [[ -n "${combined}" ]]; then
            combined="${combined} | stage_a=${ustatus} ${udetail}"
          else
            combined="stage_a=${ustatus} ${udetail}"
          fi
        fi
        if [[ -n "${combined}" ]]; then
          printf '  %-40s %-14s %s\n' "${name}" "${status}" "${combined}"
        else
          printf '  %-40s %-14s\n' "${name}" "${status}"
        fi
      done
      printf 'embed_finished_at: %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')"
    } >> "${SESSION_LOG}"
  fi
fi

# ============================================================
# Cost: cost extraction — per-workdir cost.json (token usage,
# per-phase wall-clock, on-disk footprint). Aggregated across
# cohorts by scripts/resource_costs/collect.py.
# ============================================================

EXTRACT_COST_SCRIPT="${SRC_ROOT}/scripts/agentic_embedding/extract_cost.py"
if [[ -f "${EXTRACT_COST_SCRIPT}" ]]; then
  printf 'Cost: extracting per-run cost summaries.\n'
  cost_out_dir="${REPO_ROOT}/outs/embeddings/${DATASET}"
  for workdir in "${WORKDIRS[@]}"; do
    python "${EXTRACT_COST_SCRIPT}" "${workdir}" || true
    # Mirror cost.json next to the embedding parquet under outs/embeddings/<dataset>/
    # so downstream cost/quality joins can read from outs/ without walking work/.
    # Matches the existing Embed summary.json mirror pattern.
    if [[ -s "${workdir}/cost.json" && "$(cat "${workdir}/.embed_status" 2>/dev/null)" == "OK" ]]; then
      detail=$(cat "${workdir}/.embed_detail")
      if [[ "${detail}" == *.parquet ]]; then
        stem=$(basename "${detail}" .parquet)
        mkdir -p "${cost_out_dir}"
        cp "${workdir}/cost.json" "${cost_out_dir}/${stem}_cost.json"
      fi
    fi
  done
fi

printf 'Session log: %s\n' "${SESSION_LOG}"
