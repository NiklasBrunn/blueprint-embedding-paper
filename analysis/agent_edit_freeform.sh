#!/usr/bin/env bash
# Free-form agent baseline (E3): the agent writes a single Python file
# `embed_freeform.py` from scratch — no blueprint API, no harness, no
# review step. Finalize runs the agent's script against the dataset to
# produce a parquet sibling to the blueprint embeddings, named
# freeform_<variant>_<model_slug>_<reasoning>_<dataset>_seed<N>.parquet.
#
# Two variants quantify *what* the constrained blueprint pipeline buys:
#
#   --variant constrained  (default)  uses
#       agent_prompt_freeform_constrained.md
#       Mirrors agent_prompt_blueprint_edit.md content rules verbatim
#       (IF>10 sources, DOI sanity check, no model training / pretrained
#       inference, notes-preference rule). Only the deliverable changes:
#       agent writes embed_freeform.py instead of editing build_blueprint().
#       Isolates the "constrained API + harness" effect cleanly.
#
#   --variant naive                  uses
#       agent_prompt_freeform_naive.md
#       Bare-bones prompt (~25 lines): "Read scenario, write script,
#       output parquet with cell_id + dim_*; no model training /
#       pretrained inference". No IF>10 / DOI / taxonomy rules. Worst-
#       case anchor showing what happens without prompt scaffolding.
#
#   --variant noref                  uses
#       agent_prompt_freeform_noref.md
#       Mirrors agent_prompt_blueprint_edit_noref.md: keeps the scenario
#       schema, script contract, input-format spec, dimensionality
#       guidance, and the no-training / no-pretrained-inference rule, but
#       drops the literature-grounding machinery (IF>10 sources, DOI /
#       Crossref check, resource kinds) — the agent MAY cite inline but
#       is not required to. The freeform counterpart of the blueprint
#       `_noref` arm; isolates "constrained API + harness" effect under
#       relaxed (resource-optional) content rules.
#
# Phases reduce to:
#
#   Stage A. Codex writes embed_freeform.py per the chosen
#              freeform prompt and the scenario_config.jsonc.
#   Extract (Stage B)  — Session log finalization.
#   Finalize  — Embedding. Runs each accepted workdir's embed_freeform.py
#              with --input/--output/--seed and validates the parquet
#              has cell_id + dim_*. Failures are logged with a reason.
#
# Usage / flags: run with --help (the usage() heredoc below is the single
# authoritative copy). Workdirs land under
# work/freeform_<variant>_<model_slug>_<reasoning>_<dataset>/<YYMMDD>_freeform_runNN/.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

PROMPT_TEMPLATE_DIR="${REPO_ROOT}/prompt_templates"
WORK_ROOT="${REPO_ROOT}/work"

DATASET_CHOICES=(mouse_brain human_pancreas)
REASONING_CHOICES=(low medium high xhigh)
VARIANT_CHOICES=(constrained naive noref)

usage() {
  cat <<'EOF'
Usage: agent_edit_freeform.sh --dataset NAME --model SLUG --reasoning-effort LEVEL
                              [--variant constrained|naive|noref]
                              [--no-embed] [-h|--help]

Required:
  --dataset NAME             mouse_brain | human_pancreas.
  --model SLUG               Codex model slug, e.g. gpt-5.5, gpt-5.4-mini.
  --reasoning-effort LEVEL   low | medium | high | xhigh.

Optional:
  --variant VARIANT          constrained | naive | noref. Default: constrained.
                             constrained: full content rules (IF>10, DOI,
                                          notes-preference, no pretrained
                                          models). Mirrors blueprint edit
                                          prompt minus the API.
                             naive:       bare-bones prompt; no content rules.
                             noref:       like constrained but drops the
                                          literature-grounding machinery
                                          (IF>10/DOI/Crossref/resources);
                                          mirrors the blueprint edit_noref
                                          prompt minus the API.
  --no-embed                 Skip Finalize. Default: Finalize runs.
  -h, --help                 Print this message and exit.

Output layout:
  Workdirs at  work/freeform_<variant>_<model_slug>_<reasoning>_<dataset>/<YYMMDD>_freeform_runNN/
  Session log  work/freeform_<variant>_.../sessions/<session_id>.log
  Finalize      outs/embeddings/<dataset>/freeform_<variant>_<model_slug>_<reasoning>_<dataset>_seed<N>.parquet
EOF
}

DATASET=""
MODEL=""
REASONING_EFFORT=""
VARIANT="constrained"
EMBED="on"

while (( $# > 0 )); do
  case "$1" in
    --dataset)            DATASET="$2"; shift 2 ;;
    --dataset=*)          DATASET="${1#--dataset=}"; shift ;;
    --model)              MODEL="$2"; shift 2 ;;
    --model=*)            MODEL="${1#--model=}"; shift ;;
    --reasoning-effort)   REASONING_EFFORT="$2"; shift 2 ;;
    --reasoning-effort=*) REASONING_EFFORT="${1#--reasoning-effort=}"; shift ;;
    --variant)            VARIANT="$2"; shift 2 ;;
    --variant=*)          VARIANT="${1#--variant=}"; shift ;;
    --no-embed)           EMBED="off"; shift ;;
    --embed)              EMBED="on"; shift ;;
    -h|--help)            usage; exit 0 ;;
    *) printf 'Unknown argument: %s\n' "$1" >&2; usage >&2; exit 1 ;;
  esac
done

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
valid_choice "${DATASET}"          "${DATASET_CHOICES[@]}"   || { printf -- '--dataset must be one of: %s\n' "${DATASET_CHOICES[*]}" >&2; exit 1; }
valid_choice "${REASONING_EFFORT}" "${REASONING_CHOICES[@]}" || { printf -- '--reasoning-effort must be one of: %s\n' "${REASONING_CHOICES[*]}" >&2; exit 1; }
valid_choice "${VARIANT}"          "${VARIANT_CHOICES[@]}"   || { printf -- '--variant must be one of: %s\n' "${VARIANT_CHOICES[*]}" >&2; exit 1; }

FREEFORM_INPUT="${REPO_ROOT}/data/preprocessed/${DATASET}_pp.h5ad"

SOURCE_PROMPT_FREEFORM="${PROMPT_TEMPLATE_DIR}/agent_prompt_freeform_${VARIANT}.md"
SOURCE_SCENARIO="${REPO_ROOT}/configs/scenario_config_${DATASET}.jsonc"

MODEL_SLUG="${MODEL//./_}"
MODEL_SLUG="${MODEL_SLUG//-/_}"

PARENT_NAME="freeform_${VARIANT}_${MODEL_SLUG}_${REASONING_EFFORT}_${DATASET}"
PARENT_DIR="${WORK_ROOT}/${PARENT_NAME}"
SESSION_LOG_DIR="${PARENT_DIR}/sessions"

CODEX_ARGS=(
  -c "tools.web_search=true"
  -c "sandbox_workspace_write.network_access=true"
  -m "${MODEL}"
  -c "model_reasoning_effort=\"${REASONING_EFFORT}\""
)
CODEX_ARGS+=(--disable fast_mode)

for required_file in "${SOURCE_PROMPT_FREEFORM}" "${SOURCE_SCENARIO}"; do
  if [[ ! -f "${required_file}" ]]; then
    printf 'Required file not found: %s\n' "${required_file}" >&2
    exit 1
  fi
done

printf 'How many parallel runs? [default: 1]: '
read -r N_PARALLEL_INPUT
N_PARALLEL="${N_PARALLEL_INPUT:-1}"
if ! [[ "${N_PARALLEL}" =~ ^[0-9]+$ ]] || (( N_PARALLEL < 1 )); then
  printf 'Invalid input: must be a positive integer (got: %s)\n' "${N_PARALLEL}" >&2
  exit 1
fi

mkdir -p "${PARENT_DIR}" "${SESSION_LOG_DIR}"

RUN_DATE="$(date '+%y%m%d')"
SESSION_ID="${RUN_DATE}_$(date '+%H%M%S')_$$"
SESSION_LOG="${SESSION_LOG_DIR}/${SESSION_ID}.log"
NEXT_RUN=1

shopt -s nullglob
for existing_dir in "${PARENT_DIR}/${RUN_DATE}_freeform_run"*; do
  existing_name="$(basename "${existing_dir}")"
  if [[ "${existing_name}" =~ ^${RUN_DATE}_freeform_run([0-9]{2})$ ]]; then
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

EDIT_PROMPT="$(<"${SOURCE_PROMPT_FREEFORM}")"

WORKDIRS=()
for ((i = 0; i < N_PARALLEL; i++)); do
  run_num=$((NEXT_RUN + i))
  workdir="${PARENT_DIR}/${RUN_DATE}_freeform_run$(printf '%02d' "${run_num}")"
  mkdir -p "${workdir}"
  mkdir -p "${workdir}/logs"
  cp "${SOURCE_SCENARIO}"          "${workdir}/scenario_config.jsonc"
  cp "${SOURCE_PROMPT_FREEFORM}"   "${workdir}/agent_prompt_freeform_${VARIANT}.md"
  WORKDIRS+=("${workdir}")
  printf 'Prepared working directory: %s\n' "${workdir}"
done

{
  printf 'session_id: %s\n' "${SESSION_ID}"
  printf 'started_at: %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')"
  printf 'n_parallel: %d\n' "${N_PARALLEL}"
  printf 'cohort: freeform_%s\n' "${VARIANT}"
  printf 'variant: %s\n'                "${VARIANT}"
  printf 'prompt_template: %s\n'        "$(basename "${SOURCE_PROMPT_FREEFORM}")"
  printf 'dataset: %s\n'                "${DATASET}"
  printf 'codex_model: %s\n'            "${MODEL}"
  printf 'codex_reasoning_effort: %s\n' "${REASONING_EFFORT}"
  printf 'codex_fast_mode: off\n'
  printf 'embed: %s\n'                  "${EMBED}"
  printf 'codex_web_search: always-on\n'
  printf 'codex_network_access: always-on\n'
  printf 'parent_dir: %s\n' "${PARENT_NAME}"
  printf 'workdirs:\n'
  for w in "${WORKDIRS[@]}"; do
    printf '  - %s\n' "$(basename "${w}")"
  done
} > "${SESSION_LOG}"

printf 'Session %s; shared log: %s\n' "${SESSION_ID}" "${SESSION_LOG}"

# ============================================================
# Stage A — agent writes embed_freeform.py
# ============================================================

printf 'Stage A (free-form populator) for %d workdir(s).\n' "${N_PARALLEL}"
pids=()
for workdir in "${WORKDIRS[@]}"; do
  (
    cd "${workdir}"
    log="${workdir}/logs/agent_edit.log"
    printf '[%s] [%s] stage_a starting\n' "$(date '+%H:%M:%S')" "$(basename "${workdir}")" >> "${SESSION_LOG}"
    set +e
    _t0=$(date +%s)
    codex exec ${CODEX_ARGS[@]+"${CODEX_ARGS[@]}"} "${EDIT_PROMPT}" >"${log}" 2>&1
    rc=$?
    _t1=$(date +%s)
    set -e
    echo "$((_t1 - _t0))" > "${workdir}/.stage_a_secs"
    printf '[%s] [%s] stage_a finished (exit=%d)\n' "$(date '+%H:%M:%S')" "$(basename "${workdir}")" "${rc}" >> "${SESSION_LOG}"
    echo "${rc}" > "${workdir}/.stage_a_rc"
    exit "${rc}"
  ) &
  pids+=("$!")
done
wait

SURVIVORS=()
for workdir in "${WORKDIRS[@]}"; do
  rc=$(cat "${workdir}/.stage_a_rc")
  if [[ "${rc}" == "0" && -s "${workdir}/embed_freeform.py" ]]; then
    SURVIVORS+=("${workdir}")
  else
    printf '[%s] [%s] A-FAIL (rc=%s; embed_freeform.py present=%s)\n' \
        "$(date '+%H:%M:%S')" "$(basename "${workdir}")" "${rc}" \
        "$([[ -s "${workdir}/embed_freeform.py" ]] && echo yes || echo no)" >> "${SESSION_LOG}"
  fi
done

# ============================================================
# Extract (Stage B): session log finalization
# ============================================================

{
  printf 'finished_at: %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')"
  printf 'results:\n'
  for workdir in "${WORKDIRS[@]}"; do
    name=$(basename "${workdir}")
    rc=$(cat "${workdir}/.stage_a_rc" 2>/dev/null || echo "?")
    has_script="no"
    [[ -s "${workdir}/embed_freeform.py" ]] && has_script="yes"
    printf '  %-40s rc=%s script=%s\n' "${name}" "${rc}" "${has_script}"
  done
} >> "${SESSION_LOG}"

# ============================================================
# Finalize: embedding — run each surviving embed_freeform.py
# ============================================================

embed_log_line() {
  printf '[%s] [%s] %s\n' "$(date '+%H:%M:%S')" "$1" "$2" >> "${SESSION_LOG}"
}

if [[ "${EMBED}" == "off" ]]; then
  printf 'Finalize: disabled (--no-embed).\n'
  printf 'embed: disabled (--no-embed)\n' >> "${SESSION_LOG}"
else
  EMBED_OUT_DIR="${REPO_ROOT}/outs/embeddings/${DATASET}"

  if [[ ! -f "${FREEFORM_INPUT}" ]]; then
    printf '[warn] Finalize skipped: input not found: %s\n' "${FREEFORM_INPUT}" >&2
    printf 'embed: skipped (input not found: %s)\n' "${FREEFORM_INPUT#${REPO_ROOT}/}" >> "${SESSION_LOG}"
  else
    mkdir -p "${EMBED_OUT_DIR}"
    printf 'Finalize: running each free-form script on %s\n' "${FREEFORM_INPUT#${REPO_ROOT}/}"
    printf 'embed_started_at: %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" >> "${SESSION_LOG}"
    printf 'embed_results:\n' >> "${SESSION_LOG}"

    pids=()
    for workdir in "${SURVIVORS[@]}"; do
      name=$(basename "${workdir}")
      run_nn=""
      if [[ "${name}" =~ freeform_run([0-9]{2})$ ]]; then
        run_nn="${BASH_REMATCH[1]}"
      else
        echo "FAIL" > "${workdir}/.embed_status"
        echo "cannot parse run number from ${name}" > "${workdir}/.embed_detail"
        embed_log_line "${name}" "embed FAIL (cannot parse run number)"
        continue
      fi
      seed_n=$((10#${run_nn}))
      stem="freeform_${VARIANT}_${MODEL_SLUG}_${REASONING_EFFORT}_${DATASET}_seed${seed_n}"
      parquet="${EMBED_OUT_DIR}/${stem}.parquet"
      if [[ -e "${parquet}" ]]; then
        echo "SKIP-EXISTS" > "${workdir}/.embed_status"
        echo "${parquet#${REPO_ROOT}/}" > "${workdir}/.embed_detail"
        embed_log_line "${name}" "embed SKIP-EXISTS -> ${parquet#${REPO_ROOT}/}"
        continue
      fi
      (
        log="${workdir}/logs/embed.log"
        set +e
        _t0=$(date +%s)
        python "${workdir}/embed_freeform.py" \
            --input  "${FREEFORM_INPUT}" \
            --output "${parquet}" \
            --seed   "${seed_n}" \
            >"${log}" 2>&1
        rc=$?
        _t1=$(date +%s)
        set -e
        echo "$((_t1 - _t0))" > "${workdir}/.embed_secs"
        if (( rc == 0 )) && [[ -s "${parquet}" ]]; then
          echo "OK" > "${workdir}/.embed_status"
          echo "${parquet#${REPO_ROOT}/}" > "${workdir}/.embed_detail"
        else
          echo "FAIL" > "${workdir}/.embed_status"
          echo "rc=${rc} parquet=${parquet#${REPO_ROOT}/} log=${workdir#${REPO_ROOT}/}/logs/embed.log" > "${workdir}/.embed_detail"
        fi
        echo "${rc}" > "${workdir}/.embed_rc"
        exit "${rc}"
      ) &
      pids+=("$!")
    done
    for pid in "${pids[@]}"; do
      wait "${pid}" || true
    done

    # Mirror per-run artifacts to outs/ so downstream plot scripts
    # (plot_three_way_summary.py) can render
    # from outs/ alone without walking work/. Fires for OK, FAIL, and
    # SKIP-EXISTS — failed runs still leave a status record so failure-rate
    # panels can be computed from outs/. Parallel to the summary.json /
    # cost.json mirrors in agent_edit_blueprint{,_claude}.sh.
    freeform_script_out_dir="${REPO_ROOT}/outs/freeform_scripts/${DATASET}/${PARENT_NAME}"
    mkdir -p "${freeform_script_out_dir}"
    for workdir in "${WORKDIRS[@]}"; do
      name=$(basename "${workdir}")
      if [[ "${name}" =~ freeform_run([0-9]{2})$ ]]; then
        seed_n=$((10#${BASH_REMATCH[1]}))
      else
        continue
      fi
      stem="freeform_${VARIANT}_${MODEL_SLUG}_${REASONING_EFFORT}_${DATASET}_seed${seed_n}"
      if [[ -s "${workdir}/embed_freeform.py" ]]; then
        cp "${workdir}/embed_freeform.py" "${freeform_script_out_dir}/embed_freeform_seed${seed_n}.py"
      fi
      embed_status=$(cat "${workdir}/.embed_status" 2>/dev/null || echo "SKIP-FILTERED")
      embed_detail=$(cat "${workdir}/.embed_detail" 2>/dev/null || true)
      stage_a_rc=$(cat "${workdir}/.stage_a_rc" 2>/dev/null || echo "")
      {
        printf '{\n'
        printf '  "variant": "%s",\n'        "${VARIANT}"
        printf '  "cohort": "%s",\n'         "${PARENT_NAME}"
        printf '  "seed": %d,\n'             "${seed_n}"
        printf '  "embed_status": "%s",\n'  "${embed_status}"
        printf '  "embed_detail": "%s",\n'  "${embed_detail}"
        printf '  "stage_a_rc": "%s"\n'      "${stage_a_rc}"
        printf '}\n'
      } > "${EMBED_OUT_DIR}/${stem}_status.json"
    done

    for workdir in "${WORKDIRS[@]}"; do
      name=$(basename "${workdir}")
      status=$(cat "${workdir}/.embed_status" 2>/dev/null || echo "SKIP-FILTERED")
      detail=$(cat "${workdir}/.embed_detail" 2>/dev/null || true)
      case "${status}" in
        OK)   embed_log_line "${name}" "embed OK -> ${detail}" ;;
        FAIL) embed_log_line "${name}" "embed FAIL ${detail}" ;;
      esac
    done

    {
      for workdir in "${WORKDIRS[@]}"; do
        name=$(basename "${workdir}")
        status=$(cat "${workdir}/.embed_status" 2>/dev/null || echo "SKIP-FILTERED")
        detail=$(cat "${workdir}/.embed_detail" 2>/dev/null || true)
        if [[ -n "${detail}" ]]; then
          printf '  %-40s %-14s %s\n' "${name}" "${status}" "${detail}"
        else
          printf '  %-40s %-14s\n' "${name}" "${status}"
        fi
      done
      printf 'embed_finished_at: %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')"
    } >> "${SESSION_LOG}"
  fi
fi

# ============================================================
# Embed: cost extraction — per-workdir cost.json (token usage,
# per-phase wall-clock, on-disk footprint). Aggregated across
# cohorts by scripts/resource_costs/collect.py.
# ============================================================

EXTRACT_COST_SCRIPT="${REPO_ROOT}/scripts/agentic_embedding/extract_cost.py"
if [[ -f "${EXTRACT_COST_SCRIPT}" ]]; then
  printf 'Embed: extracting per-run cost summaries.\n'
  cost_out_dir="${REPO_ROOT}/outs/embeddings/${DATASET}"
  for workdir in "${WORKDIRS[@]}"; do
    python "${EXTRACT_COST_SCRIPT}" "${workdir}" || true
    # Mirror cost.json next to the embedding parquet under outs/embeddings/<dataset>/
    # so downstream cost/quality joins can read from outs/ without walking work/.
    # Parallel to the cost.json mirror in agent_edit_blueprint{,_claude}.sh.
    if [[ -s "${workdir}/cost.json" && -s "${workdir}/.embed_detail" ]]; then
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
