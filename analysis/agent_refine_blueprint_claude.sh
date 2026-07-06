#!/usr/bin/env bash
# Claude variant of analysis/agent_refine_blueprint.sh. Takes user-selected
# Stage A workdirs from an already-completed
# work/<base>_{unrefined|noref_unrefined}_<reasoning>_<dataset>/ cohort and runs Stage B +
# the remaining downstream steps on each, producing a sibling `_refine`
# cohort (or `_noref_refine` when the source is a `_noref_unrefined` cohort) under
# work/<base>_<tag>_<reasoning>_<dataset>/. The noref contract is auto-detected
# from the source cohort name (no flag).
#
# The script mirrors agent_refine_blueprint.sh's workflow, validation,
# and output layout one-for-one — only the agent CLI changes — so the
# cohorts produced by this script and its codex sibling commingle
# correctly in metrics.py and friends.
#
# Isolation mode is *inherited* from the source cohort slug:
#   - source slug ends in `_userctx_unrefined_<r>_<ds>` → destination
#     slug ends in `_userctx_refine_<r>_<ds>`, run with project memory
#     + skills enabled (no --disable-slash-commands, no
#     CLAUDE_CODE_DISABLE_AUTO_MEMORY=1).
#   - otherwise → destination slug ends in `_refine_<r>_<ds>`, run with
#     full isolation (matches agent_edit_blueprint_claude.sh's default).
# This makes the refine step's claude context identical to the source
# Stage A run's context — a research-integrity guarantee. If you want to
# refine in the opposite mode, run agent_edit_blueprint_claude.sh
# --with-refine from scratch instead.
#
# WORKFLOW STEPS (per session; steps 2 and 6 are parallel across
# workdirs internally):
#
#   Artifact staging    — Validation + copy. Each input WORKDIR (which must have
#                .stage_a_rc == "0" and a populated embedder_blueprint.py
#                under a shared parent matching <base>_unrefined_<r>_<ds>)
#                is copied into the new
#                work/<base>_refine_<reasoning>_<dataset>/<src_run_name>/
#                destination, carrying scenario_config.jsonc, the Stage A
#                blueprint, agent_prompt_blueprint_edit.md, .stage_a_rc,
#                .stage_a_secs, and logs/agent_edit.log (so Cost's
#                cost.json can still attribute Stage A tokens).
#   Extract (Stage A)  — Early extraction. extract_run.py writes summary.json
#                from the populator-only state. Extract (Stage B) overwrites it.
#   Stage B    — iterative refinement, isolated per-iter
#                context). Snapshots embedder_blueprint.py to
#                embedder_blueprint_unrefined.py once at loop start.
#                Discriminator/Generator alternation as fresh
#                `claude -p` sessions; up to --max-iterations
#                (default 3). Each iter runs in a workdir stripped of
#                prior reports / refinement logs / carried-over Stage A
#                transcript / Stage A snapshot (moved into
#                .refine_hidden/). The Discriminator writes untagged
#                reports/report.{md,json}; regression_guard.py reverts
#                unauthorized Generator removals + re-spawns once with
#                a corrective addendum. Post-loop restores everything
#                under reports/ and logs/ with round-tagged filenames.
#                See agent_refine_blueprint.sh for the full description —
#                behavior is byte-identical (modulo log-prefix s/codex/agent/).
#   Extract (Stage B)    — Late extraction. Overwrites Extract (Stage A)'s summary.
#   Finalize    — Session-log finalization.
#   Embed    — Embedding (default on; opt out with --no-embed). The
#                Stage A sibling embed lands at the existing source-cohort
#                parquet path; the SKIP-EXISTS guard handles the duplicate.
#   Cost    — Cost extraction.
#
# HOW TO RUN (from the repo root):
#
#   bash analysis/agent_refine_blueprint_claude.sh \
#       work/claude_opus_4_7_unrefined_xhigh_human_pancreas/260525_blueprint_run01
#
#   # Refine multiple runs from the same cohort with brace expansion:
#   bash analysis/agent_refine_blueprint_claude.sh \
#       work/claude_opus_4_7_unrefined_xhigh_mouse_brain/260525_blueprint_run0{1,3,5}
#
#   # Skip the dataset-level embedding (no parquet written):
#   bash analysis/agent_refine_blueprint_claude.sh --no-embed \
#       work/claude_opus_4_7_unrefined_xhigh_human_pancreas/260525_blueprint_run01
#
#   # See all flags:
#   bash analysis/agent_refine_blueprint_claude.sh --help
#
# FLAGS: run with --help for the full flag reference (the usage() heredoc
# below is the single authoritative copy).

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

SOURCE_PROMPT_DISCRIMINATE="${REPO_ROOT}/prompt_templates/agent_prompt_blueprint_discriminate.md"
SOURCE_PROMPT_DISCRIMINATE_NOREF="${REPO_ROOT}/prompt_templates/agent_prompt_blueprint_discriminate_noref.md"
SOURCE_PROMPT_GENERATE="${REPO_ROOT}/prompt_templates/agent_prompt_blueprint_generate.md"
SOURCE_PROMPT_GENERATE_NOREF="${REPO_ROOT}/prompt_templates/agent_prompt_blueprint_generate_noref.md"
SOURCE_PROMPT_EDIT="${REPO_ROOT}/prompt_templates/agent_prompt_blueprint_edit.md"
SOURCE_EXTRACT_SCRIPT="${REPO_ROOT}/scripts/agentic_embedding/extract_run.py"
WORK_ROOT="${REPO_ROOT}/work"

DEFAULT_MAX_REFINE_ITERATIONS=3
DATASET_CHOICES=(mouse_brain human_pancreas mouse_gastrulation)
REASONING_CHOICES=(low medium high xhigh)

usage() {
  cat <<'EOF'
Usage: agent_refine_blueprint_claude.sh [OPTIONS] WORKDIR [WORKDIR ...]

Claude variant of agent_refine_blueprint.sh. Runs Stage B
(Discriminator/Generator refinement) + Phases 4-7 on a user-selected
subset of pre-existing Stage A workdirs.

Required positional:
  WORKDIR                One or more
                         work/<base>_{unrefined|noref_unrefined}_<reasoning>_<dataset>/<YYMMDD>_blueprint_runNN
                         paths. All must share a common parent cohort dir,
                         and each must have .stage_a_rc == "0" with
                         embedder_blueprint.py + scenario_config.jsonc. A
                         `_noref` source auto-selects the noref refine contract.

Optional:
  --max-iterations N         Stage B Discriminator-round cap (default: 3).
  --no-embed                 Skip Embed.
  --force                    Overwrite an existing destination workdir
                             in the _refine cohort (default: refuse).
  --model SLUG               Override auto-detected claude model.
  --reasoning-effort LEVEL   Override auto-detected reasoning effort
                             (low | medium | high | xhigh).
  -h, --help                 Print this message and exit.

Auto-detection:
  --model + --reasoning-effort are read from
  <shared_parent>/sessions/*.log (claude_model: / claude_reasoning_effort:).
  Isolation mode is inherited from the source cohort slug: if the slug
  contains _userctx, the destination keeps _userctx and runs claude with
  project memory + skills enabled; otherwise the destination runs in
  isolated mode (CLAUDE_CODE_DISABLE_AUTO_MEMORY=1,
  --disable-slash-commands --setting-sources user).

Output layout (<tag> = refine for an _unrefined source, noref_refine for a _noref source):
  Destination workdirs at  work/<base>_<tag>_<reasoning>_<dataset>/<YYMMDD>_blueprint_runNN/
  Session log at           work/<base>_<tag>_<reasoning>_<dataset>/sessions/<session_id>.log
  Embed parquets at      outs/embeddings/<dataset>/blueprint_<base>_<tag>_<reasoning>_<dataset>_seed<N>.parquet
EOF
}

INPUT_WORKDIRS=()
MAX_REFINE_ITERATIONS="${DEFAULT_MAX_REFINE_ITERATIONS}"
EMBED="on"
FORCE="off"
MODEL_OVERRIDE=""
REASONING_OVERRIDE=""

while (( $# > 0 )); do
  case "$1" in
    --max-iterations)
      [[ $# -ge 2 ]] || { printf -- '--max-iterations requires a value\n' >&2; exit 1; }
      MAX_REFINE_ITERATIONS="$2"; shift 2 ;;
    --max-iterations=*)
      MAX_REFINE_ITERATIONS="${1#--max-iterations=}"; shift ;;
    --no-embed)      EMBED="off"; shift ;;
    --embed)         EMBED="on";  shift ;;
    --force)         FORCE="on";  shift ;;
    --no-force)      FORCE="off"; shift ;;
    --model)
      [[ $# -ge 2 ]] || { printf -- '--model requires a value\n' >&2; exit 1; }
      MODEL_OVERRIDE="$2"; shift 2 ;;
    --model=*)
      MODEL_OVERRIDE="${1#--model=}"; shift ;;
    --reasoning-effort)
      [[ $# -ge 2 ]] || { printf -- '--reasoning-effort requires a value\n' >&2; exit 1; }
      REASONING_OVERRIDE="$2"; shift 2 ;;
    --reasoning-effort=*)
      REASONING_OVERRIDE="${1#--reasoning-effort=}"; shift ;;
    -h|--help)       usage; exit 0 ;;
    --*)
      printf 'Unknown flag: %s\n' "$1" >&2; usage >&2; exit 1 ;;
    *)
      INPUT_WORKDIRS+=("$1"); shift ;;
  esac
done

if (( ${#INPUT_WORKDIRS[@]} == 0 )); then
  printf 'At least one WORKDIR positional argument is required.\n' >&2
  usage >&2
  exit 1
fi

if ! [[ "${MAX_REFINE_ITERATIONS}" =~ ^[0-9]+$ ]] || (( MAX_REFINE_ITERATIONS < 1 )); then
  printf -- '--max-iterations must be a positive integer (got: %s)\n' \
    "${MAX_REFINE_ITERATIONS}" >&2
  exit 1
fi

valid_choice() {
  local needle="$1"; shift
  local item
  for item in "$@"; do [[ "${item}" == "${needle}" ]] && return 0; done
  return 1
}

# ============================================================
# Validate inputs, normalize to absolute paths, enforce single parent
# ============================================================

SOURCE_PARENT=""
NORM_WORKDIRS=()
for raw in "${INPUT_WORKDIRS[@]}"; do
  raw="${raw%/}"
  if [[ "${raw}" != /* ]]; then
    raw="$(cd "$(dirname "${raw}")" 2>/dev/null && pwd)/$(basename "${raw}")" \
      || { printf 'Cannot resolve: %s\n' "$1" >&2; exit 1; }
  fi
  if [[ ! -d "${raw}" ]]; then
    printf 'Not a directory: %s\n' "${raw}" >&2; exit 1
  fi
  if [[ ! -f "${raw}/embedder_blueprint.py" ]]; then
    printf 'Missing embedder_blueprint.py in: %s\n' "${raw}" >&2; exit 1
  fi
  if [[ ! -f "${raw}/scenario_config.jsonc" ]]; then
    printf 'Missing scenario_config.jsonc in: %s\n' "${raw}" >&2; exit 1
  fi
  rc_file="${raw}/.stage_a_rc"
  if [[ ! -f "${rc_file}" ]] || [[ "$(cat "${rc_file}")" != "0" ]]; then
    rc_value="$(cat "${rc_file}" 2>/dev/null || echo "<missing>")"
    printf 'Workdir does not have a successful Stage A (.stage_a_rc=%s): %s\n' \
      "${rc_value}" "${raw}" >&2
    exit 1
  fi
  parent="$(dirname "${raw}")"
  if [[ -z "${SOURCE_PARENT}" ]]; then
    SOURCE_PARENT="${parent}"
  elif [[ "${parent}" != "${SOURCE_PARENT}" ]]; then
    printf 'All WORKDIR inputs must share a common parent dir.\n  Got: %s\n  And: %s\n' \
      "${SOURCE_PARENT}" "${parent}" >&2
    exit 1
  fi
  NORM_WORKDIRS+=("${raw}")
done

# Parse the shared parent dir name. Stage-A source cohorts always end in
# `_unrefined_<reasoning>_<dataset>`; a noref source additionally carries a
# `_noref` grounding word just before that (`<base>_noref_unrefined_<r>_<ds>`),
# which selects the resource-unconstrained refine contract (NOREF=on). Peel the
# `_unrefined_<r>_<ds>` suffix first, then detect+strip a trailing `_noref` so
# BASE is the bare model slug (any `_userctx` marker survives in BASE for the
# isolation-mode check below).
SOURCE_PARENT_NAME="$(basename "${SOURCE_PARENT}")"
DATASET=""
REASONING=""
BASE=""
NOREF="off"
for ds in "${DATASET_CHOICES[@]}"; do
  for r in "${REASONING_CHOICES[@]}"; do
    suffix="_unrefined_${r}_${ds}"
    if [[ "${SOURCE_PARENT_NAME}" == *"${suffix}" ]]; then
      BASE="${SOURCE_PARENT_NAME%${suffix}}"
      REASONING="${r}"
      DATASET="${ds}"
      if [[ "${BASE}" == *_noref ]]; then
        NOREF="on"
        BASE="${BASE%_noref}"
      fi
      break 2
    fi
  done
done
if [[ -z "${BASE}" || -z "${REASONING}" || -z "${DATASET}" ]]; then
  printf 'Source parent dir name does not match <base>[_noref]_unrefined_<reasoning>_<dataset>: %s\n' \
    "${SOURCE_PARENT_NAME}" >&2
  printf 'This script only refines from `_unrefined` or `_noref_unrefined` source cohorts.\n' >&2
  exit 1
fi

# Under a noref source cohort, run the whole Stage B loop on the
# resource-optional contract by swapping the Discriminator/Generator prompts
# to their `_noref` siblings (the Stage-A populator contract +
# embedder_blueprint.py are already the noref variants, carried over per-run
# in Artifact staging). Mirrors agent_edit_blueprint_claude.sh's --noref --with-refine swap.
if [[ "${NOREF}" == "on" ]]; then
  SOURCE_PROMPT_DISCRIMINATE="${SOURCE_PROMPT_DISCRIMINATE_NOREF}"
  SOURCE_PROMPT_GENERATE="${SOURCE_PROMPT_GENERATE_NOREF}"
fi

# Inherit isolation mode from the source slug (presence of `_userctx`
# trailing the model identifier in BASE). This must run BEFORE
# claude/codex args are assembled.
ISOLATE="on"
if [[ "${BASE}" == *_userctx ]]; then
  ISOLATE="off"
fi

# ============================================================
# Auto-detect claude model + reasoning-effort from the source session log
# ============================================================

MODEL_AUTO=""
REASONING_AUTO=""
shopt -s nullglob
for log_path in "${SOURCE_PARENT}/sessions/"*.log; do
  while IFS= read -r line; do
    case "${line}" in
      claude_model:*)
        MODEL_AUTO="${line#claude_model:}"
        MODEL_AUTO="${MODEL_AUTO# }"
        ;;
      claude_reasoning_effort:*)
        REASONING_AUTO="${line#claude_reasoning_effort:}"
        REASONING_AUTO="${REASONING_AUTO# }"
        ;;
    esac
  done < "${log_path}"
done
shopt -u nullglob

MODEL="${MODEL_OVERRIDE:-${MODEL_AUTO}}"
REASONING_EFFORT="${REASONING_OVERRIDE:-${REASONING_AUTO:-${REASONING}}}"

if [[ -z "${MODEL}" ]]; then
  printf 'Could not auto-detect claude model from %s/sessions/*.log.\n' "${SOURCE_PARENT}" >&2
  printf 'Pass --model SLUG to override (e.g. claude-opus-4-7).\n' >&2
  exit 1
fi
if ! valid_choice "${REASONING_EFFORT}" "${REASONING_CHOICES[@]}"; then
  printf -- '--reasoning-effort must be one of: %s (got: %s)\n' \
    "${REASONING_CHOICES[*]}" "${REASONING_EFFORT}" >&2
  exit 1
fi
if [[ "${REASONING_EFFORT}" != "${REASONING}" ]]; then
  printf 'Reasoning-effort mismatch: parent dir encodes %s, resolved %s. Aborting.\n' \
    "${REASONING}" "${REASONING_EFFORT}" >&2
  printf 'Drop --reasoning-effort or fix the source parent dir name.\n' >&2
  exit 1
fi

# ============================================================
# Compute destination cohort + slugs
# ============================================================

# Output cohort tags mirror the edit script: a noref source refines into the
# `_noref_refine` cohort with its Stage-A snapshot landing in the `_noref_unrefined`
# cohort; a constrained source refines into `_refine` with an `_unrefined`
# snapshot. BASE already carries any `_userctx` marker from the source slug.
if [[ "${NOREF}" == "on" ]]; then
  MODEL_SLUG="${BASE}_noref_refine"
  UNREFINED_SLUG="${BASE}_noref_unrefined"
else
  MODEL_SLUG="${BASE}_refine"
  UNREFINED_SLUG="${BASE}_unrefined"
fi

PARENT_NAME="${MODEL_SLUG}_${REASONING_EFFORT}_${DATASET}"
PARENT_DIR="${WORK_ROOT}/${PARENT_NAME}"
SESSION_LOG_DIR="${PARENT_DIR}/sessions"
BLUEPRINT_INPUT="${REPO_ROOT}/data/preprocessed/${DATASET}_pp.h5ad"

CLAUDE_ARGS=(
  -p
  --model "${MODEL}"
  --effort "${REASONING_EFFORT}"
  --dangerously-skip-permissions
  --output-format text
)
if [[ "${ISOLATE}" == "on" ]]; then
  CLAUDE_ARGS+=(--disable-slash-commands --setting-sources user)
  export CLAUDE_CODE_DISABLE_AUTO_MEMORY=1
fi

for required_file in \
    "${SOURCE_PROMPT_DISCRIMINATE}" \
    "${SOURCE_PROMPT_GENERATE}" \
    "${SOURCE_PROMPT_EDIT}" \
    "${SOURCE_EXTRACT_SCRIPT}"; do
  if [[ ! -f "${required_file}" ]]; then
    printf 'Required upstream file not found: %s\n' "${required_file}" >&2
    exit 1
  fi
done

mkdir -p "${PARENT_DIR}" "${SESSION_LOG_DIR}"

RUN_DATE="$(date '+%y%m%d')"
SESSION_ID="${RUN_DATE}_$(date '+%H%M%S')_$$"
SESSION_LOG="${SESSION_LOG_DIR}/${SESSION_ID}.log"

# ============================================================
# Artifact staging: copy source workdirs into destination cohort
# ============================================================

WORKDIRS=()
for src in "${NORM_WORKDIRS[@]}"; do
  src_name="$(basename "${src}")"
  dst="${PARENT_DIR}/${src_name}"
  if [[ -e "${dst}" ]]; then
    if [[ "${FORCE}" != "on" ]]; then
      printf 'Destination workdir already exists: %s\n' "${dst}" >&2
      printf 'Pass --force to overwrite.\n' >&2
      exit 1
    fi
    rm -rf "${dst}"
  fi
  mkdir -p "${dst}/logs" "${dst}/reports"
  cp "${src}/embedder_blueprint.py" "${dst}/embedder_blueprint.py"
  cp "${src}/scenario_config.jsonc" "${dst}/scenario_config.jsonc"
  if [[ -f "${src}/agent_prompt_blueprint_edit.md" ]]; then
    cp "${src}/agent_prompt_blueprint_edit.md" "${dst}/agent_prompt_blueprint_edit.md"
  else
    cp "${SOURCE_PROMPT_EDIT}" "${dst}/agent_prompt_blueprint_edit.md"
  fi
  cp "${SOURCE_PROMPT_DISCRIMINATE}" "${dst}/agent_prompt_blueprint_discriminate.md"
  cp "${SOURCE_PROMPT_GENERATE}"     "${dst}/agent_prompt_blueprint_generate.md"
  cp "${SOURCE_EXTRACT_SCRIPT}" "${dst}/"
  if [[ -f "${src}/logs/agent_edit.log" ]]; then
    cp "${src}/logs/agent_edit.log" "${dst}/logs/agent_edit.log"
  fi
  cp "${src}/.stage_a_rc" "${dst}/.stage_a_rc"
  [[ -f "${src}/.stage_a_secs" ]] && cp "${src}/.stage_a_secs" "${dst}/.stage_a_secs"
  WORKDIRS+=("${dst}")
  printf 'Copied %s -> %s\n' "${src}" "${dst}"
done

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

{
  printf 'session_id: %s\n' "${SESSION_ID}"
  printf 'started_at: %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')"
  printf 'mode: refine-from-existing\n'
  printf 'source_cohort: %s\n'   "${SOURCE_PARENT_NAME}"
  printf 'n_workdirs: %d\n'      "${#WORKDIRS[@]}"
  printf 'dataset: %s\n'                 "${DATASET}"
  printf 'agent_cli: claude\n'
  printf 'claude_model: %s\n'            "${MODEL}"
  printf 'claude_reasoning_effort: %s\n' "${REASONING_EFFORT}"
  printf 'claude_isolate: %s\n'          "${ISOLATE}"
  printf 'with_refine: on\n'
  printf 'max_refine_iterations: %d\n'   "${MAX_REFINE_ITERATIONS}"
  printf 'embed: %s\n'                   "${EMBED}"
  printf 'claude_web_search: always-on (default in -p mode)\n'
  printf 'claude_network_access: always-on (--dangerously-skip-permissions)\n'
  if [[ "${ISOLATE}" == "on" ]]; then
    printf 'claude_auto_memory: disabled (CLAUDE_CODE_DISABLE_AUTO_MEMORY=1)\n'
    printf 'claude_skills: disabled (--disable-slash-commands)\n'
    printf 'claude_setting_sources: user (ignores project/local settings)\n'
  else
    printf 'claude_auto_memory: enabled (inherited from source _userctx)\n'
    printf 'claude_skills: enabled (inherited from source _userctx)\n'
    printf 'claude_setting_sources: user,project,local (default)\n'
  fi
  printf 'parent_dir: %s\n' "${PARENT_NAME}"
  printf 'workdirs:\n'
  for w in "${WORKDIRS[@]}"; do
    printf '  - %s\n' "$(basename "${w}")"
  done
} > "${SESSION_LOG}"

printf 'Session %s; shared log: %s\n' "${SESSION_ID}" "${SESSION_LOG}"

# ============================================================
# Extract (Stage A): early extraction (populator-only summary)
# ============================================================

printf 'Extract (Stage A): early extraction (from carried-over Stage A blueprint).\n'
for workdir in "${WORKDIRS[@]}"; do
  python "${workdir}/extract_run.py" "${workdir}" || true
done

# ============================================================
# Stage B (Discriminator/Generator refinement loop).
# Mirrors agent_edit_blueprint_claude.sh Stage B: isolated per-iter
# context, regression_guard.py, post-loop restore.
# ============================================================

DISCRIMINATE_PROMPT="$(<"${SOURCE_PROMPT_DISCRIMINATE}")"
GENERATE_PROMPT="$(<"${SOURCE_PROMPT_GENERATE}")"

SURVIVORS=("${WORKDIRS[@]}")
REGRESSION_GUARD="${REPO_ROOT}/scripts/agentic_embedding/regression_guard.py"

printf 'Stage B (refinement, isolated per-iter context) for %d workdir(s); max %d iter(s).\n' \
    "${#SURVIVORS[@]}" "${MAX_REFINE_ITERATIONS}"
pids=()
for workdir in "${SURVIVORS[@]}"; do
  (
    cd "${workdir}"
    cp embedder_blueprint.py embedder_blueprint_unrefined.py
    [[ -f summary.json ]] && cp summary.json summary_unrefined.json
    mkdir -p .refine_hidden/snapshot .refine_hidden/stage_a reports logs
    # Park the carried-over Stage A transcript so refinement agents never read it.
    [[ -f logs/agent_edit.log ]] && mv logs/agent_edit.log .refine_hidden/stage_a/agent_edit.log

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
      claude "${CLAUDE_ARGS[@]}" "${DISCRIMINATE_PROMPT}" \
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
        claude "${CLAUDE_ARGS[@]}" "${gen_prompt_used}" \
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

    # ---- Post-loop restoration ----
    mkdir -p reports logs
    [[ -f .refine_hidden/stage_a/agent_edit.log ]] && \
      mv .refine_hidden/stage_a/agent_edit.log logs/agent_edit.log
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
      [[ -f "${iter_dir}/discriminate.log"     ]] && cp "${iter_dir}/discriminate.log"     "logs/agent_discriminate_${nn}.log"
      [[ -f "${iter_dir}/regression_guard.log" ]] && cp "${iter_dir}/regression_guard.log" "logs/regression_guard_${nn}.log"
      for gen_attempt_log in "${iter_dir}"/generate_attempt*.log; do
        [[ -f "${gen_attempt_log}" ]] || continue
        attempt=$(basename "${gen_attempt_log}" .log | sed 's/^generate_attempt//')
        if [[ "${attempt}" == "01" ]]; then
          cp "${gen_attempt_log}" "logs/agent_generate_${nn}.log"
        else
          cp "${gen_attempt_log}" "logs/agent_generate_${nn}_attempt${attempt}.log"
        fi
      done
    done
    shopt -u nullglob
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

# ============================================================
# Extract (Stage B): late extraction
# ============================================================

printf 'Extract (Stage B): extracting blueprint summaries.\n'
for workdir in "${WORKDIRS[@]}"; do
  python "${workdir}/extract_run.py" "${workdir}" || true
done

# ============================================================
# Finalize: session log finalization
# ============================================================

{
  printf 'finished_at: %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')"
  printf 'results:\n'
  for workdir in "${WORKDIRS[@]}"; do
    name=$(basename "${workdir}")
    b_status=$(cat "${workdir}/.stage_b_status" 2>/dev/null || echo "?")
    printf '  %-40s refine=%s\n' "${name}" "${b_status}"
  done
} >> "${SESSION_LOG}"

# ============================================================
# Embed: embedding
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
    printf 'Embed: embedding refined blueprints onto %s\n' "${BLUEPRINT_INPUT#${REPO_ROOT}/}"
    printf 'embed_started_at: %s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" >> "${SESSION_LOG}"
    printf 'embed_results:\n' >> "${SESSION_LOG}"

    pids=()

    for workdir in "${WORKDIRS[@]}"; do
      name=$(basename "${workdir}")
      b_status_line=$(cat "${workdir}/.stage_b_status" 2>/dev/null || echo "")
      eligible=0
      case "${b_status_line}" in
        OK*|MAX-ITER*) eligible=1 ;;
      esac

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

      unrefined_src="${workdir}/embedder_blueprint_unrefined.py"
      unrefined_stem="blueprint_${UNREFINED_SLUG}_${REASONING_EFFORT}_${DATASET}_seed${seed_n}"
      unrefined_parquet="${EMBED_OUT_DIR}/${unrefined_stem}.parquet"
      unrefined_report="${EMBED_OUT_DIR}/${unrefined_stem}_report.md"
      emit_unrefined=0
      if [[ -f "${unrefined_src}" && ! -e "${unrefined_parquet}" ]]; then
        emit_unrefined=1
      fi

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
            blueprint_model_dir="${REPO_ROOT}/outs/blueprints/${DATASET}/${MODEL_SLUG}_${REASONING_EFFORT}"
            blueprint_dst="${blueprint_model_dir}/embedder_blueprint_seed${seed_n}.py"
            mkdir -p "${blueprint_model_dir}"
            cp "${workdir}/embedder_blueprint.py" "${blueprint_dst}"
            python "${workdir}/extract_run.py" "${workdir}" \
                --report-out "${report}" \
                --dataset "${DATASET}" \
                --model "${MODEL}" --reasoning "${REASONING_EFFORT}" \
                --input "${BLUEPRINT_INPUT#${REPO_ROOT}/}" \
                --parquet "${parquet#${REPO_ROOT}/}" \
                >/dev/null 2>&1 \
              || printf '# Blueprint embedding\n\n(report generation failed)\n' > "${report}"
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

    for workdir in "${WORKDIRS[@]}"; do
      name=$(basename "${workdir}")
      status=$(cat "${workdir}/.embed_status" 2>/dev/null || echo "SKIP-FILTERED")
      detail=$(cat "${workdir}/.embed_detail" 2>/dev/null || true)
      case "${status}" in
        OK)            embed_log_line "${name}" "embed OK -> ${detail}" ;;
        FAIL)          embed_log_line "${name}" "embed FAIL ${detail}" ;;
        SKIP-EXISTS|SKIP-FILTERED) ;;
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
# Cost: cost extraction
# ============================================================

EXTRACT_COST_SCRIPT="${REPO_ROOT}/scripts/agentic_embedding/extract_cost.py"
if [[ -f "${EXTRACT_COST_SCRIPT}" ]]; then
  printf 'Cost: extracting per-run cost summaries.\n'
  cost_out_dir="${REPO_ROOT}/outs/embeddings/${DATASET}"
  for workdir in "${WORKDIRS[@]}"; do
    python "${EXTRACT_COST_SCRIPT}" "${workdir}" || true
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
