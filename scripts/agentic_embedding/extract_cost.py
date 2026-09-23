#!/usr/bin/env python3
"""Extract a cost / wall-clock / storage summary from one agent run directory.

Reads the orchestrator's per-phase status/timing files, greps the codex/claude
transcripts for token usage, parses the embedder's `[timing]` line, and stats
the run's artifacts, writing <run_dir>/cost.json for cross-run aggregation by
scripts/resource_costs/collect.py. Works for both blueprint runs
(work/<cohort>/<date>_blueprint_runNN/) and free-form runs
(work/freeform_<cohort>/<date>_freeform_runNN/).

The orchestrators (analysis/agent_edit_*.sh) call this once per workdir as a
final phase. It is also harmless to run by hand on any existing run dir.

Usage:
    python scripts/agentic_embedding/extract_cost.py <run_dir>
"""
from __future__ import annotations

import argparse
import datetime as _dt
import json
import re
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

# Codex prints the session token total on the line after a bare "tokens used"
# line, with "." (or " ") as a thousands separator (e.g. "48.158" -> 48158).
# Claude (`--output-format json`) instead prints one JSON result with a `usage` block.
_TOKENS_LINE_RE = re.compile(r"^\s*tokens used\s*$", re.IGNORECASE)
_NUMERIC_LINE_RE = re.compile(r"^[\d.,\s]+$")
_TIMING_RE = re.compile(r"\[timing\]\s+(.*)$", re.MULTILINE)


def _read_text(p: Path) -> str | None:
    try:
        return p.read_text(errors="ignore")
    except OSError:
        return None


def _log_path(run_dir: Path, name: str) -> Path | None:
    """Locate a session log: prefer <run>/logs/<name>, fall back to <run>/<name>.

    The orchestrators (analysis/agent_edit_*.sh) write logs under <run>/logs/
    since the Stage-B refinement change. Run dirs committed before that
    convention have logs at the workdir root; the fallback keeps them readable.
    """
    new = run_dir / "logs" / name
    if new.exists():
        return new
    legacy = run_dir / name
    return legacy if legacy.exists() else None


def _log_glob(run_dir: Path, pattern: str) -> list[Path]:
    """Glob for session logs: prefer <run>/logs/<pattern>, fall back to <run>/<pattern>.

    If any matches exist under logs/, those are returned exclusively (no
    commingling with stragglers at the workdir root that would double-count).
    """
    hits = sorted((run_dir / "logs").glob(pattern))
    if hits:
        return hits
    return sorted(run_dir.glob(pattern))


def _read_head(p: Path, n: int = 4096) -> str:
    try:
        with p.open(errors="ignore") as fh:
            return fh.read(n)
    except OSError:
        return ""


def _first_line(p: Path) -> str | None:
    txt = _read_text(p)
    if not txt or not txt.strip():
        return None
    return txt.splitlines()[0].strip()


def _read_int(p: Path) -> int | None:
    s = _first_line(p)
    if s is None:
        return None
    try:
        return int(s)
    except ValueError:
        return None


def _path_bytes(p: Path) -> int:
    """Bytes of a file, recursive sum of a directory, or 0 if missing."""
    if not p.exists():
        return 0
    if p.is_file():
        return p.stat().st_size
    return sum(f.stat().st_size for f in p.rglob("*") if f.is_file())


_SESSION_LINE_RE = re.compile(
    r"^\[(?P<h>\d\d):(?P<m>\d\d):(?P<s>\d\d)\]\s+\[(?P<run>[^\]]+)\]\s+(?P<rest>.+)$"
)


def _session_phase_secs(run_dir: Path) -> dict[str, int | None]:
    """Fallback per-phase wall-clock (seconds), parsed from the cohort
    sessions/*.log when the orchestrator's `.{stage_a,embed}_secs`
    dotfiles are absent (runs predating the cost-extraction hook).

    The session log only timestamps stage_a per run, so stage_b / embed
    stay unknown for such runs."""
    out: dict[str, int | None] = {"stage_a_s": None, "stage_b_s": None, "embed_s": None}
    log_path = None
    sess = _read_text(run_dir / "session.txt")
    if sess:
        for line in sess.splitlines():
            if line.startswith("session_log:"):
                log_path = Path(line.split(":", 1)[1].strip())
                break
    txt = _read_text(log_path) if log_path else None
    if not txt:
        return out
    run_id = run_dir.name
    a_start = a_end = None
    for line in txt.splitlines():
        mm = _SESSION_LINE_RE.match(line)
        if not mm or mm.group("run") != run_id:
            continue
        secs = int(mm.group("h")) * 3600 + int(mm.group("m")) * 60 + int(mm.group("s"))
        rest = mm.group("rest")
        if rest.startswith("stage_a starting"):
            a_start = secs
        elif rest.startswith("stage_a finished"):
            a_end = secs
    if a_start is not None and a_end is not None:
        out["stage_a_s"] = a_end - a_start if a_end >= a_start else a_end - a_start + 86400
    return out


def _claude_json_tokens(lines: list[str]) -> int | None:
    """Total tokens from the last `claude -p --output-format json` result line."""
    for line in reversed(lines):
        if not line.lstrip().startswith("{"):
            continue
        try:
            usage = json.loads(line).get("usage") or {}
        except (ValueError, AttributeError):
            continue
        keys = ("input_tokens", "output_tokens",
                "cache_creation_input_tokens", "cache_read_input_tokens")
        if any(k in usage for k in keys):
            return sum(int(usage.get(k) or 0) for k in keys)
    return None


def _tokens_in_log(log_path: Path) -> int | None:
    """Sum every `tokens used\\n<N>` block in a codex transcript, else read a
    claude JSON result; None if neither is present."""
    txt = _read_text(log_path)
    if not txt:
        return None
    lines = txt.splitlines()
    total, found = 0, False
    for i, line in enumerate(lines):
        if not _TOKENS_LINE_RE.match(line):
            continue
        for j in range(i + 1, min(i + 4, len(lines))):
            cand = lines[j].strip()
            if not cand:
                continue
            if _NUMERIC_LINE_RE.match(cand):
                digits = re.sub(r"\D", "", cand)
                if digits:
                    total += int(digits)
                    found = True
            break
    return total if found else _claude_json_tokens(lines)


def _parse_timing_line(log_path: Path) -> dict:
    """Parse the embedder's `[timing] k=v k=v ...` line into a {k: value} dict
    (numbers parsed as int/float; non-numeric values — device, host — kept as
    strings)."""
    txt = _read_text(log_path)
    out: dict[str, object] = {}
    if not txt:
        return out
    m = _TIMING_RE.search(txt)
    if not m:
        return out
    for tok in m.group(1).split():
        if "=" not in tok:
            continue
        k, v = tok.split("=", 1)
        if k.endswith("_bytes"):
            try:
                out[k] = int(v)
                continue
            except ValueError:
                pass
        try:
            out[k] = float(v)            # "3.42", "0.0", "nan", "1e-3" all parse
        except ValueError:
            out[k] = v                   # device=cpu, host=<machine>, model_source=authored
    return out


def _build_cost(run_dir: Path) -> dict:
    run_id = run_dir.name
    cohort = run_dir.parent.name
    is_freeform = "freeform" in cohort or "_freeform_run" in run_id
    if "freeform_naive" in cohort:
        variant = "freeform_naive"
    elif "freeform_constrained" in cohort:
        variant = "freeform_constrained"
    elif is_freeform:
        variant = "freeform"
    else:
        variant = "blueprint"

    # Stage A transcript: codex_edit.log (blueprint+codex), agent_edit.log
    # (blueprint+claude, OR free-form+codex — that orchestrator names its codex
    # log `agent_edit.log`), or opencode_edit.log (blueprint+opencode).
    # Disambiguate the CLI by the log filename / transcript header below.
    edit_log = next(
        (p for p in (
            _log_path(run_dir, "codex_edit.log"),
            _log_path(run_dir, "agent_edit.log"),
            _log_path(run_dir, "opencode_edit.log"),
        ) if p is not None),
        None,
    )
    is_codex = bool(edit_log) and ("OpenAI Codex" in _read_head(edit_log)
                                   or (edit_log and edit_log.name.startswith("codex")))
    is_opencode = bool(edit_log) and edit_log.name.startswith("opencode")

    # Stage B (refinement) transcripts: Discriminator + Generator, per iter.
    refine_logs = (
        _log_glob(run_dir, "codex_discriminate_*.log")
        + _log_glob(run_dir, "agent_discriminate_*.log")
        + _log_glob(run_dir, "opencode_discriminate_*.log")
        + _log_glob(run_dir, "codex_generate_*.log")
        + _log_glob(run_dir, "agent_generate_*.log")
        + _log_glob(run_dir, "opencode_generate_*.log")
    )

    tokens_a = _tokens_in_log(edit_log) if edit_log else None
    tokens_b_refine_parts = [t for t in (_tokens_in_log(p) for p in refine_logs) if t is not None]
    tokens_b_refine = sum(tokens_b_refine_parts) if tokens_b_refine_parts else None
    tokens_total = None
    tokens_found = (tokens_a is not None) or (tokens_b_refine is not None)
    if tokens_found:
        tokens_total = (tokens_a or 0) + (tokens_b_refine or 0)

    embed_log = _log_path(run_dir, "embed.log")
    embed_timing = _parse_timing_line(embed_log) if embed_log is not None else {}

    # Stage A sibling embed (--with-refine workdirs only): parallel timing
    # extraction from logs/embed_unrefined.log. Yields embed_unrefined_*
    # counterparts to every embed_* key so collect.py can emit a second cost
    # row for the Stage A cohort.
    embed_unrefined_log = _log_path(run_dir, "embed_unrefined.log")
    embed_unrefined_timing = (
        _parse_timing_line(embed_unrefined_log) if embed_unrefined_log is not None else {}
    )

    # The agent's "model" is the script it produced.
    model_file = (run_dir / "embed_freeform.py") if is_freeform else (run_dir / "embedder_blueprint.py")
    model_storage_bytes = _path_bytes(model_file)
    # Stage A sibling model file (only present when refinement snapshotted it).
    unrefined_model_file = (
        None if is_freeform else (run_dir / "embedder_blueprint_unrefined.py")
    )
    unrefined_model_storage_bytes = (
        _path_bytes(unrefined_model_file) if unrefined_model_file is not None else 0
    )

    embed_detail = _first_line(run_dir / ".embed_detail")
    embedding_parquet = None
    embedding_bytes = embed_timing.get("embedding_bytes")
    if embed_detail and embed_detail.endswith(".parquet"):
        cand = REPO_ROOT / embed_detail
        if cand.exists():
            embedding_parquet = embed_detail
            embedding_bytes = cand.stat().st_size

    embed_unrefined_detail = _first_line(run_dir / ".embed_unrefined_detail")
    embed_unrefined_parquet = None
    embed_unrefined_bytes = embed_unrefined_timing.get("embedding_bytes")
    if embed_unrefined_detail and embed_unrefined_detail.endswith(".parquet"):
        cand = REPO_ROOT / embed_unrefined_detail
        if cand.exists():
            embed_unrefined_parquet = embed_unrefined_detail
            embed_unrefined_bytes = cand.stat().st_size

    # Per-phase wall-clock: prefer the orchestrator's exact dotfiles; fall back
    # to the session log for runs that predate that hook.
    sess_secs = (
        None if all((run_dir / f".{p}_secs").exists() for p in ("stage_a", "stage_b"))
        else _session_phase_secs(run_dir)
    )

    def _phase_secs(name: str) -> int | None:
        dotfile = run_dir / f".{name}_secs"
        if dotfile.exists():
            return _read_int(dotfile)
        return sess_secs.get(f"{name}_s") if sess_secs else None

    stage_b_status = _first_line(run_dir / ".stage_b_status")
    stage_b_s = _phase_secs("stage_b")
    stage_b_iters = _read_int(run_dir / ".stage_b_iters")
    stage_b_final_rating = _first_line(run_dir / ".stage_b_final_rating")

    return {
        "run_id": run_id,
        "cohort": cohort,
        "variant": variant,
        "agent_cli": (
            "opencode" if is_opencode
            else "codex" if is_codex
            else ("claude" if edit_log else None)
        ),
        "extracted_at": _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds"),
        # Status (mirrors the orchestrator's dotfiles).
        "stage_a_rc": _read_int(run_dir / ".stage_a_rc"),
        # Stage B (refinement loop) — only populated when --with-refine ran.
        "stage_b_status": stage_b_status,
        "stage_b_iters": stage_b_iters,
        "stage_b_final_rating": stage_b_final_rating,
        "embed_status": _first_line(run_dir / ".embed_status"),
        "embed_detail": embed_detail,
        # Wall-clock per phase (seconds).
        "stage_a_s": _phase_secs("stage_a"),
        "stage_b_s": stage_b_s,
        "embed_s": _phase_secs("embed"),
        # Token cost (codex transcript or claude JSON result).
        "agent_tokens_stage_a": tokens_a,
        "agent_tokens_stage_b_refine": tokens_b_refine,
        "agent_tokens_total": tokens_total,
        "agent_tokens_note": None if tokens_found
        else "no 'tokens used' line or claude JSON usage in the agent transcript",
        # Embed embedding-build timing/footprint, parsed from embed.log's
        # `[timing]` line. New runs emit the 4-bucket taxonomy (walltime = load_data
        # + obtain_model + build_embedding + other); fall back to the legacy
        # fit_s/transform_s names for embed.log files that predate it.
        "embed_walltime_s": embed_timing.get("walltime_s"),
        "embed_load_data_s": embed_timing.get("load_data_s"),
        "embed_obtain_model_s": embed_timing.get("obtain_model_s", embed_timing.get("fit_s")),
        "embed_build_embedding_s": embed_timing.get("build_embedding_s", embed_timing.get("transform_s")),
        "embed_other_s": embed_timing.get("other_s"),
        "embed_model_source": embed_timing.get("model_source", "authored"),
        "embed_model_download_s": embed_timing.get("model_download_s", 0.0),
        "embed_peak_rss_mb": embed_timing.get("peak_rss_mb"),
        "embed_rss_after_load_mb": embed_timing.get("rss_after_load_mb"),
        "embed_rss_after_obtain_mb": embed_timing.get("rss_after_obtain_mb"),
        "embed_rss_after_build_mb": embed_timing.get("rss_after_build_mb"),
        "embed_device": embed_timing.get("device"),   # the local embedding-build device (cpu)
        "embed_host": embed_timing.get("host"),
        # The Stage A / Stage B compute runs on the LLM provider's infra, not here.
        "agent_compute_device": "remote-llm",
        # On-disk footprint.
        "model_file": model_file.name if model_file.exists() else None,
        "model_storage_bytes": model_storage_bytes,
        "embedding_parquet": embedding_parquet,
        "embedding_bytes": embedding_bytes,
        "run_dir_bytes": _path_bytes(run_dir),
        # Stage A sibling — populated when embedder_blueprint_unrefined.py was
        # snapshotted by --with-refine and Embed emitted a Stage A parquet.
        # collect.py emits a second cost row keyed on these when present.
        "embed_unrefined_status": _first_line(run_dir / ".embed_unrefined_status"),
        "embed_unrefined_detail": embed_unrefined_detail,
        "embed_unrefined_s": _read_int(run_dir / ".embed_unrefined_secs"),
        "embed_unrefined_walltime_s": embed_unrefined_timing.get("walltime_s"),
        "embed_unrefined_load_data_s": embed_unrefined_timing.get("load_data_s"),
        "embed_unrefined_obtain_model_s": embed_unrefined_timing.get(
            "obtain_model_s", embed_unrefined_timing.get("fit_s")),
        "embed_unrefined_build_embedding_s": embed_unrefined_timing.get(
            "build_embedding_s", embed_unrefined_timing.get("transform_s")),
        "embed_unrefined_other_s": embed_unrefined_timing.get("other_s"),
        "embed_unrefined_model_source": embed_unrefined_timing.get("model_source", "authored"),
        "embed_unrefined_model_download_s": embed_unrefined_timing.get("model_download_s", 0.0),
        "embed_unrefined_peak_rss_mb": embed_unrefined_timing.get("peak_rss_mb"),
        "embed_unrefined_rss_after_load_mb": embed_unrefined_timing.get("rss_after_load_mb"),
        "embed_unrefined_rss_after_obtain_mb": embed_unrefined_timing.get("rss_after_obtain_mb"),
        "embed_unrefined_rss_after_build_mb": embed_unrefined_timing.get("rss_after_build_mb"),
        "embed_unrefined_device": embed_unrefined_timing.get("device"),
        "embed_unrefined_host": embed_unrefined_timing.get("host"),
        "unrefined_model_file": (
            unrefined_model_file.name
            if (unrefined_model_file is not None and unrefined_model_file.exists())
            else None
        ),
        "unrefined_model_storage_bytes": unrefined_model_storage_bytes,
        "embed_unrefined_parquet": embed_unrefined_parquet,
        "embed_unrefined_bytes": embed_unrefined_bytes,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("run_dir", type=Path, help="Path to a work/<cohort>/<run>/ directory.")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.resolve()
    if not run_dir.is_dir():
        print(f"not a directory: {run_dir}", file=sys.stderr)
        sys.exit(1)
    cost = _build_cost(run_dir)
    out_path = run_dir / "cost.json"
    out_path.write_text(json.dumps(cost, indent=2) + "\n")
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
