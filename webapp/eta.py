"""Blueprint-construction time estimates from baked timing data."""
from __future__ import annotations
import json
from functools import lru_cache
from pathlib import Path

_TABLE_PATH = Path(__file__).resolve().parent / "eta_table.json"


@lru_cache(maxsize=1)
def _load():
    try:
        return json.loads(_TABLE_PATH.read_text())
    except FileNotFoundError:
        return {}


def _model_key(model: str) -> str:
    return model.lower().replace("-", "_").replace(".", "_")


def estimate(model: str, effort: str, refine: bool):
    """Return (estimated_seconds | None, n_runs)."""
    key = f"{_model_key(model)}|{effort}|{int(bool(refine))}"
    entry = _load().get(key)
    if not entry:
        return None, 0
    seconds = entry.get("stage_a_s") or 0.0
    if refine and entry.get("stage_b_s"):
        seconds += entry["stage_b_s"]
    return seconds, entry.get("n", 0)


def format_eta(seconds, n):
    if not seconds:
        return "No timing data for this combination yet."
    minutes = seconds / 60.0
    return (f"Estimated blueprint time: ~{minutes:.0f} min "
            f"(based on {n} past run{'s' if n != 1 else ''}).")
