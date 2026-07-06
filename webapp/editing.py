"""In-browser blueprint editing.

The generated blueprint can be edited in the Generate tab; the edited copy is
saved per-session under ``RUNTIME/<token>/edited/`` and **preferred** over the
as-generated file when mapping with the in-session source — no download / edit /
re-upload round-trip.
"""
from __future__ import annotations

from pathlib import Path

from webapp.mapping import validate_blueprint
from webapp.pipeline import RUNTIME, resolve_run_dir


def generated_blueprint_path(active):
    """Path to the as-generated in-session blueprint (the driver's output)."""
    if not active or not active.get("workdir"):
        return None
    return resolve_run_dir(active["workdir"]) / "embedder_blueprint.py"


def edited_blueprint_path(token) -> Path:
    return RUNTIME / token / "edited" / "embedder_blueprint.py"


def load_generated_source(active) -> str:
    """Text of the as-generated blueprint, or '' if not available yet."""
    p = generated_blueprint_path(active)
    return p.read_text() if p and p.exists() else ""


def save_edited(token, source):
    """Write the edited source to the session's edited path and validate it.

    The file is written even when invalid (so work-in-progress persists); the
    returned validity lets the UI warn, and mapping re-validates regardless.
    Returns ``(ok: bool, errors: list[str])``."""
    dest = edited_blueprint_path(token)
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(source or "")
    ok, _summary, errors = validate_blueprint(dest)
    return ok, errors


def revert_edited(token) -> None:
    """Remove the edited copy so mapping falls back to the generated blueprint."""
    p = edited_blueprint_path(token)
    if p.exists():
        p.unlink()


def resolve_session_blueprint(active, token):
    """Blueprint to map for the in-session source: the edited copy if it exists,
    otherwise the as-generated file."""
    edited = edited_blueprint_path(token)
    if edited.exists():
        return edited
    return generated_blueprint_path(active)
