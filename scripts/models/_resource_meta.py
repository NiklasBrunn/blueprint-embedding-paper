"""Shared resource-cost bookkeeping for the model scripts (scripts/models/*.py).

One small module, imported by every model script, so the cost/wall-clock/storage
record is *standardised* across methods (and parseable by
scripts/resource_costs/collect.py). It defines:

  * `rss_mb()`            — current peak RSS high-water mark (MiB), self+children.
  * `path_bytes(p)`       — file size / recursive dir size / 0.
  * `time_split(...)`     — derive `other_s` so the four time buckets sum to
                            `walltime_s` exactly.
  * `resource_lines(...)` — the `## Resources` markdown block every report ends
                            with (timing breakdown · RSS trajectory · storage).

Time taxonomy (additive — sums to the script's wall-clock):
  walltime_s = load_data_s + obtain_model_s + build_embedding_s + other_s
    load_data_s       read the input AnnData from disk.
    obtain_model_s    get a usable model: train it (de-novo VAEs, PCA loadings)
                      or load a persisted/pretrained checkpoint + adapt it to the
                      dataset's gene space. `model_source` tags which: "trained" /
                      "loaded" / "none" (UCell — nothing acquired) /
                      "authored" (agent methods — not used here). `None`/`n/a`
                      when the upstream code loads-and-forwards atomically (UCE).
    build_embedding_s the forward pass / scoring that turns (data + model) → the
                      embedding array.
    other_s           the residual within the timed span: densify/copy/setup
                      glue + writing the parquet. Often small; carried explicitly
                      so the sum is exact.
Memory is a high-water *max*, not a sum, so it is not decomposed; instead we
record a 4-point trajectory (RSS after load / after model / after embedding /
final peak) plus the headline peak. Download time/footprint of pretrained
checkpoints is a one-time cost (recorded as `model_download_s`, 0 when cached),
never folded into the per-run total.
"""
from __future__ import annotations

import resource
import sys


def rss_mb() -> float:
    """Current peak resident-set-size high-water mark of this process tree, MiB
    (ru_maxrss is bytes on macOS, KiB on Linux/BSD; self+children so it also
    covers any subprocess/worker — children is 0 for the in-process methods)."""
    rss = (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
           + resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
    return rss / 1024.0**2 if sys.platform == "darwin" else rss / 1024.0


def path_bytes(p) -> int:
    """Bytes of a file, recursive sum over a directory, or 0 if missing."""
    from pathlib import Path
    p = Path(p)
    if not p.exists():
        return 0
    return p.stat().st_size if p.is_file() else sum(
        f.stat().st_size for f in p.rglob("*") if f.is_file())


def time_split(walltime_s: float, load_data_s: float,
               obtain_model_s: float | None, build_embedding_s: float | None) -> float:
    """`other_s` = walltime - load - obtain - build (treating an unknown
    obtain/build component as 0). Kept so the four buckets sum to walltime."""
    return float(walltime_s) - float(load_data_s) - float(obtain_model_s or 0.0) - float(build_embedding_s or 0.0)


def _fmt_s(x) -> str:
    return f"{float(x):.2f} s" if x is not None else "n/a"


def resource_lines(s: dict, *, storage_what: str, obtain_what: str, build_what: str) -> list[str]:
    """The standardised `## Resources` markdown block, from a stats dict carrying:
    walltime_s, load_data_s, obtain_model_s (or None), build_embedding_s, other_s,
    model_source, peak_rss_mb, rss_after_load_mb, rss_after_obtain_mb,
    rss_after_build_mb, model_storage_bytes, model_download_s, embedding_bytes,
    device, host."""
    g = s.get
    obtain = g("obtain_model_s")
    obtain_txt = _fmt_s(obtain) if obtain is not None else "n/a (loaded & run atomically upstream — folded into build-embedding)"
    dl = float(g("model_download_s", 0.0) or 0.0)
    return [
        "",
        "## Resources",
        f"- walltime: {float(g('walltime_s', 0.0)):.2f} s   "
        f"(= load {float(g('load_data_s', 0.0)):.2f} s + obtain-model {_fmt_s(obtain)} "
        f"+ build-embedding {_fmt_s(g('build_embedding_s'))} + other {float(g('other_s', 0.0)):.3f} s)",
        f"- load (read input data from disk): {float(g('load_data_s', 0.0)):.2f} s",
        f"- obtain model — {obtain_what} [{g('model_source', '?')}]: {obtain_txt}",
        f"- build embedding ({build_what}): {_fmt_s(g('build_embedding_s'))}",
        f"- other (data prep + parquet write + glue): {float(g('other_s', 0.0)):.3f} s",
        f"- peak RSS: {float(g('peak_rss_mb', 0.0)):.1f} MB   "
        f"(RSS trajectory: after-load {float(g('rss_after_load_mb', 0.0)):.1f} → "
        f"after-model {float(g('rss_after_obtain_mb', 0.0)):.1f} → "
        f"after-embedding {float(g('rss_after_build_mb', 0.0)):.1f} MB)",
        f"- model storage ({storage_what}): {int(g('model_storage_bytes', 0))} bytes",
        f"- model download: {dl:.2f} s ({'cold-cache fetch this run' if dl > 0 else '0 — model already cached on disk / nothing to download'})",
        f"- embedding parquet: {int(g('embedding_bytes', 0))} bytes",
        f"- device: {g('device', 'cpu')}",
        f"- host: {g('host', '?')}",
    ]
