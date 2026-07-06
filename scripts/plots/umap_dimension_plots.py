"""
Per-dimension UMAP plots for ONE selected embedding run.

Takes the same 2D UMAP layout produced by scripts/plots/umap_plots.py (the coords
parquet is reused if present, otherwise recomputed with identical settings and
cached) and recolours it by each embedding *dimension's activity value* — one
panel per dim_* column, plus a small-multiples grid of all of them.

Panel titles: for interpretable models (blueprints), the axis name + type (e.g.
"Endocrine identity (plain)") parsed from the run's <stem>_report.md companion;
for everything else, "Dimension <i>".

Inputs:
    outs/embeddings/<dataset>/<stem>.parquet                       (cell_id + dim_* columns)
    outs/plots/umap/<dataset>/<model>/<stem>_umap_coords.parquet   (reused if present)
    outs/embeddings/<dataset>/<stem>_report.md                     (optional; axis-title source)

Outputs:
    outs/plots/umap/<dataset>/<model>/<stem>_umap_coords.parquet   (written if recomputed;
                                                                    shared with umap_plots.py)
    outs/plots/umap_axes/<dataset>/<model>/<stem>_dims/<stem>_dim<NN>[__<axis-slug>].png / .pdf
    outs/plots/umap_axes/<dataset>/<model>/<stem>_dims_grid.png / .pdf   (contact sheet of all dims)

Usage:
    python scripts/plots/umap_dimension_plots.py --dataset human_pancreas \
        --model blueprint_gpt_5_5_high_count --seed 24
    python scripts/plots/umap_dimension_plots.py --dataset human_pancreas \
        --input outs/embeddings/human_pancreas/pca_human_pancreas_seed0.parquet
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import tempfile
import textwrap
import time
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import umap


REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_brain", "human_pancreas", "mouse_gastrulation", "mouse_cortex", "human_immune")

# Per-dataset scatter point size — same values umap_plots.py uses (it keeps them in
# DEFAULT_PARAMS; only point_size is needed here).
POINT_SIZE_BY_DATASET = {
    "mouse_brain": 14.0,        # ~1.8k cells: prominent dots
    "human_pancreas": 4.0,      # ~16k cells: small dots to reduce overdraw
    "mouse_gastrulation": 1.5,  # ~116k cells: tiny dots
    "mouse_cortex": 1.5,        # ~45k cells: small dots
    "human_immune": 2.0,        # ~33k cells: small dots
}

PLOT_RC = {
    "font.family": "DejaVu Sans",
    "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 10,
    "savefig.facecolor": "white", "savefig.transparent": False,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
}


# ---------------------------------------------------------------- run resolution

def parse_stem(stem: str, dataset: str) -> tuple[str, int | None]:
    """Recover (model, seed) from '<model>_<dataset>[_seed<N>]' — same logic as umap_plots.py."""
    if f"_{dataset}" not in stem:
        raise ValueError(f"Stem {stem!r} does not contain dataset {dataset!r}.")
    head, _, tail = stem.partition(f"_{dataset}")
    if not head:
        raise ValueError(f"Empty model name parsed from stem {stem!r}.")
    if not tail:
        return head, None
    if tail.startswith("_seed"):
        return head, int(tail[len("_seed"):])
    return head + tail, None


def resolve_parquet(args: argparse.Namespace, embeddings_dir: Path) -> Path:
    """Pick the single embedding parquet for the selected run."""
    if args.input is not None:
        path = args.input if args.input.is_absolute() else (REPO_ROOT / args.input)
        if not path.exists():
            raise FileNotFoundError(f"--input not found: {path}")
        return path.resolve()
    if not args.model:
        raise SystemExit("Pass either --input PATH or --model SLUG [--seed N].")
    if args.seed is not None:
        path = embeddings_dir / f"{args.model}_{args.dataset}_seed{args.seed}.parquet"
        if not path.exists():
            raise FileNotFoundError(f"No such embedding: {path}")
        return path.resolve()
    # No seed given: accept a seedless parquet, or a unique seeded one.
    seedless = embeddings_dir / f"{args.model}_{args.dataset}.parquet"
    if seedless.exists():
        return seedless.resolve()
    seeded = sorted(embeddings_dir.glob(f"{args.model}_{args.dataset}_seed*.parquet"))
    if not seeded:
        raise FileNotFoundError(
            f"No embedding for model {args.model!r} in {embeddings_dir}."
        )
    if len(seeded) > 1:
        seeds = sorted(parse_stem(p.stem, args.dataset)[1] for p in seeded)
        raise SystemExit(
            f"model {args.model!r} has multiple seeds ({', '.join(map(str, seeds))}); pass --seed N."
        )
    return seeded[0].resolve()


# ---------------------------------------------------------------- UMAP (kept in sync with umap_plots.py)

def compute_or_load_umap(parquet_path: Path, coords_path: Path, params: dict, overwrite: bool):
    """Return (cell_ids, coords, dim_cols, df_sorted, cached). Coords are bit-for-bit
    shareable with umap_plots.py: identical row order (mergesort on cell_id), identical
    umap.UMAP(...) call, identical cache schema (cell_id, umap_0, umap_1)."""
    df = pd.read_parquet(parquet_path)
    if "cell_id" not in df.columns:
        raise ValueError(f"{parquet_path} has no 'cell_id' column.")
    dim_cols = sorted(
        (c for c in df.columns if c.startswith("dim_")),
        key=lambda n: int(n.split("_")[1]),
    )
    if not dim_cols:
        raise ValueError(f"{parquet_path} has no dim_* columns.")
    df = df.sort_values("cell_id", kind="mergesort").reset_index(drop=True)
    cell_ids = df["cell_id"].astype(str)

    if coords_path.exists() and not overwrite:
        cached = pd.read_parquet(coords_path)
        if cached["cell_id"].astype(str).tolist() != cell_ids.tolist():
            raise ValueError(f"Cached coords {coords_path} disagree with {parquet_path}.")
        return cell_ids, cached[["umap_0", "umap_1"]].to_numpy(), dim_cols, df, True

    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=params["n_neighbors"], min_dist=params["min_dist"], metric=params["metric"],
        init="spectral",
        random_state=params["random_state"], transform_seed=params["random_state"],
        low_memory=True,
    )
    coords = reducer.fit_transform(df[dim_cols].to_numpy(dtype=float))
    coords_path.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"cell_id": cell_ids, "umap_0": coords[:, 0], "umap_1": coords[:, 1]}).to_parquet(
        coords_path, index=False
    )
    return cell_ids, coords, dim_cols, df, False


# ---------------------------------------------------------------- axis titles

_TABLE_HEADER = "## Embedding dimensions"


def _parse_report_axis_table(report_path: Path, n_dims: int) -> list[tuple[str, str]] | None:
    """(axis name, axis kind) for each row of the '## Embedding dimensions' table, in dim order.
    Kind is "" if the table has no Kind column."""
    if not report_path.exists():
        return None
    lines = report_path.read_text().splitlines()
    try:
        start = next(i for i, ln in enumerate(lines) if ln.strip() == _TABLE_HEADER)
    except StopIteration:
        return None
    rows: list[tuple[str, str]] = []
    seen_separator = False
    for ln in lines[start + 1:]:
        s = ln.strip()
        if not s.startswith("|"):
            if rows or seen_separator:
                break
            continue
        cells = [c.strip() for c in s.strip("|").split("|")]
        if not cells:
            continue
        if set(cells[0]) <= {"-", ":"}:  # the |---|---| separator row
            seen_separator = True
            continue
        if cells[0].lower() == "axis":     # the header row
            continue
        rows.append((cells[0], cells[1] if len(cells) > 1 else ""))
    if len(rows) != n_dims:
        return None
    return rows


def _pretty(name: str) -> str:
    s = name.replace("_", " ").strip()
    return s[:1].upper() + s[1:] if s else name


def _slug(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", name.lower()).strip("_")


def resolve_titles(
    parquet_path: Path, n_dims: int, override: str | None
) -> tuple[list[str], list[str | None]]:
    """Return (display_titles, filename_slugs). slug is None for non-interpretable dims."""
    if override:
        candidate = Path(override)
        if candidate.exists():
            titles = [ln.strip() for ln in candidate.read_text().splitlines() if ln.strip()]
        else:
            titles = [t.strip() for t in override.split(",") if t.strip()]
        if len(titles) != n_dims:
            raise SystemExit(
                f"--axis-titles has {len(titles)} entries but embedding has {n_dims} dims."
            )
        return titles, [_slug(t) for t in titles]

    report_path = parquet_path.with_name(parquet_path.stem + "_report.md")
    rows = _parse_report_axis_table(report_path, n_dims)
    if rows is not None:
        titles = [f"{_pretty(name)} ({kind})" if kind else _pretty(name) for name, kind in rows]
        return titles, [_slug(name) for name, _ in rows]

    return [f"Dimension {i}" for i in range(n_dims)], [None] * n_dims


# ---------------------------------------------------------------- styling / rendering

def style_axes(ax, title: str, *, title_fontsize: float | None = None) -> None:
    ax.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("UMAP 1", labelpad=6)
    ax.set_ylabel("UMAP 2", labelpad=6)
    if title:
        ax.set_title(title, pad=6, fontsize=title_fontsize)


def _upper_half_cmap(base_name: str):
    """The white->high half of a diverging colormap (e.g. RdBu_r -> white->red),
    for one-sided axes whose negative pole carries no meaning."""
    import matplotlib.colors as mcolors
    base = plt.get_cmap(base_name)
    return mcolors.LinearSegmentedColormap.from_list(
        f"{base_name}_upper", base(np.linspace(0.5, 1.0, 256)))


def _draw_dim(ax, coords, values: np.ndarray, params: dict, has_neg: bool = True):
    """Scatter one dimension; high-activation cells drawn last (on top). Returns the mappable.

    Axes with negative modules (contrast / gated-with-negatives) are signed and use
    the full diverging colormap centered at 0. Axes without negatives (plain,
    gated-with-only-a-gate) have no meaningful negative pole, so they use only the
    upper (white->high) half of that colormap with the scale floored at 0.
    """
    vmax = params["vmax"]
    if vmax is None:
        finite = np.abs(values[np.isfinite(values)])
        vmax = float(np.nanpercentile(finite, 99)) if finite.size else 1.0
    vmax = max(vmax, 1e-6)
    if has_neg:
        cmap, vmin = params["cmap"], -vmax
        order = np.argsort(np.abs(np.nan_to_num(values)), kind="mergesort")
    else:
        cmap, vmin = params["seq_cmap"], 0.0
        order = np.argsort(np.nan_to_num(values), kind="mergesort")  # highest on top
    sc = ax.scatter(
        coords[order, 0], coords[order, 1],
        c=values[order], cmap=cmap, vmin=vmin, vmax=vmax,
        s=params["point_size"], alpha=params["alpha"], linewidths=0, rasterized=True,
    )
    ax.set_aspect("equal", adjustable="datalim")
    ax.margins(0.03)
    return sc


def render_dim_panel(coords, values, title: str, out_png: Path, out_pdf: Path, params: dict,
                     has_neg: bool = True) -> None:
    panel = params["panel_size"]
    fig, ax = plt.subplots(figsize=(panel + 0.9, panel))
    sc = _draw_dim(ax, coords, values, params, has_neg)
    style_axes(ax, title)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.03, shrink=0.85)
    cbar.set_label("dimension activity", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    fig.tight_layout(pad=0.4)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=params["png_dpi"], bbox_inches="tight", facecolor="white")
    fig.savefig(out_pdf, dpi=params["pdf_dpi"], bbox_inches="tight", facecolor="white")
    plt.close(fig)


def render_grid(coords, dim_values: list[np.ndarray], panel_titles: list[str],
                suptitle: str, out_png: Path, out_pdf: Path, params: dict,
                has_neg: list[bool] | None = None) -> None:
    n = len(dim_values)
    if has_neg is None:
        has_neg = [True] * n
    ncols = math.ceil(math.sqrt(n))
    nrows = math.ceil(n / ncols)
    tile = params["grid_tile"]
    # Wrap long titles so they stay within their panel (max 2 lines, rest elided);
    # blueprint axis names are otherwise wide enough to spill across neighbours.
    wrapped = [textwrap.fill(t, width=26, max_lines=2, placeholder="…") for t in panel_titles]
    n_lines = max((t.count("\n") + 1 for t in wrapped), default=1)
    fig, axes = plt.subplots(
        nrows, ncols,
        figsize=(tile * ncols, (tile + 0.24 * n_lines) * nrows + 0.6),
        constrained_layout=True,
    )
    fig.get_layout_engine().set(h_pad=0.05, w_pad=0.03, hspace=0.045, wspace=0.0)
    axes = np.atleast_1d(axes).ravel()
    for i, ax in enumerate(axes):
        if i >= n:
            ax.axis("off")
            continue
        _draw_dim(ax, coords, dim_values[i], params, has_neg[i])
        ax.set_facecolor("white")
        for spine in ax.spines.values():
            spine.set_visible(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(wrapped[i], fontsize=8, pad=2.5, linespacing=0.92)
    fig.suptitle(suptitle, fontsize=10)
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=params["png_dpi"], facecolor="white")
    fig.savefig(out_pdf, dpi=params["pdf_dpi"], facecolor="white")
    plt.close(fig)


# ---------------------------------------------------------------- CLI / main

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Per-dimension UMAP plots for one embedding run.")
    p.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    p.add_argument("--input", type=Path, help="Embedding parquet for the selected run.")
    p.add_argument("--model", type=str, help="Model slug (with --seed) instead of --input.")
    p.add_argument("--seed", type=int, default=None, help="Seed for --model (omit for seedless / unique).")
    p.add_argument("--axis-titles", type=str, default=None,
                   help="Override panel titles: comma-separated list, or a path to one-title-per-line file.")
    p.add_argument("--n-neighbors", type=int, default=30)
    p.add_argument("--min-dist", type=float, default=0.15)
    p.add_argument("--metric", type=str, default="euclidean")
    p.add_argument("--random-state", type=int, default=42)
    p.add_argument("--cmap", type=str, default="RdBu_r", help="Diverging colormap, centered at 0.")
    p.add_argument("--vmax", type=float, default=None,
                   help="Fixed symmetric colour limit (+/-) across all dims; default = per-dim 99th pct of |value|.")
    p.add_argument("--point-size", type=float, default=None, help="Override per-dataset point size.")
    p.add_argument("--alpha", type=float, default=0.9)
    p.add_argument("--panel-size", type=float, default=4.0)
    p.add_argument("--grid-tile", type=float, default=2.2, help="Per-panel size in the combined grid.")
    p.add_argument("--max-grid-dims", type=int, default=64, help="Skip the grid above this many dims.")
    p.add_argument("--png-dpi", type=int, default=600)
    p.add_argument("--pdf-dpi", type=int, default=600)
    p.add_argument("--grid-only", action="store_true",
                   help="Only (re)write the combined grid, reusing the cached UMAP coords; skip per-dim panels.")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    matplotlib.rcParams.update(PLOT_RC)

    embeddings_dir = REPO_ROOT / "outs" / "embeddings" / args.dataset
    parquet_path = resolve_parquet(args, embeddings_dir)
    stem = parquet_path.stem
    model, seed = parse_stem(stem, args.dataset)

    # Coords cache lives in umap_plots.py's tree so the two scripts share it; the
    # per-dimension figures go in a separate umap_axes/ tree.
    coords_path = REPO_ROOT / "outs" / "plots" / "umap" / args.dataset / model / f"{stem}_umap_coords.parquet"
    out_dir = REPO_ROOT / "outs" / "plots" / "umap_axes" / args.dataset / model
    dims_dir = out_dir / f"{stem}_dims"
    grid_png = out_dir / f"{stem}_dims_grid.png"
    grid_pdf = out_dir / f"{stem}_dims_grid.pdf"

    existing = [p for p in (grid_png, grid_pdf) if p.exists()]
    if not args.grid_only and dims_dir.exists():
        existing += sorted(dims_dir.glob("*.png")) + sorted(dims_dir.glob("*.pdf"))
    if existing and not args.overwrite:
        raise FileExistsError(
            "Refusing to overwrite (pass --overwrite):\n  " + "\n  ".join(str(p) for p in existing)
        )
    if args.overwrite and not args.grid_only and dims_dir.exists():
        for p in list(dims_dir.glob("*.png")) + list(dims_dir.glob("*.pdf")):
            p.unlink()

    umap_params = {k: getattr(args, k) for k in ("n_neighbors", "min_dist", "metric", "random_state")}
    t0 = time.perf_counter()
    cell_ids, coords, dim_cols, df, cached = compute_or_load_umap(
        parquet_path, coords_path, umap_params, args.overwrite and not args.grid_only
    )
    n_dims = len(dim_cols)
    umap_msg = "cached" if cached else f"{time.perf_counter() - t0:.1f}s"
    seed_repr = f"seed={seed}" if seed is not None else "seed=-"
    print(f"[ok] {parquet_path.name}: model={model} {seed_repr} n_cells={len(cell_ids)} "
          f"n_dims={n_dims} umap={umap_msg}")

    titles, slugs = resolve_titles(parquet_path, n_dims, args.axis_titles)
    interpretable = any(s is not None for s in slugs)
    print(f"[ok] titles: source={'report/override' if interpretable else 'dim index'}")

    point_size = args.point_size if args.point_size is not None else POINT_SIZE_BY_DATASET[args.dataset]
    plot_params = {
        "cmap": args.cmap, "seq_cmap": _upper_half_cmap(args.cmap), "vmax": args.vmax,
        "point_size": point_size, "alpha": args.alpha,
        "panel_size": args.panel_size, "grid_tile": args.grid_tile,
        "png_dpi": args.png_dpi, "pdf_dpi": args.pdf_dpi,
    }

    # Per-dim signedness from the blueprint: axes with no negative modules use the
    # one-sided (white->high) colormap. Defaults to signed if no summary is found.
    has_neg = [True] * n_dims
    summary_path = parquet_path.with_name(stem + "_summary.json")
    if summary_path.exists():
        axes_meta = json.loads(summary_path.read_text()).get("blueprint", {}).get("axes", [])
        for i in range(min(n_dims, len(axes_meta))):
            has_neg[i] = len(axes_meta[i].get("negatives", [])) > 0

    dim_values = [df[col].to_numpy(dtype=float) for col in dim_cols]
    if not args.grid_only:
        width = max(2, len(str(n_dims - 1)))
        for i, col in enumerate(dim_cols):
            suffix = f"__{slugs[i]}" if slugs[i] else ""
            out_png = dims_dir / f"{stem}_dim{i:0{width}d}{suffix}.png"
            out_pdf = dims_dir / f"{stem}_dim{i:0{width}d}{suffix}.pdf"
            t1 = time.perf_counter()
            render_dim_panel(coords, dim_values[i], titles[i], out_png, out_pdf, plot_params,
                             has_neg[i])
            print(f"[ok]   dim {i}: '{titles[i]}' render={time.perf_counter() - t1:.1f}s "
                  f"-> {out_png.relative_to(REPO_ROOT)}")

    short_titles = [titles[i] if slugs[i] else f"dim {i}" for i in range(n_dims)]
    if n_dims > args.max_grid_dims:
        print(f"[skip] grid: n_dims={n_dims} > --max-grid-dims={args.max_grid_dims}")
    else:
        t2 = time.perf_counter()
        render_grid(coords, dim_values, short_titles, stem, grid_png, grid_pdf, plot_params, has_neg)
        print(f"[ok]   grid: {n_dims} panels render={time.perf_counter() - t2:.1f}s "
              f"-> {grid_png.relative_to(REPO_ROOT)}")


if __name__ == "__main__":
    main()
