"""
Authored-blueprint "card" for ONE run: a summary panel + minimal grounding stat.

Renders, from the run's structured summary, a compact panel that makes the
authored model concrete: headline counts (axes by kind, modules, genes,
resources, DOIs), the literature-grounding integrity number (fraction of cited
marker genes actually present in the dataset's var_names), and one fully
expanded example axis -- its algebra, its positive/gate/negative gene modules
with the actual cited gene symbols, and the resource citations + DOIs those
genes are grounded in.

Inputs:
    outs/embeddings/<dataset>/<stem>_summary.json     (axes, modules, resources)
    outs/embeddings/<dataset>/<stem>_coverage.json    (per-module genes-present)

Outputs:
    outs/plots/blueprint_card/<dataset>/<stem>_card.{png,pdf}

Usage:
    python scripts/plots/blueprint_card.py --dataset mouse_cortex
    python scripts/plots/blueprint_card.py --dataset mouse_cortex --axis myelin_maturation_program
"""
from __future__ import annotations

import argparse
import json
import os
import tempfile
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "matplotlib"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[2]

DATASET_DEFAULTS = {
    "mouse_cortex":   ("blueprint_gpt_5_5_unrefined_high", 2),
    "human_pancreas": ("blueprint_gpt_5_5_unrefined_high", 1),
    "human_immune":   ("blueprint_gpt_5_5_unrefined_high", 1),
}
# Preferred example axis per dataset: a canonical, textbook cell type (not a
# process/maturation axis), chosen if present in the run; else auto-picked.
DATASET_AXIS_PREF = {
    "mouse_cortex": ["pvalb_program", "sst_program", "vip_program",
                     "astrocyte_identity", "microglia_homeostatic_program"],
}
KIND_COLOR = {"plain": "#4c78a8", "contrast": "#f58518", "gated": "#54a24b"}
ALGEBRA = {
    "plain":    r"$a = s_{\mathrm{pos}}$",
    "contrast": r"$a = s_{\mathrm{pos}} - w\,s_{\mathrm{neg}}$",
    "gated":    r"$a = s_{\mathrm{gate}}\,(s_{\mathrm{pos}} - w\,\overline{s_{\mathrm{neg}}})$",
}


def gene_coverage(summary: dict) -> tuple[int, int, float]:
    """(unique genes present, unique declared, fraction) across all modules."""
    declared, missing = set(), set()
    for m in summary["blueprint"]["modules"]:
        declared.update(m["genes"])
        missing.update(m.get("missing_genes", []))
    present = declared - missing
    frac = len(present) / len(declared) if declared else float("nan")
    return len(present), len(declared), frac


def pick_axis(summary: dict, want: str | None, prefs: list[str]) -> dict:
    axes = summary["blueprint"]["axes"]
    by_name = {a["name"]: a for a in axes}
    if want:
        if want in by_name:
            return by_name[want]
        raise SystemExit(f"--axis {want!r} not found. Choices: {list(by_name)}")
    # Canonical-type preference list for this dataset, first match wins.
    for name in prefs:
        if name in by_name:
            return by_name[name]
    # Default: a gated axis with few negatives whose positive module carries a DOI.
    modules = {m["name"]: m for m in summary["blueprint"]["modules"]}
    res = {r["key"]: r for r in summary["blueprint"]["resources"]}

    def has_doi(mod_name: str) -> bool:
        m = modules.get(mod_name, {})
        return any(res.get(k, {}).get("doi") for k in m.get("resource_keys", []))

    gated = [a for a in axes if a["kind"] == "gated" and has_doi(a["positive"])]
    # Prefer a gated axis with 1-2 negatives (richest algebra to display), then fewest.
    gated.sort(key=lambda a: (0 if 1 <= len(a.get("negatives", [])) <= 2 else 1,
                              len(a.get("negatives", []))))
    return gated[0] if gated else axes[0]


def fmt_genes(genes: list[str], n: int = 9) -> str:
    shown = ", ".join(genes[:n])
    return shown + (f", … (+{len(genes) - n})" if len(genes) > n else "")


def render(summary: dict, axis: dict, cov: tuple, dataset: str, seed: int,
           out_stem: Path) -> None:
    bp = summary["blueprint"]
    modules = {m["name"]: m for m in bp["modules"]}
    res = {r["key"]: r for r in bp["resources"]}
    kc = bp["n_axes_by_kind"]
    n_present, n_decl, frac = cov
    W = 6.0   # tightened to the content width (less right whitespace; larger when scaled)

    def clean(name: str) -> str:
        for suf in ("_interneuron", "_program", "_identity", "_cell"):
            name = name.replace(suf, "")
        return name

    # ---- Two-pass layout: collect (x_frac, inches-from-top, text, kwargs) ----
    ops: list[tuple] = []
    lines: list[float] = []      # y-in-inches of horizontal dividers
    c = 0.30                      # cursor, inches from the top

    def T(x, s, adv, **kw):
        nonlocal c
        ops.append((x, c, s, kw)); c += adv

    def P(x, s, **kw):           # place without advancing (same line)
        ops.append((x, c, s, kw))

    T(0.04, "Authored blueprint", 0.30, fontsize=14, weight="bold")
    T(0.04, f"{dataset} · literature-grounded · GPT-5.5 high · seed{seed}", 0.26,
      fontsize=8.5, color="0.35")
    counts = (f"{bp['n_axes']} axes ({kc.get('plain',0)} plain · "
              f"{kc.get('contrast',0)} contrast · {kc.get('gated',0)} gated)  ·  "
              f"{bp['n_modules']} modules  ·  {bp['n_unique_genes']} genes  ·  "
              f"{bp['n_resources']} resources  ·  {bp['n_dois_unique']} DOIs")
    T(0.04, counts, 0.24, fontsize=8.2)
    T(0.04, f"literature grounding: {n_present}/{n_decl} cited marker genes present "
            f"in the data ({100 * frac:.1f}%)", 0.20, fontsize=8.4, color="#1a6b3a")
    lines.append(c); c += 0.18

    T(0.04, "Example axis", 0.26, fontsize=10, weight="bold")
    P(0.55, f"[{axis['kind']}]   {ALGEBRA[axis['kind']]}", fontsize=9,
      color=KIND_COLOR.get(axis["kind"], "k"))
    T(0.04, axis["name"], 0.34, fontsize=10.5, family="monospace",
      weight="bold", color=KIND_COLOR.get(axis["kind"], "k"))

    pos_m = modules[axis["positive"]]
    # `detail_is_genes` flags rows whose detail is a list of gene symbols
    # (gate/pos) vs sibling *module* names (neg). Gene lists are prefixed
    # "Genes:" and italicised per nomenclature; module lists are prefixed
    # "Modules:" and stay upright (capitalised program names, not symbols).
    def cap(s: str) -> str:
        return s[:1].upper() + s[1:]

    rows = []
    if axis.get("gate"):
        rows.append(("gate", clean(axis["gate"]), fmt_genes(modules[axis["gate"]]["genes"]), True))
    rows.append(("pos", clean(axis["positive"]), fmt_genes(pos_m["genes"]), True))
    if axis.get("negatives"):
        negs = [cap(clean(n)) for n in axis["negatives"]]
        shown = ", ".join(negs[:6]) + (f", … (+{len(negs) - 6})" if len(negs) > 6 else "")
        rows.append(("neg", f"{len(negs)} sibling programs", shown, False))
    for role, mod, detail, detail_is_genes in rows:
        P(0.06, role, fontsize=8.3, color="0.45", family="monospace")
        T(0.15, mod, 0.17, fontsize=8.3, family="monospace", weight="bold")
        T(0.15, detail, 0.21, fontsize=7.4, color="0.25",
          _prefix="Genes: " if detail_is_genes else "Modules: ",
          _italic_body=detail_is_genes)

    c += 0.04
    # Distinct resource keys used by the example axis (positive, gate, negatives).
    res_keys = list(pos_m.get("resource_keys", []))
    if axis.get("gate"):
        res_keys += modules.get(axis["gate"], {}).get("resource_keys", [])
    for nm in axis.get("negatives", []):
        res_keys += modules.get(nm, {}).get("resource_keys", [])
    seen, ordered_keys = set(), []
    for k in res_keys:
        if k not in seen:
            seen.add(k); ordered_keys.append(k)

    P(0.06, "cites", fontsize=8.3, color="0.45", family="monospace")
    for k in ordered_keys[:4]:
        r = res.get(k, {})
        tail = f"doi:{r['doi']}" if r.get("doi") else "(database identifier)"
        T(0.15, f"• {r.get('citation', k)[:62]}", 0.155, fontsize=7.3)
        T(0.165, tail, 0.20, fontsize=7.0, color="#1a4fa0", family="monospace")

    H = c + 0.18                  # total figure height in inches

    # ---- Draw ----
    fig = plt.figure(figsize=(W, H))
    ax = fig.add_axes([0, 0, 1, 1]); ax.axis("off")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.add_patch(plt.Rectangle((0.012, 0.014), 0.976, 0.972, fill=False,
                               edgecolor="0.6", lw=1.0))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    for x, ci, s, kw in ops:
        kw = dict(kw)
        prefix = kw.pop("_prefix", None)
        italic_body = kw.pop("_italic_body", False)
        y = 1 - ci / H
        if prefix is None:
            ax.text(x, y, s, va="top", **kw)
            continue
        # Upright label, then the list placed right after it (measured width).
        t = ax.text(x, y, prefix, va="top", **kw)
        bb = t.get_window_extent(renderer)
        dx = inv.transform((bb.x1, 0))[0] - inv.transform((bb.x0, 0))[0]
        body_kw = dict(kw, style="italic") if italic_body else kw
        ax.text(x + dx, y, s, va="top", **body_kw)
    for ci in lines:
        ax.plot([0.04, 0.96], [1 - ci / H] * 2, color="0.8", lw=0.8)

    out_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        fig.savefig(out_stem.with_suffix(f".{ext}"), dpi=300, facecolor="white")
    plt.close(fig)
    print(f"[ok] wrote {out_stem}.{{png,pdf}}")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Authored-blueprint summary card.")
    p.add_argument("--dataset", choices=tuple(DATASET_DEFAULTS), default="mouse_cortex")
    p.add_argument("--cohort", default=None)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--axis", default=None, help="Axis name to expand (default: auto).")
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cohort, seed = DATASET_DEFAULTS[args.dataset]
    cohort = args.cohort or cohort
    seed = args.seed if args.seed is not None else seed
    stem = f"{cohort}_{args.dataset}_seed{seed}"
    emb_dir = REPO_ROOT / "outs" / "embeddings" / args.dataset
    summary = json.loads((emb_dir / f"{stem}_summary.json").read_text())

    out_stem = REPO_ROOT / "outs" / "plots" / "blueprint_card" / args.dataset / f"{stem}_card"
    if out_stem.with_suffix(".png").exists() and not args.overwrite:
        raise FileExistsError(f"Refusing to overwrite {out_stem}.png (pass --overwrite).")

    axis = pick_axis(summary, args.axis, DATASET_AXIS_PREF.get(args.dataset, []))
    render(summary, axis, gene_coverage(summary), args.dataset, seed, out_stem)



if __name__ == "__main__":
    main()
