"""
Build PanglaoDB cell-type marker .gmt files for the UCell / expiMap baselines.

Author-faithful gene-set sourcing: the original UCell method scores
user-supplied *cell-type marker* signatures, and expiMap uses PanglaoDB
marker sets (for cell-type identification) alongside
Reactome pathways (for molecular processes). This script produces, per benchmark
dataset:

  - panglao_<dataset>.gmt              — PanglaoDB cell-type markers (UCell input)
  - reactome_panglao_<dataset>.gmt     — Reactome ∪ PanglaoDB (expiMap input)

Gene sets are taken *off the shelf*: every PanglaoDB cell type for the dataset's
organ(s), with no hand-mapping to the dataset's label vocabulary (no
cherry-picking). Organ map:
  mouse_cortex   -> Brain
  human_pancreas -> Pancreas
  human_immune   -> Immune system + Blood

PanglaoDB lists human HUGO symbols. We keep them uppercase for human datasets and
title-case them for mouse (e.g. CTRB1 -> Ctrb1), the standard mouse-ortholog
convention, so UCell's *case-sensitive* symbol match hits mouse var_names
(e.g. Pou5f1). expiMap uppercases both sides, so the merged file's case is moot
for it.

Inputs (downloaded + cached on first run):
    models/genesets/panglaodb/PanglaoDB_markers_27_Mar_2020.tsv.gz
    models/genesets/<organism>/<reactome>.symbols.gmt   (already cached / MSigDB)

Outputs:
    models/genesets/<organism>/panglao_<dataset>.gmt
    models/genesets/<organism>/reactome_panglao_<dataset>.gmt

Usage:
    python scripts/IO/build_panglao_gmt.py --dataset human_pancreas
    python scripts/IO/build_panglao_gmt.py --dataset mouse_cortex --canonical-only
"""
from __future__ import annotations

import argparse
import gzip
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]

PANGLAO_URL = ("https://panglaodb.se/markers/"
               "PanglaoDB_markers_27_Mar_2020.tsv.gz")
# PanglaoDB blocks bare user-agents (403); present a browser UA.
_UA = ("Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) "
       "AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36")

MSIGDB_URL = ("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
              "{dbver}/{category}.v{dbver}.symbols.gmt")

# dataset -> organism, PanglaoDB species code, and organ(s) to pull (off-the-shelf).
DATASET_CONFIG = {
    "mouse_cortex":   {"organism": "mouse", "species": "Mm",
                       "organs": ("Brain",),
                       "reactome": ("m2.cp.reactome", "2024.1.Mm")},
    "human_pancreas": {"organism": "human", "species": "Hs",
                       "organs": ("Pancreas",),
                       "reactome": ("c2.cp.reactome", "2024.1.Hs")},
    "human_immune":   {"organism": "human", "species": "Hs",
                       "organs": ("Immune system", "Blood"),
                       "reactome": ("c2.cp.reactome", "2024.1.Hs")},
}

# PanglaoDB column order (1-indexed in the file; 0-indexed here).
COL_SPECIES, COL_SYMBOL, COL_CELLTYPE, COL_CANONICAL, COL_ORGAN = 0, 1, 2, 7, 9


def _download_panglao(cache: Path) -> Path:
    cache.mkdir(parents=True, exist_ok=True)
    dst = cache / "PanglaoDB_markers_27_Mar_2020.tsv.gz"
    if dst.exists():
        return dst
    req = urllib.request.Request(PANGLAO_URL, headers={"User-Agent": _UA})
    try:
        with urllib.request.urlopen(req, timeout=60) as r, open(dst, "wb") as f:
            f.write(r.read())
    except Exception as e:
        if dst.exists():
            dst.unlink()
        raise RuntimeError(
            f"Failed to download PanglaoDB markers from {PANGLAO_URL}: {e}. "
            "Download it manually (browser) to "
            f"{dst} and re-run."
        ) from e
    return dst


def _reactome_gmt(cfg: dict) -> Path:
    """Locate (or download) the organism's cached Reactome .gmt."""
    cache = REPO_ROOT / "models" / "genesets" / cfg["organism"]
    category, dbver = cfg["reactome"]
    gmt = cache / f"{category}.v{dbver}.symbols.gmt"
    if gmt.exists():
        return gmt
    cache.mkdir(parents=True, exist_ok=True)
    url = MSIGDB_URL.format(dbver=dbver, category=category)
    urllib.request.urlretrieve(url, gmt)
    return gmt


def _case(sym: str, organism: str) -> str:
    """Map a PanglaoDB HUGO symbol to the dataset's symbol convention."""
    sym = sym.strip()
    return sym.upper() if organism == "human" else sym.capitalize()


def build_panglao(tsv_gz: Path, cfg: dict, canonical_only: bool) -> dict[str, list[str]]:
    organs = set(cfg["organs"])
    species = cfg["species"]
    sigs: dict[str, list[str]] = {}
    seen: dict[str, set[str]] = {}
    with gzip.open(tsv_gz, "rt") as f:
        next(f)  # header
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) <= COL_ORGAN:
                continue
            if p[COL_ORGAN] not in organs:
                continue
            if species not in p[COL_SPECIES].split():  # 'Mm Hs' applies to both
                continue
            if canonical_only and p[COL_CANONICAL].strip() != "1":
                continue
            cell_type = p[COL_CELLTYPE].strip()
            gene = _case(p[COL_SYMBOL], cfg["organism"])
            if not gene:
                continue
            bucket = seen.setdefault(cell_type, set())
            if gene not in bucket:
                bucket.add(gene)
                sigs.setdefault(cell_type, []).append(gene)
    return sigs


def write_gmt(sigs: dict[str, list[str]], path: Path, source: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    for name in sorted(sigs):
        # GMT: name \t description \t gene1 \t gene2 ...
        # Collapse whitespace in the program name to underscores: scArches'
        # add_annotations splits term names on whitespace, which would corrupt
        # and collide multi-word PanglaoDB cell-type names (e.g. "Pancreatic
        # progenitor cells" and "Pancreatic stellate cells" both -> "PANCREATIC").
        safe = "_".join(name.split())
        lines.append("\t".join([safe, source, *sigs[name]]))
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description="Build PanglaoDB cell-type marker .gmt files.")
    ap.add_argument("--dataset", choices=tuple(DATASET_CONFIG), required=True)
    ap.add_argument("--canonical-only", action="store_true",
                    help="Keep only PanglaoDB 'canonical marker' genes (col 8 == 1).")
    ap.add_argument("--min-genes", type=int, default=2,
                    help="Drop signatures with fewer than this many genes (default 2).")
    ap.add_argument("--markers", type=str,
                    help="Path to a local PanglaoDB markers .tsv.gz (skips download).")
    args = ap.parse_args()

    cfg = DATASET_CONFIG[args.dataset]
    organism = cfg["organism"]

    if args.markers:
        tsv_gz = Path(args.markers).resolve()
        if not tsv_gz.exists():
            raise FileNotFoundError(f"--markers file not found: {tsv_gz}")
    else:
        tsv_gz = _download_panglao(REPO_ROOT / "models" / "genesets" / "panglaodb")

    sigs = build_panglao(tsv_gz, cfg, canonical_only=args.canonical_only)
    sigs = {k: v for k, v in sigs.items() if len(v) >= args.min_genes}
    if not sigs:
        raise RuntimeError(
            f"No PanglaoDB signatures for organs={cfg['organs']} species={cfg['species']}. "
            "Check the markers file / column layout."
        )

    out_dir = REPO_ROOT / "models" / "genesets" / organism
    source_tag = "PanglaoDB:" + "+".join(cfg["organs"])
    panglao_gmt = out_dir / f"panglao_{args.dataset}.gmt"
    write_gmt(sigs, panglao_gmt, source_tag)

    # Merged Reactome ∪ PanglaoDB for expiMap (concatenate raw .gmt text).
    reactome_gmt = _reactome_gmt(cfg)
    merged_gmt = out_dir / f"reactome_panglao_{args.dataset}.gmt"
    merged_text = reactome_gmt.read_text().rstrip("\n") + "\n" + panglao_gmt.read_text()
    merged_gmt.write_text(merged_text)

    n_genes = len({g for v in sigs.values() for g in v})
    n_reactome = sum(1 for ln in reactome_gmt.read_text().splitlines() if ln.strip())
    print(f"[ok] {args.dataset}: {len(sigs)} PanglaoDB cell-type signatures "
          f"({n_genes} unique genes, organs={cfg['organs']}, species={cfg['species']})")
    print(f"[ok] wrote {panglao_gmt}")
    print(f"[ok] wrote {merged_gmt}  (+{n_reactome} Reactome programs)")


if __name__ == "__main__":
    main()
