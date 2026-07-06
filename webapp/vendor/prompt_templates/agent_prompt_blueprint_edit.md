# Task

Design an embedding for the cell-type/state landscape of the
biological scenario in `scenario_config.jsonc`. The embedding must be
derived from knowledge of the scientific literature for that scenario,
instantiated by editing ONLY the `build_blueprint()` function in
`embedder_blueprint.py`. Do not create or modify any other file.

## Scenario

`scenario_config.jsonc` defines the scenario. Fields:

| Field | Type | `null` means |
| --- | --- | --- |
| `species` | array of str (≥1) | — (required) |
| `tissue` | array of str | unspecified |
| `celltype` | array of str | all cell types in the tissue |
| `age` | str | unspecified |
| `sex` | str | unspecified |
| `condition` | array of str | baseline / wild-type / healthy |
| `notes` | str | unspecified |

Rules:

- Treat `null` as unspecified; do not invent defaults to fill the gap.
- The embedding must cover the cell types and states biologically
  expected within the species / tissue / age / sex / condition
  intersection. Do not include anything outside that scope.
- `notes` is a *preference within constraints*. It may shift module
  emphasis, axis-type choice, dimensionality, and the **structural
  shape** of the embedding — e.g. a time-forward gating tree when the
  notes describe developmental ordering, an axis topology that
  resolves a specific subtype group, a contrast axis aligned with a
  signaling pathway. It does NOT relax the literature-grounding rules
  below. If IF>10 literature does not support the requested focus or
  structure, omit it rather than invent ungrounded modules.

## Method

### Sources

Every gene module must be supported by ≥1 resource declared via
`add_resource(...)` and referred to by key in
`add_gene_module(..., resources=[<key>, ...])`. Axis builders do not
cite resources; literature grounding lives at the gene-module level.

Resource kinds:

- `kind="paper"` — peer-reviewed paper, journal impact factor **> 10**.
  `doi` is required (canonical handle).
- `kind="atlas"` — reference atlas (Allen Brain Atlas, Tabula Muris,
  Human Cell Atlas, …), cited via its companion IF>10 paper. `doi` is
  required.
- `kind="database"` — curated pathway / gene-set DB (KEGG, MSigDB,
  Reactome, GO, WikiPathways, …). `doi` is optional here — the
  gene-set identifier inside `citation` is the canonical handle.

When the scenario admits it, ground the embedding in the breadth of
the evidence base — papers, atlases, and databases — rather than
leaning on one kind.

**Database fallback (recommended).** When a module's program maps to
a curated gene-set entry, co-cite a `kind="database"` resource
alongside the paper(s). Database pages remain verifiable when the
cited papers don't.

### Axis design and dimensionality

Use the `EmbedderBlueprint` API as documented in the per-method
docstrings of `embedder_blueprint.py`. The number of axes you register
determines the embedding dimensionality — let the biology drive it,
do not pick a fixed number out of habit. Use `notes` (when non-null)
to inform axis-type selection **and the overall structural shape** of
the embedding (e.g. gating-tree topology for a developmental ordering,
contrast axes for opposing programs).

## Constraints

- EDIT ONLY `build_blueprint()`. No new helpers, classes, methods, or
  top-level definitions; no new imports beyond what the skeleton has.
- All `add_resource` / `add_gene_module` arguments are inline string /
  list literals — no aliasing through local variables. The harness
  parses the file statically.
- Module-level `resources=[...]` must be non-empty, and every key must
  refer to a resource registered earlier in `build_blueprint`.
- DOIs are in the format `10.xxxx/yyy`. Bare URLs and `doi:` prefixes
  are rejected. Declare resources by `citation` and `doi` only — no
  URL field, no URL in the `citation` text.
- Declare resources you have specific grounds to believe contain the
  gene lists you cite. Do not paste a citation from memory without
  that grounding.
- For `kind="database"` resources, include the gene-set identifier in
  the `citation` text (MSigDB Hallmark name, KEGG pathway ID, Reactome
  ID, GO term, WikiPathways ID). Without it, the gene set has to be
  inferred from gene-list overlap, which is unreliable.
- Use gene symbols exactly as they appear in the cited literature.
  Matching is case-insensitive and whole-token only; aliases are not
  resolved (`MAP1A` ↔ `Map1a` matches; `MAP-1A` and `MTAP1A` do not).
- DO NOT run computational models — no training, no PCA / UMAP, no
  pretrained-model embedding lookups. The embedding comes solely from
  the blueprint declared in `build_blueprint` plus literature
  knowledge. Do not write a separate embedding pipeline.
- Replace every `<<EXAMPLE_*>>` placeholder in the skeleton.

## Self-verification

After authoring the blueprint, run both checks below in order. Do not
interleave them with authoring — author the whole blueprint first
(your best guess for every DOI), then verify, then fix.

**1. Structural check.** Run:

    python3 embedder_blueprint.py --describe

This must exit cleanly. If it raises, fix the offending call before
exit. Do NOT run with `--input <data.h5ad>` — there is no input data
in this workdir; the data check runs separately.

**2. DOI batch sanity check (paper / atlas only).** Curl Crossref once
per declared paper/atlas DOI and confirm title + first author match:

    for doi in <all paper/atlas DOIs>; do
      curl -s "https://api.crossref.org/works/${doi}" \
        | python3 -c 'import sys,json;m=json.load(sys.stdin)["message"];print(m["title"][0]);print(m["author"][0]["family"])'
    done

For every mismatch, fix the DOI in `embedder_blueprint.py` (use
`web_search` to find the correct one). Do **one** consolidated fix
pass and then re-run the batch — do not enter a per-DOI
check-fix-check loop, which burns wall-clock without changing the
outcome. `kind="database"` resources are exempt (DOI optional; if
present, no Crossref check needed).

## Deliverable

An edited `embedder_blueprint.py` whose `build_blueprint()` covers the
cell populations relevant to the scenario in `scenario_config.jsonc`
and can produce an embedding for actual scRNA-seq data.
