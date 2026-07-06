# Task

Design an embedding for the cell-type/state landscape of the
biological scenario in `scenario_config.jsonc`. Instantiate it by
editing ONLY the `build_blueprint()` function in
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
  signaling pathway.

## Method

There is no literature-resource or citation mechanism. Define
gene modules by their gene lists alone — there is no `add_resource`,
and `add_gene_module` takes no `resources` argument. Draw on whatever
biological knowledge you have to choose markers.

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
- All `add_gene_module` arguments are inline string / list literals —
  no aliasing through local variables. The harness parses the file
  statically.
- Use gene symbols exactly as they appear in the literature you draw
  from. Matching is case-insensitive and whole-token only; aliases
  are not resolved (`MAP1A` ↔ `Map1a` matches; `MAP-1A` and `MTAP1A`
  do not).
- DO NOT run computational models — no training, no PCA / UMAP, no
  pretrained-model embedding lookups. The embedding comes solely from
  the blueprint declared in `build_blueprint`. Do not write a separate
  embedding pipeline.
- Replace every `<<EXAMPLE_*>>` placeholder in the skeleton.

## Self-verification

After authoring the blueprint, run the structural check below. Do not
interleave it with authoring — author the whole blueprint first, then
verify, then fix.

**Structural check.** Run:

    python3 embedder_blueprint.py --describe

This must exit cleanly. If it raises, fix the offending call before
exit. Do NOT run with `--input <data.h5ad>` — there is no input data
in this workdir; the data check runs separately.

## Deliverable

An edited `embedder_blueprint.py` whose `build_blueprint()` covers the
cell populations relevant to the scenario in `scenario_config.jsonc`
and can produce an embedding for actual scRNA-seq data.
