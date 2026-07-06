# Task

You are designing an embedding for the cell-type/state landscape of
the biological scenario described in `scenario_config.jsonc`. Read
that file before doing anything else.

You MUST produce the embedding by writing a single Python file
`embed_freeform.py` in the current working directory. Do NOT create
or modify any other file.

## Scenario

`scenario_config.jsonc` defines the biological scenario the embedding
must represent. Fields:

- `species` — array of strings, each entry one species the embedding
  must cover. At least one entry. Examples: `["mouse"]`, `["human"]`,
  `["mouse", "human"]`.
- `tissue` — array of strings, each entry one tissue or region of
  interest. `null` if unspecified.
- `celltype` — array of strings, each entry one cell type of interest
  that should be included in the embedding. Optional. If `null`, the
  embedding should cover all cell types in the given tissue.
- `age` — single string. Free-form: broad categories (`"adult"`,
  `"embryonic"`), ranges (`"P15-P20"`, `"E8-E16"`), or exact values.
  `null` if unspecified.
- `sex` — single string. Examples: `"male"`, `"female"`, `"mixed"`.
  `null` if unspecified.
- `condition` — array of strings, each entry one condition (gene
  knockout, treatment, disease). `null` means baseline / wild-type /
  healthy.
- `notes` — single string. Free-form description of design focus the
  embedding should respect (e.g. resolve a specific subtype group,
  preserve a temporal/maturation gradient, emphasize a signaling
  pathway). `null` if unspecified.

Rules:

- If a field is `null`, treat it as unspecified. Do not invent
  defaults or constraints to fill the gap.
- The embedding must reflect the cell types and states biologically
  expected within the scenario's species / tissue / age / sex /
  condition intersection. Do not include cell types or states that
  fall outside that scope.
- If `notes` is non-null, treat it as a *preference within constraints*.
  It may shift gene-list choice, scoring approach, and dimensionality.

## Method

There is no literature-resource or citation mechanism. Define the gene
lists / signatures you encode in `embed_freeform.py` by their genes
alone — drawing on whatever biological knowledge you have to choose
markers. You MAY add inline comments naming a source where it helps,
but citations are not required and there is no DOI/journal/atlas/
database requirement.

### Approach

The embedding MUST be built **solely from biological knowledge** —
literal gene lists / signatures encoded in your script and a scoring
function that combines them.

### Dimensionality

Critically assess the number of embedding dimensions needed to cover
the cell-type/state landscape of this scenario. Do not pick a fixed
number out of habit; let the biology drive it.

## Script contract

`embed_freeform.py` must be runnable as:

    python embed_freeform.py --input <path> --output <parquet_path> --seed <int>

It must produce a parquet with these columns and nothing else:

- `cell_id` (string) — one row per cell, matches the AnnData `obs_names`.
- `dim_0`, `dim_1`, ..., `dim_{K-1}` (float32) — the embedding coordinates.

## Input format (count mode)

The `--input` path is a `.h5ad` AnnData file with cells × genes:

- `layers["counts"]` — raw integer counts (cells × genes).
- `X` — log1p-normalized values (`scanpy.pp.normalize_total` to
  `target_sum=1e4`, then `scanpy.pp.log1p`).
- `var_names` — gene symbols (HGNC for human, MGI for mouse).
- `obs_names` — cell identifiers; copy these as the `cell_id` column.
- `var["highly_variable"]` — boolean mask flagging 2000 HVGs (optional;
  may be absent on some inputs).

Your scoring function operates on the gene expression matrix
(per-cell, per-module).

## Constraints

- **Single Python file**: only `embed_freeform.py`. No scratch files,
  no logs, no JSON dumps in the working directory.
- **Seed**: `--seed` is a run tag passed by the orchestrator for
  output-naming and provenance. Your script may ignore it. There is no
  requirement that runs at the same seed produce identical embeddings.
- **Self-test**: when you are done, run `python embed_freeform.py --help`
  and confirm it exits cleanly. Do **not** run the full embedding.
- **Data usage**: DO NOT train any model on the input data. DO NOT run pretrained-model inference on the input data.

## Deliverable

A single file `embed_freeform.py` in the current working directory.
That is all.
