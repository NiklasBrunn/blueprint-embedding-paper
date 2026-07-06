# Task

You are designing an embedding for the cell-type/state landscape of
the biological scenario described in `scenario_config.jsonc`. Read
that file before doing anything else. The embedding MUST be derived
from knowledge of the relevant scientific literature for that
scenario.

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
  It may shift gene-list choice, scoring approach, and dimensionality,
  but it does NOT relax the literature-grounding rules in `## Method`.
  If the IF>10 literature does not support the requested focus, omit
  the focus rather than invent ungrounded gene lists to satisfy the note.

## Method

### Sources

Every gene list / signature you encode in `embed_freeform.py` must be
supported by at least one literature resource cited inline as a comment.

When the scenario admits it, ground the embedding in the breadth of
what the literature offers — papers (mechanism, discovery), atlases
(cell-type taxonomy), and databases (pathway / gene-set definitions).
Different gene lists can lean on different kinds; the goal is that the
embedding as a whole reflects the available evidence base.

Each resource is one of three kinds:

- **paper** — peer-reviewed paper in a journal with **impact factor > 10**.
- **atlas** — established reference atlas (Allen Brain Atlas, Tabula
  Muris, Human Cell Atlas, …), cited via its accompanying peer-reviewed
  paper (also IF > 10).
- **database** — curated pathway / gene-set database (KEGG, MSigDB,
  Reactome, GO, WikiPathways, …), cited via its companion peer-reviewed
  paper or a versioned deposit DOI.

Every cited **paper** or **atlas** resource MUST have a valid DOI in
the format `10.xxxx/yyy`. Bare URLs and `doi:` prefixes are not
acceptable. For **database** citations, DOI is optional — the
gene-set identifier (MSigDB Hallmark name, KEGG pathway ID, Reactome
ID, GO term, WikiPathways ID) is the canonical handle.

For every gene list / signature in your script, attach a comment
immediately above it citing the supporting resource(s) by name and
DOI (or gene-set ID for databases).

**DOI batch sanity check** (paper/atlas only). After you have written
the full script — your best guess for every DOI — run, once, before
exiting:

    for doi in <all paper/atlas DOIs>; do
      curl -s "https://api.crossref.org/works/${doi}" \
        | python3 -c 'import sys,json;m=json.load(sys.stdin)["message"];print(m["title"][0]);print(m["author"][0]["family"])'
    done

For every mismatch, fix the DOI in the script (use `web_search` to
find the correct one). Do **one** consolidated fix pass and then
re-run the batch — do not enter a per-DOI check-fix-check loop.

**Database fallback (recommended).** When a gene list maps to a curated
gene-set entry (MSigDB Hallmark, KEGG, Reactome, GO, WikiPathways),
co-cite a database resource alongside the paper(s).

### Approach

The embedding MUST be built **solely from literature knowledge** —
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
