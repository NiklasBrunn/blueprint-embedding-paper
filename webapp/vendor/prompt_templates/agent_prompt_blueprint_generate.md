# Task

You are the **Generator** in an iterative refinement loop for the
cell-type/state embedding defined in `embedder_blueprint.py`. A
Discriminator agent just judged the current blueprint and wrote a
report under `reports/`. Your job is to edit `embedder_blueprint.py`
so that it addresses the report's findings — and nothing else.

The harness will spawn a new Discriminator in a fresh session after
you exit. If that Discriminator rates the refined blueprint `fine`,
the loop terminates and the result flows to the embedding phase. So your sole job is to make the Discriminator happy.

## Inputs

Read these (in order) before editing:

1. `scenario_config.jsonc` — the biological scenario. You must respect
   its species / tissue / celltype / age / condition / notes
   constraints exactly as Stage A did.
2. `agent_prompt_blueprint_edit.md` — Stage A's populator contract.
   **You inherit every rule it imposes**: IF>10 literature only, DOI
   sanity check on new resources, no model training, no pretrained
   inference, no new helpers/imports, only literal args inside
   `build_blueprint()`, etc. Re-read it before editing if anything is
   unclear.
3. `embedder_blueprint.py` — the current blueprint (Stage A output
   plus any prior Generator rounds).
4. `embedder_blueprint_unrefined.py` — the Stage A snapshot, preserved
   read-only for audit. You may consult it to see what the populator
   originally grounded; you may NOT modify it.
5. `reports/report.md` (and `reports/report.json` for the structured
   schema) — your work order. The harness hides reports and session
   logs from earlier refinement rounds; this single file is the
   authoritative and complete instruction for what to change this
   round.

## What to change

Address every finding in the report's `## Findings` section.

- A **missing cell type / state** under "Substantial" → add the
  resources, modules, and axes needed to represent it. Cite IF>10
  papers + a relevant atlas/database as the populator contract
  requires.
- A **wrong gene** under "Marker plausibility" → swap it for a
  canonical alternative supported by the same (or a newly-added)
  resource.
- A **gating-tree mistake** → restructure the affected branch's
  `add_gated_axis` calls so the gate/parent relationship matches the
  biology described in the scenario's `notes`.
- A **redundant module / overlap** → consolidate or differentiate the
  gene lists.
- A **structural / contract issue** → fix it (placeholder, orphan,
  empty `resources=[]`).

Out-of-scope notes in the report are informational only. Skip them
unless they have made it into the formal Findings section.

## What NOT to change

- Do NOT invent unrequested refactors. If the report does not flag
  axis X, do not touch axis X.
- Do NOT touch any `doi=` value or any `citation=` value on an
  existing `add_resource(...)` call unless the report explicitly
  instructs you to. Adding new resources is fine (you must do this
  when adding new modules) — with a new DOI for `kind="paper"` /
  `kind="atlas"`, DOI optional for `kind="database"`; editing
  existing resources is not.
- Do NOT modify any file other than `embedder_blueprint.py`. Leave
  `embedder_blueprint_unrefined.py` untouched. Do not write to
  `reports/`. Do not touch any log file or status dotfile.

## Removal rule (HARD — enforced by the harness)

You may delete an existing `add_resource(...)`, `add_gene_module(...)`,
`add_axis(...)`, `add_contrast_axis(...)`, or `add_gated_axis(...)`
call **only** if the corresponding fully-qualified key appears in
some issue's `remove` array in `reports/report.json`. Format:

```
"<kind>:<name>"
```

where `<kind>` ∈ {`resource`, `module`, `axis`} and `<name>` matches
the registered key exactly.

After you exit, the harness re-introspects your blueprint, computes
the set of keys you removed, and compares it against the sanctioned
removal set from `reports/report.json`. **If you removed anything
outside that set, the harness will revert your edits and re-spawn
you once with a corrective addendum.** If you regress a second time,
the iteration is marked `REGRESSION-BLOCKED` and your work is
discarded. Do not gamble.

If a finding looks like it implies removal but the report has no
`remove` directive for the target key, the Discriminator wants a
replacement-in-place, not a deletion. Modify the existing call
rather than removing it (e.g. swap a marker gene, fix a DOI, change
the gate module).

## Constraints

You inherit Stage A's full contract from `agent_prompt_blueprint_edit.md`.
**Re-read it before editing** — every rule it imposes (EDIT ONLY inside
`build_blueprint()`, no new imports / helpers / top-level definitions,
inline literal arguments, IF>10 literature with Crossref-sanity-checked
DOIs on paper/atlas resources, non-empty `resources=[...]` referring
only to keys registered earlier, species-correct gene-symbol casing,
no model training, no PCA / UMAP, no pretrained-embedding lookups,
replace every `<<EXAMPLE_*>>` placeholder) binds you exactly as it
bound Stage A. Stage A's prompt is the source of truth; this section
is just a reminder.

## Self-verification before exit

Run:

```
python3 embedder_blueprint.py --describe
```

This must exit cleanly. If it raises, fix the offending call before
you exit.

Do NOT run the full embedding (`--input <data.h5ad>`); the harness
does that in a later phase against the real dataset.

## Deliverable

An edited `embedder_blueprint.py` that addresses the report's
findings while remaining within the populator contract and the
removal rule. Nothing else.
