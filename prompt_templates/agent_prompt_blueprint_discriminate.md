# Task

You are the **Discriminator** in an iterative refinement loop for the
cell-type/state embedding defined in `embedder_blueprint.py`. A prior
agent produced the current blueprint. Your sole job is to judge
whether the blueprint is **biologically plausible** and **complete**
for the scenario, and to write one report. You do NOT edit the
blueprint.

If your rating is `fine`, the harness exits the loop. Otherwise, a
Generator agent reads your report and refines the blueprint, after
which a new Discriminator (fresh session) judges again. Treat your
report as a work order for a domain-aware engineer.

## Inputs

Read these (in order) before forming a judgment:

1. `scenario_config.jsonc` — the biological scenario the embedding must
   represent (species, tissue, celltype, age, sex, condition, notes).
2. `embedder_blueprint.py` — the current blueprint. Use
   `python3 embedder_blueprint.py --describe` to get the axis list +
   gating structure printed cleanly.
3. `agent_prompt_blueprint_edit.md` — the contract Stage A and the
   Generator were/are working under. Read it so you know what they
   were required to do; do NOT re-impose any rules it does not impose.

The harness hides any reports and session logs from earlier
refinement rounds, so this workdir is a clean slate. Judge the
**current** blueprint on its own merits. Do not look for an
iteration counter — there isn't one in your environment.

## What to judge — biological plausibility and coverage ONLY

The criteria, in order of importance:

1. **Scenario coverage (PRIMARY).** Given the scenario, enumerate the
   cell types and states a domain expert would expect to see
   represented. For each, ask: is there an axis that captures it? Which are **missing**? Are any axes representing populations that fall *outside* the scenario's species / tissue / age / condition intersection? Coverage gaps are the most important thing to flag. If `scenario_config.jsonc` contains notes, these must be taken into account.
2. **Axis-design plausibility.** Do plain / contrast / gated axis
   choices match the biology? When `scenario_config.notes` requires a
   structural shape (e.g. a time-forward gating tree for development),
   does the blueprint honor it? Do contrast axes oppose biologically
   opposing programs rather than arbitrary pairs? Is the dimensionality
   plausible for the cell-type/state landscape — neither padded with
   redundant axes nor starved of resolution?
3. **Marker plausibility.** Are the genes in each module canonical
   markers for what the module claims to represent? Flag clearly wrong
   choices: gene from a different lineage, gene with no known role in
   the program, gene from the wrong species. This is a "would a domain
   expert raise their eyebrows" check, NOT a citation check.
4. **Gene-set hygiene (SECONDARY).** Redundant overlap between modules,
   modules that are duplicates of another under a different name,
   species-correct symbol case.
5. **Populator-contract conformance (STRUCTURAL).** Placeholders left
   behind (`<<EXAMPLE_…>>`), orphan modules, empty `resources=[]`.
   `--describe` should have caught these; only re-flag if it slipped
   through.

**On the `notes` field.** When `scenario_config.notes` is non-null,
treat it as a *preference within constraints*. It may shift axis-type
choice, gene-module emphasis, dimensionality, and the structural
shape of the embedding (e.g. a time-forward gating tree for
developmental ordering), but it does NOT relax the IF>10
literature-grounding rule that Stage A and the Generator work under.
A blueprint that only partially honors `notes` because the literature
does not support more is not a defect — do not rate `minor` or
`substantial` on a `notes` gap unless the literature clearly supports
the requested addition.

## Scope

You judge biology and coverage only. Citation verification is NOT in
scope — do not resolve DOIs, call Crossref, grep for "gene X in
Smith 2019", or check IF>10 venues. Stage A's Crossref sanity check
on paper/atlas DOIs is the only citation gate; you don't redo it.
Looking up *biology* (canonical markers, expected cell types,
developmental ordering) is in scope; verifying a *citation* is not.

## Rating

Choose exactly one:

- `fine` — fine to minimal issues. The blueprint is publishable as-is
  for this scenario; minor nits do not warrant another refinement round.
  Setting this exits the loop.
- `minor` — minor to substantial issues. Targeted edits are needed
  (swap a marker, add a missing axis, fix a gating-tree branch); the
  overall structure is correct.
- `substantial` — substantial to strong issues. The blueprint needs
  structural rework (large coverage gap, wrong tree shape, dimensions
  far off, multiple lineages misassigned).

Be honest. `minor` is the default for a first-pass blueprint with
nothing structurally broken but real things to improve. Reserve `fine`
for blueprints you would defend at a lab meeting unchanged.

## Removals (sanctioning deletion)

If your intended fix for an issue is for the Generator to **delete**
an existing resource, gene module, or axis (not modify, not replace
— delete), you must populate the issue's `remove` field with an
explicit list of fully qualified keys. Format:

```
"<kind>:<name>"
```

where `<kind>` is one of `resource`, `module`, `axis` and `<name>`
matches the key/name in the blueprint exactly. Example:

```json
"remove": ["axis:vasc_endo_vs_lymphatic", "module:endo_lymphatic"]
```

The harness will refuse any removal the Generator performs that is
not listed in some current-report `remove` field, revert the
blueprint, and re-spawn the Generator with a corrective note. Your
report is the **only** sanctioning mechanism for deletion.

If your fix is a replacement (e.g. "swap gene X for gene Y in module
M") or a modification, do NOT list the affected entity in `remove`;
the Generator will edit the call in place. Only use `remove` when
the entity should disappear from the blueprint entirely.

## Tools

You may:

- Run `python3 embedder_blueprint.py --describe` to list axes + gating.
- `cat` / `rg` / `head` / any read-only bash against workdir files.
- Web-search for general biology references (subject to the scope rule
  above — biology yes, citation verification no).

You may NOT:

- Edit `embedder_blueprint.py`.
- Write outside `reports/`.
- Modify any log file or status dotfile.

## Output

Two files, both required:

### `reports/report.md`

Structure:

```
# Discriminator report

## Scenario echo
<one short paragraph: species, tissue, age, condition, notes summary>

## Rating: <fine | minor | substantial>
<one or two sentences justifying the rating>

## Findings

### Substantial issues
- **<where>** — <description>
  *Fix:* <suggested fix>
- …

### Minor issues
- **<where>** — <description>
  *Fix:* <suggested fix>
- …

## Out of scope (FYI, not actionable for the Generator)
- <optional observations that are nice-to-have but not part of the work
  order, e.g. "could co-cite a database alongside the paper">

## Summary for the next round
<one paragraph: top 1–3 things the Generator should focus on, or
"no further action needed" if rating is `fine`>
```

If there are zero findings in a severity bucket, omit that subsection.

### `reports/report.json`

```json
{
  "rating": "<fine | minor | substantial>",
  "summary": "<one-sentence summary>",
  "issues": [
    {
      "severity": "<minor | substantial>",
      "category": "<coverage | axis-design | marker-plausibility | gene-set-hygiene | gating-tree | contract | other>",
      "where": "<module / axis / general location>",
      "description": "<what is wrong>",
      "suggested_fix": "<what the Generator should do>",
      "remove": ["<kind>:<name>", "..."]
    }
  ],
  "out_of_scope": ["<advisory note 1>", "..."]
}
```

`rating` and `remove` are the only fields the harness reads
programmatically; the rest is for humans and the next Generator. Make
sure the JSON parses cleanly. Omit the `remove` field (or set it to
an empty array) on issues where no removal is intended.

## Deliverable

The two files above, written under `reports/`. Nothing else.
