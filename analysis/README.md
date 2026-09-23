# Analysis scripts

Bash orchestrators for the end-to-end workflow. All are run from the repo root with
`conda activate agentic-embeddings` (see `../environment.yml`); each takes `--help`.

| Script | Purpose |
| --- | --- |
| `run_data_pipeline.sh` | Download → QC → preprocess the datasets into `data/preprocessed/<dataset>_pp.h5ad`. |
| `agent_edit_blueprint.sh` | Author blueprints with the **Codex** CLI (Stage A, optional Stage B refinement), then embed onto the data. |
| `agent_edit_blueprint_claude.sh` | Same flow driven by Anthropic's **Claude** CLI (cross-engine check). |
| `agent_refine_blueprint.sh` / `_claude.sh` | Run only the Stage B refinement loop over an existing blueprint. |
| `agent_edit_freeform.sh` | Free-form prompting baseline (arbitrary script, same output contract, no blueprint API). |
| `reproduce_figures.sh` | Re-plot the quantitative/structural figures from the shipped `outs/` CSVs (no data, no baselines). |
| `reproduce_all.sh` | Full end-to-end rebuild from the shipped blueprints (data → embeddings → baselines → metrics → all figures). |

The authoring scripts write agent run dirs under `work/` (git-ignored) and, per accepted
run, the Embed-step artifacts under `outs/embeddings/<dataset>/` plus an archived blueprint
under `outs/blueprints/<dataset>/<cohort>/`. Both Codex and Claude orchestrators produce the
same workdir layout and aggregate as separate cohorts in the metric scripts.

## Authoring blueprints

```bash
# One authoring run (Codex, GPT-5.5 high), Stage A only, then embed onto the data:
bash analysis/agent_edit_blueprint.sh --dataset mouse_cortex

# Full flag surface:
bash analysis/agent_edit_blueprint.sh --help
```

### Stage B refinement (`--with-refine`)

With `--with-refine`, a *discriminator* rates the blueprint (`fine` / `minor` /
`substantial`) and a *generator* revises it, under a regression guard that reverts
unauthorized deletions, for up to `--max-iterations` (default 3). The Embed step then
emits **two** cohorts per run:

- the pre-refinement Stage A snapshot → `blueprint_<model>_unrefined_<reasoning>_<dataset>_seed<N>.parquet`
- the post-refinement blueprint → `blueprint_<model>_refine_<reasoning>_<dataset>_seed<N>.parquet`

Both appear in `outs/metrics/<dataset>/metrics.csv` and are joinable by `(method, dataset,
seed)`, so Stage A vs Stage A+B comparisons fall out of the aggregators with no extra
configuration.

### Resource-unconstrained variant (`--noref`)

`--noref` runs a comparison cohort without literature grounding: there is no citation
mechanism, and gene modules are plain gene lists. It swaps in
`scripts/agentic_embedding/embedder_blueprint_noref.py` and the `*_noref.md` prompt
templates, and tags the cohort `_noref_{unrefined|refine}`. It composes with
`--with-refine`. The runs surface as a distinct cohort in `metrics.csv` / `costs.csv`
alongside the `_unrefined` / `_refine` rows.

## Agent isolation

Codex runs use `-s workspace-write` (edits confined to the workdir, network on for web search
and the Crossref check) and hide the user's Codex plugins, apps, memories and skills
(`--disable plugins --disable apps --disable memories -c skills.include_instructions=false`).

### Claude script: isolation caveat

Like `codex`, the `claude` CLI auto-loads project- and user-specific context into every
non-interactive run:

- **Project auto-memory** at `~/.claude/projects/<slug>/memory/MEMORY.md` (and linked notes),
  a local per-machine artifact. If that memory holds anything
  relevant to blueprint authoring, for example a list of high-impact journals or reference
  atlases, Claude effectively starts each run with a partial cheat sheet for the
  literature-grounding sub-task that Codex never sees.
- **User-level skills, plugins and MCP servers**, which can shift the agent's process priors
  or add literature access Codex runs do not have.

To keep Claude cohorts apples-to-apples with Codex, `agent_edit_blueprint_claude.sh`
defaults to **isolated mode** (pass `--no-isolate` to opt out):

- `CLAUDE_CODE_DISABLE_AUTO_MEMORY=1` — disables project auto-memory.
- `--safe-mode` — disables CLAUDE.md, skills, plugins, hooks and custom agents.
- `--strict-mcp-config` — loads no MCP servers (e.g. literature connectors).

These knobs keep subscription OAuth working (no `--bare`), so default-isolated runs do not
require an `ANTHROPIC_API_KEY`. `--no-isolate` loads project memory + skills and tags the
cohort `_userctx`; use it only for diagnostic probes, never for the Codex comparison. See
`bash analysis/agent_edit_blueprint_claude.sh --help` for the full flag surface.
