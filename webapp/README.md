# Blueprint Pipeline Web App

A local-first [Plotly Dash](https://dash.plotly.com/) app that wraps the agentic
blueprint pipeline: define a scenario, generate a biology-grounded embedding
blueprint with your own API key, download/hand-edit it, drop in an `.h5ad`, map
it to embedding coordinates, and explore the result interactively as a UMAP.

Everything this app adds lives under `webapp/`. It **edits no other file** in
the repository and writes runtime artifacts only under `webapp/runtime/`
(git-ignored).

## Standalone (no repo needed)

The `webapp/` folder is self-contained: it ships its own virtual environment
(`webapp/.venv`) **and** vendored copies of the canonical pipeline sources under
`webapp/vendor/` (the blueprint templates, prompt templates, `extract_run.py`,
`regression_guard.py`, `extract_cost.py`). At runtime the app prefers
`webapp/vendor/` and only falls back to the parent repo if the vendor dir is
absent, so you can copy `webapp/` anywhere and run it without the rest of the
repository.

Caveats for a copied-out bundle: real blueprint generation still needs the
external `codex` CLI on your `PATH` plus an API key (unchanged), and the
token-free mock mode (`WEBAPP_MOCK_CODEX=1`) is repo-only since it seeds from
archived blueprints in `outs/`.

To refresh the vendored sources after the repo originals change (run in-tree):

```bash
bash webapp/sync_vendor.sh            # re-copies the canonical files into webapp/vendor/
```

`webapp/tests/test_webapp_vendor.py` fails if `vendor/` drifts from the repo
originals (and skips that check when the repo isn't present).

## System requirements

- **Operating system:** macOS or Linux. The app uses a POSIX shell for the
  pipeline driver and, for the native **Browse…** file picker, macOS `osascript`
  with a `tkinter` fallback. On Windows, run it under **WSL2**.
- **Python:** ≥ 3.10 (3.11 / 3.12 recommended). `setup.sh` auto-selects one.
- **Disk:** ~0.5 GB for the virtual environment (it includes a vendored
  Plotly.js). `webapp/runtime/` then grows with each Map/Generate run; it is
  git-ignored and safe to delete between sessions.
- **Memory:** enough to load your `.h5ad` and compute a UMAP — budget several GB
  of RAM for ~10⁵ cells.
- **For the Generate tab only** (Map and Explore need none of this):
  - **Node.js ≥ 18** — to install the Codex CLI.
  - The **OpenAI Codex CLI** on your `PATH`.
  - An **OpenAI account + API key** with access to the model you pick (e.g. a
    GPT-5.x model). API usage is billed to that account.

## Installation

**1 — The web app (Python).** One command creates an isolated virtual
environment and installs everything (no `conda activate` needed):

```bash
bash webapp/setup.sh
```

**2 — The Codex CLI (only if you will use the Generate tab).** Map/Explore work
without it.

```bash
# Requires Node.js >= 18. Install the CLI globally:
npm install -g @openai/codex
# (Alternative install methods — Homebrew, prebuilt binaries — are listed at
#  https://github.com/openai/codex)

codex --version                       # verify it is on your PATH
```

Authenticate the CLI in one of two ways:

- **API key (used by this app):** the app passes your key to `codex` via the
  `OPENAI_API_KEY` environment variable — paste it in the Generate tab, or
  `export OPENAI_API_KEY=...` in the shell before launching (see *API key
  handling* below). **or**
- **Sign-in:** run `codex login` once (ChatGPT account sign-in).

No manual "trust this directory" step is needed — the app runs `codex` with
`--skip-git-repo-check` and an explicit sandbox; see *Generation prerequisites* below.

**3 — Launch:**

```bash
webapp/.venv/bin/python -m webapp     # serves http://127.0.0.1:8050
```

## Run it

```bash
bash webapp/setup.sh                  # one-time: creates webapp/.venv, installs deps
webapp/.venv/bin/python -m webapp     # serves http://127.0.0.1:8050
```

The app ships its own virtual environment — you do **not** need to
`conda activate agentic-embeddings`. If `OPENAI_API_KEY` is set in your shell it
is picked up automatically.

Open the printed URL in a browser. Two tabs:

### Generate
Fill in the scenario (species, tissue, optional age/sex/condition), toggle
**Temporal analysis** if you want a time-forward axis tree, add free-form
**Notes / additional rules** (the temporal toggle auto-fills the temporal-tree
rules here, which you can edit), pick the model, reasoning effort, and pipeline
steps (Stage-B refinement, resource-unconstrained "noref"), paste your OpenAI API
key, and click **Generate**. A live phase tracker (Stage A → Stage B → Done) and
log tail show progress, with an ETA estimated from past runs. When it finishes
you can download `embedder_blueprint.py`, `summary.json`, and the report.

**Edit blueprint (in-browser).** After generation the blueprint loads into a
syntax-highlighted editor at the bottom of the Generate tab. Edit it (e.g. fix a
gene name, drop a module), click **Save & validate** — invalid edits are flagged
(placeholders left, orphaned modules) but still saved so work persists — and
**Use in-session blueprint** on the Explore tab will map your edited version.
**Revert to generated** restores the original. No download / re-upload needed.
(Editing executes Python locally; disable it for any hosted/untrusted instance.)

### Explore
Provide an `.h5ad` (or `.h5ad.gz`) one of two ways: **drag-and-drop** (loads the
whole file through the browser into memory — fine up to ~100 MB), or pick a
**local file** — click **Browse…** to open a native file chooser (or paste the
path) and the server reads it straight from disk (no size limit — recommended
for large single-cell files). Browse opens the dialog on the machine running the
app, so it is for local use (the browser cannot reveal a file's path itself). Choose the blueprint to map with: the
one generated this session (including your in-browser edits), **or** an uploaded
`embedder_blueprint.py`. Click **Map** to run the embedder and compute a UMAP —
an animated progress bar shows while mapping runs. Then color the UMAP by:

- **metadata** — categorical obs columns get a discrete palette; numeric columns
  a sequential scale;
- **gene** — log1p expression of any gene in a blueprint module (sequential,
  positive-only);
- **axis** — activation score of any axis (diverging/symmetric for
  contrast/gated axes, sequential for plain axes).

If the blueprint has a gate tree, a **Temporal trajectory** panel appears with a
branch selector and a time **ruler**: scrubbing interpolates a convex
combination of adjacent axes along the chosen root→leaf branch
(`(1-α)·vᵢ + α·vᵢ₊₁`), the same morphing view as the repo's interactive HTML.

Download the augmented `.h5ad` (embedding coords in `obsm['X_blueprint']`, UMAP
in `obsm['X_umap']`, per-axis scores in `obs`) and the current plot (PNG / SVG /
HTML).

## Generation prerequisites (codex CLI)

The **Generate** tab shells out to the external `codex` CLI, which must be on
your `PATH`, plus your own OpenAI API key. (Map/Explore need neither.)

The pipeline runs `codex` **non-interactively** inside per-run workdirs under
`webapp/runtime/`, so it works out of the box on any machine with **no manual
directory-trust step**. The app invokes codex with these flags (see
`CODEX_ARGS` in `webapp/pipeline_driver.sh`):

- **`--skip-git-repo-check`** — the runtime workdirs are not git repositories, and
  codex otherwise refuses to run in a non-repo / untrusted directory (failing with
  `Not in a trusted directory` in the generation log). This flag skips that guard.
- **`-s workspace-write`** — pins codex to the least-privilege *workspace-write*
  sandbox: the agent may write **only within its own workdir** (it edits
  `build_blueprint()`); reads are not restricted. Network access is on
  (`sandbox_workspace_write.network_access=true`) for web search and the DOI check.
- **`--disable plugins --disable apps --disable memories -c skills.include_instructions=false`**
  — hides your personal Codex plugins, apps, memories and skills from the agent.

These require a reasonably **recent codex** (the flags/config keys must be
supported) — check with `codex --version`. If an older codex still refuses to run,
trust the folder once as a fallback by adding to `~/.codex/config.toml` (absolute
path from `cd /path/to/webapp && pwd`):

  ```toml
  [projects."/abs/path/to/webapp"]
  trust_level = "trusted"
  ```

## API key handling

Your key is held **only in server-side memory** for the session, passed to the
codex subprocess via the `OPENAI_API_KEY` environment variable, never written to
disk, never logged, and cleared from the browser field on launch and dropped
from memory when the job stops. If `OPENAI_API_KEY` is already set in your
environment when you launch the app, you don't need to paste anything.

## Deployment notes (hosted instance)

The same app can be served behind TLS as a shared demo. For untrusted users:

- **Disable `.py` blueprint upload** (executing an uploaded blueprint is code
  execution; safe on your own machine, not for strangers).
- **Cap concurrent jobs** and add request limits.
- Never persist keys; keep the in-memory-only handling above.

## Token-free development

Set `WEBAPP_MOCK_CODEX=1` to skip the real `codex` CLI: the vendored driver
populates the blueprint from a real archived one instead, so you can exercise
the full Generate → Map → Explore flow without spending API tokens. This is what
the end-to-end test (`webapp/tests/test_webapp_e2e_mock.py`) uses.

## Tests

```bash
webapp/.venv/bin/python -m pytest webapp/tests -v
```

## License & third-party notices

This app is released under the **MIT License** — see [`LICENSE`](LICENSE).
Third-party components and their licenses, plus the external tools it relies on,
are documented in [`NOTICE`](NOTICE).

Key points for redistribution:

- **You are shipping this wrapper, not OpenAI's software.** The Codex CLI is
  *invoked*, never bundled; each user installs it themselves and runs it under
  **their own OpenAI account and API key**. This repository embeds no key and
  grants no access to OpenAI services. Generated blueprints/outputs are produced
  by, and subject to, OpenAI's [Terms of Use](https://openai.com/policies).
- **Dependencies are permissive.** The bundled `plotly.min.js` is MIT (its
  header is retained); the Python packages are MIT / BSD-3-Clause / Apache-2.0
  (see `NOTICE`). None are copyleft.
- **No data is bundled.** Any datasets you load remain under their own original
  licenses, independent of this software.
- **Before publishing publicly:** (1) confirm the **copyright holder** line in
  `LICENSE` matches your institution's policy (e.g. you, and/or
  Universitätsklinikum/Universität Freiburg) and clear the release with your
  group lead / technology-transfer office; (2) if you offer a *hosted* instance
  for untrusted users, follow the hardening in *Deployment notes* above.

*This section is informational, not legal advice.*
