# Discriminator report

## Scenario echo
The scenario is human pancreas, with all cell types in scope because `celltype` is null. Age, sex, and condition are unspecified/null, so the blueprint should cover the baseline human pancreatic cell-type landscape without assuming a disease, treatment, sex-specific, or developmental focus; there are no notes imposing a special topology.

## Rating: minor
The blueprint is broadly plausible and covers the major pancreatic compartments: acinar, ductal/centroacinar, endocrine islet cells, vascular, stromal, neural, adipocyte, and immune populations. It still needs targeted refinement for a few normal pancreas populations that are missing or merged and for some marker/axis choices that are weaker than the overall design.

## Findings

### Minor issues
- **mural_perivascular_identity / vascular_mural** — Pericytes and vascular smooth muscle are collapsed into one axis/module. `RGS5`, `PDGFRB`, `CSPG4`, and `MCAM` capture pericytes, while `ACTA2`, `TAGLN`, `MYH11`, and `MYLK` capture smooth muscle/myofibroblast-like contractile cells. In human pancreas these are related but distinct vascular support populations, and merging them reduces coverage resolution.
  *Fix:* Split into separate pericyte and vascular smooth muscle modules/axes, or add a gated/contrast refinement axis within the existing mural compartment.
- **general coverage** — Mesothelial/capsular surface cells are not represented. Human pancreas tissue preparations can contain mesothelial or capsule-associated cells, and an all-cell pancreas embedding should have a way to separate them from ductal epithelium and fibroblasts.
  *Fix:* Add a mesothelial module and plain axis with markers such as `MSLN`, `WT1`, `UPK3B`, `KRT19`, `KRT8`, `KRT18`, `PDPN`, and/or `LUM` only as appropriate for the chosen module composition.
- **ductal_centroacinar_terminal_state / ductal_epithelial** — Ductal biology is represented mostly as one broad ductal identity plus a centroacinar/terminal duct refinement. This misses useful normal ductal heterogeneity, especially ion/fluid transport ductal programs versus mucin/SPP1-high ductal states, which can otherwise blur with centroacinar or stress-like epithelial states.
  *Fix:* Add one ductal refinement axis for bicarbonate/ion-transport ducts (`CFTR`, `SLC4A4`, `SLC26A6`, `AQP1`, `KRT7`, `KRT19`) versus mucin/injury-associated ductal markers (`SPP1`, `MUC1`, `TACSTD2`, `KRT17`) or otherwise sharpen the existing ductal modules.
- **endocrine_enterochromaffin_serotonin_cell / enterochromaffin_serotonin_cell** — A separate pancreatic enterochromaffin/serotonin endocrine cell axis is over-specific for a baseline human pancreas embedding. Serotonergic markers can appear in rare pancreatic endocrine contexts, but this is not as canonical as beta, alpha, delta, PP, and ghrelin/epsilon cells and risks importing a gut enterochromaffin identity.
  *Fix:* Recast this as a rare endocrine serotonergic state rather than a full canonical cell type, keep it lower priority, or fold it into a broader rare endocrine-state axis with careful negative modules.
- **epsilon_ghrelin_cell** — The epsilon/ghrelin module contains several weak or shared markers (`MLN`, `IAPP`, `SCG2`, `PAX6`) alongside `GHRL`. `IAPP` is strongly beta-associated, and `SCG2`/`PAX6` are broad endocrine markers, so the axis may not cleanly isolate ghrelin cells.
  *Fix:* Strengthen the module around ghrelin-cell markers and remove or downweight beta/common endocrine genes; consider `GHRL`, `PCSK1N`, and other pancreas-supported ghrelin-cell markers while avoiding `IAPP`.
- **immune_b_plasma_cell / b_plasma_cell** — The module name says `b_plasma_cell`, but the markers are conventional B-cell markers (`MS4A1`, `CD79A`, `CD79B`, `BANK1`, `CD22`) and plasma cells are already represented separately. This is mostly a naming/semantic hygiene issue, but it can confuse downstream interpretation.
  *Fix:* Rename or reinterpret this as a B-cell module/axis, keeping the separate `plasma_cell` module for antibody-secreting cells.

## Out of scope (FYI, not actionable for the Generator)
- The global cycling, interferon, immediate-early stress, and ER-load axes are reasonable generic state axes for dissociation/secretory tissue context; they do not require a disease condition.
- The blueprint is somewhat high-dimensional, but 32 axes is defensible for an all-cell human pancreas embedding with immune and stromal subtyping.

## Summary for the next round
Focus on targeted coverage refinements: split pericyte from vascular smooth muscle, add mesothelial/capsular coverage, and sharpen ductal/endocrine rare-state handling. The existing high-level structure should be preserved; no wholesale redesign is needed.
