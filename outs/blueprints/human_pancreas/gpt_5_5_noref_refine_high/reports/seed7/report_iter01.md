# Discriminator report

## Scenario echo
The scenario is human pancreas with no specified age, sex, cell subtype, or disease/treatment condition, so the blueprint should cover the healthy pancreatic cell-type landscape broadly: exocrine acinar and ductal/centroacinar epithelium, endocrine islet subtypes, stromal/mural/endothelial/lymphatic compartments, nerves/glia, resident immune cells, and plausible baseline or dissociation-related states. There are no notes imposing a special topology or emphasis.

## Rating: minor
The blueprint covers the major pancreatic compartments and most canonical endocrine/exocrine/stromal/immune populations. It is not quite publishable unchanged because several axes are biologically awkward or over/under-specific for an unrestricted healthy pancreas scenario, but these are targeted edits rather than structural failure.

## Findings

### Minor issues
- **axis:lymphatic_endothelial_identity / module:vascular_endothelial** — The lymphatic axis is gated by `vascular_endothelial`, whose genes include blood-endothelial/capillary-skewed markers such as `VWF`, `KDR`, `EMCN`, and `PLVAP`. True lymphatic endothelial cells share a general endothelial core but may score weakly on this vascular module, so the gate can suppress the lymphatic signal.
  *Fix:* Split an `endothelial_core` module (`PECAM1`, `CDH5`, `ESAM`, `CLDN5`, `RAMP2` or similar) from a blood/vascular endothelial module, then gate both vascular and lymphatic refinement axes on `endothelial_core`.

- **axis:pericyte_smooth_muscle_vs_fibroblast / module:pericyte_smooth_muscle / module:stellate_myofibroblast** — Mural cells are contrasted directly against fibroblasts, and the pericyte/smooth-muscle and stellate/myofibroblast modules share several contractile markers (`ACTA2`, `TAGLN`, `MYL9`, `PDGFRB`). This makes the axes less interpretable and can blur pericytes, vascular smooth muscle, and activated pancreatic stellate cells.
  *Fix:* Separate pericyte (`RGS5`, `CSPG4`, `MCAM`, `NOTCH3`, `PDGFRB`) from smooth-muscle/contractile (`MYH11`, `ACTA2`, `TAGLN`, `MYL9`) emphasis where possible, and make mural identity a plain or stromal-gated refinement instead of a contrast against fibroblast matrix. Keep activated stellate/myofibroblast as a distinct stromal state if desired.

- **axis:endocrine_progenitor_like_state / module:endocrine_progenitor_like** — The module mixes ductal/progenitor markers (`SOX9`, `HNF1B`) with endocrine lineage markers (`NEUROG3`, `INSM1`, `RFX6`, `PAX6`, `ISL1`). For a no-age-specified healthy pancreas scenario this is ambiguous: it partially covers fetal/developmental endocrine progenitors, but the blueprint does not otherwise provide a coherent developmental pancreas branch, and in adult healthy pancreas this can instead look like ductal/endocrine contamination or stress.
  *Fix:* Retune this module to a clearly named rare developmental/immature endocrine program if all ages are meant to be covered, and consider adding a small complementary immature endocrine/beta axis. If the intended operating data are adult healthy pancreas, narrow or de-emphasize this axis rather than treating it as a core endocrine state.

- **axis:beta_maturation / module:beta_mature_secretory / module:beta_cell** — `beta_mature_secretory` is mostly a subset of `beta_cell`, so the `beta_maturation` dimension may mostly duplicate beta identity rather than cleanly resolve beta maturity.
  *Fix:* Make the positive module more specifically mature beta (`UCN3`, `MAFA`, `SLC30A8`, `G6PC2`, `ABCC8`, `KCNJ11`, plus other mature-function markers) and use a negative module that better represents immature/dedifferentiated beta or endocrine progenitor state rather than a broad mixed progenitor module.

- **module:pancreatic_epithelial / axis:epithelial_identity** — The pan-epithelial gate contains ductal-skewed genes (`KRT19`, `MUC1`, `CLDN4`) and may not behave as a neutral epithelial identity score across acinar, ductal, and endocrine cells.
  *Fix:* Use a cleaner epithelial/core gate with broadly expressed epithelial junction/cytokeratin genes, or avoid relying on this axis as the main gate for all epithelial refinements where acinar/endocrine cells may express fewer ductal keratins.

## Out of scope (FYI, not actionable for the Generator)
- The stress-state axes are generic but not inherently wrong; whether to keep all of interferon, inflammatory, ER stress, oxidative stress, hypoxia/glycolysis, and senescence depends on the expected dataset quality and experimental context.

## Summary for the next round
Focus on tightening the endothelial gating, making the stromal/mural axes biologically cleaner, and resolving the ambiguous progenitor/immature endocrine handling. The overall compartment coverage is good, so the next Generator should avoid wholesale redesign and instead make targeted module and axis refinements.
