# Discriminator report

## Scenario echo
The scenario is healthy/baseline human pancreas, with no specific age, sex, condition, or cell-type restriction and no special notes. The embedding should therefore cover the full expected pancreatic tissue landscape: exocrine acinar and ductal/centroacinar epithelium, endocrine islet subtypes, endothelial and mural vascular cells, fibroblast/stellate stroma, neural/glial cells, and resident immune populations.

## Rating: minor
The current blueprint broadly covers the major human pancreatic compartments and is structurally valid. It still needs targeted refinement for endothelial/mural resolution, immune coverage, and a few marker modules that mix canonical subtype markers with generic or questionable genes.

## Findings

### Minor issues
- **axis:vascular_lymphatic_bias / module:lymphatic_endothelial** — Lymphatic endothelium is present only as a gated positive program under `vascular_endothelial`, with no blood-endothelial negative module or independent identity axis. This makes lymphatics look like a state of vascular endothelial cells rather than a separable endothelial subtype expected in whole pancreas.
  *Fix:* Add a pan-endothelial gate if desired, keep blood vascular and lymphatic modules separate, and represent lymphatics with a plain `lymphatic_endothelial_identity` axis and/or a gated `lymphatic_vs_blood_endothelial` contrast using canonical lymphatic markers such as `PROX1`, `LYVE1`, `PDPN`, `FLT4`, and `CCL21`.

- **axis:mural_cell_bias / module:pericyte_smooth_muscle** — The module combines pericyte markers (`RGS5`, `PDGFRB`, `CSPG4`, `MCAM`, `NOTCH3`) with vascular smooth-muscle markers (`MYH11`, `CNN1`, `MYLK`, `ACTA2`, `TAGLN`), and the axis gates the module by itself while contrasting against endothelial/fibroblast modules. This is mostly redundant with `mural_cell_identity` and does not resolve the two main mural states.
  *Fix:* Split pericyte and vascular smooth-muscle programs and replace the self-gated mural bias with a biologically meaningful pericyte-vs-smooth-muscle gated or contrast axis.

- **immune subtype coverage** — The immune compartment includes macrophage/myeloid, T/NK, B/plasma, and mast axes, but dendritic/APC populations are not explicitly captured and T vs NK and B vs plasma are merged. For a whole-pancreas embedding, at least dendritic cells are a plausible resident/contaminating population that should not be forced into macrophage or generic leukocyte space.
  *Fix:* Add a dendritic/APC module and gated axis, for example separating conventional dendritic markers (`CD1C`, `FCER1A`, `CLEC10A`, `CLEC9A`, `XCR1`, `LAMP3`) from macrophage genes. If dimensionality permits, split `t_nk_lymphoid` into T-cell and NK/cytotoxic programs and split `b_plasma_cell` into B-cell and plasma-cell programs.

- **module:centroacinar_ductal** — The module mixes plausible centroacinar/ductal genes (`SOX9`, `KRT19`, `CA2`, `CLDN10`, `SLC4A4`) with genes that are more stress/injury or nonspecific ductal programs (`TFF1`, `MMP7`, `SPINK1`). In a baseline healthy pancreas scenario, this weakens the claim that the axis specifically resolves centroacinar cells.
  *Fix:* Tighten the module toward canonical ductal/centroacinar markers, for example emphasizing `CA2`, `CLDN10`, `CFTR`, `SOX9`, `KRT19`, `SLC4A4`, and related ductal transport/junction genes; reserve `TFF1`/`MMP7` for a separate injured/stress ductal axis only if intentionally modeling that state.

- **endocrine subtype marker hygiene** — Several endocrine subtype modules include generic endocrine genes already represented by the `endocrine_islet` gate (`CHGA`, `CHGB`, `PCSK1`) or less defensible subtype markers. Examples include `DLK1` in the core `beta_cell` module for an otherwise baseline pancreas embedding, `MLN` in `epsilon_ghrelin_cell`, and generic `CHGA`/`CHGB` in the PP and epsilon subtype modules. These will blur subtype-specific gated axes.
  *Fix:* Keep generic neuroendocrine genes mainly in `endocrine_islet`, and make subtype modules more specific. For beta cells, use mature beta-enriched markers unless an immature/fetal beta state is intentionally added. For PP and epsilon cells, rely more heavily on hormone and subtype-specific genes (`PPY`/`PYY` for PP; `GHRL` and other ghrelin-cell-enriched genes for epsilon) rather than generic endocrine markers.

## Out of scope (FYI, not actionable for the Generator)
- The blueprint has no placeholder or orphan-module problem; `python3 embedder_blueprint.py --describe` exits cleanly.

## Summary for the next round
Keep the overall structure, but refine three areas: make lymphatic and mural vascular axes biologically distinct rather than self-gated/coarse, add at least a dendritic/APC immune axis, and tighten ductal/centroacinar plus endocrine subtype gene modules so gated subtype dimensions are driven by specific markers rather than generic endocrine or stress genes.
