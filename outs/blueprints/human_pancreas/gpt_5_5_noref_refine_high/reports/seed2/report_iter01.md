# Discriminator report

## Scenario echo
The scenario is human pancreas with all pancreatic cell types in scope because `celltype` is null. Age and sex are unspecified, condition is baseline/healthy, and there are no notes requesting a special topology or emphasis.

## Rating: minor
The blueprint covers the major human pancreas compartments and has a reasonable gated/refinement structure. Targeted refinements are still needed because several expected pancreas populations are merged too coarsely and one acinar state axis is largely redundant with acinar identity.

## Findings

### Minor issues
- **t_nk_cell / cytotoxic_t_nk_state axes** — T cells and NK cells are merged into one identity module containing both CD3/TRAC markers and NK cytotoxic markers. This will identify lymphocytes broadly, but it does not let the embedding resolve canonical T-cell versus NK-cell populations expected in human pancreas immune content.
  *Fix:* Add separate `t_cell` and `nk_cell` modules and gated or contrast refinement axes within a lymphoid gate, using T markers such as `CD3D`, `CD3E`, `TRAC`, `IL7R` and NK markers such as `NKG7`, `GNLY`, `KLRD1`, `KLRF1`, `FCGR3A`, while keeping a cytotoxic state if desired.
- **b_plasma_cell axis** — B cells and plasma cells are combined in one module containing `MS4A1`/`CD79A` alongside `MZB1`/`JCHAIN`/immunoglobulin genes. These are distinct lymphoid states/cell types in pancreas single-cell data and should not be represented only as one blended axis.
  *Fix:* Split or refine into B-cell and plasma-cell axes, for example B-cell markers `MS4A1`, `CD79A`, `CD79B`, `BANK1`, `CD74` and plasma-cell markers `MZB1`, `JCHAIN`, `XBP1`, `SDC1`, immunoglobulin genes.
- **pericyte_smooth_muscle_identity** — Pericytes and vascular smooth muscle cells are merged into one mural-cell module. The listed genes are plausible, but pancreas vasculature commonly needs pericyte versus smooth-muscle resolution.
  *Fix:* Add gated refinement axes under the existing mural identity, separating pericyte markers such as `RGS5`, `PDGFRB`, `CSPG4`, `MCAM`, `KCNJ8` from smooth-muscle markers such as `ACTA2`, `MYH11`, `TAGLN`, `CNN1`, `MYL9`.
- **schwann_neural_identity** — The module mixes Schwann/glial markers (`SOX10`, `PLP1`, `MPZ`, `S100B`) with neuronal markers (`PRPH`) and does not separately capture intrinsic/autonomic neurons. Neural populations are rare but biologically expected in human pancreas.
  *Fix:* Keep or retune a Schwann/glial axis and add a separate neuronal refinement/identity module with markers such as `PRPH`, `TUBB3`, `ELAVL3`, `ELAVL4`, `PHOX2B`, `RET`, `CHAT` where appropriate.
- **acinar_zymogen_load** — `acinar_zymogen_high` strongly duplicates `exocrine_acinar`, reusing many of the same digestive enzyme genes. As a gated axis it mostly recapitulates acinar abundance rather than a distinct acinar state.
  *Fix:* Replace the module contents with a more distinct acinar state if intended, such as stress/regenerative/secretory-load markers that do not simply duplicate the identity module, or otherwise de-emphasize this axis.
- **dendritic_cell** — The dendritic module combines cDC2 markers (`CD1C`, `CLEC10A`, `FCER1A`) with pDC markers (`LILRA4`, `IRF7`, `TCF4`) into a single myeloid refinement. Both are plausible immune populations, but a single axis cannot distinguish them.
  *Fix:* Add separate cDC and pDC refinements within the myeloid gate if immune resolution is desired, while preserving the broader myeloid/macrophage/monocyte layout.

## Out of scope (FYI, not actionable for the Generator)
- Endocrine rare-cell coverage is somewhat generous, including epsilon and enterochromaffin-like axes; this is acceptable for an age-unspecified human pancreas scenario and does not require deletion.
- `pancreatic_progenitor_to_mature_endocrine` is acceptable given unspecified age, but it should not be expanded into a full developmental tree unless the scenario notes request developmental ordering.

## Summary for the next round
Focus on targeted subtyping rather than structural rework: split T versus NK, B versus plasma, pericyte versus smooth muscle, and glial versus neuronal programs; then make the acinar state axis less redundant with acinar identity. The existing major-compartment coverage and most marker choices are biologically plausible.
