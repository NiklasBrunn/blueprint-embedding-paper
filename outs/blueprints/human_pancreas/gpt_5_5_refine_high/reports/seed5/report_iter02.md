# Discriminator report

## Scenario echo
The scenario is healthy/baseline human pancreas, with all pancreatic cell types in scope and no specified age, sex, condition, or design notes. The embedding should therefore cover endocrine islet cells, exocrine acinar and ductal/centroacinar epithelium, stromal/vascular/neural compartments, and resident immune populations without emphasizing disease- or development-specific programs.

## Rating: minor
The blueprint covers the major expected human pancreas compartments and the overall gated structure is biologically plausible. Targeted edits are still warranted for a few marker choices, stromal/endothelial gate contrasts, and redundant immune axes, but the design does not need structural rework.

## Findings

### Minor issues
- **module:acinar_secretory** -- The acinar module is plausible but underuses canonical pancreatic acinar digestive-enzyme markers: amylase genes are absent, while REG1A/REG1B are more injury/regeneration-associated and not ideal as core acinar discriminators.
  *Fix:* Add or substitute canonical acinar markers such as AMY2A, AMY2B, CELA3A, and/or CELA3B while retaining strong existing markers like PRSS1/2, CPA1/2, CTRB1/2, PNLIP, CLPS, and CEL.

- **module:lymphatic_endothelial** -- TFF3 is not a canonical lymphatic endothelial marker and would raise concern because it is better known as an epithelial/secretory marker.
  *Fix:* Replace TFF3 with a more plausible lymphatic endothelial marker such as NRP2, FOXC2, SOX18, MAF, or ACKR4, keeping PROX1, PDPN, LYVE1, FLT4, CCL21, and MMRN1.

- **axis:lymphatic_endothelial_within_endothelium** -- The axis gates on endothelial cells but subtracts pericyte/smooth muscle, which is not the biological sibling population for lymphatic endothelial cells. This makes the contrast less interpretable than a blood-vs-lymphatic endothelial split.
  *Fix:* Either make this a state-scoped gated axis without negatives, or add a blood/vascular endothelial sibling module and use that as the negative module.

- **axis:contractile_stellate_myofibroblast** -- The contractile stellate/myofibroblast state is gated by matrix fibroblast/stellate cells but subtracts endothelial cells. Endothelium is not the relevant sibling state, and the current contrast does not cleanly separate quiescent/matrix stellate cells from activated contractile stellate/myofibroblast cells.
  *Fix:* Remove the endothelial negative and treat the contractile program as a gated state, or add a quiescent stellate/fibroblast sibling module for the negative side.

- **axis:pericyte_smooth_muscle_vs_endothelium** -- The axis name implies a contrast against endothelium, but the implementation is a plain pericyte/smooth muscle axis. The biology is fine, but the declared axis semantics are inconsistent.
  *Fix:* Rename it to `pericyte_smooth_muscle`, or convert it to an actual contrast axis against endothelial if the intended biological dimension is mural-vs-endothelial.

- **axis:t_nk_within_leukocytes and axis:b_plasma_within_leukocytes** -- The combined T/NK and B/plasma axes are redundant with the separate T, NK, B, and plasma cell axes and their negative sets mix coarse modules with fine sibling modules. This adds dimensionality without a clear biological state or lineage level.
  *Fix:* Prefer deleting the combined coarse immune modules/axes, or rework the immune section into a coherent two-level tree with parent lymphoid/B-lineage axes followed by gated fine subtypes.

## Out of scope (FYI, not actionable for the Generator)
- Citation DOI verification was not assessed here, per the discriminator scope.

## Summary for the next round
Keep the broad pancreas coverage and overall endocrine/exocrine/stromal/immune layout, but tighten the marker biology and axis semantics: add canonical acinar amylase/elastase markers, replace TFF3 in the lymphatic module, revise arbitrary endothelial/pericyte negatives in gated stromal/endothelial axes, and simplify or formalize the redundant immune hierarchy.
