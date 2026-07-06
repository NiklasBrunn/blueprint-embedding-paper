# Discriminator report

## Scenario echo
The scenario is human pancreas with no specified age, sex, disease, treatment, or restricted cell type, so the embedding should cover the expected healthy pancreatic landscape: acinar and ductal exocrine epithelium, endocrine islet cell types, vascular and stromal compartments, resident immune cells, adipose/mesothelial contributions, and innervation-associated cells. There are no notes imposing a special structural shape.

## Rating: minor
The blueprint covers the major healthy human pancreas compartments and uses mostly plausible canonical marker modules. A targeted refinement pass is still warranted because acinar identity is gated through a ductal-biased epithelial core, several vascular/mesenchymal identities are compressed, and one acinar state module mixes healthy adaptive markers with ductal/metaplastic markers.

## Findings

### Minor issues
- **axis:epithelial_acinar_identity / module:pancreatic_epithelial_core** -- The acinar identity axis is gated by `pancreatic_epithelial_core`, which contains several ductal-enriched epithelial markers (`KRT19`, `CLDN4`, `MUC1`) alongside broader epithelial genes. Human acinar cells are epithelial, but their scRNA signal is dominated by digestive enzyme and zymogen genes; this gate can unnecessarily dampen bona fide acinar cells and make ductal epithelial score the gate more strongly than acinar cells.
  *Fix:* Modify the acinar identity axis to avoid a ductal-heavy gate. Use a plain acinar axis, gate it with a more acinar/exocrine-appropriate module, or revise the epithelial gate so it is not dominated by ductal markers before using it for acinar identity.

- **axis:pericyte_smooth_muscle_vs_fibroblast / module:pericyte_smooth_muscle** -- Pericytes and vascular smooth muscle are both expected in human pancreas, but the current module and single contrast axis merge them (`RGS5`, `CSPG4`, `PDGFRB`, `NOTCH3` with `ACTA2`, `TAGLN`, `MYH11`, `MYL9`). That is biologically plausible as a mural-cell aggregate, but it loses a real vascular/stromal distinction.
  *Fix:* Split this into separate pericyte and smooth-muscle modules/axes, or add a smooth-muscle refinement axis within the mural compartment while preserving pericyte coverage.

- **module:mesenchymal_fibroblast_stellate / axis:fibroblast_stellate_identity** -- Fibroblasts and pancreatic stellate cells are represented only as one ECM-heavy module. Healthy pancreas commonly contains quiescent stellate-like cells as well as fibroblast/ECM states, and this axis does not resolve them.
  *Fix:* Add a stellate/quiescent-stellate refinement module/axis or split the current module into fibroblast/ECM and stellate-biased programs, using pancreas-appropriate markers such as `RBP1`, `DES`, `DCN`, `LUM`, `COL6A1`, `IGFBP4`, `ADIRF`, or related validated choices. Keep pericyte markers distinct from stellate markers where possible.

- **module:acinar_adaptive_reg / axis:acinar_adaptive_reg_state** -- `REG1A`, `REG1B`, `REG3A`, and `SPINK1` are reasonable acinar adaptive/regenerative markers, but `AGR2`, `MUC1`, and `TFF1` can pull the module toward ductal, mucinous, or metaplastic programs rather than a clean healthy acinar state. In a baseline healthy-pancreas scenario this state should not read like pancreatitis/ADM unless explicitly framed as an injury/dissociation-adaptive axis.
  *Fix:* Either rename and design this as an acinar injury/metaplastic/dissociation-adaptive state with an appropriate ductal/metaplastic contrast, or replace the ductal/mucinous genes with cleaner acinar stress/adaptive genes while keeping the REG/SPINK1 core.

- **ductal coverage** -- Ductal cells are represented by a broad identity axis and one bicarbonate/mucin state. This is acceptable for coarse coverage, but it under-resolves healthy ductal heterogeneity such as centroacinar/intercalated ductal programs versus larger duct/mucin or activated `SPP1`/`MMP7` states.
  *Fix:* Add one ductal refinement axis if dimensionality budget allows: for example a centroacinar/intercalated/bicarbonate transport program contrasted with a larger-duct or activated ductal program. Keep injury-associated markers from dominating the main ductal identity module.

## Out of scope (FYI, not actionable for the Generator)
- Reference context checked for expected human pancreas compartments and markers: Human Protein Atlas pancreas overview (https://v23.proteinatlas.org/humanproteome/tissue%2Bcell%2Btype/pancreas) and a recent human pancreas single-cell atlas summary noting endocrine, exocrine, stromal, vascular, immune, and ganglial populations (https://pmc.ncbi.nlm.nih.gov/articles/PMC11802906/).

## Summary for the next round
Focus on targeted refinement rather than rebuilding: adjust the acinar identity gate so true acinar cells are not suppressed by ductal epithelial markers, split or refine compressed mural/stellate stromal axes, and clean up the acinar adaptive plus ductal state modules so healthy-pancreas identity and injury/metaplastic signals are not conflated.
