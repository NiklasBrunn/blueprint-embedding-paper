# Discriminator report

## Scenario echo
The scenario is human pancreas, with cell type, age, sex, and condition unspecified/null; this means a baseline healthy pancreas embedding should cover the major endocrine, exocrine epithelial, vascular, stromal, immune, neural, mesothelial, and adipose-associated compartments without adding disease- or treatment-specific programs. There are no notes imposing a special topology or emphasis.

## Rating: minor
The blueprint is broadly biologically plausible and covers the expected pancreas landscape, including acinar, ductal/centroacinar, endocrine alpha/beta/delta/PP/epsilon, endothelial/lymphatic, stromal, immune, mast, Schwann/neural, adipocyte, mesothelial, and proliferation programs. Targeted edits would improve completeness and axis behavior, but the overall structure is appropriate.

## Findings

### Minor issues
- **vascular_endothelial / lymphatic_endothelial_identity / endothelial modules** — Blood endothelial cells are represented mostly as one broad vascular axis, with a lymphatic refinement, but healthy human pancreas datasets often contain distinguishable capillary/fenestrated/islet endothelial, arterial, and venous/postcapillary endothelial programs. The current `endothelial_blood_vessel` module already mixes pan-endothelial, capillary, and vascular markers, so these subtypes may collapse into one dimension.
  *Fix:* Add one or more gated endothelial subtype axes under a pan-endothelial gate, for example capillary/fenestrated/islet endothelial using markers such as `PLVAP`, `EMCN`, `RAMP2`, `RGCC`, `CA4`; arterial using `GJA5`, `SOX17`, `EFNB2`, `HEY1`; and/or venous/postcapillary using `ACKR1`, `VWF`, `SELE`, `VCAM1`, depending on desired dimensionality.

- **pericyte_smooth_muscle_identity** — The pericyte/smooth muscle axis is gated by `mesenchymal_fibroblast`, whose markers are mainly fibroblast/ECM genes (`COL1A1`, `COL1A2`, `DCN`, `LUM`, `PDGFRA`). True mural cells marked by `RGS5`, `PDGFRB`, `CSPG4`, `NOTCH3`, `ACTA2`, `TAGLN`, or `MYH11` may not score strongly on that fibroblast gate, so the gating could suppress the very cells this axis is meant to identify.
  *Fix:* Add a broader mural/mesenchymal gate module or make the mural program a plain axis/contrast separate from fibroblast-stellate identity. If kept gated, use a gate that includes mural-compatible markers rather than only fibroblast ECM markers.

- **stellate_matrix_remodeling_state / mesenchymal_fibroblast** — Pancreatic stellate/fibroblast biology is represented mainly by ECM-remodeling markers, which captures activated fibroblast-like or matrix-producing states but under-resolves quiescent pancreatic stellate cells and vitamin-A/lipid-associated stellate biology expected in healthy pancreas.
  *Fix:* Add or revise a stellate/quiescent module with markers such as `RGS5`, `PDGFRB`, `DES`, `LRAT`, `RBP1`, `PDPN`, and `COL15A1` as appropriate, and contrast it against the matrix-remodeling/activated fibroblast program rather than relying only on a generic ECM module.

- **ductal_bicarbonate_mucin** — The ductal axis is usable, but the module blends baseline ductal transport markers (`KRT19`, `KRT7`, `SOX9`, `CFTR`, `SLC4A4`, `SLC26A6`, `MUC1`) with markers that can emphasize mucinous/metaplastic or injury-associated ductal programs (`TFF1`, `TFF2`, `SPP1`). In a healthy baseline pancreas embedding, this may slightly bias the ductal identity axis toward reactive duct-like cells.
  *Fix:* Keep the core ductal transport/epithelial markers, but consider moving `TFF1`/`TFF2` and possibly high `SPP1` emphasis into a separate ductal reactive/mucinous state axis if such a state is desired; otherwise replace them with more baseline ductal markers such as `KRT17`, `KRT23`, `MMP7`, `ANXA4`, or `HNF1B`.

## Out of scope (FYI, not actionable for the Generator)
- Endocrine epsilon cells are rare in adult pancreas, but age is unspecified here, so retaining an epsilon axis is reasonable.
- Immune subtyping is broad but acceptable for a healthy whole-pancreas embedding; deeper T-cell, macrophage, or dendritic-state resolution would be optional rather than required.

## Summary for the next round
Focus on improving vascular and stromal resolution without changing the overall architecture: add endothelial subtype detail if dimensionality allows, fix the pericyte/smooth muscle gating so mural cells are not dependent on a fibroblast ECM gate, and separate quiescent stellate from activated/matrix-remodeling fibroblast biology. The ductal module can be lightly cleaned up to reduce reactive/mucinous bias.
