# Discriminator report

## Scenario echo
The scenario is healthy/baseline human pancreas, with no cell type, age, sex, or condition restriction and no additional notes. The embedding should therefore cover the major whole-pancreas cell landscape: exocrine acinar and ductal/centroacinar epithelium, endocrine islet cell types, pancreatic stromal/stellate and mural cells, blood and lymphatic endothelium, resident immune cells, and neural/glial elements.

## Rating: minor
The blueprint is biologically plausible overall and covers most expected human pancreatic compartments with sensible gated subtype axes. It still has several targeted gaps and marker-design issues that should be refined before calling it publishable as-is.

## Findings

### Minor issues
- **general coverage / epithelial-stromal interface** — A whole human pancreas embedding can reasonably encounter mesothelial/capsular cells, but there is no axis or module for this compartment. The current stromal and ductal axes may not separate MSLN/WT1/UPK3B/CALB2-like mesothelial cells from fibroblast/ductal signal.
  *Fix:* Add a small mesothelial/capsular module and plain identity axis, using human mesothelial markers such as MSLN, WT1, UPK3B, CALB2, KRT8/KRT18 as appropriate while avoiding excessive overlap with ductal markers.
- **vascular endothelial_identity / pan_endothelial / vascular_lymphatic_bias** — The endothelial coverage separates blood versus lymphatic endothelium, but it does not resolve the major blood endothelial states expected in pancreas, especially fenestrated/islet capillary endothelium versus arterial/venous programs. For a whole-tissue pancreas embedding, this leaves vascular resolution somewhat coarse.
  *Fix:* Add one or two gated blood-endothelial refinement axes, for example a fenestrated/capillary-islet endothelial module with PLVAP/RAMP2/RGCC/EMCN-like markers and, if dimensionality remains reasonable, arterial/venous bias modules with canonical markers. Keep these gated by a pan-endothelial module.
- **gamma_pp_cell and epsilon_ghrelin_cell modules** — These rare endocrine subtype modules are heavily padded with generic neuroendocrine genes shared across all islet endocrine cells. Because the subtype axes are gated by endocrine_islet, this shared CHGA/CHGB/SCG/SYP content weakens the ability to distinguish PP and ghrelin cells from other endocrine cells.
  *Fix:* Sharpen these modules toward subtype-enriched markers. Keep PPY/PYY for PP/gamma cells and GHRL/MBOAT4 for epsilon cells, but reduce generic endocrine genes and add additional subtype-biased human markers where defensible.
- **quiescent_pancreatic_stellate module** — The quiescent stellate module includes several mural/pericyte-associated markers such as RGS5, PDGFRB, and DES, while a separate mural/pericyte branch already exists. This risks conflating quiescent pancreatic stellate cells with pericytes rather than cleanly representing a stellate quiescence state.
  *Fix:* Rebalance the quiescent stellate module toward retinoid/lipid-storage and resting stromal markers such as LRAT, RBP1, APOD, DCN, LUM, and CXCL14-like genes, and reduce mural/pericyte-heavy markers unless they are intentionally used with stronger negative gating.
- **beta_maturity_within_islets / immature_stressed_beta_cell** — The negative module mixes immaturity, dedifferentiation/glycolytic shift, ER stress, heat-shock, and oxidative stress genes into one "immature_stressed" concept. In a healthy pancreas scenario this is plausible as a stress/dedifferentiation axis, but the label and marker set are biologically broad.
  *Fix:* Either rename/recast the module as beta stress/dedifferentiation or split maturity and stress into separate beta-gated axes if the Generator wants higher state resolution. Keep the axis scoped to beta cells.

## Out of scope (FYI, not actionable for the Generator)
- None.

## Summary for the next round
The Generator should keep the current overall structure, but add a mesothelial/capsular compartment, modestly refine blood endothelial substructure, and clean up subtype/state marker specificity for PP/epsilon endocrine cells, quiescent stellate cells, and the beta maturity/stress axis. No existing axis or module needs deletion.
