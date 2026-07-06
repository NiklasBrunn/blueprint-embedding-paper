# Discriminator report

## Scenario echo
The scenario is human pancreas with celltype, age, sex, and condition unspecified, so the embedding should cover the baseline human pancreatic cell-type landscape across endocrine islet, exocrine acinar/ductal, vascular, stromal/mural, neural, and resident immune populations. There are no notes imposing an extra topology or emphasis.

## Rating: minor
The blueprint has a sound overall structure and covers most expected pancreatic compartments. It still needs targeted biological refinements: some resident immune populations are missing or blurred, mural cells are gated in a way that can suppress them, and one endothelial specialization module is not marker-specific enough.

## Findings

### Minor issues
- **immune_leukocyte / immune subtype axes** - The resident immune compartment resolves macrophage/myeloid, T, B, and mast cells, but it does not explicitly capture dendritic/APC subsets or cytotoxic/NK-like lymphocytes that are commonly present in human pancreas atlases, even if sparse. The current `macrophage_myeloid` module is C1Q/TYROBP macrophage-skewed, and the `t_cell` module lacks cytotoxic/NK markers such as `NKG7`, `GNLY`, `PRF1`, or `KLRD1`.
  *Fix:* Add one or two gated leukocyte axes for cytotoxic/NK-like lymphocytes and conventional/plasmacytoid dendritic or antigen-presenting myeloid cells, using supported markers such as `NKG7`, `GNLY`, `PRF1`, `GZMB`, `KLRD1`, `FCGR3A` for cytotoxic/NK and `FCER1A`, `CD1C`, `CLEC10A`, `CLEC9A`, `LILRA4`, `IRF8` for DC subsets where literature support is available.
- **stromal_pericyte_smooth_muscle_identity** - The axis gates `pericyte_smooth_muscle` by `stromal_matrix`, a collagen/DCN/LUM/PDGFRA fibroblast-like module. Pericytes and vascular smooth muscle cells are mural cells and may not score strongly on that fibroblast matrix gate, so the existing gated axis can under-call an expected pancreatic mural compartment despite having a plausible mural marker module.
  *Fix:* Make mural identity a plain compartment axis or gate it with a mural-compatible core rather than a fibroblast ECM module. Keep stellate/fibroblast refinement under the stromal matrix compartment, and consider adding a separate activated stellate/fibroblast state if supported.
- **acinar_endothelium / acinar_endothelial_specialization** - The `acinar_endothelium` module mixes plausible endothelial/capillary genes with broad immediate-early and activation genes (`JUN`, `FOS`, `FOSB`, `KLF2`, `KLF4`) plus `CD74`. As written, the gated axis may capture dissociation stress, flow response, or immune-interacting endothelium more than acinar-region endothelial specialization.
  *Fix:* Replace the positive module genes with more anatomically or subtype-specific acinar endothelial markers from the pancreas endothelial source. If stress/activation is intended, split it into a separately named endothelial activation/state axis rather than using it as the acinar endothelial identity axis.
- **acinar_regenerating_substate** - The REG-family module is plausible for acinar stress/regeneration, but contrasting it against the whole ductal epithelial module makes the axis partly an acinar-vs-ductal separator rather than a state scoped cleanly within acinar cells. In a baseline healthy pancreas embedding, this should be clearly framed as a rare stress/regeneration program, not as a major lineage distinction.
  *Fix:* Rename or tune the module/axis toward acinar REG/stress/regeneration and consider leaving `negative_modules=[]` under the acinar gate, or add a separate ductal injury/SPP1-MMP7 state if the intent is to model both epithelial stress programs.
- **gamma_pp_cell / epsilon_cell** - The PP/gamma and epsilon endocrine modules are very thin, especially `epsilon_cell` with only `GHRL`. These are rare populations, but one- or two-gene axes are brittle and can be sensitive to dropout or dataset gene naming.
  *Fix:* Expand these modules with additional literature-supported human pancreatic PP and epsilon markers if available. If support is weak for adult/baseline samples, keep them but mark the axes as rare/low-resolution by using conservative sibling contrasts and avoiding overinterpretation.

## Out of scope (FYI, not actionable for the Generator)
- Citation DOI verification was not assessed here. The report only evaluates biological coverage, axis design, and marker plausibility.

## Summary for the next round
The next Generator should keep the overall pancreas-wide structure, then refine three areas: add missing immune resolution for cytotoxic/NK and DC/APC populations, change mural/pericyte identity so it is not dependent on a fibroblast ECM gate, and cleanly separate acinar endothelial identity from generic endothelial stress/activation. The REG acinar and rare endocrine axes only need lighter marker/axis tuning.
