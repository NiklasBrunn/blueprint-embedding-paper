# Discriminator report

## Scenario echo
The scenario is human pancreas with all expected cell types in scope because celltype, age, sex, condition, and notes are unspecified/null. The embedding should therefore cover the healthy baseline pancreas landscape across endocrine, exocrine, stromal, vascular, immune, neural, and generic cell states without adding disease- or perturbation-specific structure.

## Rating: minor
The blueprint is broadly plausible and covers the major pancreas compartments, including endocrine subtypes, acinar/ductal epithelium, vascular and stromal cells, immune cells, neural/glial cells, cell cycle, interferon response, and secretory stress. Targeted refinements are still needed for a few expected pancreas populations/states and for several marker choices that blur lineages or baseline subtype calls.

## Findings

### Minor issues
- **stromal / fibroblast_stellate / pericyte / adipocyte axes** — Healthy human pancreas commonly contains quiescent pancreatic stellate cells, which are not cleanly represented. The current `fibroblast_stellate` module is mostly generic ECM fibroblast/activated matrix, while quiescent stellate markers such as RGS5/ADIRF/GPX3/IGFBP4/FABP4-like lipid-storage programs are split across `pericyte` and `adipocyte`, risking misassignment of quiescent stellate cells.
  *Fix:* Add or refocus a quiescent pancreatic stellate module/axis, and keep `matrix_remodeling_within_stellate_fibroblast` as the activated/remodeling state rather than the only stellate refinement.
- **lymphatic_endothelial module** — Most markers are canonical lymphatic endothelial markers, but `TFF3` is an epithelial/secretory marker and is out of place in a lymphatic endothelial module.
  *Fix:* Replace `TFF3` with a lymphatic endothelial marker such as `RELN`, `FOXC2`, or `NRP2`.
- **venous_capillary_endothelial module** — `SELE` and `ICAM1` are endothelial activation/inflammation markers rather than baseline venous/capillary identity markers. In a healthy/null-condition pancreas blueprint, they make the vascular subtype axis partly state-driven.
  *Fix:* Replace `SELE` and `ICAM1` with more stable capillary/venous endothelial markers, e.g. `CA4`, `EDNRB`, `APLNR`, `ENG`, or `GNG11`, while retaining `ACKR1`, `PLVAP`, `RGCC`, and `RAMP3`.
- **T/NK lymphocyte axes** — The blueprint separates T cells and NK cells, but it does not explicitly capture cytotoxic T cells. The `t_cell` module is biased toward general/naive-like T markers (`IL7R`, `CCR7`, `LTB`), while cytotoxic genes are placed in `nk_cell`; CD8/cytotoxic T cells could therefore be pulled toward the NK contrast axis.
  *Fix:* Add a `cytotoxic_t_cell` gated axis under T cells using markers such as `CD8A`, `CD8B`, `GZMK`, `GZMA`, `CCL5`, `NKG7`, and/or `TRAC`, contrasted against NK markers where appropriate.
- **neural_glial / neuron_within_neural axes** — Pancreatic Schwann/glial cells and autonomic neurons are expected rare populations, but the current neuron axis is gated by `neural_glial`. True neuronal cells may score weakly if they do not coexpress enough glial markers, and `neural_glial` itself includes `NRXN1`, which is more neuronal than glial.
  *Fix:* Either make `neuron` a plain axis alongside a cleaned Schwann/glial axis, or add a broader neural gate that contains both Schwann/glial and neuronal markers before applying the neuron refinement.
- **myeloid_core / dendritic_cell modules** — `myeloid_core` is somewhat macrophage-skewed because of `C1QA`, `C1QB`, and `C1QC`, and `dendritic_cell` collapses cDC1, cDC2, and pDC markers into one module. This is acceptable for coarse immune coverage but underspecified for a full pancreas all-cell embedding.
  *Fix:* Make the core myeloid gate less macrophage-specific by emphasizing shared markers such as `LYZ`, `AIF1`, `TYROBP`, `FCER1G`, `LST1`, `HLA-DRA`, and `ITGAM`; optionally split dendritic cells into cDC1/cDC2/pDC if adding dimensions is acceptable.

## Out of scope (FYI, not actionable for the Generator)
- The current inclusion of rare endocrine epsilon cells, mesothelial cells, adipocytes, cycling cells, interferon response, and secretory/ER load is defensible for an all-cell human pancreas embedding and does not require deletion.

## Summary for the next round
Focus on targeted refinements rather than structural replacement: add/refocus quiescent pancreatic stellate coverage, add cytotoxic T-cell resolution, and clean the vascular/neural/myeloid marker sets so baseline subtype axes are not driven by epithelial or inflammatory-state genes.
