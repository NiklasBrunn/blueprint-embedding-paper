# Discriminator report

## Scenario echo
Human pancreas, with no age, sex, or disease condition specified and no special notes. Because `celltype` is null, the embedding should cover the expected healthy human pancreatic landscape: acinar and ductal epithelium, endocrine islet cell types, stromal/mural/mesothelial/adipose components, vascular and lymphatic endothelium, resident immune populations, peripheral neural/glial cells, and common assay/stress states.

## Rating: minor
The blueprint is broadly biologically plausible and covers the major human pancreatic compartments. Another targeted refinement round is warranted because one axis is a redundant mirror dimension and the mural-cell design is likely to under-represent pericytes/smooth muscle by gating them through a fibroblast-like ECM module.

## Findings

### Minor issues
- **axis:ductal_vs_acinar** — This is the exact reciprocal of `axis:exocrine_acinar_vs_ductal`, so it adds a padded dimension rather than new biological resolution. A single acinar-versus-ductal contrast is enough; if additional exocrine resolution is desired, it should capture centroacinar/ductal/tuft or acinar state biology rather than the same contrast with reversed sign.
  *Fix:* Delete `axis:ductal_vs_acinar` and keep one acinar/ductal contrast axis, or replace the redundant dimension with a non-duplicative exocrine state/refinement axis.

- **axis:stromal_pericyte_smooth_muscle_refinement** — The pericyte/smooth-muscle program is gated by `mesenchymal_matrix`, whose markers are fibroblast/ECM-heavy (`COL1A1`, `DCN`, `LUM`, `PDGFRA`, `DPT`, etc.). Human pancreatic mural cells can be `RGS5`/`PDGFRB`/`CSPG4`/`MCAM`/`NOTCH3` or contractile `ACTA2`/`MYH11` without strong fibroblast matrix identity, so this gate can suppress an expected vascular mural population.
  *Fix:* Represent mural cells with an ungated plain axis or introduce a broader stromal/mural gate that includes mural markers, then keep fibroblast/stellate and myofibroblast programs as separate refinements.

- **module:epsilon_cell / axis:endocrine_epsilon_refinement** — Epsilon/ghrelin cells are rare in healthy human pancreas, especially outside fetal/young contexts, and the module is thin. `GHRL` is appropriate, but `MLN` is more characteristic of intestinal enteroendocrine motilin cells than pancreatic epsilon identity, while `SLC18A2` and `PCSK1` are shared endocrine markers.
  *Fix:* Either make this a conservative rare endocrine/ghrelin axis with stronger pancreatic epsilon specificity and remove `MLN`, or de-emphasize it if the Generator chooses to focus the unspecified-age baseline on common adult pancreatic cell types.

- **module:t_nk_cell / axis:immune_t_nk_refinement** — The module name says T/NK, but the marker list is essentially T-cell/lymphocyte (`CD3D`, `CD3E`, `TRAC`, `TRBC1`, `TRBC2`, `IL7R`, `CCR7`, `LTB`), while NK/cytotoxic identity is already handled separately by `nk_cytotoxic_cell`.
  *Fix:* Rename or reinterpret this as a T-cell axis, or broaden it deliberately if a combined lymphoid axis is intended; keep NK/cytotoxic cells separated by the existing NK module.

## Out of scope (FYI, not actionable for the Generator)
- The blueprint does not need disease-specific axes such as acinar-to-ductal metaplasia or activated inflammatory fibrosis for this healthy/baseline scenario, though the generic stress axes are reasonable for scRNA-seq data quality/state capture.

## Summary for the next round
Focus on removing the redundant reciprocal acinar/ductal dimension, fixing mural-cell representation so pericytes/smooth muscle are not hidden behind a fibroblast ECM gate, and tightening the rare epsilon and T-cell labeling/markers. The rest of the blueprint has good coverage of the expected human pancreas cell landscape.
