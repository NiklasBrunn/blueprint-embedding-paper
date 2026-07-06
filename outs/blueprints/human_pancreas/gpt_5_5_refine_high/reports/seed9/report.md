# Discriminator report

## Scenario echo
The scenario is healthy/baseline human pancreas, with age, sex, and condition unspecified and `celltype=null`, so the embedding should cover the broad adult human pancreatic cell landscape: acinar, ductal, endocrine islet subtypes, stromal/mesenchymal, vascular and lymphatic endothelium, resident immune cells, neural/glial elements, mesothelium/capsule, and any well-supported rare resident populations. There are no notes imposing an extra structural preference.

## Rating: minor
The blueprint is broadly plausible and structurally valid: it captures the major exocrine, endocrine, stromal, vascular, immune, glial, and mesothelial compartments with mostly canonical human markers. It still needs targeted refinement for a few expected subpopulations and marker/axis naming choices before I would defend it unchanged.

## Findings

### Minor issues
- **myeloid_macrophage / immune coverage** — The only myeloid axis is a macrophage-skewed module (`LYZ`, `CD68`, `CD163`, `CSF1R`, `C1QA/B`, `MSR1`, `AIF1`). Healthy human pancreas atlases can contain dendritic/antigen-presenting myeloid populations that are not well represented by this resident macrophage signature.
  *Fix:* Add a small dendritic/APC-focused module and axis if supported by the cited pancreas atlas literature, or broaden the myeloid design to distinguish macrophage from DC-like antigen-presenting cells. Candidate high-level markers to consider include `FCER1A`, `CD1C`, `CLEC10A`, `CLEC9A`, `XCR1`, `LILRA4`, and HLA class II genes, chosen only where literature-supported.
- **t_nk_lymphocyte / nk_cell** — The `t_nk_lymphocyte` module is actually T-cell focused (`CD3D`, `CD3E`, `TRAC`, `TRBC1/2`, `IL7R`, `CCR7`) while NK cells already have a separate `nk_cell` axis. The current name implies a mixed T/NK population and may obscure interpretation.
  *Fix:* Rename/refocus `t_nk_lymphocyte` to a T-cell module/axis, or explicitly redesign it as a pan-cytotoxic/T-NK gate with a separate T-cell subtype axis. The simpler fix is to make it `t_cell` and keep `nk_cell` separate.
- **epsilon_cell** — The rare epsilon/ghrelin cell axis is reasonable to include, but the module mixes strong ghrelin-cell genes (`GHRL`, `MBOAT4`) with less canonical or potentially non-pancreas-enter endocrine/metabolic genes (`MLN`, `FFAR3`, `CD320`, `ACSL1`). For a very rare subtype, marker specificity matters.
  *Fix:* Recheck and tighten the epsilon module against human pancreas/islet evidence. Keep high-confidence ghrelin-cell markers and replace weak entries with better-supported epsilon markers, or reduce the module to fewer high-confidence genes rather than padding it.
- **stellate_fibroblast_ecm** — The stromal axis primarily scores broad collagen/ECM fibroblast activation (`COL1A1`, `COL1A2`, `COL3A1`, `MMP2`, `SPARC`, `VIM`) and may under-resolve pancreatic stellate quiescent versus activated fibroblast/myofibroblast-like states. In healthy pancreas, this is not a fatal gap, but it is a meaningful state distinction for the mesenchymal compartment.
  *Fix:* Either split the stromal design into quiescent stellate/fibroblast and activated ECM/myofibroblast axes, or adjust marker balance so the existing axis is clearly a resident stellate/fibroblast identity axis rather than mostly an ECM activation score.

## Out of scope (FYI, not actionable for the Generator)
- I did not re-verify paper or atlas DOI metadata; citation verification is outside this discriminator pass.

## Summary for the next round
Keep the current broad architecture, but refine the immune and stromal resolution and clean up rare/subtype marker specificity. The highest-priority fixes are adding or representing dendritic/APC myeloid coverage, renaming/refocusing the T-cell axis, and tightening the epsilon-cell marker list.
