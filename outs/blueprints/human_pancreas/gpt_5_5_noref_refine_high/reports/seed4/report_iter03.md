# Discriminator report

## Scenario echo
The scenario is human pancreas with celltype, age, sex, and condition unspecified; condition null means healthy/baseline, and notes are null. The embedding therefore should cover the expected human pancreatic epithelial endocrine/exocrine lineages plus stromal, vascular, immune, neural, and other resident/supporting populations without adding disease- or development-only structure as a requirement.

## Rating: minor
The blueprint covers the major human pancreas compartments well and uses mostly plausible human markers. A few all-pancreas populations/states are still too coarse or under-resolved, but the overall structure is biologically sound and does not need a wholesale redesign.

## Findings

### Minor issues
- **immune myeloid compartment** — The blueprint has a C1Q/TYROBP-biased `myeloid_macrophage` module and an MHC-II macrophage state, but it lacks a separate monocyte/inflammatory myeloid axis. Healthy pancreas datasets can contain resident macrophages plus recruited/circulating monocyte-like cells, and the current module will tend to fold FCN1/S100A8/S100A9/LST1/VCAN/CD14/CCR2-like cells into macrophage or dendritic space.
  *Fix:* Add a monocyte/inflammatory myeloid gene module and a gated immune or myeloid axis contrasting it against resident C1Q macrophage and dendritic APC modules. Keep the existing macrophage module; no deletion is needed.
- **dendritic_apc module / immune_dendritic_apc_identity axis** — `dendritic_apc` mixes cDC2 markers (`CD1C`, `FCER1A`, `CLEC10A`), pDC markers (`LILRA4`, `IRF7`), cDC1 markers (`CLEC9A`, `BATF3`), and generic antigen-presentation genes. This is acceptable as a coarse APC signal, but for an all-cell-type pancreas embedding it under-resolves distinct dendritic subsets.
  *Fix:* Either split into cDC2, cDC1, and pDC modules with separate immune-gated axes, or keep the broad APC axis and add one or two subtype-refinement axes where dimensionality must remain compact.
- **stellate/fibroblast state resolution** — The blueprint distinguishes quiescent stellate cells from ECM fibroblast signal mainly through `stellate_quiescent_vs_fibroblast_ecm`, where the ECM/activated side is represented only as the negative end of the axis. Human pancreas stromal biology commonly separates quiescent stellate, fibroblast/ECM, and activated myofibroblast-like programs; the current representation may make activated stellate/fibroblast states hard to interpret.
  *Fix:* Add a positive mesenchymal-gated ECM remodeling / activated stellate-fibroblast state axis using markers such as `FN1`, `POSTN`, `THY1`, `PDPN`, `COL1A1`, `COL3A1`, `MMP2`, `ACTA2`, and `TAGLN`, while retaining the quiescent stellate axis.
- **venous_endothelial module** — The venous module includes `VCAM1`, `SELE`, and `SELP`, which can reflect endothelial activation or post-capillary inflammatory phenotype rather than baseline venous identity. In a healthy/baseline pancreas scenario this can blur venous identity with activation state.
  *Fix:* Bias the venous identity module toward steadier venous markers such as `ACKR1`, `NR2F2`, `VWF`, `IL33`, and related venous/endothelial genes, and, if desired, move `VCAM1`/`SELE`/`SELP` into a separate endothelial activation state axis.

## Out of scope (FYI, not actionable for the Generator)
- Rare ductal tuft/chemosensory cells and pancreas-region-specific endothelial states could be added in a more exhaustive atlas-style embedding, but they are not essential enough here to force another structural expansion.

## Summary for the next round
Focus on targeted refinements: add a monocyte/inflammatory myeloid axis, improve dendritic subset resolution, and make the activated/ECM-producing stellate-fibroblast state positively represented rather than only the negative side of a quiescent axis. The current major epithelial, endocrine, exocrine, vascular, stromal, immune, neural, adipocyte, mesothelial, and cycling coverage should otherwise be preserved.
