# Discriminator report

## Scenario echo
The scenario is healthy/baseline human pancreas, with no cell type, age, sex, or condition restriction and no additional notes. The embedding should therefore cover the major epithelial, endocrine, stromal, vascular, neural, adipose/mesothelial, and resident immune compartments expected in pancreas without specializing to a disease or developmental trajectory.

## Rating: minor
The blueprint is broadly biologically plausible and covers the expected major human pancreas lineages, including endocrine islet subtypes, acinar, ductal/centroacinar, vascular, stromal/mural, neural, adipocyte/mesothelial, and immune cells. It still needs targeted refinement because the ductal/centroacinar separation is not specific, pancreatic stellate/fibroblast biology is partly conflated with mural/pericyte programs, and a few marker choices reduce specificity.

## Findings

### Minor issues
- **axis:ductal_centroacinar_terminal / module:centroacinar_terminal_duct** - The gated axis is intended to resolve centroacinar/terminal duct cells, but the module heavily overlaps with the general ductal module (CFTR, ANXA4, SLC4A4, SOX9, KRT19) and is contrasted only against acinar/endocrine modules. As written, it mostly reiterates ductal epithelial identity rather than distinguishing centroacinar/terminal duct from broader ductal cells.
  *Fix:* Make the axis explicitly discriminate terminal duct/centroacinar from general ductal cells, for example by adding `ductal_bicarbonate` as a negative module or by narrowing/renaming the module to a ductal-continuum state. Keep canonical terminal duct/centroacinar markers such as CFTR, CLDN10, AQP1/AQP5, ANXA4, and SOX9, but reduce direct duplication with the general ductal module where possible.
- **axis:stromal_fibroblast_vs_mural / modules:fibroblast_stellate, vascular_mural** - Pancreatic stellate cells are an important pancreatic stromal population, and their marker profile can overlap with pericyte/mural programs (RGS5, PDGFRB, ACTA2, DES in activated states). The current fibroblast/stromal axis is mostly extracellular matrix fibroblast markers, while the mural module contains genes that can also score stellate/perivascular stromal cells. This may force pancreatic stellate cells into a fibroblast-vs-mural contrast that is biologically too coarse.
  *Fix:* Refine stromal representation by adding or modifying a pancreatic stellate/perivascular stromal module with PSC-relevant markers such as LRAT, RBP1, PDGFRB, RGS5, COL1A1/COL3A1, DCN/LUM, and optionally ACTA2/DES for activated stellate/myofibroblast signal. Use gated stromal axes to separate fibroblast/stellate from true vascular mural cells rather than relying only on a broad contrast.
- **module:epsilon_cell** - The epsilon module includes strong ghrelin-cell markers (GHRL, MBOAT4), but several genes in the list (TM4SF4, TM4SF5, HES4, SSTR1, SSTR2) are not specific canonical epsilon markers and may pull the axis toward ductal/epithelial or general endocrine receptor programs. Epsilon cells are rare in human pancreas, so marker specificity matters.
  *Fix:* Keep GHRL and MBOAT4 and replace less-specific genes with more defensible epsilon/ghrelin-cell markers where available, or keep the module compact rather than padded. Avoid epithelial/ductal-associated genes unless there is a clear epsilon-cell rationale.
- **module:endothelial_venous_postcapillary** - ACKR1 and NR2F2 are plausible venous/postcapillary endothelial markers, but SELE and VCAM1 are also endothelial activation/inflammation markers. In a healthy baseline pancreas scenario, this could make the axis partly reflect activation rather than venous identity.
  *Fix:* Rebalance this module toward venous/postcapillary identity markers and, if activation is worth preserving, represent it as a separate gated endothelial state rather than mixing it into venous subtype identity.

## Out of scope (FYI, not actionable for the Generator)
- No citations are required or expected for this blueprint variant.

## Summary for the next round
Focus on targeted refinements rather than redesign: improve the ductal versus centroacinar/terminal duct axis, separate pancreatic stellate/perivascular stromal identity from true mural/pericyte identity, and tighten the less-specific epsilon and venous endothelial marker modules. The overall lineage coverage and gated subtype structure are otherwise suitable for healthy human pancreas.
