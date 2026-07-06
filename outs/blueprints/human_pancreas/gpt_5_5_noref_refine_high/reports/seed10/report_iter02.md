# Discriminator report

## Scenario echo
The scenario is human pancreas, with no age, sex, or disease restriction and no additional notes. Because `celltype` is null, the embedding should cover the normal pancreas cell-type landscape: endocrine islet cells, exocrine acinar and ductal/centroacinar epithelium, vascular and lymphatic endothelium, mural/stromal cells, resident immune populations, neural/Schwann cells, adipocytes, and mesothelial/surface-associated cells where present.

## Rating: minor
The blueprint is broadly biologically plausible and on-scope for human pancreas, with good coverage of the major epithelial, endocrine, stromal, vascular, immune, neural, adipocyte, and mesothelial compartments. Targeted refinement is still warranted because several expected cell populations are collapsed into mixed modules/axes, especially T versus NK cells, B versus plasma cells, and blood endothelial subtypes.

## Findings

### Minor issues
- **axis:immune_t_nk_cell / module:t_nk_cell** — T cells and NK cells are represented as one mixed module containing both TCR/CD3 genes and cytotoxic/NK genes. These are distinct expected immune populations in human pancreas scRNA-seq and should be separable rather than occupying one blended axis.
  *Fix:* Replace the mixed T/NK representation with separate immune-gated T-cell and NK-cell axes. Use T-cell markers such as CD3D, CD3E, TRAC, IL7R, CCR7, CD2 or CD247 for T cells, and NK/cytotoxic markers such as NKG7, GNLY, PRF1, GZMB, KLRD1, FCGR3A or NCR1 for NK cells.

- **axis:immune_b_plasma_cell / module:b_plasma_cell** — B cells and plasma cells are similarly collapsed into one module combining MS4A1/CD79A/CD79B with MZB1/JCHAIN/XBP1/immunoglobulin genes. A normal pancreas embedding should distinguish lymphoid B cells from antibody-secreting plasma cells when both are present.
  *Fix:* Replace the mixed B/plasma representation with separate immune-gated B-cell and plasma-cell axes. Use MS4A1, CD79A, CD79B, BANK1, CD74 or HLA-DRA for B cells, and MZB1, JCHAIN, XBP1, IGKC, IGHG1 or SDC1 for plasma cells.

- **vascular_blood_endothelium / module:endothelial_blood** — The blood endothelial module mixes capillary/general endothelial markers (PLVAP, EMCN, RGCC, APLNR), venous markers (ACKR1), and arterial markers (GJA5) into a single axis. This captures blood endothelium broadly but loses expected artery/capillary/venous resolution in pancreas.
  *Fix:* Keep a broad endothelial identity if desired, but add endothelial-core-gated subtype axes for capillary/fenestrated endothelium, arterial endothelium, and venous/postcapillary endothelium. For example, use PLVAP/RGCC/EMCN/CA4/APLNR for capillary, GJA5/EFNB2/SOX17/DLL4/HEY1 for arterial, and ACKR1/SELE/VCAM1/NR2F2 for venous where appropriate.

- **mural_pericyte_vs_smooth_muscle** — The axis encodes pericyte positive versus smooth muscle negative, so smooth muscle is only represented as the opposite pole of a pericyte axis. Vascular smooth muscle is an expected pancreas-associated mural population and deserves a directly interpretable positive axis.
  *Fix:* Add a smooth-muscle-positive gated axis under `vascular_mural`, or replace the current single contrast-style gated axis with two sibling gated axes: one pericyte axis and one smooth muscle axis.

- **module:epsilon_cell** — The epsilon/ghrelin module contains canonical GHRL and MBOAT4, but the rest of the module is dominated by broad endocrine developmental/transcriptional markers shared across endocrine lineages (ISL1, PCSK1, PAX6, RFX6). This makes an already rare cell type less specific.
  *Fix:* Retain GHRL and MBOAT4, but improve specificity with additional ghrelin/epsilon-associated markers where available and reduce reliance on pan-endocrine transcription factors that are already present in the gate or sibling modules.

## Out of scope (FYI, not actionable for the Generator)
- The blueprint does not explicitly model disease, regeneration, fibrosis, or developmental trajectories, which is appropriate because the scenario gives no condition, age, or notes requiring those states.

## Summary for the next round
Focus on splitting mixed immune axes into separate T, NK, B, and plasma-cell representations, adding clearer blood endothelial subtype resolution, and making mural and epsilon-cell axes more directly interpretable. The major pancreas compartment structure should otherwise be preserved.
