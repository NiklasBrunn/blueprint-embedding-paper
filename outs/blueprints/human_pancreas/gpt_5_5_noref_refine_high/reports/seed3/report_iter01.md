# Discriminator report

## Scenario echo
The scenario is healthy/baseline human pancreas, with no cell-type restriction, no specified age or sex, and no additional notes. The embedding should therefore cover the expected broad pancreatic compartment landscape: acinar, ductal/centroacinar, endocrine islet subtypes, vascular and lymphatic endothelium, mural cells, stromal/stellate cells, neural/Schwann-associated cells, adipocytes, and resident immune populations/states.

## Rating: minor
The blueprint is broadly plausible and structurally coherent for human pancreas, with good coverage of the major epithelial, endocrine, stromal, neural, and immune compartments. It still needs targeted refinement for a few coverage and marker/axis-design issues before I would call it publishable unchanged.

## Findings

### Minor issues
- **vascular compartment** — The blueprint has a generic blood endothelial module/axis but no lymphatic endothelial representation. Lymphatic endothelial cells are expected resident vascular cells in pancreatic tissue and are not captured by markers such as PECAM1/VWF/CDH5/PLVAP alone.
  *Fix:* Add a lymphatic endothelial module with canonical markers such as PROX1, PDPN, LYVE1, FLT4, CCL21, MMRN1, and/or RELN, plus an endothelial-gated lymphatic axis or a blood-vs-lymphatic endothelial contrast.
- **immune_common / immune gated axes** — The immune gate is partly pan-immune but is heavily weighted toward myeloid/APC genes (LST1, TYROBP, FCER1G, HLA genes). This can suppress T, NK, B, and plasma-cell axes because the gate itself will score much better in myeloid cells than in lymphocytes.
  *Fix:* Replace or rebalance immune_common toward pan-leukocyte genes such as PTPRC, CD52, LAPTM5, LCP1, CORO1A, FYB1, and similar broad immune markers; keep myeloid/APC markers in the myeloid/DC modules.
- **immune_b_plasma_cell axis/module** — B cells and plasma cells are merged into one module containing mutually divergent markers (MS4A1/CD79A/CD79B/BANK1 versus MZB1/JCHAIN/XBP1/IGHG1/IGKC). A mean-expression module will under-score each pure population and loses the biologically meaningful B-to-plasma distinction.
  *Fix:* Split into separate immune-gated B-cell and plasma-cell axes/modules, for example MS4A1/CD79A/CD79B/BANK1/CD74 for B cells and MZB1/JCHAIN/XBP1/SDC1/IGHG1/IGKC for plasma cells.
- **dendritic_antigen_presenting module** — This combines conventional DC markers (FCER1A, CLEC10A, CD1C, ITGAX) with plasmacytoid DC markers (LILRA4, CLEC4C). If both are present, a single averaged module can blur distinct immune subtypes.
  *Fix:* Either split cDC and pDC into separate immune-gated axes, or narrow the module to the main pancreas DC population intended by the axis.
- **epithelial_endocrine_vs_exocrine axis** — The axis is named as endocrine versus exocrine, but its negative side is only acinar_zymogen and does not include ductal epithelium, another major exocrine compartment. The current behavior is more accurately endocrine versus acinar.
  *Fix:* Rename the axis to reflect endocrine_vs_acinar, or redesign the exocrine contrast so ductal identity is also accounted for through an additional axis or revised topology.
- **epithelial_secretory_er_load axis** — This axis is gated only by endocrine_common but is named epithelial_secretory_er_load. It therefore does not represent all epithelial cells and is redundant in name with the separate acinar_secretory_er_load axis.
  *Fix:* Rename it to endocrine_secretory_er_load, or add a true ductal/epithelial-scoped ER-load axis if broad epithelial secretory stress is intended.
- **rare endocrine subtype sibling structure** — The epsilon_ghrelin and enterochromaffin_serotonin axes are plausible rare endocrine/state-like refinements, but they are not treated as sibling negatives of each other despite sharing endocrine/neuroendocrine genes. The enterochromaffin module also needs to be framed carefully because TPH1/serotonin programs can overlap pancreatic endocrine states rather than a canonical adult pancreas cell type.
  *Fix:* Include epsilon_ghrelin_cell and enterochromaffin_serotonin_cell in each other's endocrine subtype negative sets where appropriate, and tighten the serotonin module toward a clearly intended rare EC-like endocrine state.

## Out of scope (FYI, not actionable for the Generator)
- Healthy pancreas atlases vary in whether they recover very rare populations such as mesothelial/capsular cells or fetal-like serotonin endocrine states; absence of a note means these should not dominate dimensionality.

## Summary for the next round
Focus on adding lymphatic endothelial coverage, making the immune gate genuinely pan-leukocyte and splitting B versus plasma cells, then clean up misleading axis labels/design around endocrine-vs-exocrine and secretory ER load. The current blueprint is close enough that targeted edits should be sufficient; it does not need wholesale rework.
