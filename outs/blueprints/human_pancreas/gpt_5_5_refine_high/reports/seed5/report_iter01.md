# Discriminator report

## Scenario echo
Human pancreas with no specified age, sex, or disease condition, and no extra notes. Because `celltype` is null, the embedding should cover the healthy human pancreas cell landscape: endocrine islet cells, exocrine acinar and ductal/centroacinar epithelium, vascular and stromal cells, neural/glial cells, and resident immune populations.

## Rating: minor
The blueprint is broadly biologically plausible and covers the major pancreatic compartments, with a sensible endocrine/exocrine/stromal/immune layout. It still needs targeted refinement for expected minor populations and one immune-gating marker problem before I would call it complete.

## Findings

### Minor issues
- **Immune coverage** — Mast cells are an expected resident immune population in human pancreas atlases, but there is no mast-cell module or axis. The current leukocyte set resolves macrophage/monocyte, dendritic/APC, T/NK, and B/plasma only.
  *Fix:* Add a `mast_cell` module with canonical markers such as `KIT`, `TPSAB1`, `TPSB2`, `CPA3`, `MS4A2`, `HDC`, `GATA2`, and add a leukocyte-gated mast-cell axis.

- **Vascular coverage** — The endothelial module is blood-vascular biased (`VWF`, `KDR`, `PLVAP`, `EMCN`, `CLDN5`) and does not resolve lymphatic endothelium, which is within the healthy pancreas tissue landscape.
  *Fix:* Add a lymphatic endothelial module/axis, preferably gated by a broad endothelial module, using markers such as `PROX1`, `PDPN`, `LYVE1`, `FLT4`, `CCL21`, `MMRN1`, and `TFF3` where supported.

- **leukocyte_common / immune gated axes** — `leukocyte_common` is not truly pan-leukocyte: most genes are APC/myeloid/B biased (`LST1`, HLA class II genes, `CD74`, `TYROBP`). This gate can under-score T, NK, and mast cells and inflate myeloid/APC programs.
  *Fix:* Rework `leukocyte_common` toward pan-leukocyte markers such as `PTPRC`, `LCP1`, `CORO1A`, `LAPTM5`, `CD52`, `CD53`, `RAC2`; keep HLA/TYROBP genes in the dendritic/macrophage modules rather than the common gate.

- **pericyte_smooth_muscle_vs_endothelium axis** — Pericytes/smooth muscle and endothelial cells are adjacent vascular compartments, not opposite ends of a single biological program, and endothelial already has its own plain axis.
  *Fix:* Change this to a plain `pericyte_smooth_muscle` axis, or gate a mural-cell axis by a broader vascular/stromal context if the Generator wants to keep vascular scoping.

- **t_nk_cell and b_plasma_cell resolution** — These modules combine distinct populations. T cells and NK cells can be separated by `CD3D/CD3E/TRAC/IL7R` versus `NKG7/GNLY/KLRD1/NCR1`, and B cells and plasma cells by `MS4A1/CD79A/CD79B` versus `MZB1/JCHAIN/XBP1/IGHG*`.
  *Fix:* If dimensionality permits, split into separate T-cell, NK-cell, B-cell, and plasma-cell modules/axes. If kept aggregated, rename or document them as broad lymphoid axes and ensure the markers do not blur the intended score.

## Out of scope (FYI, not actionable for the Generator)
- I did not verify DOI correctness or whether every cited resource contains each exact marker list; that is outside this discriminator pass.

## Summary for the next round
Focus on targeted completion rather than structural rework: add mast-cell and lymphatic-endothelial coverage, make the leukocyte gate genuinely pan-leukocyte, and convert the mural-vs-endothelial contrast into a plain mural/pericyte-smooth-muscle axis. Consider splitting the aggregated lymphoid modules if the added dimensions still fit the desired resolution.
