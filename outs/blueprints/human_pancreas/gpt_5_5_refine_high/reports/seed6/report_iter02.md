# Discriminator report

## Scenario echo
The scenario is healthy/baseline human pancreas with no specified age, sex, condition, or notes, and `celltype=null`, so the embedding should cover the expected whole-pancreas landscape: endocrine islet cells, exocrine acinar and ductal epithelium, vascular and lymphatic endothelium, stromal/stellate/mural cells, resident immune cells, and neural/Schwann-associated cells.

## Rating: minor
The blueprint is broadly biologically plausible and covers the major pancreas compartments and canonical endocrine subtypes. Targeted refinement is still warranted because the immune gate is not a clean pan-leukocyte program, dendritic/APC immune coverage is thin, and one stromal state risks over-emphasizing disease-like myofibroblast activation in an otherwise healthy pancreas scenario.

## Findings

### Minor issues
- **module:immune_core / immune gated axes** — The immune gate mixes pan-leukocyte markers with myeloid/APC markers (`LST1`, `TYROBP`, `CD74`, `HLA-DRA`, `HLA-DPA1`). This will gate macrophage/B/APC-like cells well but can attenuate T and NK axes, because those lineages may only strongly contribute `PTPRC` and `CORO1A` from the current module.
  *Fix:* Rebuild `immune_core` as a cleaner pan-leukocyte module, for example emphasizing `PTPRC`, `CORO1A`, `LCP1`, `LAPTM5`, `CD52`, `CD53`, `ARHGDIB`, and similar broadly hematopoietic genes. Keep antigen-presentation genes in B-cell, macrophage/DC, or a separate APC module rather than in the universal gate.
- **immune coverage** — Resident dendritic/APC populations are not explicitly represented. A broad pancreas embedding can reasonably distinguish macrophage/monocyte from conventional dendritic cells or APC-like dendritic cells, especially because the current `macrophage_monocyte` module does not capture markers such as `FCER1A`, `CD1C`, `CLEC10A`, `CLEC9A`, `XCR1`, `LILRA4`, or `IRF8`.
  *Fix:* Add a dendritic/APC module and a gated immune axis if supported by the cited pancreas/atlas resources. A compact cDC-oriented axis is enough; plasmacytoid DCs are optional unless the source atlases clearly support them in pancreas.
- **axis:mesenchymal_stellate_myofibroblast_activation / module:stellate_myofibroblast_activation** — A stellate activation axis is biologically interpretable, but in a baseline healthy pancreas scenario it should not read as a pancreatitis/cancer fibrosis program. The current module contains several contractile mural-cell genes (`ACTA2`, `TAGLN`, `MYL9`, `TPM2`) that overlap with `pericyte_smooth_muscle`, so this axis may partly re-score smooth muscle/pericytes rather than a pancreatic stellate/fibroblast activation state.
  *Fix:* Either reframe the axis as a stromal matrix-remodeling/activated PSC state with less mural-cell overlap, or balance it with clearer quiescent stellate/fibroblast markers. Keep `POSTN`, `FN1`, `SPARC`, and collagens if justified, but reduce reliance on contractile markers that are also mural identities.
- **module:epsilon_cell / axis:endocrine_epsilon_identity** — Epsilon/ghrelin cells are rare in human pancreas and are more age-dependent than alpha, beta, delta, and PP cells. Keeping an epsilon axis is acceptable because age is unspecified, but the module should be conservative: `GHRL` is canonical, whereas `MLN` and `DEFB1` are less convincing as adult pancreatic epsilon identity markers and may pull toward enteroendocrine or ductal programs.
  *Fix:* Keep the axis only if the cited pancreas literature supports enough epsilon markers; otherwise tighten the module around canonical ghrelin-cell/endocrine markers and remove less specific genes from the module.

## Out of scope (FYI, not actionable for the Generator)
- I did not re-check DOI validity or whether each cited paper contains each exact gene list; that is outside this Discriminator pass.

## Summary for the next round
Focus on tightening immune modeling: make `immune_core` a true pan-leukocyte gate and add a small dendritic/APC axis if the pancreas/atlas literature supports it. Then reduce stromal ambiguity by making the stellate activation axis less mural-cell-like, and make the rare epsilon module more conservative.
