# Discriminator report

## Scenario echo
The scenario is human pancreas with all cell types in scope, no specified age, sex, disease, treatment, or design notes. The embedding should therefore cover the baseline pancreatic epithelial/endocrine/exocrine compartments plus stromal, vascular, neural, and resident immune populations without imposing a disease- or developmental-specific topology.

## Rating: minor
The blueprint is broadly biologically plausible and covers the main human pancreas compartments: endocrine subtypes, acinar, ductal, stellate/fibroblast, mural, blood and lymphatic endothelium, immune subsets, and Schwann/neural cells. It is not quite publishable unchanged because a few biologically distinct resident populations are collapsed into single axes, and one broad immune module includes a ubiquitous marker that weakens specificity.

## Findings

### Minor issues
- **mural_cell / mural_identity** — The mural module merges pericyte markers (RGS5, PDGFRB, CSPG4, MCAM, NOTCH3) with vascular smooth muscle markers (ACTA2, TAGLN, MYH11). Human pancreas contains both pericytes and smooth muscle-like vascular mural cells, and a single identity axis cannot separate these expected stromal populations.
  *Fix:* Split the mural representation into pericyte and smooth-muscle modules/axes, or keep a mural core and add gated subtype axes for pericyte versus vascular smooth muscle. Plausible marker sets would separate RGS5/PDGFRB/CSPG4/NOTCH3/ABCC9/KCNJ8 from ACTA2/TAGLN/MYH11/CNN1/MYL9/TPM2.
- **dendritic_cell / dendritic_cell_within_immune** — The dendritic module combines cDC2 markers (FCER1A, CLEC10A, CD1C), cDC1 markers (CLEC9A, XCR1), and pDC-associated LILRA4 into one axis. That is acceptable for coarse leukocyte annotation but under-resolves expected human pancreas antigen-presenting subsets in an all-cell tissue embedding.
  *Fix:* Either rename/reshape this as a broad dendritic/APC axis and add cDC1/cDC2/pDC gated immune subaxes, or split the current module into separate dendritic subtype modules with sibling negative modules under immune_core.
- **immune_core** — B2M is included in the immune core, but it is broadly expressed across nucleated human cells and will dilute leukocyte specificity. CXCR4 is also broad across immune states and some nonimmune contexts, though less problematic than B2M.
  *Fix:* Remove B2M from the immune_core gene list and replace it with a more leukocyte-enriched pan-immune marker such as CD53, LAPTM5, or another cited immune atlas-supported marker. Keep PTPRC/LCP1/CORO1A-type genes as the core identity signal.
- **epsilon_cell / epsilon_within_endocrine** — The epsilon module is driven by one highly specific marker, GHRL, plus broad endocrine/prohormone genes shared with other endocrine modules. Because ghrelin-positive epsilon cells are rare in human pancreas and age is unspecified, this axis is plausible but could score noisy endocrine cells rather than a clean epsilon population.
  *Fix:* If the literature resources support additional human pancreatic ghrelin-cell markers, add them and reduce reliance on CHGA/CHGB/PCSK1/PCSK2. If support is weak, keep the axis but make the module more explicitly GHRL-centered and avoid presenting it as equally resolved as beta/alpha/delta/PP cells.

## Out of scope (FYI, not actionable for the Generator)
- Citation DOI verification is out of scope for this discriminator report; the structural `--describe` run succeeded.

## Summary for the next round
Focus on targeted refinement rather than structural rebuild: split mural cells into pericyte and smooth-muscle resolution, split or clarify the mixed dendritic-cell axis, and tighten pan-immune and epsilon marker specificity. The rest of the pancreas coverage is broadly appropriate.
