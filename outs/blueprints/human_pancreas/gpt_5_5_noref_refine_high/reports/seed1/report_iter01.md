# Discriminator report

## Scenario echo
The scenario is healthy human pancreas with no specified age, sex, condition, or celltype restriction, so the embedding should broadly cover pancreatic epithelial compartments, endocrine islet cell types, acinar and ductal epithelium, stromal/mural cells, vasculature, immune cells, and neural-associated populations. There are no notes imposing a special topology or emphasis.

## Rating: minor
The blueprint is broadly biologically plausible and covers the main human pancreatic lineages, including acinar, ductal/centroacinar, endocrine subtypes, stromal/mural, endothelial, immune, mast, adipocyte, mesothelial, Schwann/neural, and cell-cycle states. It needs targeted refinement because several gates/marker modules would bias or blur important populations, especially non-myeloid immune cells and rare neural/mesothelial compartments.

## Findings

### Minor issues
- **module:leukocyte_core and axes gated by leukocyte_core** — The leukocyte gate is myeloid/APC-heavy: LST1, TYROBP, SPI1, FCER1G, and HLA-DRA are not pan-leukocyte markers and will preferentially score monocytes/macrophages/DCs. Because all T, NK, B, plasma, and mast axes are gated by this module, those populations could be artificially suppressed or distorted.
  *Fix:* Replace the leukocyte_core gene set with broader pan-immune markers such as PTPRC, LAPTM5, LCP1, CORO1A, CD53, DOCK2, ARHGAP15, and FYB1, leaving lineage-specific myeloid genes in macrophage_monocyte or dendritic modules.
- **axis:mesothelial_surface_epithelium** — Mesothelial cells are plausible for pancreatic surface/capsular tissue, and the module genes are mostly canonical, but gating them through epithelial_core may under-detect them because epithelial_core contains pancreatic epithelial/EPCAM/TACSTD2/KRT19/CDH1 features rather than a clean mesothelial gate. KRT8/KRT18 overlap helps, but the gate is still too pancreatic-epithelial-biased.
  *Fix:* Make mesothelial_surface_epithelium a plain mesothelial axis, or gate it with a broader epithelial/serosal module that does not depend on EPCAM/TACSTD2/KRT19/CDH1. Consider adding CALB2, KRT7, or CLDN15 if the module needs more mesothelial support.
- **module:neuron_neuroendocrine_fiber and axis:neuronal_fiber_within_neural_compartment** — The neural refinement mixes neuronal genes with generic neuroendocrine vesicle genes (SNAP25, SYP) that are also strong in endocrine cells, and then gates the signal by Schwann/glial markers. This may fail to represent actual enteric/intrapancreatic neurons while also relying on genes that overlap endocrine identity.
  *Fix:* Separate the neuronal program more cleanly from endocrine vesicle biology. Use a neuronal module enriched for TUBB3, ELAVL4, RBFOX3, PHOX2B, RET, ASCL1/CHAT/TH as appropriate, and either make it a plain rare-neural axis or gate it with a broad neural/glial compartment rather than only Schwann markers.
- **general endocrine/developmental modules** — The blueprint includes pancreatic_epithelial_progenitor, endocrine_progenitor_neurogenin, epsilon_ghrelin, and enterochromaffin-like endocrine axes. These are biologically plausible for an age-unspecified human pancreas scenario, but if the intended data are adult healthy pancreas they may add rare/fetal-skewed dimensions at the expense of adult endocrine state resolution.
  *Fix:* Keep them if age truly remains unspecified. If the Generator wants a more adult-leaning healthy pancreas embedding, de-emphasize progenitor/fetal endocrine axes and add adult endocrine state axes such as beta stress/dedifferentiation or alpha/beta maturity only if grounded by canonical markers.

## Out of scope (FYI, not actionable for the Generator)
- The blueprint has no placeholder or orphan-module structural problems; `python3 embedder_blueprint.py --describe` exits cleanly.
- No deletions are required for the current issues; all can be addressed by marker swaps or axis/gate modifications.

## Summary for the next round
Focus on replacing the myeloid-biased leukocyte_core with a true pan-leukocyte gate, then adjust the mesothelial and neural axes so they are not suppressed by inappropriate pancreatic epithelial or Schwann-only gates. The overall lineage coverage and dimensionality are otherwise reasonable for a healthy human pancreas embedding.
