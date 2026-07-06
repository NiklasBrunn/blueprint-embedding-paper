# Discriminator report

## Scenario echo
The scenario is healthy/baseline human pancreas with no specified age, sex, or focused notes, and `celltype` is null, so the embedding should cover the full pancreatic tissue landscape: endocrine islet cells, exocrine acinar and ductal/centroacinar epithelium, vascular and lymphatic endothelium, stromal/stellate/mural/neural support cells, and resident immune populations.

## Rating: minor
The blueprint is broadly biologically plausible and covers the major human pancreas compartments, including endocrine subtypes, acinar/ductal epithelium, vascular, stromal, neural, and immune lineages. It still needs targeted refinement for ductal/centroacinar coverage and a few marker/module-design issues, but the overall structure does not need replacement.

## Findings

### Minor issues
- **ductal/centroacinar coverage** — The blueprint has a ductal core and a mucin-like ductal substate, but it does not explicitly resolve centroacinar/terminal duct cells, which are expected in whole human pancreas and are not equivalent to a generic ductal or mucin program.
  *Fix:* Add or recast a ductal refinement axis for centroacinar/terminal duct identity using supported markers such as `SPP1`, `CLDN10`, `KRT19`, `SOX9`, `CFTR`, `ANXA4`, and related atlas-supported ductal subtype genes. Keep the existing broad ductal axis, but make the ductal substructure biologically interpretable rather than only mucin-focused.

- **endothelial_lymphatic module** — `TFF3` is not a canonical lymphatic endothelial marker and reads more like epithelial/goblet-like or ductal mucin biology. Its presence in the lymphatic module will bleed epithelial signal into the lymphatic-vs-vascular axis.
  *Fix:* Remove `TFF3` from `endothelial_lymphatic` and replace it with a better-supported lymphatic endothelial marker if available from the cited atlases, for example `ACKR4`, `RELN`, `TBX1`, `CCL21`, or another human pancreas lymphatic marker supported by the resources.

- **stellate_quiescent / stromal_fibroblast modules and axes** — `stellate_quiescent` substantially duplicates `stromal_fibroblast` (`PDGFRA`, `DCN`, `LUM`, collagens, `FAP`, `CXCL14`) and includes `FAP`, which is more consistent with activated fibroblast biology than quiescent pancreatic stellate identity. As written, `stromal_fibroblast_compartment` and `stellate_matrix_compartment` are likely redundant rather than resolving fibroblast versus quiescent stellate biology.
  *Fix:* Rework `stellate_quiescent` toward a supported pancreatic stellate/quiescent or perivascular-stellate program, dropping activation-biased markers such as `FAP` and using better-supported quiescent/stellate markers where available, such as retinoid/lipid-storage or stellate-associated genes (`RBP1`, `LRAT`, `DES`, `GFAP`, `COLEC11`, `APOD`, `RGS5`/`PDGFRB` if justified by the atlas context). Then make `stellate_activation_state` contrast or gate activated myofibroblast markers against the quiescent/stellate module rather than leaving two near-duplicate plain stromal axes.

## Out of scope (FYI, not actionable for the Generator)
- Citation verification was not repeated, per task scope.
- The immune compartment could be further subdivided into monocyte, cDC1/cDC2, pDC, CD4/CD8/Treg states, but the current broad immune axes are acceptable for a pancreas-wide embedding unless a future scenario asks for immune resolution.

## Summary for the next round
Keep the current broad pancreas structure. Focus on adding an explicit centroacinar/terminal duct refinement, cleaning the lymphatic endothelial marker set by replacing `TFF3`, and separating true quiescent/stellate biology from generic fibroblast and activated matrix programs.
