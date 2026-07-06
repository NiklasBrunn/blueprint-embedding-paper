# Discriminator report

## Scenario echo
The scenario is human pancreas, with no specified age, sex, disease, perturbation, or cell-type restriction, so the embedding should cover the healthy human pancreas cell-type landscape broadly: endocrine islet cells, acinar and ductal epithelium, stromal/vascular/lymphatic compartments, resident immune cells, neural/Schwann cells, mesothelium/capsule, and common states such as cycling.

## Rating: minor
The blueprint is broadly plausible and covers most major pancreatic compartments with reasonable human marker sets. It needs targeted refinement, mainly removing one redundant inverse exocrine contrast, tightening immune gates/labels, and adding a small number of expected tissue populations that are not explicitly represented.

## Findings

### Minor issues
- **axes: acinar_vs_ductal_exocrine and ductal_vs_acinar_exocrine** — These two axes are exact sign-reversed duplicates using the same two modules. They add dimensionality without adding biological resolution; keeping both pads the embedding and can over-weight the acinar/ductal split.
  *Fix:* Keep one acinar/ductal contrast axis and remove the reverse duplicate. If ductal positivity needs its own readable dimension, replace the duplicate with a plain or gated ductal/centroacinar program rather than the mathematical inverse.

- **immune_core / immune gated axes** — The immune gate is biased toward myeloid/APC genes (`LST1`, `TYROBP`, `FCER1G`, `AIF1`, `CD74`, `HLA-DRA`) rather than a pan-leukocyte program. Because all immune subtype axes are gated by this module, T, NK, and some B cells may be under-scored. The `t_nk_lymphocyte` module is actually T-cell markers, while NK cells are separate; `b_plasma_cell` is actually a B-cell module while plasma cells are separate.
  *Fix:* Make `immune_core` more pan-leukocyte, e.g. centered on `PTPRC`, `LCP1`, `CORO1A`, `CD53`, `LAPTM5`, `HCLS1`, `RAC2`, `SRGN` or similar human leukocyte markers. Rename or relabel the T and B modules/axes to match their biology, or adjust the marker sets if the broader labels are intended.

- **immune compartment coverage** — Dendritic cells/antigen-presenting dendritic subsets are expected among resident pancreatic immune populations, but the blueprint only resolves macrophage/monocyte, T, NK, B, plasma, and mast cells. The current macrophage/APC module will not cleanly identify cDC-like cells.
  *Fix:* Add an immune-gated dendritic-cell axis/module if supported by the same pancreas/atlas resources, using canonical human cDC/APC markers such as `CD1C`, `FCER1A`, `CLEC10A`, `ITGAX`, `IRF8`, `HLA-DRA`, with care not to duplicate macrophage markers.

- **whole-pancreas stromal coverage** — Healthy human pancreas tissue can include adipocytes or adipocyte-like fat cells, especially in adult whole-tissue samples, but no module captures an adipocyte program. This is less central than epithelial/endocrine coverage but is a gap for an all-cell-types pancreas embedding.
  *Fix:* Add a plain or stromal-scoped adipocyte axis if the cited pancreas/atlas literature supports it, with markers such as `ADIPOQ`, `PLIN1`, `LPL`, `FABP4`, `CIDEC`, and `LEP`/`LIPE` as appropriate.

## Out of scope (FYI, not actionable for the Generator)
- I did not verify DOIs or journal impact factors; citation verification is outside this discriminator pass.

## Summary for the next round
The next Generator should keep the current broad structure but remove the duplicate reverse acinar/ductal contrast, revise the immune gate to be truly pan-leukocyte with clearer T/B labels, and add dendritic-cell coverage. Add an adipocyte axis if the available human pancreas atlas/literature resources support whole-tissue adipocyte capture.
