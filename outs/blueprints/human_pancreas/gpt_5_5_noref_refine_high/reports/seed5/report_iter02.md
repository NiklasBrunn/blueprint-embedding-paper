# Discriminator report

## Scenario echo
Human pancreas, with no age, sex, or disease restriction and no additional notes. Because `celltype` is null, the embedding should cover the expected healthy human pancreatic landscape: endocrine islet cells, exocrine acinar/ductal/centroacinar epithelium, vascular cells, stromal/stellate cells, resident immune populations, mesothelium/adipose where captured, and neural/glial elements.

## Rating: minor
The blueprint has broad, biologically plausible pancreas coverage and no structural conformance problems. It needs targeted refinements to avoid a few misleading contrasts and improve some marker/state definitions, but the overall structure is usable.

## Findings

### Minor issues
- **axis:neuron_within_neural / modules:pan_neural, neural_glial, neuron** — The neural branch is biologically awkward. `pan_neural` mixes Schwann/glial genes with neuronal genes, but the gated neuron axis subtracts `endocrine_core` rather than a neural sibling module. Endocrine cells and neurons both have neuroendocrine/secretory features, but they are not sibling neural states; the current contrast may partly encode "neuron versus endocrine contamination" instead of neuron versus Schwann/glial identity.
  *Fix:* Make `neuron_within_neural` contrast `neuron` against `neural_glial` under the `pan_neural` gate, and consider adding a reciprocal gated Schwann/glial axis under `pan_neural` if dimensionality allows.

- **module:quiescent_pancreatic_stellate / axis:quiescent_stellate_within_stromal** — The quiescent stellate module contains several genes that also strongly mark pericytes or adipocytes (`RGS5`, `PDGFRB`, `FABP4`), while the axis simultaneously uses pericyte and adipocyte modules as negatives. This is not fatal, but it can blur quiescent pancreatic stellate cells with mural/adipose programs.
  *Fix:* Keep the stellate axis but rebalance the module toward quiescent PSC/fibroblast markers such as `SPARCL1`, `COL15A1`, `COL6A1`, `COL6A2`, `DCN`, `LUM`, `GPX3`, and `IGFBP4`, reducing dependence on pericyte/adipocyte-heavy markers.

- **module:ductal_secretory_mucin / axis:ductal_secretory_mucin_within_ductal** — A MUC5B/TFF/LCN2/MMP7 ductal program is plausible in human pancreas atlases, but it leans toward secretory/metaplastic or stressed ductal states for an otherwise healthy, condition-null scenario. As written it may overemphasize an injury-like ductal state relative to core ductal regional identities.
  *Fix:* Retain the ductal subtype axis, but tune the module toward canonical MUC5B-positive/secretory ductal identity and avoid making inflammatory injury genes (`LCN2`, `MMP7`) dominate. If the intended target is a stress/metaplastic state, rename the module/axis accordingly.

- **general immune coverage** — The immune compartment covers myeloid, DC subsets, T/NK, B/plasma, and mast cells, but lacks a simple pan-leukocyte anchor. This is not a major coverage gap because lineage-specific immune axes are present, but a `PTPRC`-centered immune axis would make immune versus non-immune separation cleaner.
  *Fix:* Add a compact `immune_common` module/axis with markers such as `PTPRC`, `LCP1`, `CORO1A`, `CD52`, and `ARHGAP15`, or otherwise ensure existing lineage axes are sufficient for datasets with mixed low-depth immune cells.

## Out of scope (FYI, not actionable for the Generator)
- The blueprint is relatively high dimensional for a hand-authored marker embedding, but the dimensionality is defensible given whole-pancreas coverage.

## Summary for the next round
Focus on cleaning the neural gating contrast, improving the quiescent stellate marker set so it does not fight its own negatives, and slightly de-emphasizing injury-heavy ductal secretory markers. A pan-immune anchor would be a useful small addition, but the current blueprint does not need wholesale restructuring.
