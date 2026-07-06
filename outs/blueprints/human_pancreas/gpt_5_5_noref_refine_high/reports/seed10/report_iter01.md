# Discriminator report

## Scenario echo
The scenario is healthy/baseline human pancreas with all pancreatic cell types in scope; age and sex are unspecified, and there are no notes imposing a special developmental, disease, or pathway-focused topology.

## Rating: minor
The blueprint covers the main human pancreas landscape well: endocrine islet subtypes, acinar, ductal/centroacinar, endothelial/lymphatic, stromal/mural, mesothelial, neural/Schwann, and major immune compartments. It still needs targeted edits for a few missing expected populations/substates, a redundant stromal axis, and several marker choices that are not canonical enough for the claimed modules.

## Findings

### Minor issues
- **Coverage: adipocytes and antigen-presenting dendritic cells** — Broad human pancreas tissue commonly includes adipocyte/fat-associated signal and tissue immune APC diversity beyond C1Q macrophages, especially conventional dendritic cells. The current blueprint has no adipocyte axis and collapses all myeloid/APC signal into `myeloid_macrophage`.
  *Fix:* Add an `adipocyte` module/axis with markers such as `ADIPOQ`, `PLIN1`, `FABP4`, `LPL`, `LEP`, `CIDEC`. Add a dendritic-cell/APC refinement under `immune_core` or split myeloid into macrophage and DC modules using markers such as `FCER1A`, `CD1C`, `CLEC10A`, `CLEC9A`, `XCR1`, `CCR7`, `LAMP3` as appropriate.

- **Axis design: `stromal_fibroblast_vs_mural` and `mural_contractile_vs_matrix_stroma`** — These two contrast axes are exact inverses of one another and therefore redundant dimensions. A single fibroblast-vs-mural contrast, plus a positive mural identity/subtype axis if needed, is biologically cleaner.
  *Fix:* Remove `mural_contractile_vs_matrix_stroma`. If additional mural resolution is desired, replace the lost dimension with a non-redundant pericyte/smooth-muscle refinement rather than another inverse contrast.

- **Coverage/axis design: mural compartment resolution** — `vascular_mural` combines pericytes (`RGS5`, `PDGFRB`, `CSPG4`, `NOTCH3`) and vascular smooth muscle/contractile cells (`ACTA2`, `TAGLN`, `MYH11`). Those are related but separable vascular support populations in pancreas.
  *Fix:* Consider splitting into `pericyte` and `smooth_muscle_contractile` modules and adding either a gated/contrast refinement or at least a clearer mural identity axis after removing the inverse duplicate.

- **Marker plausibility: `endothelial_lymphatic`** — `TFF3` is not a canonical lymphatic endothelial marker and reads more epithelial/ductal or mucosal than lymphatic. The rest of the module is plausible.
  *Fix:* Replace `TFF3` with a stronger lymphatic marker such as `CCL21`, `RELN`, `TBX1`, `FOXC2`, or another accepted lymphatic endothelial gene not already overused.

- **Marker plausibility: `immune_core`** — `HLA-DRA` is useful for antigen-presenting cells but is not a pan-leukocyte core marker and can bias the general immune axis toward macrophage/DC/B-cell populations while underrepresenting T/NK cells.
  *Fix:* Keep `PTPRC`, `LCP1`, and `CORO1A`, but replace `HLA-DRA` in the core with a broader leukocyte marker such as `RAC2`, `ARHGAP15`, `FYB1`, or `ITGAL`; reserve `HLA-DRA` for APC/DC/myeloid modules.

- **Gene-set hygiene: `RBP4` overlap in endocrine subtype modules** — `RBP4` appears in both `beta_cell` and `delta_cell`, which weakens subtype discrimination. `RBP4` is not one of the strongest canonical beta-cell markers compared with the other beta genes already present.
  *Fix:* Remove or replace `RBP4` in `beta_cell` with a more beta-selective marker such as `UCN3`, `DLK1`, `NPTX2`, `HADH`, or `ABCC8`; keep the delta module focused on `SST`, `HHEX`, `RGS4`, `LEPR`, `BCHE`, etc.

- **Marker plausibility: `epsilon_cell`** — `GHSR` is a ghrelin receptor rather than the usual defining secretory product marker, and `MLN` is not a core pancreatic epsilon-cell marker. Because epsilon cells are rare, weak markers can make this axis noisy.
  *Fix:* Keep `GHRL` and endocrine-context genes, but replace `GHSR`/`MLN` with stronger epsilon/ghrelin-cell-associated markers where available, or reduce the module to a conservative rare-cell signature rather than padding it.

## Out of scope (FYI, not actionable for the Generator)
- No citations are required by the blueprint contract, so literature-source validation was not assessed.

## Summary for the next round
Focus on targeted cleanup rather than restructuring: add missing adipocyte and DC/APC coverage, remove the redundant inverse stromal contrast, optionally split mural cells into pericyte versus smooth muscle, and tighten the few questionable markers in lymphatic endothelial, immune core, beta/delta, and epsilon modules.
