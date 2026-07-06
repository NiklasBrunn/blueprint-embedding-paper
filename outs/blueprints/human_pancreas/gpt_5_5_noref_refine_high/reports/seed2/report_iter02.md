# Discriminator report

## Scenario echo
The scenario is human pancreas, with all pancreatic cell types in scope because `celltype` is null; age and sex are unspecified, and condition is null/healthy baseline. There are no extra notes requesting a specific developmental, disease, or pathway-focused structure.

## Rating: minor
The blueprint has broad, plausible coverage of major pancreatic epithelial, endocrine, exocrine, stromal, vascular, immune, neural, and cycling compartments. It needs targeted biological cleanup, mainly because a few state axes are mislabeled or skew toward stress/metaplastic programs that are not a clean healthy-pancreas baseline representation.

## Findings

### Minor issues
- **axis:acinar_zymogen_load / module:acinar_zymogen_high** - The module is not actually a zymogen-load program. Most genes are unfolded-protein-response/ER-stress markers (`XBP1`, `HSPA5`, `HSP90B1`, `HERPUD1`, `DDIT3`, `ATF4`, `DNAJB9`, `PDIA4`) rather than acinar digestive enzyme/zymogen markers. Because the general acinar identity module already contains the digestive enzymes, the current axis is biologically mislabeled and may encode stressed acinar cells instead.
  *Fix:* Either rename/reframe this as an acinar secretory-stress/UPR state, or replace the gene list and axis name with a nonredundant healthy acinar subtype/state if one is intended. If keeping a zymogen/load state, use acinar digestive enzyme/secretory genes and avoid duplicating the core acinar module too heavily.
- **axis:ductal_transport_mucin_state / module:ductal_transport_mucin** - This mixes healthy ductal ion/fluid transport markers (`CFTR`, `SLC4A4`, `AQP1`) with mucin/gastric metaplasia-associated markers (`TFF1`, `TFF2`, partly `MUC13`). In a healthy baseline pancreas scenario, the TFF-heavy mucinous component is not a clean expected state and could pull the embedding toward pathology/metaplasia.
  *Fix:* Retune this module toward healthy ductal/centroacinar transport and secretion, e.g. keep `CFTR`, `SLC4A4`, `AQP1` and use additional ductal markers such as `KRT19`, `KRT7`, `ANXA4`, `SPP1`, `LCN2` only if justified by baseline datasets; avoid TFF/metaplastic markers unless the axis is explicitly renamed as a metaplastic/stress state.
- **axis:pdc_within_myeloid / module:plasmacytoid_dendritic_cell** - The pDC subtype is biologically in scope, but gating it with the broad `myeloid_core` module (`LYZ`, `LST1`, `AIF1`, `FCER1G`, `TYROBP`, `CSF1R`, `HLA-DRA`) can suppress canonical pDCs, which are better captured by `LILRA4`, `CLEC4C`, `IL3RA`, `TCF4`, `IRF7`, `PLD4` and may not score strongly on monocyte/macrophage-associated myeloid genes.
  *Fix:* Use a dendritic/APC gate that includes HLA/FLT3/IRF8/CST3 or add a separate pDC identity axis, and keep the pDC contrast against cDC/monocyte/macrophage modules scoped to an appropriate antigen-presenting immune compartment.
- **axis:pancreatic_progenitor_to_mature_endocrine / module:pancreatic_progenitor_bias** - The axis is ungated and uses genes that overlap adult ductal/centroacinar and beta-cell programs (`SOX9`, `HNF1B`, `ONECUT1`, `PDX1`, `NKX6-1`, `HES1`, `PROM1`). With no developmental note and healthy baseline condition, this can behave more like a ductal/progenitor-like or immature epithelial signal than a true pancreas progenitor-to-mature-endocrine axis.
  *Fix:* If the intent is broad age coverage, gate the axis to pancreatic epithelium and make the contrast explicitly immature/progenitor-like versus mature endocrine/acinar/ductal programs. If the intent is adult healthy pancreas, retune it into a ductal/centroacinar progenitor-like state or drop the developmental contrast in favor of mature compartment resolution.

## Out of scope (FYI, not actionable for the Generator)
- Adipocytes, neutrophils, and endocrine rare subtypes are not heavily represented, but their omission is acceptable for a baseline pancreas blueprint unless the downstream dataset is known to include substantial peripancreatic fat or blood contamination.
- Some modules include overlapping shared markers by design; this is not a structural issue as long as the targeted cleanup above reduces misleading state labels.

## Summary for the next round
Keep the overall 43-axis structure and broad compartment coverage, but clean up the biology of the state axes: reframe or retune acinar UPR versus zymogen load, make the ductal state healthy-baseline rather than TFF/mucinous-metaplastic, and adjust pDC gating so canonical pDCs are not filtered by monocyte/macrophage-biased myeloid markers.
