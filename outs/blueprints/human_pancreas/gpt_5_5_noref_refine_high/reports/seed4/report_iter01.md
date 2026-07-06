# Discriminator report

## Scenario echo
Human pancreas, with no cell-type, age, sex, condition, or notes restrictions. The embedding should therefore cover the main healthy pancreatic epithelial, endocrine, stromal, vascular, neural, and immune compartments without emphasizing disease- or development-specific programs unless they are plausible baseline/rare states.

## Rating: minor
The blueprint has a broadly sensible pancreatic scaffold and captures the major exocrine, endocrine, stromal, vascular, immune, and neural compartments. It still needs targeted biological cleanup: expected vascular and immune substructure is under-resolved, one developmental progenitor axis is out of scope for an unqualified healthy pancreas scenario, and several modules contain markers that would confound neighboring lineages.

## Findings

### Minor issues
- **vascular axes** — Blood endothelium is represented as one coarse compartment plus lymphatic, but a domain expert would expect at least arterial, venous, and capillary endothelial distinctions in a whole-pancreas embedding.
  *Fix:* Add gated endothelial subtype axes/modules for arterial endothelium (e.g. GJA5, EFNB2, SOX17, HEY1/HEY2), venous endothelium (e.g. ACKR1, NR2F2, VWF, VCAM1 where appropriate), and capillary endothelium (e.g. RGCC, PLVAP, EMCN, CA4/RAMP2), while retaining lymphatic as a separate endothelial branch.
- **immune compartment** — The immune coverage misses dendritic-cell/APC subtypes, and `immune_core` is biased toward myeloid/APC genes (`LST1`, HLA genes, `CD74`), which will weaken T/NK gating. The `t_cell` module is also skewed toward naive/IL7R T cells and may miss cytotoxic T cells.
  *Fix:* Broaden `immune_core` to pan-leukocyte markers such as `PTPRC`, `LAPTM5`, `LCP1`, `CORO1A`, `CD53`, `ARHGAP15`, `RAC2`; add a dendritic/APC axis using markers such as `CD1C`, `FCER1A`, `CLEC10A`, `ITGAX`, `LILRA4`, `IRF7` as appropriate; and expand T-cell markers with `CD2`, `CD247`, `CD8A`, `CD8B` or split helper/cytotoxic T if adding one extra axis is acceptable.
- **epithelial_pancreatic_progenitor_state / pancreatic_progenitor** — A fetal/developmental pancreatic progenitor program (`NEUROG3`, `PTF1A`, `NKX2-2` mixed with ductal and beta genes) is not expected as a generic healthy human pancreas axis when age and developmental notes are unspecified. It also risks conflating normal ductal SOX9/HNF1B signal with true developmental progenitors.
  *Fix:* Remove the developmental progenitor axis/module, or replace it with an in-scope mature ductal/centroacinar refinement if additional epithelial resolution is desired.
- **epithelial_endocrine_vs_exocrine** — The axis name says endocrine versus exocrine, but the negative side is only `acinar_cell`; ductal epithelium is also exocrine and is not represented in the contrast.
  *Fix:* Either rename/reframe this as endocrine versus acinar, or add/use an exocrine module that includes both acinar and ductal programs as the negative side.
- **marker hygiene across several modules** — `lymphatic_endothelial` includes `TFF3`, which is epithelial/goblet-like rather than a canonical lymphatic marker; `ductal_bicarbonate_mucin` includes less canonical or confounding choices (`SLC26A3`, `MUC5B`, `TFF1`, `AQP1`) for normal pancreatic ducts; `mesothelial` includes fibroblast-associated `LUM`; and `autonomic_neuronal` includes developmental `ASCL1`.
  *Fix:* Replace `TFF3` with lymphatic markers such as `NRP2`, `RELN`, `TBX1`, or `CCL21`-adjacent lymphatic genes; sharpen ductal transport markers toward `CFTR`, `SLC4A4`, `SLC26A6`, `CA2`, `KRT7`, `KRT19`, `SOX9`; swap mesothelial `LUM` for `CALB2`, `ITLN1`, `C3`, or another mesothelial marker; and replace neuronal `ASCL1` with adult autonomic markers such as `TH`, `DBH`, `SLC18A3`, `VIP`, or `NOS1`.
- **acinar mature/adaptive modules** — `acinar_zymogen_mature` is largely a duplicate of `acinar_cell`, while `acinar_cell` itself includes `REG1A`, which is also used in the adaptive/stress module. This blurs identity versus state.
  *Fix:* Keep the acinar identity module focused on constitutive acinar enzymes/transcriptional identity, and make the mature/zymogen state more distinct with additional secretory/acinar maturation markers such as `RBPJL`, `SYCN`, `CEL`, `GP2`, `CELA2A`, `CPB1`; reserve `REG1A/REG1B/REG3A/SPINK1` for the adaptive/stress axis.

## Out of scope (FYI, not actionable for the Generator)
- No citation/resource audit was performed; this report only judges biological coverage and marker plausibility.

## Summary for the next round
Focus on targeted refinements: add endothelial arterial/venous/capillary resolution, repair immune gating and add dendritic/APC coverage, and remove or repurpose the developmental pancreatic progenitor axis. Then clean the most confounding markers in lymphatic, ductal, mesothelial, autonomic neuronal, and acinar state modules.
