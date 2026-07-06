# Discriminator report

## Scenario echo
The scenario is human pancreas, with no specific age, sex, condition, or cell type restriction and no extra notes. The embedding should therefore cover the normal pancreatic epithelial, endocrine, stromal/vascular, immune, neural, adipose, and generic proliferative/stress-state landscape without adding disease- or perturbation-specific programs.

## Rating: minor
The blueprint is broadly plausible and structurally complete for a healthy human pancreas: it covers endocrine islet cells, acinar/ductal epithelium, vascular and mural cells, stellate/fibroblast cells, adipocytes, immune lineages, Schwann/neural cells, and cycling state. Targeted edits are still warranted because several immune modules blur distinct lineages, and one rare endocrine population expected in human pancreas atlases is missing.

## Findings

### Minor issues
- **Endocrine coverage** — The endocrine branch includes beta, alpha, delta, PP/gamma, and epsilon cells, but it lacks the rare enterochromaffin/serotonergic endocrine population reported in human pancreatic/islet datasets. For a broad human pancreas embedding, this is a small but real coverage gap alongside the other rare endocrine axes.
  *Fix:* Add an enterochromaffin/serotonergic endocrine module with canonical markers such as `TPH1`, `DDC`, `SLC18A1`, `TAC1`, `LMX1A`, and endocrine markers as appropriate, then add a gated axis under `endocrine_common` contrasted against the sibling endocrine modules.
- **t_nk_cell and nk_cytotoxic_cell / axes `t_nk_within_leukocyte` and `nk_cytotoxic_within_leukocyte`** — The T/NK module contains strong NK/cytotoxic markers (`NKG7`, `KLRD1`, `GNLY`) while a separate NK/cytotoxic module uses nearly the same program. The gated T and NK axes also do not contrast against each other, so cytotoxic NK cells and cytotoxic T cells may project similarly.
  *Fix:* Make the T-cell module T-lineage specific (`CD3D`, `CD3E`, `TRAC`, `CD2`, `CD247`, plus markers such as `IL7R`, `CCR7`, `CD4`/`CD8A` if supported), keep cytotoxic genes in the NK/cytotoxic module, and include T and NK modules as mutual negatives in their gated leukocyte axes.
- **b_plasma_cell and plasma_cell / axes `b_plasma_within_leukocyte` and `plasma_within_leukocyte`** — `b_plasma_cell` mixes B-cell markers (`MS4A1`, `CD79A`, `CD79B`, `BANK1`, `CD37`) with plasma-cell markers (`MZB1`, `JCHAIN`, `XBP1`), while `plasma_cell` separately captures plasma cells. This makes the B-cell axis partially redundant with the plasma axis and weakens discrimination of naive/memory B cells from antibody-secreting cells.
  *Fix:* Convert `b_plasma_cell` into a cleaner B-cell module by removing plasma genes and adding B-cell markers such as `CD19`, `CD22`, `CD74`, `HLA-DRA`, `TCL1A`, or `MS4A1`-centered markers as supported; keep `plasma_cell` as the plasma-cell module and make the two axes mutual negatives.
- **dendritic_cell** — The dendritic module combines cDC2 (`CD1C`, `FCER1A`, `CLEC10A`), cDC1 (`CLEC9A`), and pDC (`LILRA4`) markers into one axis. A broad DC axis is acceptable, but in a 24-dimensional pancreas embedding this loses biologically meaningful immune resolution while other rarer immune populations get separate axes.
  *Fix:* Either narrow this module to a broad antigen-presenting DC program without subtype-exclusive markers, or split it into cDC and pDC axes if the supporting atlas resources justify the added dimensionality.
- **endothelial / lymphatic_vs_blood_endothelial** — The `endothelial` module is pan-endothelial/blood-endothelial leaning, but the lymphatic contrast uses it directly as the negative even though lymphatic endothelial cells also express pan-endothelial genes such as `PECAM1` and `CDH5`. This can make the contrast less clean than a vascular-common gate with blood-vs-lymphatic sibling programs.
  *Fix:* Refactor to a vascular endothelial common gate plus separate blood endothelial (`VWF`, `EMCN`, `RAMP2`, `KDR`, `ESAM`, etc.) and lymphatic (`PROX1`, `LYVE1`, `FLT4`, `PDPN`, `CCL21`) modules, or at least make the negative module explicitly blood-endothelial rather than pan-endothelial.

## Out of scope (FYI, not actionable for the Generator)
- Citation and DOI validation were not rechecked; this review only judged biological coverage, marker plausibility, and structural fit.

## Summary for the next round
Keep the overall pancreatic compartment structure. The next Generator should focus on cleaning lineage specificity in the leukocyte branch, adding a rare enterochromaffin/serotonergic endocrine axis if the cited human pancreas/islet literature supports it, and making the endothelial lymphatic contrast use a blood-endothelial rather than pan-endothelial negative program.
