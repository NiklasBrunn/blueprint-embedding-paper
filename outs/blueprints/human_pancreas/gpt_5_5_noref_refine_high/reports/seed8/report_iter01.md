# Discriminator report

## Scenario echo
The scenario is human pancreas, with no specified age, sex, disease, treatment, or requested focus; `celltype` is null, so the embedding should cover the expected whole-pancreas cell-type/state landscape rather than a single lineage.

## Rating: minor
The blueprint covers the dominant human pancreatic compartments well: acinar, ductal/centroacinar, endocrine alpha/beta/delta/PP/epsilon, endothelial/lymphatic, stromal/mural, immune, Schwann/neural, and proliferation states. It still needs targeted refinement for a few missing whole-tissue populations and several marker/axis choices that would be questioned by a pancreas domain expert.

## Findings

### Minor issues
- **general coverage** — Whole-pancreas coverage is missing mesothelial/serosal capsule cells and adipocyte/lipid-rich stromal cells. These are plausible constituents of human pancreas tissue, especially when age is unspecified and `celltype` is null. The current stromal axis does not cleanly capture either population.
  *Fix:* Add modules and axes for mesothelial cells, e.g. `MSLN`, `WT1`, `UPK3B`, `CALB2`, `KRT8`, `KRT18`, `PDPN`, and for adipocytes/lipid stromal cells, e.g. `ADIPOQ`, `PLIN1`, `FABP4`, `LPL`, `CIDEA`, `LEP`, `CFD`, `PPARG`.

- **axis:endocrine_enterochromaffin_like_cell / module:enterochromaffin_serotonin** — A distinct enterochromaffin-like endocrine cell axis is not a standard expected population for healthy human pancreas. Serotonergic genes can mark rare endocrine states or off-scope gut-like cells, but this axis currently presents them as a core pancreatic endocrine cell type alongside alpha/beta/delta/PP/epsilon.
  *Fix:* Delete this module and axis unless the Generator can recast it as a clearly scoped rare serotonergic endocrine state. Prefer using the freed dimension for missing whole-tissue coverage or stronger immune/stromal resolution.

- **axis:lymphatic_endothelial_identity / module:lymphatic_endothelial** — The lymphatic axis is gated by `endothelial_blood_vessel`, whose marker set is blood vascular biased (`VWF`, `KDR`, `EMCN`, `PLVAP`). That can suppress bona fide lymphatic endothelial cells. The lymphatic module also includes `TFF3`, which is primarily epithelial/secretory rather than a lymphatic endothelial marker.
  *Fix:* Split the endothelial gate into a pan-endothelial module or broaden the current endothelial gate with pan-endothelial genes (`PECAM1`, `CDH5`, `CLDN5`, `ESAM`, `RAMP2`) and keep blood-vessel markers in a separate vascular module if needed. Replace `TFF3` with a lymphatic marker such as `TBX1`, `FOXC2`, `CCL21`, or another accepted human lymphatic endothelial gene.

- **module:ductal_bicarbonate_mucin** — `SLC26A3` is more characteristic of intestinal/colonic chloride-bicarbonate exchange than pancreatic duct identity; pancreatic duct bicarbonate transport is better represented by `CFTR`, `SLC4A4`, and `SLC26A6`. The module otherwise has reasonable ductal genes.
  *Fix:* Replace `SLC26A3` with `SLC26A6` or another pancreatic ductal marker such as `KRT7`, `KRT19`, `SOX9`, `HNF1B`, `MUC1`, `ANXA4`, or `MMP7` if not already represented.

- **immune subtype axes** — `t_cell_nk_cell` and `b_cell_plasma_cell` each merge biologically distinct immune populations into one gated axis. For a whole-pancreas embedding this is acceptable as a coarse first pass, but it leaves expected immune structure underresolved.
  *Fix:* If adding only a few dimensions, split T vs NK with `CD3D/CD3E/TRAC/IL7R/CCR7` versus `NKG7/GNLY/KLRD1/PRF1/GZMB`, and split B vs plasma with `MS4A1/CD79A/CD79B/BANK1` versus `MZB1/JCHAIN/SDC1/IGHG1/IGKC`.

## Out of scope (FYI, not actionable for the Generator)
- The developmental/progenitor axes are not automatically wrong because age is unspecified, but they should not dominate the design unless a future scenario explicitly requests fetal/developmental ordering.

## Summary for the next round
Keep the overall pancreas structure, but add the missing mesothelial and adipocyte/lipid stromal coverage, delete or strongly recast the enterochromaffin-like axis, and clean up the lymphatic/ductal marker choices. If dimensions are added, consider splitting T/NK and B/plasma immune axes rather than adding more broad contrasts.
