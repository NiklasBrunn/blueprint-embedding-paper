# Discriminator report

## Scenario echo
The scenario is human pancreas with no specified age, sex, disease, or focused notes, and `celltype: null`, so the embedding should cover the normal pancreatic cell-type/state landscape across endocrine islet, exocrine epithelial, vascular, stromal, neural, and resident immune compartments without adding disease-only programs.

## Rating: minor
The blueprint is broadly biologically plausible and covers the major human pancreas compartments: acinar, ductal, endocrine alpha/beta/delta/PP/epsilon, endothelial subtypes, stromal/mural cells, immune cells, mast cells, and neural/Schwann cells. It still has a few coverage and marker-resolution gaps that are worth a targeted refinement round, especially plasma/antibody-secreting immune cells and activated stellate/myofibroblast state separation.

## Findings

### Minor issues
- **immune compartment** — The immune tree includes macrophage/myeloid, T, cytotoxic/NK, dendritic/APC, B, and mast axes, but it lacks a plasma-cell or antibody-secreting B-lineage axis. Adult human pancreas datasets can contain rare plasma/Ig-high cells, and the current `b_cell` module (`MS4A1`, `CD79A/B`, `BANK1`) will preferentially capture conventional B cells rather than `MZB1`/`JCHAIN`/immunoglobulin-high plasma cells.
  *Fix:* Add a `plasma_cell` or `antibody_secreting_cell` module with markers such as `MZB1`, `JCHAIN`, `XBP1`, `SDC1`, `IGHG1`, `IGKC`, `IGLC2` as appropriate for the cited atlas evidence, add an immune-gated plasma axis, and include this module among sibling negatives for the other immune subtype axes if retained.

- **stromal compartment** — The stromal branch covers broad ECM fibroblast/stromal signal, mural/pericyte/smooth muscle, and quiescent stellate cells, but it does not explicitly represent activated pancreatic stellate/myofibroblast state. That state is biologically important in pancreas and can appear even in otherwise healthy atlas material as tissue handling or baseline remodeling signal. The current `pericyte_smooth_muscle` module contains contractile genes (`ACTA2`, `TAGLN`, `MYL9`) that can blur activated stellate/myofibroblast signal with mural identity.
  *Fix:* Add an activated stellate/myofibroblast module gated by `stromal_matrix` using markers such as `ACTA2`, `TAGLN`, `COL1A1`, `COL3A1`, `POSTN`, `THY1`, `MMP2`, or `COL11A1` where supported. Consider tightening the mural module toward pericyte/smooth-muscle identity (`RGS5`, `PDGFRB`, `CSPG4`, `MCAM`, `NOTCH3`, plus contractile markers as needed) and using contrast/gating so activated stellate and mural cells are not collapsed.

- **epsilon_cell and minor endocrine modules** — The common endocrine identities are present, but some rare-endocrine marker choices are noisy. `epsilon_cell` relies on `GHRL` plus `MLN`, `ACSL1`, and `DEFB1`; `DEFB1` is more epithelial/duct-associated than canonical ghrelin-cell identity, and `MLN` marks motilin endocrine biology rather than a pancreas-specific epsilon signature. `gamma_pp_cell` and `alpha_cell` also share `TTR`, which may reduce subtype specificity.
  *Fix:* Tighten rare endocrine modules around the most pancreas-supported markers. For epsilon cells, keep `GHRL` and replace weak epithelial/off-lineage genes with better-supported ghrelin-cell markers if the cited human pancreas atlases support them; otherwise keep the epsilon axis very small and conservative. Reduce unnecessary overlap such as repeated `TTR` unless the literature source specifically uses it for both populations.

- **lymphatic_endothelium** — The lymphatic module has canonical anchors (`PROX1`, `PDPN`, `LYVE1`) but also includes less canonical pancreas-atlas-specific genes (`NTS`, `NRG3`, `SEMA3D`) while omitting common lymphatic markers such as `FLT4`, `CCL21`, `MMRN1`, or `TFF3` depending on source support.
  *Fix:* Keep the canonical lymphatic anchors and either add or substitute additional well-supported lymphatic endothelial markers from the cited pancreas endothelial resource so the axis is less dependent on noncanonical genes.

## Out of scope (FYI, not actionable for the Generator)
- Citation verification was not assessed; this report only judges biological plausibility, marker fit, and coverage.

## Summary for the next round
The next Generator should preserve the overall structure, which is largely correct for human pancreas, and make targeted additions/tightening: add a plasma-cell immune axis, add or separate activated stellate/myofibroblast biology from mural cells, and clean up noisy rare endocrine and lymphatic marker choices. No existing axis or module needs deletion solely on biological grounds.
