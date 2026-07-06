# Discriminator report

## Scenario echo
The scenario is healthy/baseline human pancreas with no specified age, sex, condition, or requested design focus. Because `celltype` is null, the embedding should cover the expected whole-pancreas cell landscape: endocrine islet subtypes, acinar and ductal epithelium, vascular and lymphatic endothelium, mural cells, fibroblast/stellate compartments, resident immune populations, neural/Schwann elements, and broadly applicable states such as proliferation.

## Rating: minor
The blueprint is broadly plausible and structurally valid, with good coverage of the major epithelial, endocrine, stromal, immune, neural, and cycling programs. It is not quite publishable unchanged because several expected pancreas populations are absent or collapsed into mixed modules that would blur biologically distinct cell types.

## Findings

### Minor issues
- **vascular compartment / endothelial module and axes** — The blueprint has a generic blood endothelial axis, but no lymphatic endothelial representation. Lymphatic endothelial cells are an expected minor vascular population in human pancreas atlases and should not be forced into the PECAM1/VWF/KDR blood endothelial program.
  *Fix:* Add a lymphatic endothelial module and axis using markers such as `PROX1`, `LYVE1`, `PDPN`, `FLT4`, `CCL21`, `MMRN1`, and/or `TFF3`, ideally contrasted or gated against the existing blood endothelial module.
- **pericyte_smooth_muscle / mural axes** — The module mixes pericyte markers (`RGS5`, `PDGFRB`, `CSPG4`, `NOTCH3`) with vascular smooth muscle markers (`ACTA2`, `TAGLN`, `MYH11`). These are related mural populations but are separable in single-cell pancreas data; one contrast against endothelium cannot resolve pericytes versus vascular smooth muscle.
  *Fix:* Split the current mural program into pericyte and smooth-muscle/contractile modules, and add an axis that distinguishes them while preserving a broad mural-vs-endothelial signal if desired.
- **immune subtype coverage** — The leukocyte branch covers macrophage, T/NK, B/plasma, and mast cells, but it lacks dendritic cell representation. Resident/classical dendritic and antigen-presenting myeloid cells are expected in a whole human pancreas landscape and are not well captured by the macrophage-heavy `LST1`/`C1Q` module.
  *Fix:* Add a dendritic-cell module and gated leukocyte axis, for example with markers such as `CD1C`, `FCER1A`, `CLEC10A`, `CLEC9A`, `LILRA4`, `IRF8`, `ITGAX`, and `HLA-DRA`, choosing the subset supported by the cited pancreas/immune atlas resources.
- **t_nk_cell and b_plasma_cell modules / leukocyte gated axes** — The current T/NK module combines T-cell markers (`CD3D`, `CD3E`, `TRAC`, `IL7R`) with cytotoxic NK markers (`NKG7`, `GNLY`, `KLRD1`, `GZMB`), and the B/plasma module similarly combines naive/B-cell markers (`MS4A1`, `CD79A`, `BANK1`) with plasma markers (`MZB1`, `JCHAIN`, immunoglobulins). These are biologically distinct leukocyte populations, so the gated axes will blur cell identity.
  *Fix:* Split T from NK/cytotoxic and B from plasma into separate modules and gated axes under `leukocyte_common`; retain broad immune axes only if they serve a separate coarse-grained purpose.
- **epsilon_cell marker specificity** — The epsilon-cell module is very small in unique ghrelin-cell signal and includes broad endocrine markers plus `MLN`, which is more classically enteroendocrine/intestinal than a canonical pancreatic epsilon marker. Because epsilon cells are rare and age-dependent, this axis risks weak or off-target scoring.
  *Fix:* Keep epsilon only if the cited human pancreas resources support it for the age-unspecified scenario, but make the module more pancreas-specific around `GHRL` and supported co-markers; replace or drop questionable `MLN` if unsupported.

## Out of scope (FYI, not actionable for the Generator)
- I did not verify DOI/resource correctness; this report only judges biological coverage and marker plausibility.

## Summary for the next round
Focus the next refinement on adding lymphatic endothelial and dendritic-cell coverage, then split the over-broad mural and immune subtype modules so the embedding resolves pericyte versus smooth muscle, T versus NK, and B versus plasma cells. After that, tighten the epsilon-cell marker set so rare endocrine scoring is not driven mainly by generic endocrine genes.
