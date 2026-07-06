# Discriminator report

## Scenario echo
The scenario is healthy/baseline human pancreas with no age, sex, condition, or note-specific emphasis. Because `celltype` is null, the embedding should cover the normal pancreas cell-type landscape: endocrine islet cells, exocrine acinar and ductal epithelium, vascular and lymphatic endothelium, mesenchymal/stromal and mural cells, resident immune populations, and peripheral neural/Schwann-like cells when present in atlas data.

## Rating: minor
The blueprint is broadly biologically plausible and covers the major human pancreas compartments, including endocrine subtypes and key exocrine/stromal/vascular populations. It needs targeted refinement because some expected immune populations are collapsed into mixed modules/axes, and one stromal axis is a little too coarse for pancreatic fibroblast/stellate biology.

## Findings

### Minor issues
- **immune_t_nk_identity / t_nk_cell and immune_b_plasma_identity / b_plasma_cell** — The blueprint merges biologically distinct lymphoid populations into two mixed modules: T and NK cells are combined, and B and plasma cells are combined. For an all-cell healthy human pancreas embedding, resident or contaminating T cells, NK/cytotoxic lymphocytes, B cells, and antibody-secreting plasma cells are expected enough that the embedding should not force them onto the same subtype axes. The current marker sets are individually plausible but internally mixed: `CD3D/CD3E/TRAC/IL7R` are T-cell markers, `NKG7/KLRD1/GNLY/PRF1` are NK/cytotoxic markers, `MS4A1/CD79A/CD79B/BANK1/IGHM` are B-cell markers, and `JCHAIN/MZB1/SDC1` are plasma-cell markers.
  *Fix:* Replace these two broad lymphoid modules/axes with separate gated immune axes for T cells, NK/cytotoxic cells, B cells, and plasma cells, using the current markers split into coherent modules plus canonical additions if needed, such as `XBP1` or `PRDM1` for plasma cells.
- **mesenchymal_fibroblast_stellate_identity / fibroblast_stellate** — The module mostly captures ECM-producing fibroblast/activated stellate signal (`COL1A1`, `COL1A2`, `COL3A1`, `MMP2`, `DCN`, `LUM`) and does not clearly distinguish pancreatic fibroblasts from pancreatic stellate/activation state. Pancreatic stellate cells are a standard stromal population in pancreas datasets, and conflating them entirely with generic fibroblast ECM may reduce biological resolution.
  *Fix:* Add or revise mesenchymal gated axes so stromal scoring can separate fibroblast/ECM identity from stellate or myofibroblast/activated-stellate state. Keep the existing mural/pericyte axis, but consider an additional stellate/activation module using canonical pancreatic stellate or myofibroblast markers where supported by the declared literature.
- **gamma_pp_cell and epsilon_cell marker hygiene** — These endocrine rare-cell modules include several broadly endocrine or overlapping markers (`IAPP`, `PCSK1`, `ISL1`, `RBP4`) that may weaken subtype specificity against beta, delta, and endocrine-core axes. `PPY` and `GHRL` anchor the intended populations, so this is not a wrong-lineage problem, but the rare endocrine axes would be cleaner with fewer shared endocrine genes.
  *Fix:* Trim or replace the least specific overlapping genes in the PP/gamma and epsilon modules with more subtype-selective markers from the cited pancreas atlases, while retaining enough genes for robust scoring.

## Out of scope (FYI, not actionable for the Generator)
- Citation/DOI verification was not reassessed here, per task scope.
- The extra `epithelial_ductal_bias` axis is somewhat redundant with the ductal plain axis and acinar-to-ductal contrast, but it is not biologically wrong and does not by itself justify another refinement round.

## Summary for the next round
Focus on splitting the mixed lymphoid modules into coherent T, NK/cytotoxic, B, and plasma-cell axes. Secondarily, improve stromal resolution for fibroblast versus stellate/activation biology and clean up rare endocrine marker overlap. The core pancreas compartment structure is otherwise sound.
