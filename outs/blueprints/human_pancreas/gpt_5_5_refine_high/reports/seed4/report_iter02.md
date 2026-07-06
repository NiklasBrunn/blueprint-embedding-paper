# Discriminator report

## Scenario echo
The scenario is healthy/baseline human pancreas with no specified age, sex, cell-type subset, or design notes. The embedding should therefore cover the major pancreatic epithelial, endocrine, vascular, stromal, immune, neural, and common state compartments expected in human pancreas without adding disease- or development-specific programs.

## Rating: minor
The blueprint has broad and biologically plausible coverage of the main human pancreas lineages, including endocrine subtypes, acinar and ductal epithelium, vascular/mural/stromal compartments, immune subsets, neural/Schwann cells, adipocytes, and cycling state. It needs targeted cleanup of a few markers, labels, and one contrast axis, but the overall structure is sound.

## Findings

### Minor issues
- **module:lymphatic_endothelial** - The module is mostly canonical for lymphatic endothelium, but `TFF3` is an epithelial/ductal trefoil marker rather than a lymphatic endothelial marker. In this scenario it risks pulling ductal signal into the lymphatic axis.
  *Fix:* Replace `TFF3` with a more canonical lymphatic/endothelial gene such as `NRP2`, `RELN`, `TBX1`, `FOXC2`, or another literature-supported lymphatic marker.

- **module:leukocyte_common** - `HLA-DRA` is not a pan-leukocyte marker; it biases the common leukocyte gate toward antigen-presenting cells and overlaps with the dendritic/B-cell modules. The gate should remain broad for T, NK, myeloid, B/plasma, mast, and DC refinement.
  *Fix:* Swap `HLA-DRA` for a broader leukocyte marker such as `CD53`, `LAPTM5`, `DOCK2`, `RAC2`, or another supported pan-immune gene.

- **module:t_nk_cell / axis:t_nk_within_leukocyte and module:b_plasma_cell / axis:b_plasma_within_leukocyte** - The biology is mostly represented, but the names are misleading. `t_nk_cell` is a T-cell module, while NK cells are already represented separately by `nk_cytotoxic_cell`; `b_plasma_cell` is a B-cell module, while plasma cells are already represented separately by `plasma_cell`.
  *Fix:* Rename or relabel these in place as T-cell and B-cell modules/axes, or adjust the gene lists if the intent is genuinely combined T/NK or B/plasma programs. No deletion is needed.

- **axis:mural_contractile_vs_endothelial** - Smooth muscle/contractile mural cells and endothelial cells are distinct vascular-associated compartments, but they are not a clear biological opposition in the same sense as acinar-vs-ductal or subtype-within-gate axes. Since an endothelial axis already exists, this contrast is slightly arbitrary.
  *Fix:* Convert this in place to a plain smooth muscle/contractile mural axis, or gate it within a vascular/stromal context if supported by existing modules.

## Out of scope (FYI, not actionable for the Generator)
- Epsilon cells and adipocytes are rare or sampling-dependent in many adult human pancreas datasets, but they are still within the broad human pancreas scenario and do not require removal.

## Summary for the next round
Focus on targeted biological cleanup rather than structural redesign: replace the non-lymphatic `TFF3` marker, make the leukocyte gate more pan-immune, correct the misleading T/NK and B/plasma labels or compositions, and make the smooth muscle/mural axis a direct compartment axis instead of contrasting it against endothelial cells.
