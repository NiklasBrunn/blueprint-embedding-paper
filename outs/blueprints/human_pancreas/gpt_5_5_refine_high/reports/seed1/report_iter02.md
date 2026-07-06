# Discriminator report

## Scenario echo
Human pancreas, with no requested restriction by age, sex, cell type, or disease condition; notes are null. The embedding therefore should represent the broad healthy human pancreatic tissue landscape: endocrine, exocrine epithelial, vascular, stromal/mesenchymal, immune, and neural/ganglial resident populations without imposing a developmental or disease-specific structure.

## Rating: minor
The blueprint captures the major pancreatic compartments and most canonical endocrine, exocrine, endothelial, immune, and neural/stromal classes. It still under-resolves a few expected resident populations and conflates some biologically distinct subtypes, so targeted refinement is warranted.

## Findings

### Minor issues
- **mesenchymal/stromal modules and axes** — The stromal compartment is dominated by fibroblast/ECM markers, and the `mural_pericyte_smooth_muscle` module combines pericytes and vascular smooth muscle in one axis. Pancreatic stellate/fibroblast-like cells, pericytes, and smooth muscle cells are expected but not cleanly separable as represented.
  *Fix:* Refine the stromal branch with a broad mesenchymal gate and distinct axes/modules for pancreatic stellate/fibroblast/ECM cells, pericytes, and vascular smooth muscle where supported. Reduce the heavy overlap between `stromal_core` and `fibroblast_matrix` so the gate is not essentially a duplicate fibroblast program.

- **neural_ganglial module / neural_ganglial_compartment axis** — The neural module is mostly Schwann/glial (`SOX10`, `PLP1`, `S100B`, `MPZ`) with only weak neuronal/ganglionic representation. Human pancreas can contain autonomic neuronal/ganglionic elements in addition to Schwann/glial cells, so the current single axis may miss neurons or collapse them into a glial score.
  *Fix:* Either rename/refocus the current module as Schwann/glial and add a separate neuronal/ganglionic module/axis using canonical neuronal/autonomic markers, or add a gated neural contrast that separates glial from neuronal/ganglionic identity.

- **immune myeloid coverage** — `macrophage` is a resident/tissue macrophage program and `dendritic_apc` is a broad composite of cDC1, cDC2, and pDC markers. A healthy human pancreas atlas can also contain monocytes and finer dendritic/APC subsets; these are not represented separately.
  *Fix:* Add at least a monocyte/myeloid-inflammatory axis under `immune_core` if supported by the cited pancreas/atlas literature. If dimensionality allows, consider splitting `dendritic_apc` into cDC and pDC modules rather than one mixed APC axis.

- **ductal epithelial resolution** — `ductal_cell` mixes canonical duct markers with remodeling/activated duct-associated genes such as `SPP1` and `MMP7`, and there is no ductal sub-axis for centroacinar/terminal duct or activated/secretory duct states. For a whole-pancreas embedding, ductal heterogeneity is modest but biologically expected.
  *Fix:* Keep the broad ductal axis, but add one gated ductal refinement if literature support is available, such as canonical/CFTR-high duct versus SPP1/MMP7-high remodeling duct, or a centroacinar/terminal duct program.

## Out of scope (FYI, not actionable for the Generator)
- I did not verify paper, atlas, or database citation handles; this review is limited to biological coverage and marker plausibility.

## Summary for the next round
The next Generator should keep the overall compartmental structure but refine under-resolved edge populations: split the stromal/mesenchymal branch into stellate/fibroblast, pericyte, and smooth muscle programs; separate Schwann/glial from neuronal/ganglionic identity; and add a monocyte/APC or ductal-substate refinement if supported by the existing literature base.
