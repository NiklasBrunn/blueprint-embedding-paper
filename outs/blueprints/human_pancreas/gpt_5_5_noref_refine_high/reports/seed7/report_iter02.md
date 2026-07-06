# Discriminator report

## Scenario echo
The scenario is human pancreas with no specified age, sex, disease, treatment, or narrower cell type. With `celltype` and `notes` both null, the blueprint should cover the broad healthy human pancreatic tissue landscape: exocrine acinar/ductal/centroacinar epithelium, endocrine islet cell types, stromal/mural/vascular/lymphatic compartments, neural/glial and mesothelial/adipose neighbors, and resident or contaminating immune populations and common capture/stress states.

## Rating: minor
The blueprint is broadly biologically plausible and covers the major pancreatic compartments with mostly canonical human markers. It needs targeted refinement rather than structural rework, chiefly to improve immune resolution and to add a general endothelial identity axis.

## Findings

### Minor issues
- **immune axes (`t_nk_cell`, `b_plasma_cell`, `t_nk_within_immune`, `b_plasma_within_immune`)** — T cells and NK cells are collapsed into one module/axis, and B cells and plasma cells are collapsed into one module/axis. For a full human pancreas embedding, these are expected immune populations and can carry distinct biology; a single mixed axis will blur them.
  *Fix:* Split or supplement these with separate gated immune axes for T cells (`CD3D`, `CD3E`, `TRAC`, `CD2`, `IL7R`), NK cells (`NKG7`, `GNLY`, `PRF1`, `KLRD1`, `FCGR3A`), B cells (`MS4A1`, `CD79A`, `CD79B`, `BANK1`, `CD74`), and plasma cells (`MZB1`, `XBP1`, `JCHAIN`, `SDC1`, `IGKC`). A separate pDC axis is optional if dimensionality allows.
- **endothelial compartment (`endothelial_core`, vascular/lymphatic axes)** — Endothelial cells are represented only through gated vascular and lymphatic subtype axes. This captures major subtypes but lacks a plain pan-endothelial identity dimension, so endothelial cells with weak subtype marker expression may be underrepresented in the embedding.
  *Fix:* Add a plain `endothelial_identity` axis using `endothelial_core`, then keep vascular and lymphatic gated refinement axes. If adding resolution, distinguish arterial/venous/capillary or islet/acinar endothelial programs only with canonical markers and without over-padding.
- **endocrine developmental/plasticity states (`endocrine_progenitor_like_state`, `immature_beta_like`)** — The NEUROG3/PAX4/SOX9/HES1-centered modules are biologically plausible for developmental or dedifferentiation-like programs, but the scenario has no developmental or disease note. As currently framed, these axes may over-emphasize rare or context-dependent progenitor-like states relative to the baseline human pancreas.
  *Fix:* Keep this as a scoped low-priority endocrine plasticity/dedifferentiation state, or retune it toward adult endocrine immaturity/stress markers. Do not let it substitute for core adult islet subtype coverage.
- **mural/stromal axis naming and design (`pericyte_smooth_muscle_vs_fibroblast`)** — The axis name implies a contrast against fibroblasts, but the implementation is a plain pericyte/smooth muscle module. The biology is mostly plausible, but the naming and axis type are inconsistent and may confuse the intended interpretation.
  *Fix:* Either rename it to a plain pericyte/mural identity axis or make it an actual contrast/gated refinement against fibroblast/stromal modules. Keep contractile smooth muscle as a separate mural state.

## Out of scope (FYI, not actionable for the Generator)
- The broad stress axes are acceptable for scRNA-seq pancreas data, but they should remain secondary to lineage coverage in downstream interpretation.
- Some marker overlap is expected in pancreas modules, especially KRT8/KRT18 in epithelial and mesothelial contexts and matrix genes across stromal states.

## Summary for the next round
Focus on targeted refinement: split the mixed immune axes into separate T, NK, B, and plasma dimensions; add a plain pan-endothelial identity axis before subtype refinements; and clarify the endocrine progenitor/plasticity and mural axis intent. The overall cell-type coverage and marker choices are otherwise defensible.
