# Discriminator report

## Scenario echo
The scenario is human pancreas with no specified age, sex, disease, treatment, or focal cell type, and no additional notes. The embedding should therefore cover the normal pancreatic epithelial/exocrine, endocrine islet, stromal, vascular, neural/glial, and resident immune landscape without assuming a specific disease or developmental program.

## Rating: minor
The blueprint has the right broad structure and mostly plausible human pancreatic marker modules. It is not quite publishable unchanged because it misses/under-resolves a few expected normal support compartments and has a couple of marker/design choices that should be cleaned up.

## Findings

### Minor issues
- **vascular compartment / endothelial axis** — The endothelial module is essentially blood vascular endothelium and does not capture lymphatic endothelium, which is a normal pancreatic support compartment and is expected in broad tissue coverage. Markers such as PECAM1/CDH5/VWF/KDR/CLDN5/RAMP2 are appropriate for vascular endothelium but do not resolve LYVE1/PDPN/PROX1/FLT4/CCL21-type lymphatic identity.
  *Fix:* Add a lymphatic endothelial module and either a separate plain axis or an endothelial-gated contrast/refinement axis for blood vascular versus lymphatic endothelium.

- **mural_pericyte_smooth_muscle** — Pericytes and vascular smooth muscle cells are merged into one axis. The genes are plausible mural markers, but RGS5/PDGFRB/MCAM/CSPG4/NOTCH3 and ACTA2/TAGLN/MYH11 represent related but distinct mural programs. A single combined score will blur pericytes with smooth muscle in a broad pancreas embedding.
  *Fix:* Split or refine the mural compartment into pericyte-enriched and smooth-muscle-enriched modules/axes, for example preserving a broad mural axis only if paired with a gated pericyte-vs-smooth-muscle refinement.

- **stellate_fibroblast** — The blueprint has one generic stellate/fibroblast axis, but healthy human pancreas commonly contains pancreatic stellate/fibroblast populations with quiescent ECM-producing and activated/myofibroblast-like programs. The current module is dominated by collagen/ECM genes plus CYGB, while ACTA2/TAGLN are assigned to the merged mural module, so activated stellate/myofibroblast-like state may be conflated with vascular smooth muscle.
  *Fix:* Add a stellate/fibroblast refinement axis separating quiescent stellate/fibroblast identity from activated myofibroblast-like activation, scoped to the stromal/stellate compartment.

- **epsilon_cell** — The epsilon-cell axis uses only GHRL. GHRL is canonical, but a one-gene module is fragile and makes this very rare endocrine population overly dependent on one detected transcript.
  *Fix:* Add additional human pancreas ghrelin/epsilon-enriched markers if literature-supported, or keep GHRL but treat epsilon as a very rare endocrine refinement rather than an equal-resolution subtype.

- **schwann_glial** — Most markers are appropriate Schwann/glial genes, but NRXN1 is more neuronal/synaptic than a canonical Schwann marker and is not needed for this module.
  *Fix:* Replace NRXN1 with a stronger peripheral glial/Schwann marker supported by the cited atlas resources, such as NGFR, ERBB3, MBP, GFRA3, or another atlas-supported Schwann/glial marker.

## Out of scope (FYI, not actionable for the Generator)
- I did not re-check DOI correctness or IF thresholds; citation verification is outside this Discriminator pass.
- The REG-positive acinar state is acceptable as a healthy-atlas acinar state/refinement and does not need removal.

## Summary for the next round
Keep the overall blueprint structure, but refine the normal support compartments: add lymphatic endothelial coverage, split or refine pericyte versus vascular smooth muscle, and add a stromal stellate/fibroblast activation refinement. Then clean up the fragile epsilon module and the Schwann/glial marker choice without otherwise overhauling the endocrine or immune layout.
