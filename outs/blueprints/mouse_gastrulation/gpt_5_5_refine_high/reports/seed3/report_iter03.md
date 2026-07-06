# Discriminator report

## Scenario echo
Mouse gastrulating embryo, mixed sex, wild-type/baseline, E6.5-E8.5, with all cell types in the tissue in scope. The notes require a time-forward hierarchical embedding: earliest/root programs are plain axes, later states are gated by their parent module, and branch-point siblings should be used as negatives where appropriate while preserving maturation gradients.

## Rating: minor
The blueprint is broadly plausible and covers the main E6.5-E8.5 mouse gastrulation landscape: epiblast, primitive streak, germ-layer derivatives, extraembryonic endoderm/ectoderm, major mesoderm branches, early blood/endothelium, neural derivatives, PGCs, and cycling. It needs targeted refinement for extraembryonic trophoblast resolution, endoderm trajectory specificity, and a few marker/negative-set choices, but the overall tree shape is not broken.

## Findings

### Minor issues
- **trophoblast branch: `extraembryonicEctoderm`, `chorionTrophoblast`, `ectoplacentalConeTrophoblast`** - The trophoblast branch misses trophoblast giant cells, which are expected by E6.5-E8.5 and around the ectoplacental cone/early placental compartment. The current `ectoplacentalConeTrophoblast` module also mixes EPC/progenitor markers with differentiated TGC-like prolactin/cathepsin markers, so one axis is carrying two distinct states.
  *Fix:* Add a trophoblast giant-cell module and gated child axis, preferably under a refined EPC/progenitor parent if supported. Keep `extraembryonicEctoderm -> chorion/EPC` as the upstream split, but move strongly differentiated TGC markers out of the EPC program and into the new child module.

- **endoderm branch: `posteriorVisceralGutEndoderm`, `gutTubeEndoderm`, `foregutEndoderm`, `hindgutEndoderm`** - The endoderm coverage recognizes both visceral endoderm and definitive endoderm/gut, but `posteriorVisceralGutEndoderm` is biologically muddled: its genes mix general visceral endoderm/yolk-sac endoderm, posterior/hindgut patterning, and definitive endoderm markers. Since the scenario window includes convergence/intercalation of VE-derived and primitive-streak-derived endoderm into gut, this branch should preserve that trajectory more explicitly.
  *Fix:* Rework `posteriorVisceralGutEndoderm` in place into a cleaner posterior/extraembryonic visceral endoderm or VE-to-hindgut-transition module. Keep the definitive-endoderm-derived `gutTubeEndoderm -> foregut/hindgut` path separate, and use sibling negatives only where they discriminate true branch alternatives.

- **blood branch: `blood_erythroid`, `blood_myeloid`, `blood_megakaryocyteEmp`** - These axes are correctly gated by `bloodProgenitor`, but they include `endothelium` as a negative even though endothelium is a sibling of `bloodProgenitor` one level up, not a sibling of erythroid/myeloid/megakaryocyte fate under blood progenitor. This weakens the time-forward tree semantics requested in the notes.
  *Fix:* Remove `endothelium` from the negative sets on the blood-progenitor child axes. If additional resolution is desired, add a hemogenic-endothelium/EMP transition axis between `hematoEndothelial` and differentiated blood fates rather than contrasting late blood children directly against endothelium.

- **parietal endoderm module: `parietalEndoderm`** - The module contains several plausible matrix/parietal endoderm genes, but also mesenchymal/general EMT genes such as `Snai1`, `Vim`, and shared `Pdgfra`, making it less specific than it should be for parietal endoderm in this window.
  *Fix:* Replace the least specific EMT/mesenchymal genes with more canonical parietal endoderm markers where supported, such as `Gata4`, `Sox7`, `Pth1r`, `Tinagl1`, `Srgn`, or other atlas-supported parietal endoderm-enriched genes.

- **coarse extraembryonic mesoderm: `amnionYolkSacMesoderm`** - Amnion and yolk-sac/extraembryonic mesoderm are merged into one module and axis. This is acceptable as a broad compartment, but it underserves the "full cell landscape" instruction because these are distinct extraembryonic mesoderm states across E7.5-E8.5.
  *Fix:* If literature support in the existing atlas resources is clear, split this into separate amnion and yolk-sac/extraembryonic mesoderm axes under `nascentMesoderm`, with allantois kept as its own sibling.

## Out of scope (FYI, not actionable for the Generator)
- I did not verify resource DOIs or journal impact factors; that is outside this discriminator pass.

## Summary for the next round
Keep the current time-forward tree, but refine the extraembryonic side of the landscape: add/separate trophoblast giant-cell biology from EPC, clean the visceral-endoderm-to-gut trajectory, and tighten blood-child negatives plus a few marker modules. These are targeted edits rather than a structural rebuild.
