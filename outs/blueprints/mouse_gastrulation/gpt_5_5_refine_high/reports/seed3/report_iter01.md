# Discriminator report

## Scenario echo
Mouse gastrulating embryo, mixed sex, wild-type/baseline, spanning E6.5-E8.5. The notes explicitly ask for a time-forward hierarchical gating tree in which non-root axes are gated by their parent program, sibling modules are used as negatives at branch points, and the embedding covers the full developmental landscape while preserving maturation gradients.

## Rating: minor
The blueprint is broadly plausible and covers many expected E6.5-E8.5 embryonic lineages: epiblast/primitive streak, major germ-layer derivatives, early neural, endoderm, mesoderm, blood/endothelium, trophoblast/extraembryonic ectoderm, visceral/parietal endoderm, PGCs, and cycling. It still needs targeted additions and local branch/marker cleanup before I would defend it unchanged for the full gastrulating embryo landscape.

## Findings

### Minor issues
- **extraembryonic and trophoblast coverage** — The extraembryonic side is under-resolved for E6.5-E8.5. There is an extraembryonic ectoderm root and one `chorionTrophoblast` child, but no separate ectoplacental cone / trophoblast giant cell / spongiotrophoblast-like branch, and the `chorionTrophoblast` module mixes chorion genes with `Tpbpa`/`Prl3d1`-type trophoblast differentiation markers. The extraembryonic mesoderm branch also only resolves allantois and lacks amnion/yolk-sac mesoderm or mesothelium-like states expected in this window.
  *Fix:* Split or refine the trophoblast branch under `extraembryonicEctoderm` into biologically meaningful sibling modules such as chorion and ectoplacental cone/TGC/spongiotrophoblast-like trophoblast, with sibling negatives. Add an extraembryonic mesoderm continuation under `nascentMesoderm` or `allantoisExtraMesoderm` for amnion/yolk-sac mesoderm/mesothelium where supported, rather than collapsing all extraembryonic mesoderm into allantois.

- **primitive/parietal/visceral endoderm coverage** — `visceralEndoderm` and `parietalEndoderm` roots are useful, but the visceral endoderm branch is too coarse for early gastrulation. E6.5-E7.5 datasets commonly distinguish anterior/distal visceral endoderm patterning from more general visceral endoderm, and this matters because it is one of the earliest axis-patterning states in the scenario.
  *Fix:* Add a gated or plain early visceral-endoderm refinement for AVE/DVE/anterior visceral endoderm markers if supported by the cited literature. Keep it within the visceral endoderm lineage rather than treating it as an epiblast-derived branch.

- **neural and caudal patterning coverage** — The ectoderm branch captures neural ectoderm, neural tube, surface ectoderm, neural crest, and NMP, but it does not preserve much of the E8.0-E8.5 anterior-posterior neural maturation gradient. `neuralEctoderm` is heavily anterior (`Otx2`, `Hesx1`, `Six3`), while `neuralTube` mixes generic neural genes with later dorsoventral spinal markers (`Nkx6-1`, `Dbx1`, `Olig2`) and has no forebrain/midbrain/hindbrain/posterior-neural separation.
  *Fix:* Add one or more time-forward neural patterning axes under `neuralEctoderm`/`neuralTube`, for example anterior neural versus posterior/hindbrain-spinal neural programs, using sibling negatives. Consider replacing very ventral/later neural tube markers in the broad `neuralTube` module with earlier generic neural tube markers, reserving DV markers for a later/refined branch only if appropriate.

- **paraxial/NMP branch design** — `somite_somiticMesoderm` is gated by `paraxialMesoderm` but uses `caudalNeuralNmp` as the negative. NMP is an upstream/adjacent bipotent state, not the natural sibling of somitic mesoderm within paraxial mesoderm. The branch also lacks a presomitic mesoderm axis, which is a key E7.5-E8.5 maturation state between nascent/paraxial mesoderm and somites.
  *Fix:* Add a presomitic/paraxial maturation module and contrast presomitic versus somitic states under `paraxialMesoderm`. Keep NMP as its own posterior progenitor branch, or use it as an upstream parent only where the tree semantics are biologically justified.

- **non-sibling negatives in gated axes** — A few gated axes violate the note's sibling-negative intent. `troph_chorionTrophoblast` contrasts chorion/trophoblast against `visceralEndoderm` and `allantoisExtraMesoderm`, which are not sibling states under extraembryonic ectoderm. `cardiac_cardiomyocyte` contrasts cardiomyocytes against endothelium and blood progenitors, which are not children of `cardiacMesoderm`.
  *Fix:* For these axes, either use no negatives for a parent-scoped maturation signal or introduce true sibling modules at the same branch point. For cardiac mesoderm, possible siblings include cardiomyocyte versus pharyngeal/splanchnic/second-heart-field-like mesoderm if the literature resources support them.

- **marker plausibility in late/differentiated modules** — Several marker sets lean too mature or too mixed for the E6.5-E8.5 window. `somiticMesoderm` includes `Pax7`, `MyoD1`, and `Myf5`, which bias toward later skeletal muscle/dermomyotome rather than broad early somite. `myeloid` includes granulocyte/adult-like markers such as `Mpo`, `Itgam`, and `Ctsg`; primitive yolk-sac myeloid/macrophage markers would be a better fit for E8.5. `axialMesendoderm` includes `Tbx6`, which is more paraxial/presomitic than axial node/notochord.
  *Fix:* Replace these with earlier stage-appropriate markers: early somite/presomitic markers for the paraxial branch, primitive macrophage/myeloid markers for the blood branch, and cleaner node/notochord/axial mesendoderm markers for `axialMesendoderm`.

## Out of scope (FYI, not actionable for the Generator)
- I did not verify DOIs or whether each cited paper directly supports every gene list; that is outside the Discriminator scope.

## Summary for the next round
Keep the overall time-forward tree, but refine the missing E6.5-E8.5 coverage: extraembryonic trophoblast/mesoderm sub-branches, AVE/visceral endoderm patterning, neural AP maturation, and presomitic-to-somitic progression. Also adjust non-sibling negatives and swap a few mature or misplaced markers so the gated axes better match biological parent-child and sibling relationships.
