# Discriminator report

## Scenario echo
Human pancreas with no age, sex, or disease restriction and no special design notes. Because `celltype` is null, the embedding should cover the major epithelial endocrine/exocrine compartments plus vascular, stromal, neural/glial, immune, and common resident/rare pancreatic cell states expected in healthy pancreas tissue.

## Rating: minor
The blueprint is broadly plausible and covers most major human pancreas lineages, including acinar, ductal, endocrine alpha/beta/delta/PP/epsilon, stromal, endothelial, immune, mast, neural/glial, and generic stress/cycling programs. It still needs targeted coverage and marker/module cleanup before I would defend it unchanged.

## Findings

### Minor issues
- **coverage: mesothelial and rare duct/tuft-like epithelial cells** — Whole human pancreas preparations can include mesothelial surface cells, and modern pancreas atlases also resolve rare tuft-like ductal/epithelial populations. The current blueprint has no mesothelial axis, and `MSLN` is instead embedded in the ductal module where it can blur ductal and mesothelial identity. Tuft-like epithelial biology is also absent.
  *Fix:* Add a mesothelial module/axis using markers such as `WT1`, `UPK3B`, `MSLN`, `CALB2`, `KRT8`, `KRT18`, `PRG4`, and optionally add a gated duct/epithelial tuft-like axis using markers such as `POU2F3`, `TRPM5`, `DCLK1`, `CHAT`, `AVIL`, `GFI1B`, `HPGDS`.
- **ductal_epithelial / centroacinar_ductal / ductal_centroacinar_refinement** — The ductal and centroacinar modules are nearly the same marker set, so the gated centroacinar axis does not really resolve centroacinar cells from ordinary ductal epithelium. `MSLN` and `TFF1` in `ductal_epithelial` are also more suggestive of mesothelial or mucinous/stress programs than a clean baseline duct identity.
  *Fix:* Make `ductal_epithelial` a cleaner broad duct module (`KRT19`, `KRT7`, `SOX9`, `HNF1B`, `CFTR`, `SLC4A4`, `MUC1`, `ANXA4`, `CLDN10/CLDN4` as appropriate), move `MSLN` to mesothelial if that axis is added, and either strengthen the centroacinar/ductal-subtype markers or redefine that refinement as a ductal state/subtype axis with a distinct marker program.
- **lymphatic_endothelium** — Most lymphatic markers are plausible, but `TFF3` is an epithelial/secretory marker and is not a canonical lymphatic endothelial marker. It risks contaminating the lymphatic score with ductal/epithelial expression.
  *Fix:* Replace `TFF3` with a lymphatic endothelial marker such as `CCL21`, `TFF3`-free `PROX1`/`LYVE1`/`PDPN`/`FLT4`/`MMRN1` support, or additional markers like `RELN`, `TBX1`, `ACKR4`, or `MRC1` depending on the intended lymphatic subtype.
- **immune_core and immune gated axes** — `immune_core` includes several myeloid/APC-biased genes (`LST1`, `HLA-DRA`, `HLA-DPA1`, `TYROBP`) and is used as the gate for T, NK, B, plasma, mast, macrophage, and dendritic axes. That gate can under-score lymphoid cells and over-emphasize myeloid/APC biology before the subtype refinements are applied.
  *Fix:* Revise `immune_core` toward pan-leukocyte markers such as `PTPRC`, `LAPTM5`, `LCP1`, `CORO1A`, `CD52`, `ARHGAP15`, `RAC2`, `FYB1`, while keeping myeloid/APC genes in the macrophage/DC modules.
- **exocrine_acinar_vs_ductal axis** — This axis name implies a contrast, but it is a plain acinar identity axis. The actual ductal contrast is handled separately by `ductal_vs_acinar`, creating an asymmetric and slightly confusing exocrine design.
  *Fix:* Rename the plain axis to `acinar_identity`, or make it a true acinar-vs-ductal contrast if the goal is a bidirectional exocrine epithelial dimension.

## Out of scope (FYI, not actionable for the Generator)
- The blueprint includes epsilon cells even though they are rare in adult pancreas; because age is unspecified, this is acceptable.
- More detailed vascular zonation, such as arterial, venous, capillary, and islet/acinar endothelial specialization, would be useful in a vascular-focused embedding but is not required here.

## Summary for the next round
Keep the overall structure, but add mesothelial coverage, consider a rare ductal/tuft-like epithelial refinement, clean up ductal/centroacinar marker separation, replace the non-lymphatic `TFF3` marker, and make the immune gate more pan-leukocyte so downstream immune subtype axes behave biologically.
