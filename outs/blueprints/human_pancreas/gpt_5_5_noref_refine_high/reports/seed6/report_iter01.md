# Discriminator report

## Scenario echo
The scenario is human pancreas with no specified age, sex, condition, or cell-type restriction, so the embedding should cover the normal whole-pancreas landscape: endocrine islet cells, acinar and ductal/centroacinar epithelium, vascular and lymphatic endothelium, stromal/perivascular cells, resident immune populations, neural/glial cells, and common technical or physiological states. There are no notes requiring a special topology.

## Rating: minor
The blueprint is biologically plausible overall and covers most major expected human pancreatic compartments. It still needs targeted refinement because it has a redundant acinar/ductal axis pair, lacks a clean acinar identity dimension, and merges several expected immune subtypes too coarsely for a whole-tissue pancreas embedding.

## Findings

### Minor issues
- **axes: exocrine_acinar_vs_ductal and ductal_vs_acinar** — These two contrast axes are exact reciprocal versions of the same acinar-ductal signal, so they add dimensionality without adding biological resolution. They also leave acinar identity represented only as a contrast against ductal cells rather than as its own strong exocrine compartment axis.
  *Fix:* Replace the reciprocal pair with non-redundant epithelial axes, for example a plain `acinar_identity` axis using `exocrine_acinar` and a plain or contrast `ductal_identity`/`ductal_vs_acinar` axis. Do not keep both exact inverse contrasts.
- **immune_t_nk_refinement and module:t_nk_cell** — T cells and NK cells are collapsed into one module containing both CD3/TRAC/IL7R and cytotoxic NK markers. For all-cell healthy pancreas, resident T cells and NK/cytotoxic lymphocytes are expected distinguishable immune populations; a combined axis will blur CD3+ T cells, NK cells, and cytotoxic T states.
  *Fix:* Split into separate T-cell and NK/cytotoxic modules/axes, or at least add one additional gated refinement that separates CD3D/CD3E/TRAC/IL7R from NKG7/GNLY/PRF1/GZMB/KLRD1.
- **immune_b_plasma_refinement and module:b_plasma_cell** — B cells and plasma cells are similarly collapsed into one module with MS4A1/CD79A/CD79B/BANK1 plus MZB1/JCHAIN/XBP1/IGKC. These are biologically related but distinct enough that a whole-pancreas embedding should separate naive/memory B cells from antibody-secreting plasma cells when present.
  *Fix:* Split `b_plasma_cell` into B-cell and plasma-cell modules and add separate immune-gated refinement axes.
- **coverage: adipocyte/lipid stromal population** — Whole human pancreas samples can include adipocytes or adipose-associated stromal cells, especially in adult tissue, but the blueprint has no adipocyte/lipid-storage module. This is not as central as endocrine/acinar/ductal coverage, but it is a realistic pancreas compartment missing from the all-cell scenario.
  *Fix:* Add an adipocyte/lipid module and plain or stromal-gated axis using markers such as `ADIPOQ`, `PLIN1`, `FABP4`, `LPL`, `CIDEA`, `LIPE`, and `DGAT2`, if the Generator wants fuller whole-tissue coverage.
- **module:centroacinar_ductal** — The module mixes canonical ductal/centroacinar genes with injury/activated ductal markers such as `SPP1`, `MMP7`, and `REG3A`. That may capture pancreatitis-like ductal activation more than a clean healthy centroacinar/terminal duct identity.
  *Fix:* Retune this module toward cleaner centroacinar/duct markers, such as `KRT7`, `KRT19`, `SOX9`, `CFTR`, `HNF1B`, `SLC4A4`, `ANXA4`, and `MUC1`; use a separate stress/activated ductal axis only if intentionally representing that state.

## Out of scope (FYI, not actionable for the Generator)
- The rare endocrine epsilon axis is acceptable because age is unspecified; it would be less important for a strictly adult-only healthy pancreas scenario.
- Generic cycling, interferon, ER stress, and immediate-early axes are reasonable as cross-lineage state dimensions for scRNA-seq data.

## Summary for the next round
Focus on reducing acinar/ductal redundancy, adding a clean acinar identity dimension, and improving immune resolution by separating T vs NK and B vs plasma cells. Consider adding adipocyte coverage and cleaning the centroacinar/ductal module so it is not dominated by injury-like ductal activation markers.
