# Discriminator report

## Scenario echo
Human pancreas, with celltype, age, sex, and condition unspecified; this should cover the normal whole-pancreas cell-type and resident-state landscape without adding disease-, fetal-, or sex-specific programs. There are no additional notes to bias the topology.

## Rating: minor
The blueprint is broadly biologically plausible and covers the main pancreatic epithelial, endocrine, vascular, stromal, immune, and neural compartments. It needs targeted refinement for a missing adipocyte/adipose compartment and for a few gates/modules whose marker composition could under-score expected resident populations.

## Findings

### Minor issues
- **coverage: adipocyte/adipose compartment** — Whole human pancreas commonly includes intrapancreatic adipose/adipocytes, especially in adult tissue, but there is no adipocyte lipid-storage/metabolic module or axis. This is a coverage gap for an otherwise whole-tissue embedding.
  *Fix:* Add an adipocyte/adipose axis with canonical human adipocyte markers such as ADIPOQ, PLIN1, FABP4, LPL, LEP, CFD, CIDEA, PPARG, and/or DGAT2, grounded in a pancreas atlas or broad human cell atlas/resource.
- **axis-design: immune_leukocyte gate** — The leukocyte gate is biased toward myeloid/APC and inflammatory expression (LST1, TYROBP, CD74, HLA-DRA/HLA-DPA1, CXCL8). Multiplying all immune subtype axes by this gate can suppress T, NK, and some B cells, so the gated lymphocyte axes are biologically less robust than the labels imply.
  *Fix:* Refactor immune_leukocyte into a broader pan-leukocyte gate using markers such as PTPRC, CORO1A, LAPTM5, LCP1, CD52, IL2RG, and/or RAC2, and keep HLA/TYROBP/LST1/CXCL8 in myeloid/APC/inflammatory modules where appropriate.
- **axis-design: stromal_matrix gate for pericyte/smooth muscle** — Pericytes and vascular smooth muscle are gated by a fibroblast/ECM-heavy stromal_matrix module containing COL1A1/COL1A2/COL3A1/DCN/LUM/PDGFRA. This may under-score PDGFRB/RGS5/MCAM/NOTCH3/ACTA2/MYH11 vascular mural cells that do not strongly express fibroblast ECM or PDGFRA.
  *Fix:* Use a broader mesenchymal/mural-compatible gate, add a separate mesenchymal_core gate, or make the pericyte/smooth-muscle identity a plain axis while retaining stellate/myofibroblast refinements under the ECM/stellate gate.
- **coverage / axis-design: NK versus cytotoxic T resolution** — The cytotoxic_nk_cell module merges NK-cell markers with generic cytotoxic effector genes, while the T-cell module contains CD8A. This captures cytotoxic lymphocytes but does not cleanly separate NK cells from cytotoxic CD8 T cells, both of which are plausible resident/vascular-associated immune populations in human pancreas.
  *Fix:* Split the current signal into an NK-cell identity module/axis (for example KLRD1, NCR1, NKG7, GNLY, PRF1, KLRF1, FCGR3A) and a T-cell/cytotoxic T state or CD8 T refinement gated by T-cell markers (CD3D, CD3E, TRAC, CD2, CD8A/CD8B, GZMB/PRF1 as state markers).
- **marker-plausibility: epsilon_cell** — The epsilon-cell axis is based on GHRL alone. GHRL is canonical, but a one-gene module is fragile for real datasets and may become a near-dropout axis for a rare pancreatic endocrine population.
  *Fix:* Add additional literature-supported epsilon/ghrelin-cell markers if available for human pancreas, or explicitly keep the module as a rare single-marker axis only if the available high-quality literature does not support a broader marker set.

## Out of scope (FYI, not actionable for the Generator)
- I did not verify DOI/title matching or IF>10 compliance; that citation gate is outside this Discriminator scope.
- Some modules cite very broad resources or databases, but the report only flags biological coverage and marker/gating plausibility.

## Summary for the next round
Focus on adding adipocyte coverage, broadening the immune gate so it does not behave like a myeloid/APC gate, and revising the stromal/mural gating so pericytes and smooth muscle are not forced through a fibroblast ECM axis. If there is room, separate NK identity from cytotoxic T-cell state and strengthen the epsilon-cell marker set.
