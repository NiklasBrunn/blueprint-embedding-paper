# Blueprint variability — gpt_5_5_refine_high_human_pancreas  (10 runs)

**Scenario:** species=human, tissue=pancreas, celltype=None

**Stability composite:** **0.752**  _(equal-weighted; components in `cohort_summary.json`)_

## Blueprint size (mean ± SD [min, max], CV)

- **n_modules**: 26.70 ± 3.56  [21.00, 34.00], CV=0.133
- **n_axes**: 26.20 ± 3.65  [21.00, 34.00], CV=0.139
- **n_unique_genes**: 181.00 ± 22.49  [157.00, 234.00], CV=0.124
- **n_resources**: 11.90 ± 3.84  [6.00, 18.00], CV=0.323
- **n_dois_unique**: 5.90 ± 1.66  [4.00, 10.00], CV=0.282
- **gene_redundancy**: 1.13 ± 0.09  [1.02, 1.34], CV=0.075
- **mean_genes_per_module**: 7.68 ± 0.40  [7.04, 8.32], CV=0.052

## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)

- **genes**:        0.541  (95% CI [0.531, 0.551])
- **DOIs**:         0.472  (95% CI [0.439, 0.505])
- **module names**: 0.567  (95% CI [0.547, 0.588])

## Structural consistency (per-name mean agreement, ∈ [0,1])

- mean axis-kind agreement:        0.986
- mean axis positive agreement:    0.991
- mean axis gate agreement:        0.991
- mean axis negatives Jaccard:     0.740
- mean module role-set agreement:  0.856
- mean within-module gene Jaccard: 0.574

## Recurrence (universal / majority / singleton / unique-total; frac in ≥80 %; Gini)

- **modules**: 12 / 21 / 21 / 60; ≥80 %: 28.3 %; Gini: 0.451
- **axes**: 6 / 16 / 64 / 106; ≥80 %: 9.4 %; Gini: 0.461
- **genes**: 78 / 156 / 136 / 415; ≥80 %: 28.0 %; Gini: 0.443
- **resources**: 2 / 7 / 35 / 54; ≥80 %: 7.4 %; Gini: 0.439

## Saturation slope at N (absolute new-item gain when adding the N-th run)

- genes:    13.85
- modules:  2.20
- DOIs:     0.75

## Cohort coverage (totals across runs)

- unique modules:   **60**
- unique axes:      **106**
- unique genes:     **415**
- unique resources (DOI): **54**

## Resource-kind composition (mean fraction per run ± SD)

- paper:    0.287 ± 0.148
- atlas:    0.233 ± 0.135
- database: 0.479 ± 0.122

## Top hubs (most-reused modules across all runs × axes)

### As positive (axis identity carrier)
- `acinar_cell` — used in 14 (run × axis) slots
- `pericyte` — used in 12 (run × axis) slots
- `fibroblast_stellate` — used in 11 (run × axis) slots
- `endothelial` — used in 11 (run × axis) slots
- `beta_cell` — used in 11 (run × axis) slots
### As gate (axis activation switch)
- `immune_core` — used in 54 (run × axis) slots
- `endocrine_core` — used in 50 (run × axis) slots
- `endothelial` — used in 17 (run × axis) slots
- `fibroblast_stellate` — used in 12 (run × axis) slots
- `ductal_cell` — used in 5 (run × axis) slots
### As negative (axis contrast partner)
- `plasma_cell` — used in 47 (run × axis) slots
- `mast_cell` — used in 47 (run × axis) slots
- `macrophage` — used in 47 (run × axis) slots
- `dendritic_cell` — used in 44 (run × axis) slots
- `alpha_cell` — used in 41 (run × axis) slots

## Top resources (most-cited DOIs)

- Segerstolpe 2016 Cell Metabolism  (`10.1016/j.cmet.2016.08.020`, paper) — 10/10 runs (100 %)
- Tabula Sapiens Consortium 2022 Science  (`10.1126/science.abl4896`, atlas) — 10/10 runs (100 %)
- Baron 2016 Cell Systems  (`10.1016/j.cels.2016.08.011`, atlas) — 8/10 runs (80 %)
- Muraro 2016 Cell Systems  (`10.1016/j.cels.2016.09.002`, atlas) — 8/10 runs (80 %)
- KEGG hsa04972 Pancreatic secretion  (``, database) — 7/10 runs (70 %)

## Sanity

- runs with placeholders:                0.0 %
- mean orphaned modules / run:           0.00
- mean modules missing resources / run:  0.00
- mean unused resources / run:           0.00
- runs with clip override:               0.0 %
- canonicalization merges (concepts with ≥ 2 distinct raw spellings): modules=28, axes=37, genes=0
- post-canonical near-duplicate pairs (edit ≤ 2): 8; of which 0 never co-occur (higher-signal drift candidates).
- DOI-kind agreement:                    0.971

## Top axes by structural drift (low = unstable)

| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |
|---|---|---|---|---|
| `pericyte_in_pericyte` | 2 | 0.50 | 1.00 | 1.00 |
| `pericyte_vs_endothelial` | 2 | 0.50 | 1.00 | 0.00 |
| `alpha_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 1.00 |
| `beta_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 1.00 |
| `delta_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 1.00 |
| `endothelial` | 10 | 1.00 | 1.00 | 1.00 |
| `epsilon_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 1.00 |
| `pp_gamma_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 1.00 |
| `endocrine_core` | 9 | 1.00 | 1.00 | 1.00 |
| `schwann_cell` | 9 | 1.00 | 1.00 | 1.00 |
| `ductal_cell` | 8 | 1.00 | 1.00 | 1.00 |
| `fibroblast_stellate` | 8 | 1.00 | 1.00 | 1.00 |
| `acinar_cell` | 7 | 1.00 | 1.00 | 1.00 |
| `immune_core` | 7 | 1.00 | 1.00 | 1.00 |
| `macrophage_in_immune_core` | 7 | 1.00 | 1.00 | 0.64 |

## Least-stable modules (within-name gene Jaccard)

- `ductal_bicarbonate`: mean Jaccard 0.214 (seen in 2 runs)
- `epsilon_cell`: mean Jaccard 0.246 (seen in 10 runs)
- `cell_cycle`: mean Jaccard 0.300 (seen in 2 runs)
- `pp_gamma_cell`: mean Jaccard 0.309 (seen in 10 runs)
- `stellate_activated`: mean Jaccard 0.316 (seen in 3 runs)
- `acinar_endothelial`: mean Jaccard 0.388 (seen in 3 runs)
- `immune_core`: mean Jaccard 0.419 (seen in 7 runs)
- `t_nk_cell`: mean Jaccard 0.420 (seen in 4 runs)
- `stellate_quiescent`: mean Jaccard 0.455 (seen in 2 runs)
- `schwann_cell`: mean Jaccard 0.473 (seen in 10 runs)
