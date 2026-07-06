# Blueprint variability — gpt_5_5_noref_refine_high_human_pancreas  (10 runs)

**Scenario:** species=human, tissue=pancreas, celltype=None

**Stability composite:** **0.749**  _(equal-weighted; components in `cohort_summary.json`)_

## Blueprint size (mean ± SD [min, max], CV)

- **n_modules**: 37.10 ± 5.97  [28.00, 44.00], CV=0.161
- **n_axes**: 36.20 ± 5.14  [30.00, 45.00], CV=0.142
- **n_unique_genes**: 247.40 ± 34.59  [194.00, 296.00], CV=0.140
- **n_resources**: 0.00 ± 0.00  [0.00, 0.00], CV=nan
- **n_dois_unique**: 0.00 ± 0.00  [0.00, 0.00], CV=nan
- **gene_redundancy**: 1.22 ± 0.09  [1.12, 1.37], CV=0.076
- **mean_genes_per_module**: 8.20 ± 0.81  [7.03, 9.21], CV=0.099

## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)

- **genes**:        0.522  (95% CI [0.514, 0.529])
- **DOIs**:         0.000  (95% CI [0.000, 0.000])
- **module names**: 0.604  (95% CI [0.591, 0.617])

## Structural consistency (per-name mean agreement, ∈ [0,1])

- mean axis-kind agreement:        0.986
- mean axis positive agreement:    0.997
- mean axis gate agreement:        0.982
- mean axis negatives Jaccard:     0.672
- mean module role-set agreement:  0.845
- mean within-module gene Jaccard: 0.538

## Recurrence (universal / majority / singleton / unique-total; frac in ≥80 %; Gini)

- **modules**: 17 / 33 / 19 / 70; ≥80 %: 37.1 %; Gini: 0.403
- **axes**: 7 / 27 / 37 / 108; ≥80 %: 15.7 %; Gini: 0.433
- **genes**: 100 / 222 / 185 / 563; ≥80 %: 26.8 %; Gini: 0.434
- **resources**: 0 / 0 / 0 / 0; ≥80 %: 0.0 %; Gini: nan

## Saturation slope at N (absolute new-item gain when adding the N-th run)

- genes:    18.29
- modules:  1.91
- DOIs:     0.00

## Cohort coverage (totals across runs)

- unique modules:   **70**
- unique axes:      **108**
- unique genes:     **563**
- unique resources (DOI): **0**

## Resource-kind composition (mean fraction per run ± SD)

- paper:    0.000 ± 0.000
- atlas:    0.000 ± 0.000
- database: 0.000 ± 0.000

## Top hubs (most-reused modules across all runs × axes)

### As positive (axis identity carrier)
- `endothelial` — used in 14 (run × axis) slots
- `pericyte` — used in 13 (run × axis) slots
- `fibroblast_stellate` — used in 12 (run × axis) slots
- `endocrine_core` — used in 11 (run × axis) slots
- `acinar_cell` — used in 11 (run × axis) slots
### As gate (axis activation switch)
- `immune_core` — used in 59 (run × axis) slots
- `endocrine_core` — used in 55 (run × axis) slots
- `endothelial` — used in 27 (run × axis) slots
- `fibroblast_stellate` — used in 16 (run × axis) slots
- `ductal_cell` — used in 13 (run × axis) slots
### As negative (axis contrast partner)
- `macrophage` — used in 60 (run × axis) slots
- `alpha_cell` — used in 47 (run × axis) slots
- `mast_cell` — used in 47 (run × axis) slots
- `beta_cell` — used in 46 (run × axis) slots
- `dendritic_cell` — used in 46 (run × axis) slots

## Sanity

- runs with placeholders:                0.0 %
- mean orphaned modules / run:           0.00
- mean modules missing resources / run:  37.10
- mean unused resources / run:           0.00
- runs with clip override:               0.0 %
- canonicalization merges (concepts with ≥ 2 distinct raw spellings): modules=44, axes=67, genes=0
- post-canonical near-duplicate pairs (edit ≤ 2): 9; of which 0 never co-occur (higher-signal drift candidates).
- DOI-kind agreement:                    nan

## Top axes by structural drift (low = unstable)

| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |
|---|---|---|---|---|
| `ductal_cell_in_pancreatic_epithelial` | 3 | 0.67 | 1.00 | 1.00 |
| `endothelial` | 12 | 1.00 | 0.83 | 0.70 |
| `lymphatic_endothelial_in_endothelial` | 11 | 1.00 | 0.82 | 0.30 |
| `alpha_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 0.93 |
| `beta_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 0.93 |
| `delta_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 0.93 |
| `epsilon_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 0.91 |
| `pp_gamma_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 0.93 |
| `immune_core` | 9 | 1.00 | 1.00 | 1.00 |
| `adipocyte` | 8 | 1.00 | 1.00 | 1.00 |
| `cell_cycle` | 8 | 1.00 | 1.00 | 1.00 |
| `centroacinar_cell_in_ductal_cell` | 8 | 1.00 | 1.00 | 0.17 |
| `dendritic_cell_in_immune_core` | 8 | 1.00 | 1.00 | 0.46 |
| `fibroblast_stellate` | 8 | 1.00 | 1.00 | 1.00 |
| `macrophage_in_immune_core` | 8 | 1.00 | 1.00 | 0.51 |

## Least-stable modules (within-name gene Jaccard)

- `acinar_zymogen`: mean Jaccard 0.099 (seen in 3 runs)
- `stellate_quiescent`: mean Jaccard 0.197 (seen in 3 runs)
- `epsilon_cell`: mean Jaccard 0.237 (seen in 10 runs)
- `centroacinar_cell`: mean Jaccard 0.253 (seen in 10 runs)
- `beta_immature`: mean Jaccard 0.286 (seen in 2 runs)
- `pp_gamma_cell`: mean Jaccard 0.287 (seen in 10 runs)
- `matrix_remodeling`: mean Jaccard 0.312 (seen in 2 runs)
- `stellate_activated`: mean Jaccard 0.328 (seen in 4 runs)
- `ductal_mucinous`: mean Jaccard 0.333 (seen in 2 runs)
- `t_nk_cell`: mean Jaccard 0.348 (seen in 5 runs)
