# Blueprint variability — gpt_5_5_noref_unrefined_high_human_pancreas  (10 runs)

**Scenario:** species=human, tissue=pancreas, celltype=None

**Stability composite:** **0.744**  _(equal-weighted; components in `cohort_summary.json`)_

## Blueprint size (mean ± SD [min, max], CV)

- **n_modules**: 29.90 ± 5.90  [21.00, 37.00], CV=0.197
- **n_axes**: 30.30 ± 4.83  [23.00, 37.00], CV=0.159
- **n_unique_genes**: 210.90 ± 36.96  [153.00, 266.00], CV=0.175
- **n_resources**: 0.00 ± 0.00  [0.00, 0.00], CV=nan
- **n_dois_unique**: 0.00 ± 0.00  [0.00, 0.00], CV=nan
- **gene_redundancy**: 1.17 ± 0.06  [1.07, 1.24], CV=0.048
- **mean_genes_per_module**: 8.29 ± 0.86  [7.09, 9.35], CV=0.104

## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)

- **genes**:        0.531  (95% CI [0.522, 0.540])
- **DOIs**:         0.000  (95% CI [0.000, 0.000])
- **module names**: 0.544  (95% CI [0.526, 0.564])

## Structural consistency (per-name mean agreement, ∈ [0,1])

- mean axis-kind agreement:        0.975
- mean axis positive agreement:    0.997
- mean axis gate agreement:        0.964
- mean axis negatives Jaccard:     0.651
- mean module role-set agreement:  0.837
- mean within-module gene Jaccard: 0.579

## Recurrence (universal / majority / singleton / unique-total; frac in ≥80 %; Gini)

- **modules**: 12 / 25 / 15 / 62; ≥80 %: 30.6 %; Gini: 0.406
- **axes**: 6 / 24 / 33 / 96; ≥80 %: 12.5 %; Gini: 0.430
- **genes**: 92 / 189 / 131 / 462; ≥80 %: 27.9 %; Gini: 0.421
- **resources**: 0 / 0 / 0 / 0; ≥80 %: 0.0 %; Gini: nan

## Saturation slope at N (absolute new-item gain when adding the N-th run)

- genes:    13.00
- modules:  1.52
- DOIs:     0.00

## Cohort coverage (totals across runs)

- unique modules:   **62**
- unique axes:      **96**
- unique genes:     **462**
- unique resources (DOI): **0**

## Resource-kind composition (mean fraction per run ± SD)

- paper:    0.000 ± 0.000
- atlas:    0.000 ± 0.000
- database: 0.000 ± 0.000

## Top hubs (most-reused modules across all runs × axes)

### As positive (axis identity carrier)
- `endothelial` — used in 13 (run × axis) slots
- `fibroblast_stellate` — used in 12 (run × axis) slots
- `endocrine_core` — used in 11 (run × axis) slots
- `acinar_cell` — used in 11 (run × axis) slots
- `pericyte` — used in 11 (run × axis) slots
### As gate (axis activation switch)
- `endocrine_core` — used in 57 (run × axis) slots
- `immune_core` — used in 44 (run × axis) slots
- `endothelial` — used in 19 (run × axis) slots
- `pancreatic_epithelial` — used in 15 (run × axis) slots
- `fibroblast_stellate` — used in 13 (run × axis) slots
### As negative (axis contrast partner)
- `alpha_cell` — used in 49 (run × axis) slots
- `beta_cell` — used in 46 (run × axis) slots
- `delta_cell` — used in 45 (run × axis) slots
- `pp_gamma_cell` — used in 43 (run × axis) slots
- `epsilon_cell` — used in 41 (run × axis) slots

## Sanity

- runs with placeholders:                0.0 %
- mean orphaned modules / run:           0.00
- mean modules missing resources / run:  29.90
- mean unused resources / run:           0.00
- runs with clip override:               0.0 %
- canonicalization merges (concepts with ≥ 2 distinct raw spellings): modules=38, axes=58, genes=0
- post-canonical near-duplicate pairs (edit ≤ 2): 7; of which 3 never co-occur (higher-signal drift candidates).
- DOI-kind agreement:                    nan

### Post-canonical near-duplicate names (review — possible typos)

- `b_plasma_cell` ↔ `plasma_cell` (canonical `bplasmacell` / `plasmacell`, edit=1, runs_a=7, runs_b=3)
- `nk_cell` ↔ `t_nk_cell` (canonical `nkcell` / `tnkcell`, edit=1, runs_a=4, runs_b=6)
- `t_cell` ↔ `t_nk_cell` (canonical `tcell` / `tnkcell`, edit=2, runs_a=4, runs_b=6)

## Top axes by structural drift (low = unstable)

| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |
|---|---|---|---|---|
| `ductal_cell_in_pancreatic_epithelial` | 3 | 0.67 | 1.00 | 0.78 |
| `endothelial` | 11 | 1.00 | 0.91 | 0.82 |
| `alpha_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 0.96 |
| `beta_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 0.96 |
| `delta_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 0.96 |
| `epsilon_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 0.96 |
| `pp_gamma_cell_in_endocrine_core` | 10 | 1.00 | 1.00 | 0.96 |
| `lymphatic_endothelial_in_endothelial` | 9 | 1.00 | 0.89 | 0.12 |
| `fibroblast_stellate` | 8 | 1.00 | 1.00 | 1.00 |
| `immune_core` | 8 | 1.00 | 1.00 | 1.00 |
| `macrophage_in_immune_core` | 8 | 1.00 | 1.00 | 0.49 |
| `mast_cell_in_immune_core` | 8 | 1.00 | 1.00 | 0.48 |
| `schwann_cell` | 8 | 1.00 | 1.00 | 1.00 |
| `cell_cycle` | 7 | 1.00 | 1.00 | 1.00 |
| `centroacinar_cell_in_ductal_cell` | 7 | 1.00 | 1.00 | 0.22 |

## Least-stable modules (within-name gene Jaccard)

- `centroacinar_cell`: mean Jaccard 0.233 (seen in 9 runs)
- `pp_gamma_cell`: mean Jaccard 0.266 (seen in 10 runs)
- `matrix_remodeling`: mean Jaccard 0.267 (seen in 2 runs)
- `epsilon_cell`: mean Jaccard 0.321 (seen in 10 runs)
- `acinar_regenerative`: mean Jaccard 0.333 (seen in 2 runs)
- `immune_core`: mean Jaccard 0.342 (seen in 8 runs)
- `myeloid_core`: mean Jaccard 0.364 (seen in 2 runs)
- `monocyte`: mean Jaccard 0.400 (seen in 2 runs)
- `smooth_muscle`: mean Jaccard 0.400 (seen in 2 runs)
- `neuron`: mean Jaccard 0.413 (seen in 4 runs)
