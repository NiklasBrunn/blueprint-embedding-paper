# Blueprint variability — gpt_5_5_unrefined_high_human_pancreas  (10 runs)

**Scenario:** species=human, tissue=pancreas, celltype=None

**Stability composite:** **0.761**  _(equal-weighted; components in `cohort_summary.json`)_

## Blueprint size (mean ± SD [min, max], CV)

- **n_modules**: 20.30 ± 2.67  [16.00, 24.00], CV=0.131
- **n_axes**: 20.90 ± 2.69  [16.00, 24.00], CV=0.128
- **n_unique_genes**: 142.80 ± 15.82  [122.00, 164.00], CV=0.111
- **n_resources**: 10.00 ± 2.31  [6.00, 13.00], CV=0.231
- **n_dois_unique**: 5.40 ± 0.97  [4.00, 7.00], CV=0.179
- **gene_redundancy**: 1.08 ± 0.05  [1.02, 1.14], CV=0.042
- **mean_genes_per_module**: 7.65 ± 0.41  [6.95, 8.20], CV=0.054

## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)

- **genes**:        0.537  (95% CI [0.522, 0.551])
- **DOIs**:         0.492  (95% CI [0.453, 0.531])
- **module names**: 0.628  (95% CI [0.606, 0.654])

## Structural consistency (per-name mean agreement, ∈ [0,1])

- mean axis-kind agreement:        0.992
- mean axis positive agreement:    1.000
- mean axis gate agreement:        0.984
- mean axis negatives Jaccard:     0.784
- mean module role-set agreement:  0.834
- mean within-module gene Jaccard: 0.574

## Recurrence (universal / majority / singleton / unique-total; frac in ≥80 %; Gini)

- **modules**: 9 / 18 / 8 / 38; ≥80 %: 36.8 %; Gini: 0.389
- **axes**: 6 / 18 / 22 / 58; ≥80 %: 17.2 %; Gini: 0.455
- **genes**: 54 / 123 / 105 / 326; ≥80 %: 30.1 %; Gini: 0.439
- **resources**: 1 / 7 / 20 / 38; ≥80 %: 10.5 %; Gini: 0.452

## Saturation slope at N (absolute new-item gain when adding the N-th run)

- genes:    10.62
- modules:  0.68
- DOIs:     0.56

## Cohort coverage (totals across runs)

- unique modules:   **38**
- unique axes:      **58**
- unique genes:     **326**
- unique resources (DOI): **38**

## Resource-kind composition (mean fraction per run ± SD)

- paper:    0.290 ± 0.139
- atlas:    0.262 ± 0.139
- database: 0.447 ± 0.087

## Top hubs (most-reused modules across all runs × axes)

### As positive (axis identity carrier)
- `acinar_cell` — used in 15 (run × axis) slots
- `endothelial` — used in 12 (run × axis) slots
- `fibroblast_stellate` — used in 11 (run × axis) slots
- `beta_cell` — used in 11 (run × axis) slots
- `endocrine_core` — used in 10 (run × axis) slots
### As gate (axis activation switch)
- `endocrine_core` — used in 50 (run × axis) slots
- `immune_core` — used in 30 (run × axis) slots
- `endothelial` — used in 15 (run × axis) slots
- `fibroblast_stellate` — used in 7 (run × axis) slots
- `acinar_cell` — used in 4 (run × axis) slots
### As negative (axis contrast partner)
- `alpha_cell` — used in 41 (run × axis) slots
- `delta_cell` — used in 40 (run × axis) slots
- `pp_gamma_cell` — used in 40 (run × axis) slots
- `epsilon_cell` — used in 40 (run × axis) slots
- `beta_cell` — used in 40 (run × axis) slots

## Top resources (most-cited DOIs)

- Segerstolpe 2016 Cell Metabolism  (`10.1016/j.cmet.2016.08.020`, paper) — 10/10 runs (100 %)
- Tabula Sapiens Consortium 2022 Science  (`10.1126/science.abl4896`, atlas) — 9/10 runs (90 %)
- Baron 2016 Cell Systems  (`10.1016/j.cels.2016.08.011`, atlas) — 8/10 runs (80 %)
- Muraro 2016 Cell Systems  (`10.1016/j.cels.2016.09.002`, atlas) — 8/10 runs (80 %)
- KEGG hsa04972 Pancreatic secretion  (``, database) — 7/10 runs (70 %)

## Sanity

- runs with placeholders:                0.0 %
- mean orphaned modules / run:           0.00
- mean modules missing resources / run:  0.00
- mean unused resources / run:           0.00
- runs with clip override:               0.0 %
- canonicalization merges (concepts with ≥ 2 distinct raw spellings): modules=23, axes=32, genes=0
- post-canonical near-duplicate pairs (edit ≤ 2): 7; of which 3 never co-occur (higher-signal drift candidates).
- DOI-kind agreement:                    0.964

### Post-canonical near-duplicate names (review — possible typos)

- `b_plasma_cell` ↔ `plasma_cell` (canonical `bplasmacell` / `plasmacell`, edit=1, runs_a=7, runs_b=1)
- `nk_cell` ↔ `t_nk_cell` (canonical `nkcell` / `tnkcell`, edit=1, runs_a=1, runs_b=7)
- `t_cell` ↔ `t_nk_cell` (canonical `tcell` / `tnkcell`, edit=2, runs_a=3, runs_b=7)

## Top axes by structural drift (low = unstable)

| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |
|---|---|---|---|---|
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
| `acinar_cell_vs_ductal_cell` | 7 | 1.00 | 1.00 | 0.71 |
| `immune_core` | 7 | 1.00 | 1.00 | 1.00 |
| `macrophage_in_immune_core` | 7 | 1.00 | 1.00 | 0.51 |
| `b_plasma_cell_in_immune_core` | 5 | 1.00 | 1.00 | 0.87 |

## Least-stable modules (within-name gene Jaccard)

- `stellate_quiescent`: mean Jaccard 0.143 (seen in 2 runs)
- `ductal_bicarbonate`: mean Jaccard 0.214 (seen in 2 runs)
- `pp_gamma_cell`: mean Jaccard 0.297 (seen in 10 runs)
- `cell_cycle`: mean Jaccard 0.300 (seen in 2 runs)
- `epsilon_cell`: mean Jaccard 0.390 (seen in 10 runs)
- `stellate_activated`: mean Jaccard 0.400 (seen in 2 runs)
- `immune_core`: mean Jaccard 0.448 (seen in 7 runs)
- `t_cell`: mean Jaccard 0.487 (seen in 3 runs)
- `schwann_cell`: mean Jaccard 0.492 (seen in 9 runs)
- `fibroblast_stellate`: mean Jaccard 0.527 (seen in 11 runs)
