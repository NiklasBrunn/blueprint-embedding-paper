# Blueprint variability — gpt_5_4_mini_unrefined_high_human_pancreas  (10 runs)

**Scenario:** species=human, tissue=pancreas, celltype=None

**Stability composite:** **0.648**  _(equal-weighted; components in `cohort_summary.json`)_

## Blueprint size (mean ± SD [min, max], CV)

- **n_modules**: 15.60 ± 6.00  [0.00, 21.00], CV=0.385
- **n_axes**: 14.70 ± 5.77  [0.00, 20.00], CV=0.393
- **n_unique_genes**: 138.10 ± 54.85  [0.00, 207.00], CV=0.397
- **n_resources**: 6.80 ± 2.86  [0.00, 10.00], CV=0.421
- **n_dois_unique**: 3.50 ± 2.01  [0.00, 7.00], CV=0.575
- **gene_redundancy**: 1.20 ± 0.15  [1.02, 1.47], CV=0.127
- **mean_genes_per_module**: 9.55 ± 3.59  [0.00, 13.84], CV=0.376

## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)

- **genes**:        0.488  (95% CI [0.477, 0.499])
- **DOIs**:         0.314  (95% CI [0.277, 0.348])
- **module names**: 0.101  (95% CI [0.074, 0.130])

## Structural consistency (per-name mean agreement, ∈ [0,1])

- mean axis-kind agreement:        0.965
- mean axis positive agreement:    0.884
- mean axis gate agreement:        0.960
- mean axis negatives Jaccard:     0.519
- mean module role-set agreement:  0.926
- mean within-module gene Jaccard: 0.523

## Recurrence (universal / majority / singleton / unique-total; frac in ≥80 %; Gini)

- **modules**: 0 / 1 / 63 / 99; ≥80 %: 1.0 %; Gini: 0.271
- **axes**: 0 / 3 / 57 / 93; ≥80 %: 0.0 %; Gini: 0.265
- **genes**: 0 / 125 / 136 / 362; ≥80 %: 21.5 %; Gini: 0.431
- **resources**: 0 / 5 / 24 / 35; ≥80 %: 2.9 %; Gini: 0.381

## Saturation slope at N (absolute new-item gain when adding the N-th run)

- genes:    15.97
- modules:  7.01
- DOIs:     0.59

## Cohort coverage (totals across runs)

- unique modules:   **99**
- unique axes:      **93**
- unique genes:     **362**
- unique resources (DOI): **35**

## Resource-kind composition (mean fraction per run ± SD)

- paper:    0.257 ± 0.210
- atlas:    0.220 ± 0.197
- database: 0.423 ± 0.247

## Top hubs (most-reused modules across all runs × axes)

### As positive (axis identity carrier)
- `endocrine_core` — used in 11 (run × axis) slots
- `acinar_cell` — used in 4 (run × axis) slots
- `endothelial` — used in 4 (run × axis) slots
- `ductal_cell` — used in 3 (run × axis) slots
- `beta_cells` — used in 3 (run × axis) slots
### As gate (axis activation switch)
- `endocrine_core` — used in 33 (run × axis) slots
- `immune_gate` — used in 5 (run × axis) slots
- `exocrine_core` — used in 2 (run × axis) slots
- `vascular_gate` — used in 2 (run × axis) slots
- `mesenchymal_gate` — used in 2 (run × axis) slots
### As negative (axis contrast partner)
- `alpha_cells` — used in 12 (run × axis) slots
- `delta_cells` — used in 12 (run × axis) slots
- `pp_cells` — used in 12 (run × axis) slots
- `epsilon_cells` — used in 12 (run × axis) slots
- `beta_cells` — used in 12 (run × axis) slots

## Top resources (most-cited DOIs)

- Segerstolpe et al. 2016 Cell Metabolism  (`10.1016/j.cmet.2016.08.020`, atlas) — 8/10 runs (80 %)
- Muraro et al. 2016 Cell Systems  (`10.1016/j.cels.2016.09.002`, atlas) — 5/10 runs (50 %)
- Arda et al. 2016 Cell Metabolism  (`10.1016/j.cmet.2016.04.002`, paper) — 5/10 runs (50 %)
- The Tabula Sapiens: A multiple-organ, single-cell transcriptomic atlas of humans  (`10.1126/science.abl4896`, atlas) — 5/10 runs (50 %)
- MSigDB HALLMARK_G2M_CHECKPOINT  (``, database) — 5/10 runs (50 %)

## Sanity

- runs with placeholders:                0.0 %
- mean orphaned modules / run:           0.00
- mean modules missing resources / run:  0.00
- mean unused resources / run:           0.00
- runs with clip override:               0.0 %
- canonicalization merges (concepts with ≥ 2 distinct raw spellings): modules=4, axes=7, genes=0
- post-canonical near-duplicate pairs (edit ≤ 2): 15; of which 7 never co-occur (higher-signal drift candidates).
- DOI-kind agreement:                    0.970

### Post-canonical near-duplicate names (review — possible typos)

- `acinar_cell` ↔ `acinar_cells` (canonical `acinarcell` / `acinarcells`, edit=1, runs_a=4, runs_b=2)
- `alpha_cell` ↔ `alpha_cells` (canonical `alphacell` / `alphacells`, edit=1, runs_a=1, runs_b=3)
- `beta_cell` ↔ `beta_cells` (canonical `betacell` / `betacells`, edit=1, runs_a=1, runs_b=3)
- `b_plasma` ↔ `plasma` (canonical `bplasma` / `plasma`, edit=1, runs_a=1, runs_b=1)
- `ductal_cell` ↔ `ductal_cells` (canonical `ductalcell` / `ductalcells`, edit=1, runs_a=4, runs_b=2)
- `pericyte_smc` ↔ `pericyte_vsmc` (canonical `pericytesmc` / `pericytevsmc`, edit=1, runs_a=1, runs_b=1)
- `schwann_cell` ↔ `schwann_cells` (canonical `schwanncell` / `schwanncells`, edit=1, runs_a=2, runs_b=1)

## Top axes by structural drift (low = unstable)

| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |
|---|---|---|---|---|
| `beta_vs_alpha` | 3 | 0.33 | 1.00 | 0.00 |
| `acinar` | 2 | 0.50 | 0.50 | 0.00 |
| `acinar_vs_ductal` | 2 | 0.50 | 1.00 | 0.00 |
| `alpha` | 2 | 0.50 | 1.00 | 0.00 |
| `alpha_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `b_cell` | 2 | 0.50 | 0.50 | 0.00 |
| `beta` | 2 | 0.50 | 1.00 | 0.00 |
| `beta_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `beta_maturity` | 2 | 0.50 | 0.50 | 1.00 |
| `cycling` | 2 | 0.50 | 1.00 | 1.00 |
| `delta` | 2 | 0.50 | 1.00 | 0.00 |
| `delta_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `ductal` | 2 | 0.50 | 0.50 | 0.00 |
| `endothelial_program` | 2 | 0.50 | 1.00 | 1.00 |
| `pp` | 2 | 0.50 | 1.00 | 0.00 |

## Least-stable modules (within-name gene Jaccard)

- `pp_program`: mean Jaccard 0.143 (seen in 2 runs)
- `delta_core`: mean Jaccard 0.235 (seen in 2 runs)
- `delta_cells`: mean Jaccard 0.255 (seen in 3 runs)
- `pp_core`: mean Jaccard 0.267 (seen in 2 runs)
- `pp_cells`: mean Jaccard 0.270 (seen in 3 runs)
- `epsilon_cells`: mean Jaccard 0.300 (seen in 3 runs)
- `myeloid_program`: mean Jaccard 0.333 (seen in 2 runs)
- `cell_cycle`: mean Jaccard 0.334 (seen in 3 runs)
- `myeloid`: mean Jaccard 0.375 (seen in 2 runs)
- `delta_program`: mean Jaccard 0.417 (seen in 2 runs)
