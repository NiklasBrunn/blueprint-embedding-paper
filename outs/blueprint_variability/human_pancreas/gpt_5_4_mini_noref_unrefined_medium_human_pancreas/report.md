# Blueprint variability — gpt_5_4_mini_noref_unrefined_medium_human_pancreas  (10 runs)

**Scenario:** species=human, tissue=pancreas, celltype=None

**Stability composite:** **0.687**  _(equal-weighted; components in `cohort_summary.json`)_

## Blueprint size (mean ± SD [min, max], CV)

- **n_modules**: 18.70 ± 3.65  [15.00, 27.00], CV=0.195
- **n_axes**: 17.80 ± 4.66  [10.00, 27.00], CV=0.262
- **n_unique_genes**: 161.90 ± 29.40  [122.00, 211.00], CV=0.182
- **n_resources**: 0.00 ± 0.00  [0.00, 0.00], CV=nan
- **n_dois_unique**: 0.00 ± 0.00  [0.00, 0.00], CV=nan
- **gene_redundancy**: 1.12 ± 0.05  [1.05, 1.20], CV=0.048
- **mean_genes_per_module**: 9.74 ± 1.10  [7.96, 11.61], CV=0.113

## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)

- **genes**:        0.517  (95% CI [0.506, 0.529])
- **DOIs**:         0.000  (95% CI [0.000, 0.000])
- **module names**: 0.151  (95% CI [0.124, 0.177])

## Structural consistency (per-name mean agreement, ∈ [0,1])

- mean axis-kind agreement:        0.990
- mean axis positive agreement:    0.938
- mean axis gate agreement:        0.966
- mean axis negatives Jaccard:     0.743
- mean module role-set agreement:  0.972
- mean within-module gene Jaccard: 0.554

## Recurrence (universal / majority / singleton / unique-total; frac in ≥80 %; Gini)

- **modules**: 0 / 7 / 75 / 106; ≥80 %: 1.9 %; Gini: 0.353
- **axes**: 0 / 1 / 87 / 124; ≥80 %: 0.8 %; Gini: 0.236
- **genes**: 60 / 150 / 131 / 370; ≥80 %: 26.2 %; Gini: 0.435
- **resources**: 0 / 0 / 0 / 0; ≥80 %: 0.0 %; Gini: nan

## Saturation slope at N (absolute new-item gain when adding the N-th run)

- genes:    14.25
- modules:  7.22
- DOIs:     0.00

## Cohort coverage (totals across runs)

- unique modules:   **106**
- unique axes:      **124**
- unique genes:     **370**
- unique resources (DOI): **0**

## Resource-kind composition (mean fraction per run ± SD)

- paper:    0.000 ± 0.000
- atlas:    0.000 ± 0.000
- database: 0.000 ± 0.000

## Top hubs (most-reused modules across all runs × axes)

### As positive (axis identity carrier)
- `endocrine_core` — used in 8 (run × axis) slots
- `endothelial` — used in 8 (run × axis) slots
- `acinar_cell` — used in 7 (run × axis) slots
- `schwann_cell` — used in 7 (run × axis) slots
- `fibroblast_stellate` — used in 6 (run × axis) slots
### As gate (axis activation switch)
- `endocrine_core` — used in 38 (run × axis) slots
- `immune_core` — used in 22 (run × axis) slots
- `endocrine_general` — used in 5 (run × axis) slots
- `endocrine` — used in 5 (run × axis) slots
- `immune` — used in 5 (run × axis) slots
### As negative (axis contrast partner)
- `alpha` — used in 16 (run × axis) slots
- `delta` — used in 16 (run × axis) slots
- `epsilon` — used in 16 (run × axis) slots
- `beta` — used in 16 (run × axis) slots
- `pp` — used in 12 (run × axis) slots

## Sanity

- runs with placeholders:                0.0 %
- mean orphaned modules / run:           0.00
- mean modules missing resources / run:  18.70
- mean unused resources / run:           0.00
- runs with clip override:               0.0 %
- canonicalization merges (concepts with ≥ 2 distinct raw spellings): modules=8, axes=8, genes=0
- post-canonical near-duplicate pairs (edit ≤ 2): 25; of which 9 never co-occur (higher-signal drift candidates).
- DOI-kind agreement:                    nan

### Post-canonical near-duplicate names (review — possible typos)

- `acinar_cell` ↔ `acinar_cells` (canonical `acinarcell` / `acinarcells`, edit=1, runs_a=7, runs_b=2)
- `b_cell` ↔ `b_cells` (canonical `bcell` / `bcells`, edit=1, runs_a=2, runs_b=1)
- `b_plasma` ↔ `plasma` (canonical `bplasma` / `plasma`, edit=1, runs_a=1, runs_b=3)
- `ductal_cell` ↔ `ductal_cells` (canonical `ductalcell` / `ductalcells`, edit=1, runs_a=7, runs_b=2)
- `nk_cell` ↔ `nk_cells` (canonical `nkcell` / `nkcells`, edit=1, runs_a=1, runs_b=1)
- `t_cell` ↔ `t_cells` (canonical `tcell` / `tcells`, edit=1, runs_a=1, runs_b=1)
- `b_cell` ↔ `t_cells` (canonical `bcell` / `tcells`, edit=2, runs_a=2, runs_b=1)
- `b_cells` ↔ `t_cell` (canonical `bcells` / `tcell`, edit=2, runs_a=1, runs_b=1)
- `endocrine` ↔ `endocrine_pp` (canonical `endocrine` / `endocrinepp`, edit=2, runs_a=1, runs_b=1)

## Top axes by structural drift (low = unstable)

| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |
|---|---|---|---|---|
| `acinar_program` | 2 | 0.50 | 1.00 | 1.00 |
| `alpha_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `beta_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `delta_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `ductal_program` | 2 | 0.50 | 1.00 | 1.00 |
| `endothelial_program` | 2 | 0.50 | 1.00 | 1.00 |
| `epsilon_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `fibroblast_program` | 2 | 0.50 | 1.00 | 1.00 |
| `pp_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `schwann_identity` | 2 | 0.50 | 1.00 | 1.00 |
| `acinar` | 3 | 0.67 | 1.00 | 1.00 |
| `acinar_cell` | 3 | 0.67 | 1.00 | 1.00 |
| `cell_cycle` | 3 | 0.67 | 1.00 | 1.00 |
| `ductal` | 3 | 0.67 | 1.00 | 1.00 |
| `ductal_cell` | 3 | 0.67 | 1.00 | 1.00 |

## Least-stable modules (within-name gene Jaccard)

- `pp`: mean Jaccard 0.178 (seen in 3 runs)
- `immune_myeloid`: mean Jaccard 0.250 (seen in 2 runs)
- `epsilon_cells`: mean Jaccard 0.286 (seen in 2 runs)
- `epsilon`: mean Jaccard 0.289 (seen in 4 runs)
- `delta`: mean Jaccard 0.297 (seen in 4 runs)
- `delta_cells`: mean Jaccard 0.333 (seen in 2 runs)
- `fibro_stellate`: mean Jaccard 0.353 (seen in 2 runs)
- `b_cell`: mean Jaccard 0.385 (seen in 2 runs)
- `fibroblast_stellate`: mean Jaccard 0.436 (seen in 6 runs)
- `pp_cells`: mean Jaccard 0.444 (seen in 2 runs)
