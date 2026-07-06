# Blueprint variability — gpt_5_4_mini_unrefined_medium_human_pancreas  (10 runs)

**Scenario:** species=human, tissue=pancreas, celltype=None

**Stability composite:** **0.659**  _(equal-weighted; components in `cohort_summary.json`)_

## Blueprint size (mean ± SD [min, max], CV)

- **n_modules**: 14.80 ± 2.10  [12.00, 18.00], CV=0.142
- **n_axes**: 13.70 ± 2.36  [10.00, 18.00], CV=0.172
- **n_unique_genes**: 118.30 ± 19.66  [96.00, 150.00], CV=0.166
- **n_resources**: 7.20 ± 2.15  [4.00, 10.00], CV=0.299
- **n_dois_unique**: 3.40 ± 0.84  [2.00, 4.00], CV=0.248
- **gene_redundancy**: 1.18 ± 0.10  [1.07, 1.32], CV=0.084
- **mean_genes_per_module**: 9.40 ± 1.02  [8.00, 10.69], CV=0.108

## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)

- **genes**:        0.506  (95% CI [0.494, 0.518])
- **DOIs**:         0.327  (95% CI [0.281, 0.374])
- **module names**: 0.109  (95% CI [0.082, 0.136])

## Structural consistency (per-name mean agreement, ∈ [0,1])

- mean axis-kind agreement:        0.990
- mean axis positive agreement:    0.885
- mean axis gate agreement:        0.972
- mean axis negatives Jaccard:     0.672
- mean module role-set agreement:  0.960
- mean within-module gene Jaccard: 0.506

## Recurrence (universal / majority / singleton / unique-total; frac in ≥80 %; Gini)

- **modules**: 0 / 3 / 74 / 96; ≥80 %: 1.0 %; Gini: 0.297
- **axes**: 0 / 2 / 64 / 90; ≥80 %: 0.0 %; Gini: 0.273
- **genes**: 49 / 98 / 96 / 287; ≥80 %: 24.4 %; Gini: 0.447
- **resources**: 0 / 2 / 28 / 40; ≥80 %: 2.5 %; Gini: 0.365

## Saturation slope at N (absolute new-item gain when adding the N-th run)

- genes:    9.02
- modules:  7.41
- DOIs:     0.72

## Cohort coverage (totals across runs)

- unique modules:   **96**
- unique axes:      **90**
- unique genes:     **287**
- unique resources (DOI): **40**

## Resource-kind composition (mean fraction per run ± SD)

- paper:    0.345 ± 0.238
- atlas:    0.186 ± 0.101
- database: 0.470 ± 0.267

## Top hubs (most-reused modules across all runs × axes)

### As positive (axis identity carrier)
- `endocrine_core` — used in 9 (run × axis) slots
- `acinar_cell` — used in 7 (run × axis) slots
- `beta_program` — used in 4 (run × axis) slots
- `alpha_program` — used in 4 (run × axis) slots
- `delta_program` — used in 4 (run × axis) slots
### As gate (axis activation switch)
- `endocrine_core` — used in 34 (run × axis) slots
- `immune_core` — used in 6 (run × axis) slots
- `islet_core` — used in 5 (run × axis) slots
- `endocrine_gate` — used in 5 (run × axis) slots
- `exocrine_core` — used in 3 (run × axis) slots
### As negative (axis contrast partner)
- `alpha_program` — used in 14 (run × axis) slots
- `delta_program` — used in 14 (run × axis) slots
- `pp_program` — used in 14 (run × axis) slots
- `beta_program` — used in 14 (run × axis) slots
- `epsilon_program` — used in 8 (run × axis) slots

## Top resources (most-cited DOIs)

- Segerstolpe et al. 2016 Cell Metabolism  (`10.1016/j.cmet.2016.08.020`, paper) — 9/10 runs (90 %)
- The Tabula Sapiens Consortium et al. 2022 Science  (`10.1126/science.abl4896`, atlas) — 7/10 runs (70 %)
- Arda et al. 2016 Cell Metabolism  (`10.1016/j.cmet.2016.04.002`, paper) — 4/10 runs (40 %)
- GO:0006955 immune response  (``, database) — 4/10 runs (40 %)
- KEGG hsa04972 Pancreatic secretion  (``, database) — 4/10 runs (40 %)

## Sanity

- runs with placeholders:                0.0 %
- mean orphaned modules / run:           0.00
- mean modules missing resources / run:  0.00
- mean unused resources / run:           0.00
- runs with clip override:               0.0 %
- canonicalization merges (concepts with ≥ 2 distinct raw spellings): modules=3, axes=5, genes=0
- post-canonical near-duplicate pairs (edit ≤ 2): 13; of which 7 never co-occur (higher-signal drift candidates).
- DOI-kind agreement:                    1.000

### Post-canonical near-duplicate names (review — possible typos)

- `acinar_cell` ↔ `acinar_cells` (canonical `acinarcell` / `acinarcells`, edit=1, runs_a=5, runs_b=1)
- `alpha_cell` ↔ `alpha_cells` (canonical `alphacell` / `alphacells`, edit=1, runs_a=1, runs_b=1)
- `beta_cell` ↔ `beta_cells` (canonical `betacell` / `betacells`, edit=1, runs_a=1, runs_b=1)
- `b_plasma_cell` ↔ `b_plasma_cells` (canonical `bplasmacell` / `bplasmacells`, edit=1, runs_a=1, runs_b=1)
- `ductal_cell` ↔ `ductal_cells` (canonical `ductalcell` / `ductalcells`, edit=1, runs_a=5, runs_b=1)
- `endothelial_cell` ↔ `endothelial_cells` (canonical `endothelialcell` / `endothelialcells`, edit=1, runs_a=1, runs_b=1)
- `nk_cells` ↔ `t_nk_cell` (canonical `nkcells` / `tnkcell`, edit=2, runs_a=1, runs_b=1)

## Top axes by structural drift (low = unstable)

| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |
|---|---|---|---|---|
| `immune_core` | 4 | 0.25 | 1.00 | 1.00 |
| `endothelial` | 6 | 0.33 | 1.00 | 1.00 |
| `acinar_cell` | 3 | 0.33 | 0.67 | 0.33 |
| `ductal_cell` | 3 | 0.33 | 0.67 | 0.33 |
| `acinar` | 2 | 0.50 | 1.00 | 1.00 |
| `alpha` | 2 | 0.50 | 1.00 | 0.00 |
| `beta` | 2 | 0.50 | 1.00 | 0.00 |
| `cell_cycle` | 2 | 0.50 | 1.00 | 1.00 |
| `cycling` | 2 | 0.50 | 1.00 | 1.00 |
| `delta` | 2 | 0.50 | 1.00 | 0.00 |
| `ductal` | 2 | 0.50 | 1.00 | 1.00 |
| `fibroblast_stellate` | 2 | 0.50 | 1.00 | 1.00 |
| `pericyte` | 2 | 0.50 | 1.00 | 1.00 |
| `pericyte_in_pericyte` | 2 | 0.50 | 1.00 | 1.00 |
| `pp` | 2 | 0.50 | 1.00 | 0.00 |

## Least-stable modules (within-name gene Jaccard)

- `epsilon_program`: mean Jaccard 0.100 (seen in 2 runs)
- `pp_core`: mean Jaccard 0.231 (seen in 2 runs)
- `cycling_program`: mean Jaccard 0.333 (seen in 2 runs)
- `delta_core`: mean Jaccard 0.333 (seen in 2 runs)
- `pp_program`: mean Jaccard 0.356 (seen in 4 runs)
- `immune_program`: mean Jaccard 0.375 (seen in 2 runs)
- `fibroblast_core`: mean Jaccard 0.429 (seen in 2 runs)
- `exocrine_core`: mean Jaccard 0.436 (seen in 3 runs)
- `delta_program`: mean Jaccard 0.456 (seen in 4 runs)
- `alpha_program`: mean Jaccard 0.460 (seen in 4 runs)
