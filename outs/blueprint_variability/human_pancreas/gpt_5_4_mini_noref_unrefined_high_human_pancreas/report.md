# Blueprint variability — gpt_5_4_mini_noref_unrefined_high_human_pancreas  (9 runs)

**Scenario:** species=human, tissue=pancreas, celltype=None

**Stability composite:** **0.654**  _(equal-weighted; components in `cohort_summary.json`)_

## Blueprint size (mean ± SD [min, max], CV)

- **n_modules**: 22.00 ± 4.47  [16.00, 29.00], CV=0.203
- **n_axes**: 21.89 ± 4.83  [16.00, 29.00], CV=0.221
- **n_unique_genes**: 183.00 ± 31.06  [154.00, 235.00], CV=0.170
- **n_resources**: 0.00 ± 0.00  [0.00, 0.00], CV=nan
- **n_dois_unique**: 0.00 ± 0.00  [0.00, 0.00], CV=nan
- **gene_redundancy**: 1.17 ± 0.10  [1.05, 1.30], CV=0.083
- **mean_genes_per_module**: 9.83 ± 1.10  [7.86, 11.77], CV=0.112

## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)

- **genes**:        0.505  (95% CI [0.495, 0.515])
- **DOIs**:         0.000  (95% CI [0.000, 0.000])
- **module names**: 0.075  (95% CI [0.059, 0.092])

## Structural consistency (per-name mean agreement, ∈ [0,1])

- mean axis-kind agreement:        0.972
- mean axis positive agreement:    0.922
- mean axis gate agreement:        0.970
- mean axis negatives Jaccard:     0.331
- mean module role-set agreement:  0.911
- mean within-module gene Jaccard: 0.537

## Recurrence (universal / majority / singleton / unique-total; frac in ≥80 %; Gini)

- **modules**: 0 / 4 / 94 / 134; ≥80 %: 0.0 %; Gini: 0.246
- **axes**: 0 / 4 / 123 / 150; ≥80 %: 0.0 %; Gini: 0.209
- **genes**: 73 / 151 / 159 / 424; ≥80 %: 23.1 %; Gini: 0.433
- **resources**: 0 / 0 / 0 / 0; ≥80 %: 0.0 %; Gini: nan

## Saturation slope at N (absolute new-item gain when adding the N-th run)

- genes:    17.42
- modules:  10.14
- DOIs:     0.00

## Cohort coverage (totals across runs)

- unique modules:   **134**
- unique axes:      **150**
- unique genes:     **424**
- unique resources (DOI): **0**

## Resource-kind composition (mean fraction per run ± SD)

- paper:    0.000 ± 0.000
- atlas:    0.000 ± 0.000
- database: 0.000 ± 0.000

## Top hubs (most-reused modules across all runs × axes)

### As positive (axis identity carrier)
- `endocrine_core` — used in 7 (run × axis) slots
- `fibroblast_stellate` — used in 6 (run × axis) slots
- `immune_core` — used in 5 (run × axis) slots
- `lymphatic_endothelial` — used in 5 (run × axis) slots
- `endothelial` — used in 4 (run × axis) slots
### As gate (axis activation switch)
- `endocrine_core` — used in 33 (run × axis) slots
- `immune_core` — used in 33 (run × axis) slots
- `endocrine_compartment` — used in 5 (run × axis) slots
- `fibroblast_stellate` — used in 4 (run × axis) slots
- `vascular_core` — used in 4 (run × axis) slots
### As negative (axis contrast partner)
- `beta_cells` — used in 10 (run × axis) slots
- `alpha_cells` — used in 10 (run × axis) slots
- `delta_cells` — used in 8 (run × axis) slots
- `pp_cells` — used in 8 (run × axis) slots
- `epsilon_cells` — used in 8 (run × axis) slots

## Sanity

- runs with placeholders:                0.0 %
- mean orphaned modules / run:           0.00
- mean modules missing resources / run:  22.00
- mean unused resources / run:           0.00
- runs with clip override:               0.0 %
- canonicalization merges (concepts with ≥ 2 distinct raw spellings): modules=7, axes=7, genes=0
- post-canonical near-duplicate pairs (edit ≤ 2): 31; of which 19 never co-occur (higher-signal drift candidates).
- DOI-kind agreement:                    nan

### Post-canonical near-duplicate names (review — possible typos)

- `acinar_cell` ↔ `acinar_cells` (canonical `acinarcell` / `acinarcells`, edit=1, runs_a=3, runs_b=2)
- `alpha_cell` ↔ `alpha_cells` (canonical `alphacell` / `alphacells`, edit=1, runs_a=1, runs_b=3)
- `b_cell` ↔ `b_cells` (canonical `bcell` / `bcells`, edit=1, runs_a=1, runs_b=1)
- `beta_cell` ↔ `beta_cells` (canonical `betacell` / `betacells`, edit=1, runs_a=1, runs_b=3)
- `bplasma_program` ↔ `plasma_program` (canonical `bplasmaprogram` / `plasmaprogram`, edit=1, runs_a=1, runs_b=2)
- `delta_cell` ↔ `delta_cells` (canonical `deltacell` / `deltacells`, edit=1, runs_a=1, runs_b=2)
- `ductal_cell` ↔ `ductal_cells` (canonical `ductalcell` / `ductalcells`, edit=1, runs_a=3, runs_b=2)
- `epsilon_cell` ↔ `epsilon_cells` (canonical `epsiloncell` / `epsiloncells`, edit=1, runs_a=1, runs_b=2)
- `lymphoid_cell` ↔ `lymphoid_cells` (canonical `lymphoidcell` / `lymphoidcells`, edit=1, runs_a=1, runs_b=1)
- `mast_cell` ↔ `mast_cells` (canonical `mastcell` / `mastcells`, edit=1, runs_a=2, runs_b=2)
- `myeloid_cell` ↔ `myeloid_cells` (canonical `myeloidcell` / `myeloidcells`, edit=1, runs_a=1, runs_b=2)
- `nk_program` ↔ `tnk_program` (canonical `nkprogram` / `tnkprogram`, edit=1, runs_a=2, runs_b=1)
- `plasma_cell` ↔ `plasma_cells` (canonical `plasmacell` / `plasmacells`, edit=1, runs_a=1, runs_b=1)
- `schwann_cell` ↔ `schwann_cells` (canonical `schwanncell` / `schwanncells`, edit=1, runs_a=4, runs_b=1)
- `b_cell` ↔ `t_cells` (canonical `bcell` / `tcells`, edit=2, runs_a=1, runs_b=1)

## Top axes by structural drift (low = unstable)

| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |
|---|---|---|---|---|
| `endothelial` | 4 | 0.25 | 1.00 | 1.00 |
| `acinar_vs_ductal` | 3 | 0.33 | 1.00 | 0.00 |
| `mast_cell` | 3 | 0.33 | 0.67 | 0.00 |
| `acinar_cell` | 5 | 0.40 | 0.60 | 0.30 |
| `ductal_cell` | 5 | 0.40 | 0.60 | 0.30 |
| `alpha_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `alpha_program` | 2 | 0.50 | 0.50 | 0.00 |
| `b_cell_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `beta_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `beta_program` | 2 | 0.50 | 0.50 | 0.00 |
| `delta_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `delta_program` | 2 | 0.50 | 0.50 | 0.00 |
| `dendritic_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `epsilon_identity` | 2 | 0.50 | 1.00 | 0.00 |
| `epsilon_program` | 2 | 0.50 | 0.50 | 0.00 |

## Least-stable modules (within-name gene Jaccard)

- `epsilon_program`: mean Jaccard 0.200 (seen in 2 runs)
- `epsilon_cells`: mean Jaccard 0.222 (seen in 2 runs)
- `pp_program`: mean Jaccard 0.273 (seen in 2 runs)
- `fibro_stellate`: mean Jaccard 0.312 (seen in 2 runs)
- `delta_cells`: mean Jaccard 0.333 (seen in 2 runs)
- `dendritic`: mean Jaccard 0.333 (seen in 2 runs)
- `pp_cells`: mean Jaccard 0.333 (seen in 2 runs)
- `exocrine_core`: mean Jaccard 0.368 (seen in 2 runs)
- `cell_cycle`: mean Jaccard 0.412 (seen in 2 runs)
- `acinar_cell`: mean Jaccard 0.443 (seen in 3 runs)
