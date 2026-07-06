# Blueprint variability — gpt_5_5_noref_unrefined_high_mouse_cortex  (10 runs)

**Scenario:** species=mouse, tissue=isocortex, celltype=None

**Stability composite:** **0.745**  _(equal-weighted; components in `cohort_summary.json`)_

## Blueprint size (mean ± SD [min, max], CV)

- **n_modules**: 30.10 ± 3.87  [22.00, 36.00], CV=0.129
- **n_axes**: 28.30 ± 4.03  [21.00, 35.00], CV=0.142
- **n_unique_genes**: 187.80 ± 28.36  [122.00, 221.00], CV=0.151
- **n_resources**: 0.00 ± 0.00  [0.00, 0.00], CV=nan
- **n_dois_unique**: 0.00 ± 0.00  [0.00, 0.00], CV=nan
- **gene_redundancy**: 1.23 ± 0.07  [1.11, 1.32], CV=0.061
- **mean_genes_per_module**: 7.65 ± 0.69  [6.14, 8.50], CV=0.090

## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)

- **genes**:        0.538  (95% CI [0.530, 0.546])
- **DOIs**:         0.000  (95% CI [0.000, 0.000])
- **module names**: 0.728  (95% CI [0.714, 0.742])

## Structural consistency (per-name mean agreement, ∈ [0,1])

- mean axis-kind agreement:        0.897
- mean axis positive agreement:    0.991
- mean axis gate agreement:        0.855
- mean axis negatives Jaccard:     0.665
- mean module role-set agreement:  0.844
- mean within-module gene Jaccard: 0.475

## Recurrence (universal / majority / singleton / unique-total; frac in ≥80 %; Gini)

- **modules**: 19 / 28 / 5 / 43; ≥80 %: 55.8 %; Gini: 0.285
- **axes**: 15 / 27 / 10 / 49; ≥80 %: 42.9 %; Gini: 0.372
- **genes**: 85 / 171 / 132 / 415; ≥80 %: 27.2 %; Gini: 0.430
- **resources**: 0 / 0 / 0 / 0; ≥80 %: 0.0 %; Gini: nan

## Saturation slope at N (absolute new-item gain when adding the N-th run)

- genes:    13.57
- modules:  0.52
- DOIs:     0.00

## Cohort coverage (totals across runs)

- unique modules:   **43**
- unique axes:      **49**
- unique genes:     **415**
- unique resources (DOI): **0**

## Resource-kind composition (mean fraction per run ± SD)

- paper:    0.000 ± 0.000
- atlas:    0.000 ± 0.000
- database: 0.000 ± 0.000

## Top hubs (most-reused modules across all runs × axes)

### As positive (axis identity carrier)
- `neuronal_activity_ieg` — used in 14 (run × axis) slots
- `pvalb` — used in 13 (run × axis) slots
- `sncg` — used in 12 (run × axis) slots
- `astrocyte` — used in 12 (run × axis) slots
- `mature_oligodendrocyte` — used in 11 (run × axis) slots
### As gate (axis activation switch)
- `glutamatergic_neuron` — used in 69 (run × axis) slots
- `gabaergic_neuron` — used in 60 (run × axis) slots
- `non_neuronal` — used in 41 (run × axis) slots
- `endothelial` — used in 15 (run × axis) slots
- `pan_neuronal` — used in 12 (run × axis) slots
### As negative (axis contrast partner)
- `sncg` — used in 60 (run × axis) slots
- `pvalb` — used in 58 (run × axis) slots
- `l4_it` — used in 57 (run × axis) slots
- `l5_et` — used in 57 (run × axis) slots
- `l6_ct` — used in 57 (run × axis) slots

## Sanity

- runs with placeholders:                0.0 %
- mean orphaned modules / run:           0.10
- mean modules missing resources / run:  30.10
- mean unused resources / run:           0.00
- runs with clip override:               0.0 %
- canonicalization merges (concepts with ≥ 2 distinct raw spellings): modules=35, axes=38, genes=0
- post-canonical near-duplicate pairs (edit ≤ 2): 16; of which 0 never co-occur (higher-signal drift candidates).
- DOI-kind agreement:                    nan

## Top axes by structural drift (low = unstable)

| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |
|---|---|---|---|---|
| `oligodendrocyte_maturation` | 4 | 0.75 | 0.50 | 0.42 |
| `opc_cell_cycle_in_opc` | 5 | 0.80 | 0.80 | 1.00 |
| `neuronal_activity_ieg` | 14 | 1.00 | 0.71 | 1.00 |
| `pvalb` | 13 | 1.00 | 1.00 | 0.84 |
| `astrocyte` | 12 | 1.00 | 0.50 | 0.48 |
| `sncg` | 12 | 1.00 | 1.00 | 0.83 |
| `endothelial` | 10 | 1.00 | 0.50 | 0.36 |
| `glutamatergic_neuron_vs_gabaergic_neuron` | 10 | 1.00 | 0.90 | 1.00 |
| `l4_it` | 10 | 1.00 | 1.00 | 0.70 |
| `l5_et` | 10 | 1.00 | 1.00 | 0.70 |
| `l6b` | 10 | 1.00 | 1.00 | 0.70 |
| `l6_ct` | 10 | 1.00 | 1.00 | 0.70 |
| `lamp5` | 10 | 1.00 | 1.00 | 0.91 |
| `pan_neuronal_vs_non_neuronal` | 10 | 1.00 | 1.00 | 0.64 |
| `sst` | 10 | 1.00 | 1.00 | 0.91 |

## Least-stable modules (within-name gene Jaccard)

- `l5_6_np`: mean Jaccard 0.067 (seen in 2 runs)
- `l6_it_car3`: mean Jaccard 0.140 (seen in 4 runs)
- `opc_cell_cycle`: mean Jaccard 0.154 (seen in 2 runs)
- `capillary_endothelial`: mean Jaccard 0.239 (seen in 3 runs)
- `sncg`: mean Jaccard 0.243 (seen in 12 runs)
- `l5_it`: mean Jaccard 0.244 (seen in 9 runs)
- `myeloid`: mean Jaccard 0.301 (seen in 3 runs)
- `glia`: mean Jaccard 0.333 (seen in 2 runs)
- `sst`: mean Jaccard 0.336 (seen in 10 runs)
- `l4_it`: mean Jaccard 0.355 (seen in 10 runs)
