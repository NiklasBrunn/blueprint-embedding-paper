# Blueprint variability — gpt_5_5_unrefined_high_mouse_cortex  (10 runs)

**Scenario:** species=mouse, tissue=isocortex, celltype=None

**Stability composite:** **0.788**  _(equal-weighted; components in `cohort_summary.json`)_

## Blueprint size (mean ± SD [min, max], CV)

- **n_modules**: 25.40 ± 2.95  [22.00, 30.00], CV=0.116
- **n_axes**: 23.80 ± 2.86  [19.00, 28.00], CV=0.120
- **n_unique_genes**: 145.50 ± 12.09  [128.00, 165.00], CV=0.083
- **n_resources**: 10.00 ± 1.89  [8.00, 14.00], CV=0.189
- **n_dois_unique**: 7.40 ± 1.43  [5.00, 9.00], CV=0.193
- **gene_redundancy**: 1.12 ± 0.07  [1.02, 1.25], CV=0.064
- **mean_genes_per_module**: 6.46 ± 0.40  [5.89, 7.09], CV=0.062

## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)

- **genes**:        0.568  (95% CI [0.562, 0.573])
- **DOIs**:         0.576  (95% CI [0.544, 0.611])
- **module names**: 0.786  (95% CI [0.767, 0.805])

## Structural consistency (per-name mean agreement, ∈ [0,1])

- mean axis-kind agreement:        0.938
- mean axis positive agreement:    1.000
- mean axis gate agreement:        0.938
- mean axis negatives Jaccard:     0.765
- mean module role-set agreement:  0.891
- mean within-module gene Jaccard: 0.544

## Recurrence (universal / majority / singleton / unique-total; frac in ≥80 %; Gini)

- **modules**: 20 / 25 / 7 / 37; ≥80 %: 56.8 %; Gini: 0.286
- **axes**: 3 / 21 / 22 / 54; ≥80 %: 33.3 %; Gini: 0.436
- **genes**: 73 / 127 / 109 / 315; ≥80 %: 30.2 %; Gini: 0.433
- **resources**: 2 / 9 / 18 / 32; ≥80 %: 15.6 %; Gini: 0.484

## Saturation slope at N (absolute new-item gain when adding the N-th run)

- genes:    10.69
- modules:  0.77
- DOIs:     0.37

## Cohort coverage (totals across runs)

- unique modules:   **37**
- unique axes:      **54**
- unique genes:     **315**
- unique resources (DOI): **32**

## Resource-kind composition (mean fraction per run ± SD)

- paper:    0.328 ± 0.156
- atlas:    0.419 ± 0.141
- database: 0.253 ± 0.125

## Top hubs (most-reused modules across all runs × axes)

### As positive (axis identity carrier)
- `mature_oligodendrocyte` — used in 13 (run × axis) slots
- `pan_neuronal` — used in 10 (run × axis) slots
- `glutamatergic_neuron` — used in 10 (run × axis) slots
- `l23_it` — used in 10 (run × axis) slots
- `l4_it` — used in 10 (run × axis) slots
### As gate (axis activation switch)
- `glutamatergic_neuron` — used in 67 (run × axis) slots
- `gabaergic_neuron` — used in 53 (run × axis) slots
- `pan_neuronal` — used in 10 (run × axis) slots
- `non_neuronal` — used in 7 (run × axis) slots
- `vascular` — used in 4 (run × axis) slots
### As negative (axis contrast partner)
- `l4_it` — used in 57 (run × axis) slots
- `l5_it` — used in 57 (run × axis) slots
- `l5_et` — used in 57 (run × axis) slots
- `l6_ct` — used in 57 (run × axis) slots
- `l6b` — used in 57 (run × axis) slots

## Top resources (most-cited DOIs)

- Saunders 2018 Cell molecular diversity and specializations among cells of the adult mouse brain  (`10.1016/j.cell.2018.07.028`, atlas) — 10/10 runs (100 %)
- Tasic 2018 Nature shared and distinct transcriptomic cell types across mouse neocortical areas  (`10.1038/s41586-018-0654-5`, atlas) — 10/10 runs (100 %)
- Yao 2021 Cell taxonomy of transcriptomic cell types across adult mouse isocortex and hippocampal formation  (`10.1016/j.cell.2021.04.021`, atlas) — 9/10 runs (90 %)
- Tasic 2016 Nature Neuroscience adult mouse cortical cell taxonomy  (`10.1038/nn.4216`, atlas) — 8/10 runs (80 %)
- Hrvatin 2018 Nature Neuroscience experience-dependent transcriptomic states in mouse visual cortex  (`10.1038/s41593-017-0029-5`, paper) — 8/10 runs (80 %)

## Sanity

- runs with placeholders:                0.0 %
- mean orphaned modules / run:           0.00
- mean modules missing resources / run:  0.00
- mean unused resources / run:           0.00
- runs with clip override:               0.0 %
- canonicalization merges (concepts with ≥ 2 distinct raw spellings): modules=25, axes=32, genes=0
- post-canonical near-duplicate pairs (edit ≤ 2): 10; of which 0 never co-occur (higher-signal drift candidates).
- DOI-kind agreement:                    0.935

## Top axes by structural drift (low = unstable)

| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |
|---|---|---|---|---|
| `glutamatergic_neuron_vs_gabaergic_neuron` | 10 | 1.00 | 1.00 | 1.00 |
| `neuronal_activity_ieg` | 10 | 1.00 | 1.00 | 1.00 |
| `vlmc` | 10 | 1.00 | 0.80 | 0.62 |
| `astrocyte` | 9 | 1.00 | 0.78 | 0.60 |
| `endothelial` | 9 | 1.00 | 0.78 | 0.59 |
| `l23_it` | 9 | 1.00 | 1.00 | 0.84 |
| `l4_it` | 9 | 1.00 | 1.00 | 0.84 |
| `l5_et` | 9 | 1.00 | 1.00 | 0.84 |
| `l5_it` | 9 | 1.00 | 1.00 | 0.84 |
| `l6b` | 9 | 1.00 | 1.00 | 0.84 |
| `l6_ct` | 9 | 1.00 | 1.00 | 0.84 |
| `lamp5` | 9 | 1.00 | 1.00 | 0.85 |
| `microglia` | 9 | 1.00 | 0.67 | 0.43 |
| `opc` | 9 | 1.00 | 0.67 | 0.44 |
| `pvalb` | 9 | 1.00 | 1.00 | 0.85 |

## Least-stable modules (within-name gene Jaccard)

- `l6_it_car3`: mean Jaccard 0.087 (seen in 4 runs)
- `sncg`: mean Jaccard 0.236 (seen in 10 runs)
- `l5_it`: mean Jaccard 0.272 (seen in 10 runs)
- `l5_6_np`: mean Jaccard 0.310 (seen in 3 runs)
- `l6b`: mean Jaccard 0.330 (seen in 10 runs)
- `non_neuronal`: mean Jaccard 0.369 (seen in 5 runs)
- `l23_it`: mean Jaccard 0.416 (seen in 10 runs)
- `l6_ct`: mean Jaccard 0.442 (seen in 10 runs)
- `pvalb`: mean Jaccard 0.450 (seen in 10 runs)
- `l4_it`: mean Jaccard 0.455 (seen in 10 runs)
