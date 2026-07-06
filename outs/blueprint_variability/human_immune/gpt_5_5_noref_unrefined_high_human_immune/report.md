# Blueprint variability — gpt_5_5_noref_unrefined_high_human_immune  (10 runs)

**Scenario:** species=human, tissue=peripheral blood|bone marrow, celltype=None

**Stability composite:** **0.733**  _(equal-weighted; components in `cohort_summary.json`)_

## Blueprint size (mean ± SD [min, max], CV)

- **n_modules**: 39.20 ± 4.37  [35.00, 49.00], CV=0.111
- **n_axes**: 38.10 ± 4.84  [32.00, 45.00], CV=0.127
- **n_unique_genes**: 249.10 ± 15.07  [221.00, 267.00], CV=0.060
- **n_resources**: 0.00 ± 0.00  [0.00, 0.00], CV=nan
- **n_dois_unique**: 0.00 ± 0.00  [0.00, 0.00], CV=nan
- **gene_redundancy**: 1.24 ± 0.08  [1.15, 1.39], CV=0.068
- **mean_genes_per_module**: 7.90 ± 0.57  [7.17, 8.89], CV=0.072

## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)

- **genes**:        0.640  (95% CI [0.634, 0.648])
- **DOIs**:         0.000  (95% CI [0.000, 0.000])
- **module names**: 0.298  (95% CI [0.286, 0.309])

## Structural consistency (per-name mean agreement, ∈ [0,1])

- mean axis-kind agreement:        0.963
- mean axis positive agreement:    0.972
- mean axis gate agreement:        0.919
- mean axis negatives Jaccard:     0.492
- mean module role-set agreement:  0.862
- mean within-module gene Jaccard: 0.661

## Recurrence (universal / majority / singleton / unique-total; frac in ≥80 %; Gini)

- **modules**: 8 / 28 / 67 / 140; ≥80 %: 8.6 %; Gini: 0.440
- **axes**: 0 / 3 / 213 / 278; ≥80 %: 0.0 %; Gini: 0.223
- **genes**: 123 / 251 / 85 / 431; ≥80 %: 42.5 %; Gini: 0.350
- **resources**: 0 / 0 / 0 / 0; ≥80 %: 0.0 %; Gini: nan

## Saturation slope at N (absolute new-item gain when adding the N-th run)

- genes:    8.81
- modules:  6.45
- DOIs:     0.00

## Cohort coverage (totals across runs)

- unique modules:   **140**
- unique axes:      **278**
- unique genes:     **431**
- unique resources (DOI): **0**

## Resource-kind composition (mean fraction per run ± SD)

- paper:    0.000 ± 0.000
- atlas:    0.000 ± 0.000
- database: 0.000 ± 0.000

## Top hubs (most-reused modules across all runs × axes)

### As positive (axis identity carrier)
- `plasma_cell` — used in 11 (run × axis) slots
- `megakaryocyte_platelet` — used in 10 (run × axis) slots
- `nk_cell` — used in 10 (run × axis) slots
- `classical_monocyte` — used in 10 (run × axis) slots
- `eosinophil` — used in 10 (run × axis) slots
### As gate (axis activation switch)
- `myeloid_core` — used in 26 (run × axis) slots
- `t_cell` — used in 21 (run × axis) slots
- `hematopoietic_core` — used in 14 (run × axis) slots
- `hematopoietic_stem_progenitor` — used in 10 (run × axis) slots
- `t_cell_core` — used in 10 (run × axis) slots
### As negative (axis contrast partner)
- `cDC2` — used in 19 (run × axis) slots
- `cDC1` — used in 18 (run × axis) slots
- `classical_monocyte` — used in 17 (run × axis) slots
- `plasma_cell` — used in 16 (run × axis) slots
- `mature_neutrophil` — used in 15 (run × axis) slots

## Sanity

- runs with placeholders:                0.0 %
- mean orphaned modules / run:           0.00
- mean modules missing resources / run:  39.20
- mean unused resources / run:           0.00
- runs with clip override:               0.0 %
- canonicalization merges (concepts with ≥ 2 distinct raw spellings): modules=3, axes=3, genes=0
- post-canonical near-duplicate pairs (edit ≤ 2): 19; of which 2 never co-occur (higher-signal drift candidates).
- DOI-kind agreement:                    nan

### Post-canonical near-duplicate names (review — possible typos)

- `b_lymphoid_progenitor` ↔ `lymphoid_progenitor` (canonical `blymphoidprogenitor` / `lymphoidprogenitor`, edit=1, runs_a=1, runs_b=5)
- `mast_cell` ↔ `pan_t_cell` (canonical `mastcell` / `pantcell`, edit=2, runs_a=5, runs_b=1)

## Top axes by structural drift (low = unstable)

| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |
|---|---|---|---|---|
| `erythroid_maturation` | 3 | 0.33 | 0.67 | 0.00 |
| `t_cell_vs_b_cell` | 3 | 0.33 | 1.00 | 0.33 |
| `primitive_hspc` | 4 | 0.50 | 1.00 | 1.00 |
| `cd8_cytotoxic_t_cell` | 2 | 0.50 | 0.50 | 0.50 |
| `erythroid_lineage` | 2 | 0.50 | 1.00 | 1.00 |
| `erythroid_maturation_late_vs_early` | 2 | 0.50 | 0.50 | 0.00 |
| `erythroid_vs_leukocyte` | 2 | 0.50 | 1.00 | 0.00 |
| `erythroid_vs_megakaryocyte` | 2 | 0.50 | 1.00 | 1.00 |
| `hspc_lymphoid_priming` | 2 | 0.50 | 0.50 | 0.00 |
| `mature_b_cell` | 2 | 0.50 | 1.00 | 1.00 |
| `lymphoid_vs_myeloid` | 5 | 0.60 | 1.00 | 0.60 |
| `hematopoietic_vs_marrow_niche` | 3 | 0.67 | 1.00 | 0.00 |
| `marrow_macrophage` | 3 | 0.67 | 0.67 | 0.33 |
| `myeloid_inflammatory_activation` | 3 | 0.67 | 0.33 | 1.00 |
| `myeloid_inflammatory_state` | 3 | 0.67 | 0.67 | 1.00 |

## Least-stable modules (within-name gene Jaccard)

- `erythroid_early`: mean Jaccard 0.214 (seen in 2 runs)
- `leukocyte_core`: mean Jaccard 0.300 (seen in 2 runs)
- `memory_b_cell`: mean Jaccard 0.300 (seen in 2 runs)
- `gamma_delta_t`: mean Jaccard 0.333 (seen in 2 runs)
- `neutrophil_granulocyte`: mean Jaccard 0.357 (seen in 2 runs)
- `pan_leukocyte`: mean Jaccard 0.360 (seen in 4 runs)
- `naive_b_cell`: mean Jaccard 0.400 (seen in 2 runs)
- `hematopoietic_core`: mean Jaccard 0.439 (seen in 3 runs)
- `multipotent_progenitor`: mean Jaccard 0.455 (seen in 2 runs)
- `macrophage`: mean Jaccard 0.482 (seen in 5 runs)
