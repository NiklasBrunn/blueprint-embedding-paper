# Blueprint variability — gpt_5_5_unrefined_high_human_immune  (10 runs)

**Scenario:** species=human, tissue=peripheral blood|bone marrow, celltype=None

**Stability composite:** **0.703**  _(equal-weighted; components in `cohort_summary.json`)_

## Blueprint size (mean ± SD [min, max], CV)

- **n_modules**: 31.00 ± 2.87  [27.00, 36.00], CV=0.092
- **n_axes**: 29.50 ± 3.31  [25.00, 34.00], CV=0.112
- **n_unique_genes**: 200.90 ± 19.91  [173.00, 224.00], CV=0.099
- **n_resources**: 10.50 ± 1.84  [9.00, 14.00], CV=0.175
- **n_dois_unique**: 6.00 ± 1.25  [4.00, 8.00], CV=0.208
- **gene_redundancy**: 1.14 ± 0.06  [1.04, 1.24], CV=0.050
- **mean_genes_per_module**: 7.37 ± 0.57  [6.58, 8.33], CV=0.077

## Pairwise Jaccard (mean off-diag, bootstrap 95 % CI)

- **genes**:        0.557  (95% CI [0.545, 0.568])
- **DOIs**:         0.321  (95% CI [0.284, 0.358])
- **module names**: 0.198  (95% CI [0.181, 0.215])

## Structural consistency (per-name mean agreement, ∈ [0,1])

- mean axis-kind agreement:        0.972
- mean axis positive agreement:    0.978
- mean axis gate agreement:        0.957
- mean axis negatives Jaccard:     0.504
- mean module role-set agreement:  0.899
- mean within-module gene Jaccard: 0.616

## Recurrence (universal / majority / singleton / unique-total; frac in ≥80 %; Gini)

- **modules**: 1 / 14 / 96 / 153; ≥80 %: 4.6 %; Gini: 0.398
- **axes**: 0 / 0 / 206 / 245; ≥80 %: 0.0 %; Gini: 0.147
- **genes**: 74 / 200 / 101 / 398; ≥80 %: 31.2 %; Gini: 0.387
- **resources**: 0 / 8 / 25 / 44; ≥80 %: 2.3 %; Gini: 0.404

## Saturation slope at N (absolute new-item gain when adding the N-th run)

- genes:    11.09
- modules:  9.47
- DOIs:     0.61

## Cohort coverage (totals across runs)

- unique modules:   **153**
- unique axes:      **245**
- unique genes:     **398**
- unique resources (DOI): **44**

## Resource-kind composition (mean fraction per run ± SD)

- paper:    0.278 ± 0.168
- atlas:    0.298 ± 0.129
- database: 0.424 ± 0.107

## Top hubs (most-reused modules across all runs × axes)

### As positive (axis identity carrier)
- `megakaryocyte_platelet` — used in 10 (run × axis) slots
- `plasma_cell` — used in 10 (run × axis) slots
- `classical_monocyte` — used in 9 (run × axis) slots
- `interferon_response` — used in 8 (run × axis) slots
- `eosinophil` — used in 8 (run × axis) slots
### As gate (axis activation switch)
- `hematopoietic_cd45` — used in 17 (run × axis) slots
- `t_cell` — used in 17 (run × axis) slots
- `myeloid_lineage` — used in 9 (run × axis) slots
- `t_cell_core` — used in 7 (run × axis) slots
- `immune_leukocyte` — used in 7 (run × axis) slots
### As negative (axis contrast partner)
- `classical_monocyte` — used in 25 (run × axis) slots
- `cDC2` — used in 17 (run × axis) slots
- `megakaryocyte_platelet` — used in 14 (run × axis) slots
- `nk_cell` — used in 11 (run × axis) slots
- `cDC1` — used in 11 (run × axis) slots

## Top resources (most-cited DOIs)

- Bandyopadhyay 2024 Cell  (`10.1016/j.cell.2024.04.013`, atlas) — 8/10 runs (80 %)
- Zheng et al. 2017 Nature Communications  (`10.1038/ncomms14049`, paper) — 7/10 runs (70 %)
- Velten 2017 Nature Cell Biology  (`10.1038/ncb3493`, paper) — 6/10 runs (60 %)
- Triana 2021 Nature Immunology  (`10.1038/s41590-021-01059-0`, atlas) — 6/10 runs (60 %)
- Villani et al. 2017 Science  (`10.1126/science.aah4573`, paper) — 6/10 runs (60 %)

## Sanity

- runs with placeholders:                0.0 %
- mean orphaned modules / run:           0.00
- mean modules missing resources / run:  0.00
- mean unused resources / run:           0.00
- runs with clip override:               0.0 %
- canonicalization merges (concepts with ≥ 2 distinct raw spellings): modules=9, axes=2, genes=0
- post-canonical near-duplicate pairs (edit ≤ 2): 20; of which 11 never co-occur (higher-signal drift candidates).
- DOI-kind agreement:                    0.992

### Post-canonical near-duplicate names (review — possible typos)

- `b_lymphoid_progenitor` ↔ `lymphoid_progenitor` (canonical `blymphoidprogenitor` / `lymphoidprogenitor`, edit=1, runs_a=2, runs_b=2)
- `cd4_naive_memory` ↔ `cd4_naive_memory_t` (canonical `cd4naivememory` / `cd4naivememoryt`, edit=1, runs_a=1, runs_b=4)
- `mesenchymal_stroma` ↔ `mesenchymal_stromal` (canonical `mesenchymalstroma` / `mesenchymalstromal`, edit=1, runs_a=1, runs_b=4)
- `adipocyte` ↔ `adipocyte_bm` (canonical `adipocyte` / `adipocytebm`, edit=2, runs_a=1, runs_b=1)
- `cytotoxic_T` ↔ `t_cytotoxic` (canonical `cytotoxict` / `tcytotoxic`, edit=2, runs_a=1, runs_b=1)
- `endothelial` ↔ `endothelial_bm` (canonical `endothelial` / `endothelialbm`, edit=2, runs_a=6, runs_b=1)
- `macrophage` ↔ `macrophage_bm` (canonical `macrophage` / `macrophagebm`, edit=2, runs_a=3, runs_b=1)
- `mast_cell` ↔ `pan_t_cell` (canonical `mastcell` / `pantcell`, edit=2, runs_a=1, runs_b=2)
- `mesenchymal_stromal` ↔ `mesenchymal_stromal_bm` (canonical `mesenchymalstromal` / `mesenchymalstromalbm`, edit=2, runs_a=4, runs_b=1)
- `naive_memory_t` ↔ `t_naive_memory` (canonical `naivememoryt` / `tnaivememory`, edit=2, runs_a=1, runs_b=1)
- `osteolineage` ↔ `osteolineage_bm` (canonical `osteolineage` / `osteolineagebm`, edit=2, runs_a=6, runs_b=1)

## Top axes by structural drift (low = unstable)

| axis | n_runs | positive_agree | gate_agree | mean_neg_jacc |
|---|---|---|---|---|
| `adipocyte_within_stroma` | 2 | 0.50 | 0.50 | 0.00 |
| `cdc2_identity` | 2 | 0.50 | 0.50 | 0.00 |
| `erythroid_identity` | 2 | 0.50 | 1.00 | 1.00 |
| `hsc_mpp_stemness` | 2 | 0.50 | 1.00 | 1.00 |
| `hspc_lymphoid_priming` | 2 | 0.50 | 0.50 | 0.00 |
| `mature_b_cell_identity` | 2 | 0.50 | 1.00 | 1.00 |
| `osteolineage_within_stroma` | 2 | 0.50 | 0.50 | 0.00 |
| `regulatory_t` | 2 | 0.50 | 0.50 | 0.00 |
| `b_cell_vs_t_cell` | 3 | 0.67 | 1.00 | 0.33 |
| `cdc1_identity` | 3 | 0.67 | 0.67 | 0.33 |
| `erythroid_vs_megakaryocyte` | 3 | 0.67 | 1.00 | 1.00 |
| `immune_interferon_state` | 3 | 0.67 | 0.33 | 1.00 |
| `classical_vs_nonclassical_monocyte` | 4 | 1.00 | 1.00 | 1.00 |
| `eosinophil_identity` | 3 | 1.00 | 0.67 | 0.33 |
| `marrow_macrophage_identity` | 3 | 1.00 | 0.67 | 0.33 |

## Least-stable modules (within-name gene Jaccard)

- `sinusoidal_endothelial`: mean Jaccard 0.000 (seen in 2 runs)
- `myeloid_core`: mean Jaccard 0.143 (seen in 2 runs)
- `lymphoid_lineage`: mean Jaccard 0.273 (seen in 2 runs)
- `arterial_endothelial`: mean Jaccard 0.333 (seen in 2 runs)
- `lymphoid_progenitor`: mean Jaccard 0.333 (seen in 2 runs)
- `mesenchymal_stromal`: mean Jaccard 0.390 (seen in 4 runs)
- `pan_t_cell`: mean Jaccard 0.400 (seen in 2 runs)
- `marrow_macrophage`: mean Jaccard 0.447 (seen in 3 runs)
- `conventional_dendritic`: mean Jaccard 0.455 (seen in 2 runs)
- `plasmacytoid_dendritic`: mean Jaccard 0.455 (seen in 2 runs)
