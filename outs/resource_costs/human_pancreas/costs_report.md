# Cost / wall-clock / storage — human_pancreas

Per-method medians over seeds (IQR in `costs_summary.csv`). **Time taxonomy** — the four buckets sum to `walltime` for every script we own: **`load`** = read the input AnnData from disk; **`obtain-model`** = get a usable model — train it (de-novo VAEs, PCA loadings) or load a persisted/pretrained checkpoint + adapt it to the dataset's gene space (the `source` column tags which: *trained* / *loaded* / *none* — diffmap & UCell acquire nothing — / *authored* — the LLM wrote the script); **`n/a`** when load-and-forward aren't separable upstream (UCE's eval subprocess); **`build-emb`** = the forward pass / scoring; **`other`** = the residual (data prep, parquet write, glue). For the **agent methods** `walltime` is just the Phase-5 cell-scoring run — the *model-acquisition* cost is the LLM authoring, **`Stage A`** (+ the optional **`Stage B`** refinement loop), a separate phase, also reported in tokens. **`model dl`** = one-time cold-cache fetch of a pretrained checkpoint (0 when already on disk; never part of `walltime` in steady state). **Memory** (`peak RSS`) is a high-water max, not a sum — see `costs.csv` for the after-load → after-model → after-embedding trajectory. **`model storage`** = the persisted artifact on disk (scVI/expiMap checkpoint, UCE/Census `.pt`, gene-set `.gmt`, or — for the agent methods — the populated `.py` script). **`device`** = compute device (cpu/cuda/mps; cpu for the agents' Phase-5 build). Blank cells = the producing script predates this instrumentation; re-run it to populate.

| method | family | device | source | n | walltime s | load s | obtain-model s | build-emb s | other s | model dl s | peak RSS MB | model storage | emb size | dim | Stage A s | Stage B refine s | refine iters | tokens A | tokens B refine |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| `blueprint_gpt_5_4_mini_noref_unrefined_high` | blueprint | cpu | authored | 10 | 1.8 | 0.7 | 0.0 | 1.1 | 0.0 | 0.0 | 1973.5 | 31.6 KB | 2.0 MB | 22 | 260 | — | — | 43128 | — |
| `blueprint_gpt_5_4_mini_noref_unrefined_medium` | blueprint | cpu | authored | 10 | 5.9 | 3.6 | 0.0 | 2.2 | 0.0 | 0.0 | 1960.1 | 31.4 KB | 1.9 MB | 17 | 170 | — | — | 32168 | — |
| `blueprint_gpt_5_4_mini_unrefined_high` | blueprint | cpu | authored | 10 | 5.3 | 3.6 | 0.0 | 1.6 | 0.0 | 0.0 | 1963.4 | 37.0 KB | 1.8 MB | 17 | 402 | — | — | 114514 | — |
| `blueprint_gpt_5_4_mini_unrefined_medium` | blueprint | cpu | authored | 10 | 1.8 | 0.6 | 0.0 | 1.1 | 0.0 | 0.0 | 1951.8 | 35.4 KB | 1.5 MB | 14 | 320 | — | — | 106826 | — |
| `blueprint_gpt_5_5_noref_refine_high` | blueprint | cpu | authored | 10 | 5.7 | 4.2 | 0.0 | 1.5 | 0.0 | 0.0 | 1997.1 | 37.4 KB | 3.0 MB | 35 | 200 | 709 | 3 | 52356 | 260960 |
| `blueprint_gpt_5_5_noref_unrefined_high` | blueprint | cpu | authored | 10 | 6.7 | 4.7 | 0.0 | 2.0 | 0.0 | 0.0 | 1983.1 | 35.1 KB | 2.9 MB | 30 | 200 | — | — | 52356 | — |
| `blueprint_gpt_5_5_refine_high` | blueprint | cpu | authored | 10 | 2.9 | 1.6 | 0.0 | 1.3 | 0.0 | 0.0 | 1972.4 | 43.4 KB | 2.1 MB | 26 | 310 | 742 | 3 | 123594 | 318416 |
| `blueprint_gpt_5_5_unrefined_high` | blueprint | cpu | authored | 10 | 2.1 | 0.8 | 0.0 | 1.3 | 0.0 | 0.0 | 1963.2 | 40.3 KB | 1.9 MB | 21 | 310 | — | — | 123594 | — |
| `pca` | conventional | cpu | trained | 1 | 0.9 | 0.6 | 0.2 | 0.0 | 0.1 | 0.0 | 2347.1 | 399.6 KB | 4.3 MB | 50 | — | — | — | — | — |
| `scvi` | conventional | cpu | trained | 10 | 719.1 | 0.4 | 714.5 | 0.1 | 4.2 | 0.0 | 2740.4 | 5.2 MB | 935.9 KB | 10 | — | — | — | — | — |
| `scvi_nobatch` | conventional | cpu | trained | 10 | 644.6 | 0.4 | 640.1 | 0.1 | 4.0 | 0.0 | 2711.1 | 5.2 MB | 936.6 KB | 10 | — | — | — | — | — |
| `freeform_constrained_gpt_5_5_high` | freeform | cpu | authored | 10 | 6.0 | — | — | — | — | 0.0 | — | 10.8 KB | 1.1 MB | 14 | 208 | — | — | 56940 | — |
| `freeform_naive_gpt_5_5_high` | freeform | cpu | authored | 10 | 6.0 | — | — | — | — | 0.0 | — | 11.9 KB | 2.8 MB | 40 | 158 | — | — | 26426 | — |
| `freeform_noref_gpt_5_5_high` | freeform | cpu | authored | 10 | 4.0 | — | — | — | — | 0.0 | — | 12.6 KB | 1.7 MB | 22 | 168 | — | — | 32159 | — |
| `expimap` | geneset | cpu | trained | 10 | 439.7 | 0.4 | 434.8 | 0.2 | 4.4 | 0.0 | 2949.4 | 6.2 MB | 12.7 MB | 152 | — | — | — | — | — |
| `ucell` | geneset | cpu | none | 1 | 6.5 | 0.4 | 0.0 | 6.2 | 0.0 | 0.0 | 3833.8 | 2.5 KB | 552.8 KB | 10 | — | — | — | — | — |
| `scvi_census` | pretrained | cpu | loaded | 1 | 18.1 | 0.4 | 13.6 | 0.6 | 3.6 | 0.0 | 9959.4 | 888.3 MB | 4.3 MB | 50 | — | — | — | — | — |
| `uce` | pretrained | cpu | loaded | 1 | 4582.0 | 0.5 | — | 4580.7 | 0.8 | 0.0 | 23626.5 | 5.3 GB | 108.2 MB | 1280 | — | — | — | — | — |
