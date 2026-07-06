# Embedding script — naive (count mode)

Read `scenario_config.jsonc` in the current working directory. Then
write a single Python file `embed_freeform.py` in this directory that
embeds cells from a scRNA-seq dataset using literature knowledge for
the scenario.

`embed_freeform.py` must run as:

    python embed_freeform.py --input <path> --output <parquet> --seed <int>

The output parquet must have columns `cell_id` (one row per cell) and
`dim_0`, `dim_1`, ..., `dim_{K-1}` (the embedding coordinates). `K` is
up to you.

The `--input` is a `.h5ad` AnnData file: counts in `layers["counts"]`,
log1p-normalized values in `X`, gene symbols in `var_names`, cell IDs
in `obs_names`.

The embedding must be built solely from literature knowledge. Do not
train models on the data and do not run pretrained-model inference.

`--seed` is a run tag for output naming. Your script may ignore it;
runs at the same seed need not produce identical embeddings.

Deliverable: a single file `embed_freeform.py`. Nothing else.
