"""
Held-out reference->query transfer evaluation.

Holds one study/technology out of a preprocessed dataset as an unseen *query*,
trains the conventional baselines on the *reference* only, and maps the query in
two ways that isolate the role of query-side training:

  - <method>_frozen   : query mapped through the reference model with NO query
                        fine-tune (scArches load_query_data, but train skipped).
  - <method>_scarches : query mapped WITH a lightweight, unsupervised scArches
                        surgery fine-tune (reference weights frozen).

The literature-authored *blueprint* needs neither — it scores the query zero-shot
(pure forward evaluation of cited gene programs). It is the no-training anchor.

Every method-variant's joint ref+query embedding is scored with the *full* scIB
suite (the same `run_benchmarker` used in scripts/embedding_metrics/metrics.py)
plus a transfer-specific kNN label-transfer classifier (train on the reference
embedding + labels, predict the query). One CSV row per method-variant.

No data leakage: the blueprint scores ref+query jointly in one call, so the
per-axis z-score is a single shared, label-free standardization.

Inputs:
    data/preprocessed/<dataset>_pp.h5ad
    a populated blueprint  (outs/blueprints/<dataset>/<cohort>/embedder_blueprint_seed<N>.py)

Outputs:
    outs/embeddings/<dataset>/transfer/<variant>_<dataset>_{query,joint}_seed<N>.parquet
    outs/embeddings/<dataset>/transfer/_joint_input_seed<N>.h5ad
    outs/metrics/<dataset>/transfer_metrics.csv          (one row per method-variant)
    models/<dataset>/{scvi,expimap}_transfer_ref/seed<N>/

Usage:
    python scripts/transfer/measure.py --dataset mouse_cortex \
        --blueprint-run outs/blueprints/mouse_cortex/gpt_5_5_unrefined_high/embedder_blueprint_seed1.py
    python scripts/transfer/measure.py --dataset mouse_cortex --methods expimap \
        --subsample 2000 --surgery-epochs 2 --ref-epochs 2 --seed 1
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import time
from datetime import datetime
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp

# scArches 0.6.1 imports `from anndata import read`, removed in anndata 0.10.
if not hasattr(ad, "read"):
    ad.read = ad.read_h5ad

REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET_CHOICES = ("mouse_cortex", "human_pancreas", "human_immune")

# query_key/values = what is held out as the unseen query; batch_key = the scIB
# batch (== query_key here, so iLISI/kBET measure query<->reference mixing);
# label_col = the shared label for transfer. expiMap gene-set knobs mirror
# scripts/models/expimap.py DEFAULT_PARAMS.
DEFAULT_PARAMS = {
    "mouse_cortex": {
        "query_key": "study", "query_values": ["tasic2016"], "batch_key": "study",
        "label_col": "celltype_harmonized",
        "organism": "mouse", "msigdb_dbver": "2024.1.Mm", "category": "m2.cp.reactome",
    },
    "human_pancreas": {
        "query_key": "tech", "query_values": ["smartseq2"], "batch_key": "tech",
        "label_col": "celltype",
        "organism": "human", "msigdb_dbver": "2024.1.Hs", "category": "c2.cp.reactome",
    },
    # Held out one bone-marrow donor (Oetjen_U) + one PBMC sample (Sun_sample4_TC):
    # both tissues, both from multi-batch studies so the reference keeps the same
    # studies' other batches -> all 16 final_annotation types stay present in the
    # reference (0 unseen). Override with --query-values for a different holdout.
    "human_immune": {
        "query_key": "batch", "query_values": ["Oetjen_U", "Sun_sample4_TC"], "batch_key": "batch",
        "label_col": "final_annotation",
        "organism": "human", "msigdb_dbver": "2024.1.Hs", "category": "c2.cp.reactome",
    },
}

N_LATENT = 10               # scVI latent (matches scripts/models/scvi.py)
EXPIMAP_TRAIN = dict(min_genes=12, n_epochs=400, alpha=0.7, alpha_kl=0.5,
                     hidden_layer_sizes=[256, 256, 256])
KNN_K = 15


# --------------------------------------------------------------------
# small helpers (kept inline to avoid importing scripts/models/*, which
# shadow the installed `scvi` package and need sys.path surgery)
# --------------------------------------------------------------------
def _round_counts(X):
    if sp.issparse(X):
        X.data = np.rint(X.data).astype(X.data.dtype)
        return X
    return np.rint(X).astype(X.dtype)


def _resolve_accelerator(req: str) -> str:
    req = (req or "auto").lower()
    if req in ("gpu", "cuda"):
        try:
            import torch
            if torch.cuda.is_available():
                return "cuda"
        except Exception:
            pass
        return "auto"
    if req == "mps":
        try:
            import torch
            if torch.backends.mps.is_available():
                return "mps"
        except Exception:
            pass
        return "auto"
    return req if req in ("auto", "cpu") else "auto"


def _dense(X) -> np.ndarray:
    return X.toarray() if sp.issparse(X) else np.asarray(X)


def write_parquet(emb: np.ndarray, cell_ids: list[str], path: Path, kv: dict) -> None:
    cols = {f"dim_{i}": emb[:, i].astype(np.float32) for i in range(emb.shape[1])}
    table = pa.Table.from_pydict({"cell_id": list(cell_ids), **cols})
    md = {**(table.schema.metadata or {}),
          **{k.encode(): str(v).encode() for k, v in kv.items()}}
    path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table.replace_schema_metadata(md), path, compression="zstd")


def _load_module(path: Path, name: str):
    """Import a .py file under a unique module name. Temporarily put its parent
    dir on sys.path so sibling imports (e.g. `from _resource_meta import ...`)
    resolve during exec."""
    parent = str(path.parent)
    sys.path.insert(0, parent)
    try:
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        sys.modules[name] = module
        spec.loader.exec_module(module)
        return module
    finally:
        if parent in sys.path:
            sys.path.remove(parent)


# reuse the project's scIB engine wholesale (lazy scib_metrics import inside it)
_METRICS = _load_module(REPO_ROOT / "scripts" / "embedding_metrics" / "metrics.py",
                        "_transfer_emetrics")
# reuse expiMap's MSigDB gene-set download/parse + recipe constants
_EXPIMAP = _load_module(REPO_ROOT / "scripts" / "models" / "expimap.py",
                        "_transfer_expimap")


def split_ref_query(adata: ad.AnnData, key: str, values: list[str]):
    levels = set(map(str, adata.obs[key].astype(str).unique()))
    missing = [v for v in values if v not in levels]
    if missing:
        raise ValueError(
            f"query value(s) {missing} not found in obs[{key!r}]; available: {sorted(levels)}")
    mask = adata.obs[key].astype(str).isin([str(v) for v in values]).to_numpy()
    ref, query = adata[~mask].copy(), adata[mask].copy()
    # drop the now-empty query level so scArches sees a genuinely new batch
    for a in (ref, query):
        a.obs[key] = a.obs[key].astype(str).astype("category")
    return ref, query


def _hvg_counts(adata: ad.AnnData) -> ad.AnnData:
    if "highly_variable" not in adata.var:
        raise RuntimeError("var['highly_variable'] missing; rerun preprocessing.")
    if "counts" not in adata.layers:
        raise RuntimeError("layers['counts'] missing; rerun preprocessing.")
    sub = adata[:, adata.var["highly_variable"].to_numpy()].copy()
    sub.X = _round_counts(sub.layers["counts"].copy())
    return sub


# --------------------------------------------------------------------
# scVI: frozen (no fine-tune) and scArches surgery, from one reference model
# --------------------------------------------------------------------
def run_scvi(ref, query, batch_key, seed, surgery_epochs, ref_epochs, accel, ckpt_dir):
    import scvi
    scvi.settings.seed = seed
    ref_sub, q_sub = _hvg_counts(ref), _hvg_counts(query)

    scvi.model.SCVI.setup_anndata(ref_sub, batch_key=batch_key)
    # scArches-recommended reference-mapping architecture: encode_covariates=True
    # makes the ENCODER (hence the latent z we read out) condition on the batch,
    # so query-side surgery actually adapts the embedding. With scVI's default
    # encode_covariates=False the encoder ignores batch and z is batch-invariant,
    # which would make the frozen-vs-surgery query latents bit-identical.
    ref_model = scvi.model.SCVI(
        ref_sub, n_latent=N_LATENT, n_layers=2, encode_covariates=True,
        deeply_inject_covariates=False, use_layer_norm="both", use_batch_norm="none")
    acc = _resolve_accelerator(accel)
    tkw = {"accelerator": acc, **({"devices": 1} if acc != "auto" else {})}
    if ref_epochs is not None:
        tkw["max_epochs"] = ref_epochs
    ref_model.train(**tkw)
    ckpt_dir.parent.mkdir(parents=True, exist_ok=True)
    ref_model.save(str(ckpt_dir), overwrite=True)
    ref_emb = np.asarray(ref_model.get_latent_representation(ref_sub), dtype=np.float32)

    scvi.model.SCVI.prepare_query_anndata(q_sub, ref_model)
    q = scvi.model.SCVI.load_query_data(q_sub, ref_model)
    q.is_trained = True   # read the un-fine-tuned (init query-param) latent; cf. scvi_census.py
    frozen_q = np.asarray(q.get_latent_representation(q_sub), dtype=np.float32)
    t0 = time.perf_counter()
    q.train(max_epochs=surgery_epochs, plan_kwargs={"weight_decay": 0.0})
    surgery_secs = time.perf_counter() - t0
    scarches_q = np.asarray(q.get_latent_representation(q_sub), dtype=np.float32)

    order = list(ref_sub.obs_names) + list(q_sub.obs_names)
    return {
        "scvi_frozen":   dict(emb=np.vstack([ref_emb, frozen_q]), order=order,
                              surgery="frozen", surgery_secs=0.0, surgery_epochs=0),
        "scvi_scarches": dict(emb=np.vstack([ref_emb, scarches_q]), order=order,
                              surgery="scarches", surgery_secs=surgery_secs,
                              surgery_epochs=surgery_epochs),
    }


# --------------------------------------------------------------------
# expiMap: frozen and scArches surgery, both via load_query_data (parallels scVI)
# --------------------------------------------------------------------
def run_expimap(ref, query, batch_key, seed, surgery_epochs, ref_epochs, gmt_path, ckpt_dir):
    import scarches as sca
    ref_sub, q_sub = _hvg_counts(ref), _hvg_counts(query)
    for a in (ref_sub, q_sub):
        sca.utils.add_annotations(a, str(gmt_path), min_genes=EXPIMAP_TRAIN["min_genes"],
                                  max_genes=None, clean=False)
    if int(ref_sub.varm.get("I", np.zeros((1, 0))).shape[1]) == 0:
        raise RuntimeError("No gene programs survived min_genes on the reference.")

    ref_model = sca.models.EXPIMAP(
        adata=ref_sub, condition_key=batch_key,
        hidden_layer_sizes=EXPIMAP_TRAIN["hidden_layer_sizes"], recon_loss="nb", n_ext=0)
    ref_model.train(n_epochs=(ref_epochs or EXPIMAP_TRAIN["n_epochs"]),
                    alpha=EXPIMAP_TRAIN["alpha"], alpha_kl=EXPIMAP_TRAIN["alpha_kl"],
                    alpha_epoch_anneal=100, weight_decay=0.0,
                    early_stopping_kwargs=_EXPIMAP.EARLY_STOPPING_KWARGS,
                    use_early_stopping=True, monitor_only_val=False, seed=seed)
    ckpt_dir.parent.mkdir(parents=True, exist_ok=True)
    ref_model.save(str(ckpt_dir), overwrite=True, save_anndata=False)

    # load_query_data adds the query condition node (init); train = surgery.
    q = sca.models.EXPIMAP.load_query_data(q_sub, ref_model)

    def latent(model, a):
        return np.asarray(model.get_latent(x=_dense(a.X), c=a.obs[batch_key].to_numpy(),
                                           mean=True, only_active=True), dtype=np.float32)

    order = list(ref_sub.obs_names) + list(q_sub.obs_names)
    ref_f, qry_f = latent(q, ref_sub), latent(q, q_sub)          # frozen (pre-surgery)
    out = {"expimap_frozen": dict(emb=np.vstack([ref_f, qry_f]), order=order,
                                  surgery="frozen", surgery_secs=0.0, surgery_epochs=0)}
    # surgery is the API-fragile path; keep the frozen variant even if it drifts
    try:
        t0 = time.perf_counter()
        q.train(n_epochs=surgery_epochs, alpha=EXPIMAP_TRAIN["alpha"],
                alpha_kl=EXPIMAP_TRAIN["alpha_kl"], weight_decay=0.0,
                use_early_stopping=False, seed=seed)
        surgery_secs = time.perf_counter() - t0
        ref_s, qry_s = latent(q, ref_sub), latent(q, q_sub)      # scArches (post-surgery)
        out["expimap_scarches"] = dict(emb=np.vstack([ref_s, qry_s]), order=order,
                                       surgery="scarches", surgery_secs=surgery_secs,
                                       surgery_epochs=surgery_epochs)
    except Exception as exc:
        print(f"[warn] expiMap surgery train failed ({exc!r}); reporting frozen only.")
    return out


# --------------------------------------------------------------------
# blueprint: zero-shot, scored jointly over ref+query (shared z-score), then split
# --------------------------------------------------------------------
def run_blueprint(variant, adata_joint, ref_ids, query_ids, blueprint_py, joint_h5ad, joint_parquet):
    if not joint_h5ad.exists():                 # write the shared joint input once
        joint_h5ad.parent.mkdir(parents=True, exist_ok=True)
        adata_joint.write_h5ad(joint_h5ad)
    mod = _load_module(Path(blueprint_py).resolve(), f"_transfer_bp_{variant}")
    bp = mod.build_blueprint("transfer")
    bp.build_embedding(joint_h5ad, joint_parquet, layer=None)   # adata.X holds log1p (set during preprocessing)
    X, ids = _METRICS.load_embedding(joint_parquet)
    pos = {c: i for i, c in enumerate(ids)}
    order = list(ref_ids) + list(query_ids)
    emb = X[np.array([pos[c] for c in order], dtype=int)]
    return {variant: dict(emb=emb, order=order, surgery="none",
                          surgery_secs=0.0, surgery_epochs=0)}


def knn_label_transfer(joint_emb, order, ref_ids, query_ids, labels_by_id):
    from sklearn.neighbors import KNeighborsClassifier
    from sklearn.metrics import accuracy_score, f1_score
    pos = {c: i for i, c in enumerate(order)}
    ref_idx = np.array([pos[c] for c in ref_ids], dtype=int)
    qry_idx = np.array([pos[c] for c in query_ids], dtype=int)
    y_ref = np.array([labels_by_id[c] for c in ref_ids])
    y_qry = np.array([labels_by_id[c] for c in query_ids])
    seen = set(np.unique(y_ref).tolist())
    keep = np.array([y in seen for y in y_qry])           # query labels absent from ref are unscorable
    clf = KNeighborsClassifier(n_neighbors=KNN_K, weights="distance")
    clf.fit(joint_emb[ref_idx], y_ref)
    pred = clf.predict(joint_emb[qry_idx][keep])
    return {
        "label_transfer_acc": float(accuracy_score(y_qry[keep], pred)),
        "label_transfer_macro_f1": float(f1_score(y_qry[keep], pred, average="macro")),
        "n_query_labels_unseen": int((~keep).sum()),
    }


def joint_pca(adata_joint, order, n_comp=50):
    """Unintegrated PCA reference (HVG log1p) for the scIB PCR metric."""
    try:
        from sklearn.decomposition import PCA
        sub = adata_joint[:, adata_joint.var["highly_variable"].to_numpy()]
        Xl = _dense(sub.layers["log1p"]).astype(np.float32)
        pcs = PCA(n_components=min(n_comp, min(Xl.shape) - 1)).fit_transform(Xl)
        pos = {c: i for i, c in enumerate(adata_joint.obs_names.astype(str))}
        return pcs[np.array([pos[c] for c in order], dtype=int)]
    except Exception as exc:
        print(f"[warn] joint PCA for PCR failed ({exc}); PCR will be NaN.")
        return None


def scib_row(variant_emb, order, obs_df, label_col, batch_key, unint, n_jobs):
    """Full scIB suite on one joint embedding via the project Benchmarker."""
    bm = _METRICS.run_benchmarker({variant_emb["_name"]: variant_emb["emb"]}, order,
                                  obs_df, label_col, batch_key, unint,
                                  n_neighbors=15, n_jobs=n_jobs)
    res = bm.get_results(min_max_scale=False, clean_names=True)
    stem = variant_emb["_name"]
    out = {}
    for src, dst in _METRICS.RESULT_COL_MAP.items():
        if src in res.columns:
            try:
                out[dst] = float(res.at[stem, src])
            except (TypeError, ValueError, KeyError):
                pass
    return out


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Held-out reference->query transfer evaluation.")
    p.add_argument("--dataset", choices=DATASET_CHOICES, required=True)
    p.add_argument("--blueprint-run", type=str,
                   help="Resourced (literature-cited) populated embedder_blueprint_seed<N>.py "
                        "-> 'blueprint' variant.")
    p.add_argument("--blueprint-noref-run", type=str,
                   help="Noref (citations-optional) populated embedder_blueprint_seed<N>.py "
                        "-> 'blueprint_noref' variant. Cheap to add (zero-shot, no training).")
    p.add_argument("--methods", type=str, default="blueprint,scvi,expimap",
                   help="Comma list from {blueprint,scvi,expimap}.")
    p.add_argument("--query-key", type=str)
    p.add_argument("--query-values", type=str, help="Comma-separated obs values to hold out.")
    p.add_argument("--surgery-epochs", type=int, default=100)
    p.add_argument("--ref-epochs", type=int, default=None,
                   help="Cap reference-training epochs (default: model heuristic / 400).")
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--accelerator", choices=["auto", "cpu", "mps", "cuda", "gpu"], default="auto")
    p.add_argument("--gene-sets", type=str, help="Path to a .gmt for expiMap (else MSigDB download).")
    p.add_argument("--subsample", type=int, default=0, help="Debug: cap ref+query cell counts.")
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--overwrite", action="store_true")
    return p.parse_args()


def main() -> None:
    args = parse_args()
    cfg = DEFAULT_PARAMS[args.dataset]
    query_key = args.query_key or cfg["query_key"]
    query_values = ([v.strip() for v in args.query_values.split(",")]
                    if args.query_values else cfg["query_values"])
    batch_key, label_col = cfg["batch_key"], cfg["label_col"]
    methods = [m.strip() for m in args.methods.split(",") if m.strip()]
    seed_tag = f"seed{args.seed}"

    out_emb = REPO_ROOT / "outs" / "embeddings" / args.dataset / "transfer"
    out_csv = REPO_ROOT / "outs" / "metrics" / args.dataset / "transfer_metrics.csv"

    adata = ad.read_h5ad(REPO_ROOT / "data" / "preprocessed" / f"{args.dataset}_pp.h5ad")
    ref, query = split_ref_query(adata, query_key, query_values)
    if args.subsample:
        ref = ref[:min(args.subsample, ref.n_obs)].copy()
        query = query[:min(max(args.subsample // 4, 50), query.n_obs)].copy()
    print(f"[ok] {args.dataset}: ref={ref.n_obs} query={query.n_obs} "
          f"(held out {query_key}={query_values}); label={label_col} batch={batch_key}")

    ref_ids, query_ids = list(ref.obs_names), list(query.obs_names)
    # Build the joint AnnData by indexing the original (preserves var order, HVG
    # flags, and all layers exactly — avoids ad.concat var-merge surprises).
    pos = {c: i for i, c in enumerate(adata.obs_names.astype(str))}
    joint = adata[np.array([pos[c] for c in (ref_ids + query_ids)], dtype=int)].copy()
    order_all = list(joint.obs_names.astype(str))

    variants: dict[str, dict] = {}
    if "blueprint" in methods:
        bp_specs = []
        if args.blueprint_run:
            bp_specs.append(("blueprint", args.blueprint_run))
        if args.blueprint_noref_run:
            bp_specs.append(("blueprint_noref", args.blueprint_noref_run))
        if not bp_specs:
            raise SystemExit("--blueprint-run and/or --blueprint-noref-run required "
                             "when 'blueprint' is in --methods")
        joint_h5ad = out_emb / f"_joint_input_{seed_tag}.h5ad"
        for variant, path in bp_specs:
            variants.update(run_blueprint(
                variant, joint, ref_ids, query_ids, path, joint_h5ad,
                out_emb / f"{variant}_{args.dataset}_joint_{seed_tag}.parquet"))
    if "scvi" in methods:
        variants.update(run_scvi(
            ref, query, batch_key, args.seed, args.surgery_epochs, args.ref_epochs,
            args.accelerator, REPO_ROOT / "models" / args.dataset / "scvi_transfer_ref" / seed_tag))
    if "expimap" in methods:
        if args.gene_sets:
            gmt = Path(args.gene_sets).resolve()
        else:
            cache = REPO_ROOT / "models" / "genesets" / cfg["organism"]
            gmt = _EXPIMAP._download_msigdb_gmt(cfg["category"], cfg["msigdb_dbver"], cache)
        try:
            variants.update(run_expimap(
                ref, query, batch_key, args.seed, args.surgery_epochs, args.ref_epochs, gmt,
                REPO_ROOT / "models" / args.dataset / "expimap_transfer_ref" / seed_tag))
        except Exception as exc:
            print(f"[warn] expiMap transfer failed ({exc!r}); skipping its variants.")

    if not variants:
        raise SystemExit("No method-variants produced an embedding.")

    # obs + labels shared across variants (unintegrated PCA is per-variant, in the loop below)
    obs_df = pd.DataFrame(
        {label_col: joint.obs[label_col].astype(str).to_numpy(),
         batch_key: joint.obs[batch_key].astype(str).to_numpy()},
        index=pd.Index(order_all, name="cell_id"))
    labels_by_id = dict(zip(order_all, obs_df[label_col].to_numpy()))

    rows = []
    for name, v in variants.items():
        v["_name"] = name
        # persist the joint embedding (+ query slice) for audits / UMAPs
        write_parquet(v["emb"], v["order"], out_emb / f"{name}_{args.dataset}_joint_{seed_tag}.parquet",
                      {"dataset": args.dataset, "variant": name, "seed": args.seed})
        qpos = {c: i for i, c in enumerate(v["order"])}
        qidx = np.array([qpos[c] for c in query_ids], dtype=int)
        write_parquet(v["emb"][qidx], query_ids,
                      out_emb / f"{name}_{args.dataset}_query_{seed_tag}.parquet",
                      {"dataset": args.dataset, "variant": name, "seed": args.seed, "split": "query"})

        unint = joint_pca(joint, v["order"])
        obs_v = obs_df.loc[v["order"], [label_col, batch_key]]
        scib = scib_row(v, v["order"], obs_v, label_col, batch_key, unint, args.workers)
        lt = knn_label_transfer(v["emb"], v["order"], ref_ids, query_ids, labels_by_id)
        row = {"method": name.split("_")[0], "variant": name,
               "surgery": v["surgery"], "surgery_epochs": v["surgery_epochs"],
               "surgery_secs": round(v["surgery_secs"], 3), "embedding_dim": int(v["emb"].shape[1]),
               "dataset": args.dataset, "seed": args.seed,
               "query_key": query_key, "query_values": ",".join(query_values),
               **lt, **scib}
        rows.append(row)
        print(f"[ok] {name}: acc={lt['label_transfer_acc']:.3f} "
              f"f1={lt['label_transfer_macro_f1']:.3f} total={scib.get('composite_scib', float('nan')):.3f} "
              f"surgery_s={v['surgery_secs']:.1f}")

    out_csv.parent.mkdir(parents=True, exist_ok=True)
    new = pd.DataFrame(rows)
    if out_csv.exists() and not args.overwrite:
        old = pd.read_csv(out_csv)
        key = ["variant", "dataset", "seed"]
        old = old[~old.set_index(key).index.isin(new.set_index(key).index)]
        new = pd.concat([old, new], ignore_index=True)
    new.sort_values(["dataset", "method", "variant", "seed"]).to_csv(out_csv, index=False)
    print(f"[ok] wrote {out_csv} ({len(rows)} new rows; {len(new)} total)")


if __name__ == "__main__":
    main()
