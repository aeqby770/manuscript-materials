#!/usr/bin/env python
"""scTenifoldXct L-R candidate scoring across ALL ordered cell-type pairs (full matrix)."""
import os
import re
import sys
import time
import random
import pickle
import torch
import gc
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

CONFIG = {
    "paths": {
        "vendor": Path(os.environ["LLR_WORK_ROOT"]).resolve() / "tools/competitors/optimized/scTenifoldXct",
    },
    "threads": int(os.environ.get("LLR_NATIVE_THREADS","2")),
    "seed": 12345,
    "plot": {},
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
torch.set_num_threads(CONFIG["threads"])
VENDOR = CONFIG["paths"]["vendor"]
if str(VENDOR) not in sys.path:
    sys.path.insert(0, str(VENDOR))
sys.path.insert(0,str(Path(os.environ["LLR_WORK_ROOT"])/"tools/competitors/dependencies/scTenifoldXct"))
from scTenifoldXct import scTenifoldXct, set_seed  # noqa: E402

set_seed(CONFIG["seed"])


def pick_group(adata):
    for c in ("paul15_clusters", "cell_type", "cluster", "clusters", "celltype"):
        if c in adata.obs.columns and adata.obs[c].nunique() >= 2:
            return c
    sc.pp.neighbors(adata, n_neighbors=min(15, adata.n_obs - 1))
    sc.tl.leiden(adata, resolution=0.5, key_added="sctenifold_leiden")
    return "sctenifold_leiden"


def main():
    if len(sys.argv) < 3:
        print("Usage: run_sctenifoldxct.py <h5ad> <output.csv>", file=sys.stderr)
        sys.exit(2)
    inpath, outpath = sys.argv[1], sys.argv[2]
    if os.path.exists(outpath) and os.path.getsize(outpath) > 100:
        print("CACHED")
        return

    adata = sc.read_h5ad(inpath)
    # scTenifoldXct builds gene-gene GRNs + manifold alignment (~O(n_genes^2)).
    n_global_cells,n_input_genes=adata.shape
    original_genes=adata.var_names.copy()
    if "lognorm" in adata.layers:
        x = adata.layers["lognorm"]
        adata.X = x.toarray() if sparse.issparse(x) else np.asarray(x)
    else:
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)

    gcol = pick_group(adata)
    groups = sorted(adata.obs[gcol].astype(str).unique().tolist())
    if len(groups) < 2:
        # Fewer than two groups does not satisfy the method's input contract.
        print("FAIL: fewer than 2 cell groups; scTenifoldXct needs >=2 to compare", file=sys.stderr)
        sys.exit(1)

    cache = os.path.join(os.path.dirname(outpath), "grn_cache")
    pair_dir = os.path.join(os.path.dirname(outpath), "sctenifoldxct_pairs")
    os.makedirs(pair_dir, exist_ok=True)
    all_results = []
    coverage_records=[]
    pd.DataFrame({"original_gene_id":original_genes,"native_gene_id":original_genes.str.upper()}).to_csv(Path(outpath).parent/"native_gene_id_binding.csv",index=False,encoding="utf-8")
    print("SCTENIFOLD_FULL_SCOPE",adata.n_obs,adata.n_vars,len(groups),len(groups)*(len(groups)-1),gcol,flush=True)
    # Full cell-type x cell-type matrix (no top-2 restriction). Each ordered pair
    # runs the official scTenifoldXct end to end: build GRNs, train the manifold-
    # aligned NN (get_embeds), then null_test() for the significance test — no step
    # skipped.
    for src in groups:
        for tgt in groups:
            if src == tgt:
                continue
            tag = re.sub(r"[^0-9A-Za-z_.-]", "_", f"{src}__{tgt}")
            pair_csv = os.path.join(pair_dir, f"{tag}.csv")
            # Per-pair disk cache (resumable, identical result).
            bundle_path=Path(pair_dir)/(tag+".pkl")
            if bundle_path.is_file():
                bundle=pickle.loads(bundle_path.read_bytes())
                np.random.set_state(bundle["numpy_rng_state"])
                random.setstate(bundle["python_rng_state"])
                torch.set_rng_state(bundle["torch_rng_state"])
                all_results.append(bundle["significant_calls"])
                coverage_records.append(bundle["coverage"])
                continue
            pair_started=time.perf_counter()
            print("SCTENIFOLD_PAIR_BEGIN",src,tgt,flush=True)
            xct = scTenifoldXct(
                data=adata,
                source_celltype=str(src),
                target_celltype=str(tgt),
                obs_label=gcol,
                GRN_file_dir=cache,
                rebuild_GRN=True,
                query_DB="pairs",
                verbose=False,
                nComp=5,
                n_cpus=1,
            )
            projections=xct.get_embeds(train=True,n_steps=1000,lr=0.01)
            df=xct.null_test(filter_zeros=True,pval=0.05,plot_result=False).copy().reset_index()
            aligned_candidates=xct._aligned_result.loc[xct._aligned_result.index.isin(xct._candidates)].copy().reset_index()
            assert np.isfinite(aligned_candidates["dist"].to_numpy(dtype=float)).all()
            assert np.isfinite(df["p_val"].to_numpy(dtype=float)).all()
            df["source_celltype"]=src
            df["target_celltype"]=tgt
            aligned_candidates["source_celltype"]=src
            aligned_candidates["target_celltype"]=tgt
            df.to_csv(pair_csv,index=False,encoding="utf-8")
            aligned_candidates.to_csv(Path(pair_dir)/(tag+"_native_candidate_distances.csv.gz"),index=False,encoding="utf-8")
            record=dict(dataset=Path(outpath).parent.name,method="scTenifoldXct",task="ordered_population_pair_LR_nomination",
                source_celltype=src,target_celltype=tgt,n_input_cells=n_global_cells,n_input_genes=n_input_genes,
                n_source_cells=int(adata.obs[gcol].astype(str).eq(src).sum()),n_target_cells=int(adata.obs[gcol].astype(str).eq(tgt).sum()),
                n_native_gene_vertices=len(xct._genes[src])+len(xct._genes[tgt]),n_native_candidates=len(aligned_candidates),
                n_significant_native_calls=len(df),n_training_steps=1000,learning_rate=0.01,nComp=5,n_embedding_dimensions=3,
                project_seed=CONFIG["seed"],native_pcnet_random_state=0,native_pcCoefficients_SVD_random_state=0,
                native_pvalue_cutoff=0.05,native_null_filter_zeros=True,score_direction="native_aligned_distance_lower_is_stronger",
                native_seconds=time.perf_counter()-pair_started,result_status="completed_native_ordered_pair",
                independent_molecular_truth_available=False)
            bundle=dict(significant_calls=df,aligned_candidate_distances=aligned_candidates,learned_embeddings=projections,
                source_gene_names=xct._genes[src],target_gene_names=xct._genes[tgt],
                source_cell_ids=xct._cell_data_dic[src].obs_names.to_numpy(),target_cell_ids=xct._cell_data_dic[tgt].obs_names.to_numpy(),
                coverage=record,numpy_rng_state=np.random.get_state(),python_rng_state=random.getstate(),torch_rng_state=torch.get_rng_state())
            temporary=bundle_path.with_suffix(".pkl.tmp")
            temporary.write_bytes(pickle.dumps(bundle,protocol=5))
            temporary.replace(bundle_path)
            all_results.append(df)
            coverage_records.append(record)
            pd.DataFrame(coverage_records).to_csv(Path(outpath).parent/"coverage.csv",index=False,encoding="utf-8")
            print("SCTENIFOLD_PAIR_COMPLETE",src,tgt,len(aligned_candidates),len(df),record["native_seconds"],flush=True)
            del xct,projections,aligned_candidates,bundle
            gc.collect()

    if not all_results:
        print("FAIL: no cell-type pair produced results")
        sys.exit(1)
    out = pd.concat(all_results, ignore_index=True)
    out["dataset"] = Path(outpath).parent.name
    out.to_csv(outpath, index=False,encoding="utf-8")
    print(f"OK: {len(out)} significant LR rows across {len(groups)} cell types -> {outpath}")


if __name__ == "__main__":
    main()
