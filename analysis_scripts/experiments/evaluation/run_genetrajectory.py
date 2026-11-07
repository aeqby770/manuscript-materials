#!/usr/bin/env python
"""GeneTrajectory adapter for canonical h5ad input.

Usage: python run_genetrajectory.py <canonical.h5ad> <output_gene_trajectory.csv>
"""
import os
import sys
import argparse
import time
import importlib.util
import pandas as pd
from pathlib import Path

import numpy as np
import scanpy as sc
from scipy import sparse

parser=argparse.ArgumentParser()
parser.add_argument("input_h5ad")
parser.add_argument("output_csv")
parser.add_argument("--verify-existing-result",default="")
parser.add_argument("--verify-current-graph",action="store_true")
parser.add_argument("--ot-representation",choices=("dense","graph-flow"),default="dense")
args=parser.parse_args()

CONFIG = {
    "paths": {
        "vendor": Path(os.environ["LLR_WORK_ROOT"]).resolve() / "tools/competitors/optimized/GeneTrajectory",
    },
    "threads": int(os.environ.get("LLR_NATIVE_THREADS","4")),
    "seed": 12345,
    "n_pcs": 30,
    "graph_k": 10,
    "graph_dims": 5,
    "n_genes": 2000,
    "embedding_k": 10,
    "embedding_eigs": 30,
    "n_trajectories": 2,
    "diffusion_steps": 3,
}

os.environ["OMP_NUM_THREADS"] = str(CONFIG["threads"])
np.random.seed(CONFIG["seed"])
sys.path.insert(0, str(CONFIG["paths"]["vendor"]))
sys.path.insert(0,str(Path(os.environ["LLR_WORK_ROOT"])/"tools/competitors/dependencies/wot"))

from gene_trajectory.extract_gene_trajectory import extract_gene_trajectory, get_gene_embedding  # noqa: E402
from gene_trajectory.gene_distance_shared import cal_ot_mat  # noqa: E402
from gene_trajectory.get_graph_distance import get_graph_distance  # noqa: E402


def main():
    inpath, outpath = args.input_h5ad, args.output_csv
    out=Path(outpath).parent
    out.mkdir(parents=True,exist_ok=True)
    started=time.perf_counter()
    adata = sc.read_h5ad(inpath)
    n_input_cells,n_input_genes=adata.shape

    # No subsampling: GeneTrajectory runs on the full cell set. The cell-cell graph
    # distance and the gene-gene EMD (cal_ot_mat) are parallelised across
    # CONFIG["threads"] rather than reduced in size.

    if "lognorm" in adata.layers:
        adata.X = adata.layers["lognorm"].copy()
    elif "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    # else: use adata.X as-is

    # NaN filter on sparse data BEFORE any dense conversion
    if sparse.issparse(adata.X):
        assert np.isfinite(adata.X.data).all() and (adata.X.data>=0).all()
    else:
        assert np.isfinite(adata.X).all() and (adata.X>=0).all()

    # Rank genes by sparse variance without densifying the full cell-by-gene matrix.
    # This avoids non-unique expression bins in the simulated datasets.
    n_top = min(3000, adata.n_vars - 1)
    if sparse.issparse(adata.X):
        # Var(X) = E[X^2] - (E[X])^2, works on sparse without per-column .data.var()
        x_csc = adata.X.tocsc()
        n = adata.n_obs
        mean_sq = np.array(x_csc.multiply(x_csc).sum(axis=0)).ravel() / n
        mean_x = np.array(x_csc.sum(axis=0)).ravel() / n
        gene_var = mean_sq - mean_x ** 2
        gene_var = np.nan_to_num(gene_var, nan=0.0)
    else:
        gene_var = np.asarray(adata.X.var(axis=0))
    top_var_idx = np.argsort(gene_var)[::-1][:n_top]
    adata = adata[:, top_var_idx].copy()
    adata.X = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)

    n_pcs = min(CONFIG["n_pcs"], adata.n_obs - 1, adata.n_vars - 1)
    sc.pp.pca(adata, n_comps=n_pcs, random_state=CONFIG["seed"])

    x = adata.X
    gene_sum = x.sum(axis=0)
    gene_var = x.var(axis=0)
    valid_idx = np.where(gene_sum > 0)[0]
    if valid_idx.size < 3:
        print("GeneTrajectory requires at least 3 expressed genes", file=sys.stderr)
        sys.exit(1)

    ranked_valid = valid_idx[np.argsort(gene_var[valid_idx])[::-1]]
    feature_idx = ranked_valid[: min(CONFIG["n_genes"], ranked_valid.size)]
    feature_names = adata.var_names.astype(str).to_numpy()[feature_idx]
    if args.verify_existing_result:
        feature_names=pd.read_csv(args.verify_existing_result)["gene_id"].astype(str).to_numpy()
        feature_idx=adata.var_names.get_indexer(feature_names)
        assert (feature_idx>=0).all()
    gene_expr = x[:, feature_idx].astype(np.float64)

    # Ensure the kNN graph is connected (hard requirement of get_graph_distance):
    # pick the smallest k (from CONFIG k, doubling) that yields a single connected
    # component, decided by an EXPLICIT connected-components check rather than an
    # exception-driven retry.
    from sklearn.neighbors import kneighbors_graph
    from scipy.sparse.csgraph import connected_components
    dims_use = min(CONFIG["graph_dims"], n_pcs)
    pca_emb = adata.obsm["X_pca"][:, :dims_use]
    k = min(CONFIG["graph_k"], adata.n_obs - 1)
    while k < adata.n_obs - 1:
        knn = kneighbors_graph(pca_emb, n_neighbors=k, mode="connectivity")
        if connected_components(knn, directed=False, return_labels=False) == 1:
            break
        k = min(k * 2, adata.n_obs - 1)
    graph_dist = get_graph_distance(
        adata,
        reduction="X_pca",
        k=k,
        dims=dims_use,
    )
    if args.verify_current_graph:
        import ot
        from itertools import starmap
        from functools import partial
        root=Path(os.environ["LLR_WORK_ROOT"])
        spec=importlib.util.spec_from_file_location("gene_trajectory_reference_kernel",
            root/"tools/competitors/original/GeneTrajectory/gene_trajectory/gene_distance_shared.py")
        reference=importlib.util.module_from_spec(spec);spec.loader.exec_module(reference)
        graph_i,graph_j=np.nonzero(graph_dist==1)
        n=adata.n_obs
        M=sparse.coo_matrix((np.r_[np.ones(len(graph_i)),np.zeros(n)],
            (np.r_[graph_i,np.arange(n)],np.r_[graph_j,np.arange(n)])),shape=(n,n))
        probability=gene_expr/np.asarray(tuple(map(sum,gene_expr.T)))
        lifted=probability+1.0
        t0=time.perf_counter()
        dense=list(starmap(partial(reference._cal_ot,graph_dist,gene_expr,50000),[(0,1),(0,2)]))
        dense_seconds=time.perf_counter()-t0
        t0=time.perf_counter()
        lifted_results=list(starmap(partial(ot.emd2,M=M,numItermax=50000,log=True),
            [(lifted[:,0],lifted[:,1]),(lifted[:,0],lifted[:,2])]))
        sparse_seconds=time.perf_counter()-t0
        dense_values=np.asarray(tuple(map(lambda x:x[0],dense)))
        sparse_values=np.asarray(tuple(map(lambda x:float(x[0]),lifted_results)))
        warnings=tuple(map(lambda x:x[1]["warning"],lifted_results))
        error=float(np.max(np.abs(dense_values-sparse_values)))
        pd.DataFrame([dict(n_full_input_cells=n,n_full_model_genes=len(feature_names),
            validation_unit="first_two_native_gene_pair_transport_problems_on_complete_input",
            n_complete_pair_operators=2,graph_k_used=k,n_graph_arcs=M.nnz,emd_max_iter=50000,
            dense_seconds=dense_seconds,sparse_seconds=sparse_seconds,max_abs_error=error,
            n_nonoptimal_sparse_solves=sum(map(lambda x:x is not None,warnings)))]).to_csv(outpath,index=False,encoding="utf-8")
        np.savez_compressed(out/"native_full_graph_preprocessing.npz",graph_cost=graph_dist,gene_expression=gene_expr,
            gene_names=feature_names.astype(str),cell_ids=adata.obs_names.to_numpy().astype(str))
        print("CURRENT_FULL_GRAPH_OT_EQUIVALENCE",n,k,M.nnz,error,dense_seconds,sparse_seconds,warnings,flush=True)
        assert error<=1e-9 and all(map(lambda x:x is None,warnings))
        return
    if args.verify_existing_result:
        root=Path(os.environ["LLR_WORK_ROOT"])
        original=root/"tools/competitors/original/GeneTrajectory/gene_trajectory/gene_distance_shared.py"
        spec=importlib.util.spec_from_file_location("gene_trajectory_reference_kernel",original)
        reference=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(reference)
        t0=time.perf_counter()
        native=reference.cal_ot_mat(ot_cost=graph_dist,gene_expr=gene_expr,
            show_progress_bar=False,processes=float(CONFIG["threads"]))
        native_seconds=time.perf_counter()-t0
        t0=time.perf_counter()
        fast=cal_ot_mat(ot_cost=graph_dist,gene_expr=gene_expr,
            show_progress_bar=False,processes=float(CONFIG["threads"]))
        fast_seconds=time.perf_counter()-t0
        error=float(np.max(np.abs(native-fast)))
        ranks_equal=np.array_equal(np.argsort(native,axis=1,kind="stable"),np.argsort(fast,axis=1,kind="stable"))
        assert error<=1e-12 and ranks_equal
        pd.DataFrame([dict(input_h5ad=inpath,existing_native_gene_result=args.verify_existing_result,
            n_full_input_cells=n_input_cells,n_existing_native_genes=len(feature_names),
            n_complete_gene_pairs=len(feature_names)*(len(feature_names)-1)//2,native_emd_max_iter=50000,
            native_seconds=native_seconds,optimized_seconds=fast_seconds,max_absolute_error=error,
            same_gene_distance_ranks=ranks_equal,verification_scope="complete_existing_native_16_gene_object;all_original_cells;no_new_sampling")]).to_csv(outpath,index=False,encoding="utf-8")
        import ot
        from itertools import starmap
        from functools import partial
        ii,jj=np.nonzero(graph_dist==1)
        nn=graph_dist.shape[0]
        sparse_cost=sparse.coo_matrix((np.r_[np.ones(len(ii)),np.zeros(nn)],
            (np.r_[ii,np.arange(nn)],np.r_[jj,np.arange(nn)])),shape=(nn,nn))
        normalized=gene_expr/np.asarray(tuple(map(sum,gene_expr.T)))
        lifted=normalized+1.0
        gi,gj=np.triu_indices(len(feature_names),1)
        t0=time.perf_counter()
        solved=list(starmap(partial(ot.emd2,M=sparse_cost,numItermax=50000,log=True),
            zip(lifted[:,gi].T,lifted[:,gj].T)))
        sparse_seconds=time.perf_counter()-t0
        sparse_values=np.asarray(tuple(map(lambda v:float(v[0]),solved)))
        sparse_dist=np.zeros_like(native)
        sparse_dist[gi,gj]=sparse_values
        sparse_dist[gj,gi]=sparse_values
        sparse_error=float(np.max(np.abs(native-sparse_dist)))
        sparse_rank_equal=np.array_equal(np.argsort(native,axis=1,kind="stable"),np.argsort(sparse_dist,axis=1,kind="stable"))
        solver_warnings=tuple(map(lambda v:v[1]["warning"],solved))
        pd.DataFrame([dict(n_full_input_cells=nn,n_existing_native_genes=len(feature_names),
            n_complete_gene_pairs=len(gi),n_sparse_solver_arcs=sparse_cost.nnz,
            emd_max_iter=50000,n_nonoptimal_native_sparse_solves=sum(map(lambda v:v is not None,solver_warnings)),
            dense_native_seconds=native_seconds,shared_dense_seconds=fast_seconds,sparse_graph_seconds=sparse_seconds,
            max_absolute_cost_error=sparse_error,same_gene_distance_ranks=sparse_rank_equal,
            representation="same_unweighted_cell_graph_W1;common_mass_1_bipartite_lift;explicit_zero_cost_self_arcs",
            kernel_scope="all_original_cells;all_pairs_of_existing_16_gene_object")]).to_csv(Path(outpath).with_name("genetrajectory_sparse_graph_equivalence.csv"),index=False,encoding="utf-8")
        print("GENETRAJECTORY_SPARSE_GRAPH_EQUIVALENCE",sparse_error,sparse_rank_equal,sparse_seconds,solver_warnings[:2],flush=True)
        assert sparse_error<=1e-10 and sparse_rank_equal and all(map(lambda v:v is None,solver_warnings))
        print("GENETRAJECTORY_KERNEL_EQUIVALENCE",error,ranks_equal,native_seconds,fast_seconds,flush=True)
        return
    cache=out/"native_gene_transport_distances.npz"
    if cache.is_file():
        saved=np.load(cache,allow_pickle=False)
        np.testing.assert_array_equal(saved["gene_names"],feature_names)
        np.testing.assert_array_equal(saved["cell_ids"],adata.obs_names.to_numpy())
        gene_dist=saved["gene_distance"]
    else:
        print("GENETRAJECTORY_STAGE complete_gene_OT",adata.n_obs,len(feature_names),flush=True)
        gene_dist=cal_ot_mat(ot_cost=graph_dist,gene_expr=gene_expr,
            show_progress_bar=True,processes=float(CONFIG["threads"]),
            metric_representation=args.ot_representation,
            checkpoint_path=str(out/"native_gene_transport_partial.npy"),
            status_path=str(out/"native_gene_transport_solver_status.npy"))
        np.savez_compressed(cache,gene_names=feature_names.astype(str),cell_ids=adata.obs_names.to_numpy().astype(str),gene_distance=gene_dist)
    assert np.isfinite(gene_dist).all()
    gene_embedding, eigenvalues = get_gene_embedding(
        gene_dist,
        k=min(CONFIG["embedding_k"], gene_dist.shape[0] - 1),
        n_ev=min(CONFIG["embedding_eigs"], gene_dist.shape[0] - 1),
    )

    result = extract_gene_trajectory(
        gene_embedding=gene_embedding,
        dist_mat=gene_dist,
        gene_names=feature_names.tolist(),
        t_list=np.repeat(CONFIG["diffusion_steps"], CONFIG["n_trajectories"]),
        n=CONFIG["n_trajectories"],
        dims=min(5, gene_embedding.shape[1]),
        k=min(CONFIG["embedding_k"], gene_dist.shape[0] - 1),
    )
    result.insert(0, "gene_id", result.index.astype(str))
    result.insert(1, "gene_subset_n", len(feature_names))
    result.to_csv(outpath, index=False,encoding="utf-8")
    np.savez_compressed(out/"native_gene_embedding.npz",gene_names=feature_names.astype(str),embedding=gene_embedding,eigenvalues=eigenvalues)
    pd.DataFrame([dict(dataset=out.name,method="GeneTrajectory",task="gene_ordering",
        n_input_cells=n_input_cells,n_input_genes=n_input_genes,n_graph_genes=adata.n_vars,
        n_native_genes=len(feature_names),n_native_gene_pairs=len(feature_names)*(len(feature_names)-1)//2,
        ot_representation=args.ot_representation,
        n_nonoptimal_OT_pairs=int((np.load(out/"native_gene_transport_solver_status.npy",mmap_mode="r")[np.triu_indices(len(feature_names),1)]!=1).sum()),
        n_returned_genes=len(result),n_requested_trajectories=CONFIG["n_trajectories"],graph_k=k,
        graph_dimensions=dims_use,emd_max_iter=50000,seed=CONFIG["seed"],
        native_seconds=time.perf_counter()-started,result_status="completed_native_gene_trajectory",
        independent_gene_onset_truth_available=False,source_file=str(outpath))]).to_csv(out/"coverage.csv",index=False,encoding="utf-8")

    print(f"GeneTrajectory OK: {len(result)} genes -> {outpath}")


if __name__ == "__main__":
    main()
