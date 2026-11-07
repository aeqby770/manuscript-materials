import os
import sys
import argparse
from pathlib import Path

parser=argparse.ArgumentParser()
parser.add_argument("input_h5ad")
parser.add_argument("out_dir")
parser.add_argument("--resume",choices=("none","hodge","tree","full"),default="none")
parser.add_argument("--group-column",default="paul15_clusters")
args=parser.parse_args()
CONFIG={"paths":{"root":Path(os.environ["LLR_WORK_ROOT"]),"input":Path(args.input_h5ad),"out":Path(args.out_dir)},"threads":int(os.environ.get("LLR_NATIVE_THREADS","2")),"plot":{},"seed":12345,"n_top_genes":2000,"n_pcs":40,"n_diffusion_components":30,"knn":7,"diffusion_steps":1,"n_trajectories":10000,"trajectory_cluster_eps":0.5,"min_tree_bins":20,"cut_threshold":1,"clock_field":"S0_pseudotime","solver_beta":1e-6}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"]=str(CONFIG["threads"])
import pickle
import random
import time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG["plot"])
import numpy as np
import pandas as pd
import scanpy as sc
import networkx as nx
import scipy
from scipy import sparse
from sklearn.decomposition import PCA
sys.path[:0]=[str(CONFIG["paths"]["root"]/"tools/competitors/optimized/PHLOWER"),str(CONFIG["paths"]["root"]/"tools/competitors/dependencies/PHLOWER")]
import phlower as ph
from phlower.util import test_cholesky
from phlower.tools.hodgedecomp import div

np.random.seed(CONFIG["seed"])
random.seed(CONFIG["seed"])
out=CONFIG["paths"]["out"]
out.mkdir(parents=True,exist_ok=True)
nx.nx_pydot.graphviz_layout=nx.nx_agraph.graphviz_layout
started=time.perf_counter()
match args.resume:
    case "none":
        a=sc.read_h5ad(CONFIG["paths"]["input"])
        original_cells=a.obs_names.copy()
        original_genes=a.var_names.copy()
        a.X=sparse.csr_matrix(a.layers["lognorm"],dtype=np.float64)
        assert np.isfinite(a.X.data).all() and (a.X.data>=0).all()
        sc.pp.highly_variable_genes(a,n_top_genes=CONFIG["n_top_genes"])
        work=a[:,a.var.highly_variable].copy()
        work.X=np.asarray(work.X.toarray(),dtype=np.float64)
        pca=PCA(n_components=min(CONFIG["n_pcs"],work.n_obs-1,work.n_vars-1))
        work.obsm["X_pca"]=pca.fit_transform(work.X)
        work.obs["group"]=work.obs[args.group_column].astype(str)
        pc1=pd.Series(work.obsm["X_pca"][:,0],index=work.obs_names)
        root_group=pc1.groupby(work.obs["group"],observed=True).median().idxmin()
        roots=work.obs["group"].eq(root_group).to_numpy()
        work.uns["benchmark_input"]={"n_input_cells":len(original_cells),"n_input_genes":len(original_genes),"input_genes":original_genes.to_numpy(),"source_input":str(CONFIG["paths"]["input"]),"root_group":str(root_group),"root_input_mode":"lowest_median_PC1_provided_group;no_reference_orientation","linear_solver":"SuperLU_on_identical_author_A_plus_beta_I","solver_beta":CONFIG["solver_beta"],"layout_engine":"native_Graphviz_neato_via_pygraphviz_C_library","seed":CONFIG["seed"]}
        print("PHLOWER_STAGE full_ddhodge",work.shape,root_group,flush=True)
        ph.ext.ddhodge(work,basis="X_pca",roots=roots,k=CONFIG["knn"],npc=work.obsm["X_pca"].shape[1],ndc=min(CONFIG["n_diffusion_components"],work.obsm["X_pca"].shape[1]),s=CONFIG["diffusion_steps"],lstsq_method="cholesky",layout=None,verbose=True)
        assert work.obs_names.equals(original_cells)
        g=work.uns[work.uns["graph_basis"]]
        L=sparse.csc_matrix(nx.laplacian_matrix(g.to_undirected(),weight=None))
        A=L+CONFIG["solver_beta"]*sparse.eye(L.shape[0],format="csc")
        rhs=-np.asarray(div(g)).ravel()
        reference=scipy.linalg.cho_solve(scipy.linalg.cho_factor(A.toarray(),lower=True,check_finite=False),rhs,check_finite=False)
        direct=test_cholesky(L,beta=CONFIG["solver_beta"])(rhs)
        relative_error=float(np.max(np.abs(reference-direct))/max(1,float(np.max(np.abs(reference)))))
        relative_residual=float(np.max(np.abs(A@direct-rhs))/max(1,float(np.max(np.abs(rhs)))))
        assert relative_error<=1e-7 and relative_residual<=1e-7
        pd.DataFrame([{"n_full_input_cells":work.n_obs,"n_operator_rows":L.shape[0],"beta":CONFIG["solver_beta"],"operator":"same_native_cell_Laplacian_plus_beta_I","reference":"dense_Cholesky","backend":"sparse_SuperLU","relative_solution_error":relative_error,"relative_residual":relative_residual}]).to_csv(out/"solver_equivalence.csv",index=False,encoding="utf-8")
        checkpoint={"model":work,"numpy_rng_state":np.random.get_state(),"python_rng_state":random.getstate(),"hodge_seconds":time.perf_counter()-started}
        (out/"phlower_hodge_checkpoint.pkl").write_bytes(pickle.dumps(checkpoint,protocol=5))
    case "hodge":
        checkpoint=pickle.loads((out/"phlower_hodge_checkpoint.pkl").read_bytes())
        work=checkpoint["model"]
        np.random.set_state(checkpoint["numpy_rng_state"])
        random.setstate(checkpoint["python_rng_state"])
    case "tree":
        work=pickle.loads((out/"phlower_tree_checkpoint.pkl").read_bytes())
    case "full":
        work=pickle.loads((out/"phlower_native_result.pkl").read_bytes())
work.uns["benchmark_input"]["container_compatibility"]="NumPy_object_array_for_variable_length_native_roots"
tree_started=time.perf_counter()
match args.resume:
    case "none" | "hodge":
        g=work.uns[work.uns["graph_basis"]]
        assert set(g.nodes())==set(np.arange(work.n_obs))
        layout=nx.nx_agraph.graphviz_layout(g,prog="neato")
        work.obsm[work.uns["graph_basis"]]=pd.DataFrame.from_dict(layout,orient="index").reindex(np.arange(work.n_obs)).to_numpy(dtype=float)
        print("PHLOWER_STAGE author_Delaunay_and_Hodge_decomposition",flush=True)
        ph.tl.construct_delaunay(work,cluster_name="group",node_attr="u",start_n=10,end_n=10,circle_quant=0.1,calc_layout=True,random_seed=CONFIG["seed"])
        ph.tl.L1Norm_decomp(work)
        ph.tl.knee_eigen(work)
        print("PHLOWER_STAGE author_10000_trajectories",flush=True)
        ph.tl.random_climb_knn(work,n=CONFIG["n_trajectories"],seeds=CONFIG["seed"])
        ph.tl.trajs_matrix(work)
        ph.tl.trajs_clustering(work,eps=CONFIG["trajectory_cluster_eps"])
        print("PHLOWER_STAGE author_harmonic_tree_and_cell_projection",flush=True)
        ph.tl.harmonic_stream_tree(work,trajs_use=CONFIG["n_trajectories"],trajs_clusters="trajs_clusters",min_bin_number=CONFIG["min_tree_bins"],cut_threshold=CONFIG["cut_threshold"],layout_name="cumsum",eigen_n=-1,random_seed=CONFIG["seed"],verbose=True)
        repaired_count=len(work.uns["fate_tree"].graph.get("empty_suffix_sentinels_omitted",[]))
        work.uns["benchmark_input"]["n_empty_suffix_sentinels"]=repaired_count
        work.uns["benchmark_input"]["implementation_repair"]=np.where(repaired_count>0,"omit_unobserved_suffix_sentinel_leaves;all_supported_native_nodes_edges_counts_retained","none").item()
        work.uns["benchmark_input"]["native_seconds"]=checkpoint["hodge_seconds"]+time.perf_counter()-tree_started
        work.uns["benchmark_input"]["runtime_scope"]="sum_of_native_stage_elapsed_times_including_solver_equivalence_verification;imports_excluded"
        (out/"phlower_native_result.pkl").write_bytes(pickle.dumps(work,protocol=5))
    case "tree":
        from phlower.tools.harmonic_pseudo_tree import pairwise_hbranching_dict, create_branching_tree, create_detail_tree, add_node_pca, create_bstream_tree, add_pca_to_stream
        native_tree_before=work.uns["fate_tree"].copy()
        ddf=work.uns["ddf"]
        pairwise_bdict=pairwise_hbranching_dict(ddf,bottom_up=True,cut_threshold=CONFIG["cut_threshold"])
        htree,root=create_branching_tree(pairwise_bdict,keys=None)
        work.uns["fate_tree"]=create_detail_tree(work,htree,root,ddf,trim_end=False,layout_name="cumsum",cluster="group")
        supported={n for n,d in native_tree_before.nodes(data=True) if len(d["ecount"])>0}
        original_nodes={native_tree_before.nodes[n]["original"] for n in supported}
        assert original_nodes=={d["original"] for n,d in work.uns["fate_tree"].nodes(data=True)}
        original_edges={(native_tree_before.nodes[u]["original"],native_tree_before.nodes[v]["original"]) for u,v in native_tree_before.edges() if u in supported and v in supported}
        assert original_edges=={(work.uns["fate_tree"].nodes[u]["original"],work.uns["fate_tree"].nodes[v]["original"]) for u,v in work.uns["fate_tree"].edges()}
        counts_before={d["original"]:d["ecount"] for n,d in native_tree_before.nodes(data=True) if n in supported}
        assert all(d["ecount"]==counts_before[d["original"]] for n,d in work.uns["fate_tree"].nodes(data=True))
        pd.DataFrame([{"native_node":n,"trajectory_groups":str(d["original"][0]),"sentinel_time_bin":d["original"][1],"edge_weight_sum":sum(v for k,v in d["ecount"]),"native_parent":str(list(native_tree_before.predecessors(n))),"repair":"omit_only_unobserved_suffix_sentinel_leaf;all_supported_nodes_edges_and_counts_identical"} for n,d in native_tree_before.nodes(data=True) if n not in supported]).to_csv(out/"phlower_empty_suffix_repair.csv",index=False,encoding="utf-8")
        from phlower.tools.cumsum_utils import add_root_cumsum
        add_root_cumsum(work)
        add_node_pca(work,pca_name="X_pca")
        add_node_pca(work,pca_name=work.uns["graph_basis"])
        create_bstream_tree(work,layout_name="cumsum",iscopy=False)
        add_pca_to_stream(work,attr="cumsum")
        work.uns["benchmark_input"]["native_seconds"]=np.nan
        work.uns["benchmark_input"]["runtime_scope"]="tree_resumed_from_saved_native_failure;no_uninterrupted_runtime_claim"
        work.uns["benchmark_input"]["implementation_repair"]="omit_unobserved_suffix_sentinel_leaves;all_supported_native_nodes_edges_counts_retained"
        work.uns["benchmark_input"]["n_empty_suffix_sentinels"]=len(native_tree_before)-len(supported)
        (out/"phlower_native_result.pkl").write_bytes(pickle.dumps(work,protocol=5))
    case "full":
        pass
clock=work.obs[CONFIG["clock_field"]].to_numpy(dtype=float)
pred=work.obs.copy()
pred.insert(0,"cell_id",work.obs_names)
pred["pseudotime"]=clock
pred.to_csv(out/"phlower_cell_predictions.csv",index=False,encoding="utf-8")
pd.DataFrame({"cell_id":work.obs_names,"pseudotime":clock,"native_ddhodge_u":work.obs["u"].to_numpy(dtype=float)}).to_csv(out/"phlower_pseudotime.csv",index=False,encoding="utf-8")
tree=work.uns["fate_tree"]
nx.to_pandas_edgelist(tree).to_csv(out/"phlower_native_tree.csv",index=False,encoding="utf-8")
meta=work.uns["benchmark_input"]
pd.DataFrame([{"dataset":out.name,"method":"PHLOWER","task":"full_native_harmonic_tree_and_cell_projection","n_input_cells":meta["n_input_cells"],"n_input_genes":meta["n_input_genes"],"n_native_genes":work.n_vars,"n_finite_pseudotime":int(np.isfinite(clock).sum()),"n_native_tree_nodes":tree.number_of_nodes(),"n_native_tree_edges":tree.number_of_edges(),"n_native_tree_leaves":int(pd.Series(dict(tree.out_degree())).eq(0).sum()),"n_trajectories":len(work.uns["knn_trajs"]),"seed":CONFIG["seed"],"clock_field":CONFIG["clock_field"],"root_input_mode":meta["root_input_mode"],"direction_policy":"as_returned_no_truth_flip","native_seconds":meta["native_seconds"],"implementation_repair":meta.get("implementation_repair","none"),"n_empty_suffix_sentinels":meta.get("n_empty_suffix_sentinels",0),"runtime_scope":meta.get("runtime_scope",""),"result_status":"completed_full_native_workflow","linear_solver":meta["linear_solver"],"solver_beta":meta["solver_beta"],"container_compatibility":meta["container_compatibility"]}]).to_csv(out/"coverage.csv",index=False,encoding="utf-8")
check=pickle.loads((out/"phlower_native_result.pkl").read_bytes())
assert check.obs_names.equals(work.obs_names)
assert list(check.uns["fate_tree"].edges())==list(tree.edges())
assert np.array_equal(check.obs[CONFIG["clock_field"]].to_numpy(),clock,equal_nan=True)
print("PHLOWER_FULL_NATIVE_COMPLETE",work.n_obs,tree.number_of_nodes(),tree.number_of_edges(),flush=True)
