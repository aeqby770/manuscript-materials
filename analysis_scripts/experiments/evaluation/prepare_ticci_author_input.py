import os,sys
CONFIG={"paths":{"root":sys.argv[1],"manifest":"experiments/evaluation/ticci_author_protocol_inputs.csv","output":"results/by_algorithm/TICCI_author_functional_entropy","metadata":"inputs/comparison/recovered_native/input_metadata"},"threads":int(os.environ.get("LLR_NATIVE_THREADS",os.cpu_count())),"plot":{}}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"]=str(CONFIG["threads"])
from pathlib import Path
from itertools import repeat
import pandas as pd,numpy as np,scanpy as sc,h5py,shutil
from scipy import sparse
root=Path(CONFIG["paths"]["root"]);dataset=sys.argv[2]
job=pd.read_csv(root/CONFIG["paths"]["manifest"]).set_index("dataset").loc[dataset]
out=root/CONFIG["paths"]["output"]/dataset
out.mkdir(parents=True,exist_ok=True)
adata=sc.read_h5ad(root/job.input_h5ad)
original_cells=pd.Index((root/CONFIG["paths"]["metadata"]/dataset/"cells.txt").read_text(encoding="utf-8-sig").splitlines())
original_genes=pd.Index((root/CONFIG["paths"]["metadata"]/dataset/"genes.txt").read_text(encoding="utf-8-sig").splitlines())
assert adata.obs_names.is_unique and adata.var_names.is_unique
assert set(adata.obs_names)==set(original_cells)
assert set(adata.var_names)==set(original_genes)
matrix=sparse.csr_matrix(adata.X,dtype=np.float64)
matrix.sum_duplicates();matrix.sort_indices()
assert np.isfinite(matrix.data).all() and (matrix.data>=0).all()
with h5py.File(out/"entropy_input.h5","w") as h:
 h.create_dataset("data",data=matrix.data,compression="gzip",compression_opts=1)
 h.create_dataset("indices",data=matrix.indices,compression="gzip",compression_opts=1)
 h.create_dataset("indptr",data=matrix.indptr,compression="gzip",compression_opts=1)
 h.create_dataset("shape",data=np.asarray(matrix.shape,dtype=np.int64))
pd.Series(adata.var_names).to_csv(out/"entropy_genes.txt",index=False,header=False,encoding="utf-8")
pd.Series(adata.obs_names).to_csv(out/"entropy_cells.txt",index=False,header=False,encoding="utf-8")
group_candidates=pd.Index(["clusters_enlarged","clusters","celltype","cell_type","celltype.mapped","Cell_type","paul15_clusters","ct","ct_subtype","state_info"])
group_columns=group_candidates[group_candidates.isin(adata.obs.columns)]
provided_groups=adata.obs.loc[:,group_columns].astype(str).rename_axis("cell_id").reset_index().melt(id_vars="cell_id",value_name="Classification").drop_duplicates("cell_id",keep="first")
provided_groups={"single_unannotated_RNA_population":pd.DataFrame({"cell_id":adata.obs_names,"Classification":"RNA"})}.get(job.grouping_mode,provided_groups)
provided_groups[["cell_id","Classification"]].to_csv(out/"ticci_cellchat_groups.csv",index=False,encoding="utf-8")
binding=pd.read_csv(root/"results/evaluation/final065/saved_native/native_TICCI_implementation_binding.csv")
old_paths=binding.loc[binding.dataset.eq(dataset)&(job.cellchat_mode=="reuse_saved_calls"),"source_file"].dropna().map(Path).map(lambda p:p.parent).drop_duplicates()
list(map(shutil.copy2,old_paths.map(lambda p:p/"ticci_cellchat_interactions.csv"),repeat(out/"ticci_cellchat_interactions.csv")))
meta=pd.DataFrame([{"dataset":dataset,"species":job.species,"input_h5ad":str(root/job.input_h5ad),"n_cells":adata.n_obs,"n_genes":adata.n_vars,"original_cell_ids_exact":True,"original_gene_ids_exact":True,"entropy_input":"original_benchmark_X;full_cells;all_original_genes","matrix_format":"HDF5_CSR_transposed_as_R_CSC;float64_lossless","CellChat_stage":str(job.cellchat_mode),"classification_source":str(job.grouping_mode),"upstream_comparison":"same_fixed_CCI_weight_1;resolution_1;n_neighbors_10;2000_variable_genes","old_result_directory":";".join(old_paths.astype(str)),"bootstrap_B":100,"bootstrap_genes":1000,"seed":12345}])
meta["entropy_mode"]=job.get("entropy_mode","go_kappa")
meta.to_csv(out/"input_binding.csv",index=False,encoding="utf-8")
pd.DataFrame([{"dataset":dataset,"species":job.species,"input_h5ad":job.input_h5ad,"status":"prepared_full_input","started":pd.Timestamp.now(tz="Asia/Shanghai").strftime("%Y-%m-%d %H:%M:%S")}]).to_csv(out/"run_status.csv",index=False,encoding="utf-8")
print(meta.to_string(index=False),flush=True)
