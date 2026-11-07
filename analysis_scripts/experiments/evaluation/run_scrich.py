import os
import sys
CONFIG={"paths":{"root":"D:/biotrade/LineaLinkR_0.6.5","dataset":sys.argv[1],"input":sys.argv[3],"similarity_scratch":".similarity_working.dat"},"threads":os.cpu_count(),"plot":{},"seed":12345,"species":sys.argv[2],"nsim":5,"kmin":3,"kmax":6,"n_pathways":10,"target":False,"layout":sys.argv[4]}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"]=str(CONFIG["threads"])
from pathlib import Path
from functools import partial
import importlib.util
import types
import time
import pickle
import ctypes
import struct
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from scipy.io import mmread
from scipy.sparse.linalg import eigsh
from scipy.sparse import csr_matrix
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG["plot"])
root=Path(CONFIG["paths"]["root"])
dataset=CONFIG["paths"]["dataset"]
source=Path(CONFIG["paths"]["input"])
out=root/"results/by_algorithm/scRICH"/dataset
out.mkdir(parents=True,exist_ok=True)
vendor=root/"tools/competitors/optimized/scRICH"
package=types.ModuleType("scrich_native_tools")
package.__path__=[str(vendor/"scrich/tools")]
sys.modules[package.__name__]=package
spec=importlib.util.spec_from_file_location("scrich_native_tools.mode_clustering_func",vendor/"scrich/tools/mode_clustering_func.py")
native=importlib.util.module_from_spec(spec)
sys.modules[spec.name]=native
spec.loader.exec_module(native)
sign_spec=importlib.util.spec_from_file_location("scrich_native_tools.signaling_func",vendor/"scrich/tools/signaling_func.py")
signaling=importlib.util.module_from_spec(sign_spec)
sys.modules[sign_spec.name]=signaling
sign_spec.loader.exec_module(signaling)
filenames={"prefixed":{"genes":dataset+"_genes.txt","cells":dataset+"_cells.txt","meta":dataset+"_meta.csv","matrix":dataset+"_matrix.mtx"},"plain":{"genes":"genes.txt","cells":"cells.txt","meta":"meta.csv","matrix":"expr.mtx"}}[CONFIG["layout"]]
genes=(source/filenames["genes"]).read_text(encoding="utf-8-sig").splitlines()
cells=(source/filenames["cells"]).read_text(encoding="utf-8-sig").splitlines()
obs=pd.read_csv(source/filenames["meta"],index_col=0)
obs=obs.loc[cells]
group_candidates=pd.Index(["cell_type","clusters","celltype","cluster","paul15_clusters","ct","state_info"])
group_column=group_candidates.intersection(obs.columns,sort=False)[0]
obs["scRICH_group"]=obs[group_column].astype(str)
expr=csr_matrix(mmread(source/filenames["matrix"]).T)
a=ad.AnnData(X=expr,obs=obs,var=pd.DataFrame(index=genes))
is_counts=bool(np.all(a.X.data>=0) and np.all(a.X.data==np.floor(a.X.data)))
scale=np.where(is_counts,10000/np.maximum(np.asarray(a.X.sum(axis=1)).ravel(),1),1.0)
a.X=a.X.multiply(scale[:,None]).tocsr()
a.X.data=np.where(is_counts,np.log1p(a.X.data),a.X.data)
normalization={True:"library_size_10000_log1p",False:"provided_expression"}[is_counts]
np.random.seed(CONFIG["seed"])
solver_start=np.random.default_rng(CONFIG["seed"]).standard_normal(a.n_obs)
native.eigsh=partial(eigsh,v0=solver_start)
os.chdir(vendor)
before=time.perf_counter()
native.pathways_overview(a,human=CONFIG["species"]=="human")
overview=pd.DataFrame.from_dict(a.uns["paths_overview"]["gene_number"],orient="index",columns=["n_receptors","n_ligands"])
overview.to_csv(out/"pathway_input_coverage.csv",index_label="pathway",encoding="utf-8")
memory_status=ctypes.create_string_buffer(struct.calcsize("<II7Q"))
struct.pack_into("<I",memory_status,0,len(memory_status))
assert ctypes.windll.kernel32.GlobalMemoryStatusEx(memory_status)
memory_values=struct.unpack("<II7Q",memory_status.raw)
available_matrix_bytes=min(memory_values[3],memory_values[5])
count_bytes=np.min_scalar_type(CONFIG["nsim"]*(CONFIG["kmax"]-CONFIG["kmin"]+1)).itemsize
matrix_peak_bytes=(np.dtype(np.float64).itemsize+count_bytes)*a.n_obs**2
memory_map={True:str(out/CONFIG["paths"]["similarity_scratch"]),False:None}[matrix_peak_bytes>available_matrix_bytes]
similarity_storage={True:"file_backed_float64",False:"memory"}[memory_map is not None]
print({"dataset":dataset,"n_cells":a.n_obs,"n_genes":a.n_vars,"normalization":normalization,"group_column":group_column,"eligible_pathways":len(overview),"similarity_storage":similarity_storage,"matrix_peak_bytes":matrix_peak_bytes,"available_matrix_bytes":available_matrix_bytes},flush=True)
native.all_path_sim(a,annotations="scRICH_group",nsim=CONFIG["nsim"],kmin=CONFIG["kmin"],kmax=CONFIG["kmax"],target=CONFIG["target"],human=CONFIG["species"]=="human",pathways_to_use=CONFIG["n_pathways"],seed=CONFIG["seed"],print_info=False,memory_map=memory_map)
list(map(partial(signaling.compute_ccc_matrix,a,key="scRICH_group",include_target=False,use_velocity=False,perm_test=False),a.uns["selected_pathways"]))
seconds=time.perf_counter()-before
a.obs.to_csv(out/"cell_modes.csv.gz",encoding="utf-8")
spectrum=pd.DataFrame.from_dict(a.uns["SymNMF"],orient="index").rename_axis("pathway").reset_index()
path_results=pd.DataFrame.from_dict(a.uns["pathways"],orient="index").rename_axis("pathway").reset_index()
summary=spectrum.merge(path_results,on="pathway",how="left")
summary.insert(0,"dataset",dataset)
summary.to_csv(out/"pathway_heterogeneity.csv",index=False,encoding="utf-8")
with (out/"scrich_native.pkl").open("wb") as stream:
    pickle.dump({"obs":a.obs,"uns":a.uns,"config":CONFIG},stream,protocol=5)
coverage=pd.DataFrame([{"dataset":dataset,"method":"scRICH_LR","n_input_cells":a.n_obs,"n_input_genes":a.n_vars,"n_eligible_pathways":len(overview),"n_analyzed_pathways":len(a.uns["selected_pathways"]),"pathway_selection":"native_top10_mean_expression","downstream_target":False,"normalization":normalization,"group_column":group_column,"core_seconds":seconds,"similarity_storage":similarity_storage,"matrix_peak_bytes":matrix_peak_bytes,"result_status":"scored","ground_truth_accuracy_available":False}])
coverage.to_csv(out/"coverage.csv",index=False,encoding="utf-8")
print(coverage.to_string(index=False),flush=True)

