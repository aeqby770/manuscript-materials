import os
import sys
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("dataset")
parser.add_argument("--mode",choices=("train","reuse"),default="train")
parser.add_argument("--training-seconds",type=float,default=0.0)
args=parser.parse_args()
CONFIG={"paths":{"root":"D:/biotrade/LineaLinkR_0.6.5","inputs":"inputs/benchmark/curated","input_template":"{dataset}_benchmark.h5ad","source":"tools/competitors/optimized/CellAgentChat/src","dependencies":"tools/competitors/dependencies/CellAgentChat","outputs":"results/by_algorithm/CellAgentChat","lr_prior":"mouse_lr_pair.tsv","model":"native_model.pt","rates":"conversion_rates.txt","background":"background_lr_scores.pkl","distribution":"distribution.csv","native":"cellagentchat_native.pkl","scores":"lr_scores.csv.gz","coverage":"coverage.csv"},"threads":os.cpu_count(),"plot":{},"dataset":args.dataset,"mode":args.mode,"completed_training_seconds":args.training_seconds,"species":"mouse","group_column":"clusters","normalization_total":10000,"seed":12345,"epochs":50,"batch":256,"learning_rate":0.1,"n_shuffles":10,"background_scores":10000,"max_steps":1,"significance_threshold":0.05,"device":"cpu"}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"]=str(CONFIG["threads"])
from pathlib import Path
import random
import time
import numpy as np
import pandas as pd
import scanpy as sc
import torch
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG["plot"])
root=Path(CONFIG["paths"]["root"])
sys.path[:0]=[str(root/CONFIG["paths"]["dependencies"]),str(root/CONFIG["paths"]["source"])]
import model_setup as native_setup
from preprocessor import setup_adata
from permutations import permutation_test
from bckground_distribution import get_distribution, save_distribution
from Communication import CCI, receiving_score
dataset=CONFIG["dataset"]
out=root/CONFIG["paths"]["outputs"]/dataset
out.mkdir(parents=True,exist_ok=True)
source=root/CONFIG["paths"]["inputs"]/CONFIG["paths"]["input_template"].format(dataset=dataset)
a=sc.read_h5ad(source)
sc.pp.normalize_total(a,target_sum=CONFIG["normalization_total"])
sc.pp.log1p(a)
a=setup_adata(a,coordinates_key=None,batch_label=None,cell_type_label=CONFIG["group_column"],copy=False)
random.seed(CONFIG["seed"])
np.random.seed(CONFIG["seed"])
torch.manual_seed(CONFIG["seed"])
torch.set_num_threads(CONFIG["threads"])
lig_uni,rec_uni,lr_pairs=native_setup.load_db(a,file=CONFIG["paths"]["lr_prior"],sep="\t")
tf_uni,rec_tf_uni=native_setup.load_tf_db(CONFIG["species"],a,rec_uni)
coverage=pd.DataFrame([{"dataset":dataset,"method":"CellAgentChat","n_input_cells":a.n_obs,"n_input_genes":a.n_vars,"species":CONFIG["species"],"n_cell_types":a.obs.cell_type.nunique(),"n_native_lr_pairs":len(lr_pairs),"n_native_receptors":len(rec_uni),"n_native_tfs":len(tf_uni),"normalization":"library_size_10000_log1p","epochs":CONFIG["epochs"],"n_feature_shuffles":CONFIG["n_shuffles"],"background_score_target":CONFIG["background_scores"],"distance_mode":False,"native_fit_started":True,"resumed_from_saved_model":CONFIG["mode"]=="reuse","result_status":"running_native","ground_truth_accuracy_available":False}])
coverage.to_csv(out/CONFIG["paths"]["coverage"],index=False,encoding="utf-8")
print(coverage.to_string(index=False),flush=True)
before=time.perf_counter()
print({"dataset":dataset,"stage":{"train":"native_neural_network","reuse":"reuse_saved_neural_network"}[CONFIG["mode"]]},flush=True)
prepare,prepare_args,prepare_kwargs={
    "train":(native_setup.train,(a,lig_uni,rec_uni,tf_uni,rec_tf_uni,lr_pairs),{"path":str(out/CONFIG["paths"]["model"]),"epochs":CONFIG["epochs"],"lr":CONFIG["learning_rate"],"batch":CONFIG["batch"],"device":CONFIG["device"]}),
    "reuse":(native_setup._build_io_matrices,(a,lig_uni,rec_uni),{})
}[CONFIG["mode"]]
mat,C=prepare(*prepare_args,**prepare_kwargs)[:2]
training_seconds={"train":time.perf_counter()-before,"reuse":CONFIG["completed_training_seconds"]}[CONFIG["mode"]]
model_dl=native_setup.load_model(str(out/CONFIG["paths"]["model"]),device=CONFIG["device"])
print({"dataset":dataset,"stage":"native_receptor_feature_selection","training_seconds":training_seconds},flush=True)
feature_start=time.perf_counter()
conversion_rates=native_setup.feature_selection(model_dl,mat,C,rec_uni,n_shuffles=CONFIG["n_shuffles"],batch=CONFIG["batch"],seed=CONFIG["seed"],device=CONFIG["device"])
assert np.isfinite(conversion_rates).all()
native_setup.save_conversion_rate(conversion_rates,file=str(out/CONFIG["paths"]["rates"]))
rates=native_setup.add_rates(conversion_rates,rec_uni)
feature_seconds=time.perf_counter()-feature_start
del mat,C
print({"dataset":dataset,"stage":"native_background_distribution","feature_seconds":feature_seconds,"min_conversion_rate":float(np.min(conversion_rates)),"max_conversion_rate":float(np.max(conversion_rates))},flush=True)
background_start=time.perf_counter()
test_scores,background_model,distance=permutation_test(CONFIG["background_scores"],a.n_obs,a,lig_uni,rec_uni,rates,dist=False,seed=CONFIG["seed"])
pd.to_pickle(test_scores,out/CONFIG["paths"]["background"],protocol=5)
del background_model
params=get_distribution(test_scores,dist=distance,scaled=False)
save_distribution(params,path=str(out/CONFIG["paths"]["distribution"]))
background_seconds=time.perf_counter()-background_start
print({"dataset":dataset,"stage":"native_CCI","background_seconds":background_seconds},flush=True)
network_start=time.perf_counter()
clusters=a.obs.cell_type.astype(str).unique().tolist()
model=CCI(a.n_obs,a,lig_uni,rec_uni,rates,params,clusters,dist=False,max_steps=CONFIG["max_steps"],plot=False,path=str(out),threshold=CONFIG["significance_threshold"],net=model_dl)
receiving=receiving_score(model,path=str(out))
native_results={"obs":a.obs,"results":model.results,"results2":model.results2,"receiving":receiving,"lig_uni":lig_uni,"rec_uni":rec_uni,"conversion_rates":rates,"distribution":params,"config":CONFIG}
pd.to_pickle(native_results,out/CONFIG["paths"]["native"],protocol=5)
scores=pd.DataFrame(model.results2).rename_axis("lr").rename_axis("cell_pair",axis=1).stack().rename("score").reset_index()
cell_types=pd.DataFrame({"sender":clusters})
cell_pairs=cell_types.merge(cell_types.rename(columns={"sender":"receiver"}),how="cross")
cell_pairs["cell_pair"]=cell_pairs.sender+"_"+cell_pairs.receiver
lr_map=pd.Series(lig_uni).explode().rename("receptor").rename_axis("ligand").reset_index().drop_duplicates()
lr_map["lr"]=lr_map.ligand+"-"+lr_map.receptor
scores=scores.merge(cell_pairs,on="cell_pair",how="left",validate="many_to_one").merge(lr_map,on="lr",how="left",validate="many_to_one")
scores=scores[["sender","receiver","ligand","receptor","score"]]
assert scores[["sender","receiver","ligand","receptor"]].notna().all().all()
assert np.isfinite(scores.score).all()
scores.to_csv(out/CONFIG["paths"]["scores"],index=False,encoding="utf-8",compression="gzip")
coverage["training_seconds"]=training_seconds
coverage["feature_seconds"]=feature_seconds
coverage["background_seconds"]=background_seconds
coverage["network_seconds"]=time.perf_counter()-network_start
coverage["core_seconds"]=time.perf_counter()-before+{"train":0.0,"reuse":training_seconds}[CONFIG["mode"]]
coverage["n_scored_lr_keys"]=len(scores)
coverage["result_status"]="scored"
coverage.to_csv(out/CONFIG["paths"]["coverage"],index=False,encoding="utf-8")
print(coverage.to_string(index=False),flush=True)


