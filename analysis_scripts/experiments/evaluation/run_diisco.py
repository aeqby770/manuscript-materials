import os
import sys
import argparse
parser=argparse.ArgumentParser()
parser.add_argument("dataset")
parser.add_argument("metadata")
parser.add_argument("time_column")
parser.add_argument("sample_column")
parser.add_argument("group_column")
parser.add_argument("--elbo",choices=("native","analytic-kl"),default="native")
parser.add_argument("--unit",default="")
parser.add_argument("--stratum-column",default="")
parser.add_argument("--stratum-value",default="")
parser.add_argument("--exclude-stratum",action="store_true")
parser.add_argument("--initial-status",default="failed_unstable_optimization")
args=parser.parse_args()
CONFIG={"paths":{"root":os.environ["LLR_WORK_ROOT"],"dataset":args.dataset,"metadata":args.metadata},"threads":int(os.environ.get("LLR_NATIVE_THREADS","1")),"plot":{},"seed":12345,"time_column":args.time_column,"sample_column":args.sample_column,"group_column":args.group_column,"n_iter":10000,"learning_rate":0.00005,"guide":"MultivariateNormalFactorized","prior":"uninformative_off_diagonal","elbo_estimator":args.elbo,"unit":args.unit,"stratum_column":args.stratum_column,"stratum_value":args.stratum_value,"exclude_stratum":args.exclude_stratum,"initial_status":args.initial_status,"implementation":{"native":"original/DIISCO","analytic-kl":"optimized/DIISCO"}[args.elbo]}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"]=str(CONFIG["threads"])
from pathlib import Path
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG["plot"])
root=Path(CONFIG["paths"]["root"])
base_out=root/"results/by_algorithm/DIISCO"/CONFIG["paths"]["dataset"]
out=base_out/CONFIG["unit"]
out.mkdir(parents=True,exist_ok=True)
meta=pd.read_csv(CONFIG["paths"]["metadata"],index_col=0)
n_source_cells=len(meta)
stratum=meta.reindex(columns=[CONFIG["stratum_column"]]).iloc[:,0].astype("string")
selection=stratum.eq(CONFIG["stratum_value"]).fillna(False)^CONFIG["exclude_stratum"]
meta=meta.loc[selection|(CONFIG["stratum_column"]=="")]
time_values=pd.to_numeric(meta[CONFIG["time_column"]].astype(str).str.extract(r"([0-9]+(?:\.[0-9]+)?)",expand=False))
groups=meta[CONFIG["group_column"]].astype("string")
samples=meta[CONFIG["sample_column"]].astype("string")
assert time_values.notna().all() and groups.notna().all() and samples.notna().all()
sample_time=pd.DataFrame({"sample":samples,"time":time_values}).drop_duplicates()
assert not sample_time["sample"].duplicated().any()
sample_time=sample_time.sort_values(["time","sample"]).set_index("sample")
counts=pd.crosstab(samples,groups).reindex(sample_time.index)
assert int(counts.to_numpy().sum())==len(meta)
proportions=counts.div(counts.sum(axis=1),axis=0)
inputs=proportions.join(sample_time).assign(n_cells=counts.sum(axis=1))
inputs.to_csv(out/"input_sample_proportions.csv",encoding="utf-8")
cell_types=proportions.columns.to_numpy()
sys.path.insert(0,str(root/"tools/competitors/dependencies/DIISCO"))
sys.path.insert(0,str(root/"tools/competitors"/CONFIG["implementation"]))
import torch
import pyro
from diisco.diisco import DIISCO
torch.set_default_dtype(torch.float64)
torch.set_num_threads(CONFIG["threads"])
pyro.set_rng_seed(CONFIG["seed"])
times=torch.as_tensor(sample_time["time"].to_numpy()).reshape(-1,1)
values=torch.as_tensor(proportions.to_numpy())
prior=torch.ones((len(cell_types),len(cell_types)))-torch.eye(len(cell_types))
model=DIISCO(lambda_matrix=prior,hypers_init_vals={},verbose=True,verbose_freq=100)
before=time.perf_counter()
model.fit(times,values,n_iter=CONFIG["n_iter"],lr=CONFIG["learning_rate"],guide=CONFIG["guide"])
seconds=time.perf_counter()-before
losses=np.asarray(model.losses)
interaction=(model.guide.W_mean.detach()*prior.unsqueeze(0)).numpy()
latent=model.guide.F_mean.detach().squeeze(-1).numpy()
assert np.isfinite(losses).all() and np.isfinite(interaction).all() and np.isfinite(latent).all()
index=pd.MultiIndex.from_product([sample_time.index,cell_types,cell_types],names=["sample","receiver","sender"])
net=pd.DataFrame({"posterior_mean":interaction.reshape(-1)},index=index).reset_index()
net=net.merge(sample_time.reset_index(),on="sample",how="left")
net.to_csv(out/"interaction_posterior_mean.csv.gz",index=False,encoding="utf-8")
pd.DataFrame(latent,index=sample_time.index,columns=cell_types).to_csv(out/"latent_proportions.csv",encoding="utf-8")
pd.DataFrame({"iteration":np.arange(1,len(losses)+1),"loss":losses}).to_csv(out/"loss.csv",index=False,encoding="utf-8")
torch.save({"parameters":pyro.get_param_store().get_state(),"config":CONFIG,"hyperparameters":model.hypers_init_vals,"times":times,"proportions":values,"lambda_matrix":prior,"cell_types":cell_types.tolist(),"samples":sample_time.index.tolist()},out/"diisco_native.pt")
coverage=pd.DataFrame([{"dataset":CONFIG["paths"]["dataset"],"method":"DIISCO","n_input_cells":len(meta),"n_samples":len(sample_time),"n_timepoints":sample_time["time"].nunique(),"n_cell_types":len(cell_types),"guide":CONFIG["guide"],"prior":CONFIG["prior"],"n_iter":len(losses),"learning_rate":CONFIG["learning_rate"],"first_loss":losses[0],"last_loss":losses[-1],"core_seconds":seconds,"result_status":str(np.where(losses[-1]<=losses[0],"scored","failed_unstable_optimization")),"task":"known_time_cell_type_interaction","ground_truth_accuracy_available":False,"initial_result_status":CONFIG["initial_status"],"initial_coverage_file":str(base_out/"failed_initial"/"coverage.csv"),"repair_detail":"native tutorial optimizer: lr=0.00005, n_iter=10000; same full input and seed"}])
coverage["n_source_cells"]=n_source_cells
coverage["unit"]=CONFIG["unit"]
coverage["stratum_column"]=CONFIG["stratum_column"]
coverage["stratum_value"]=CONFIG["stratum_value"]
coverage["exclude_stratum"]=CONFIG["exclude_stratum"]
coverage["elbo_estimator"]=CONFIG["elbo_estimator"]
coverage["repair_detail"]=coverage["repair_detail"]+"; ELBO="+CONFIG["elbo_estimator"]
coverage.to_csv(out/"coverage.csv",index=False,encoding="utf-8")
print(coverage.to_string(index=False),flush=True)

