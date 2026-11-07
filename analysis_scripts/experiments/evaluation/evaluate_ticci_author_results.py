import os
from pathlib import Path
CONFIG={"paths":{"root":str(Path(__file__).resolve().parents[2]),"native":"results/by_algorithm/TICCI_author_functional_entropy","output":"results/evaluation/final065/saved_native","reference":"inputs/comparison/results/paper_tables/reference"},"threads":os.cpu_count(),"plot":{},"component_dataset":"dentategyrus"}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"]=str(CONFIG["threads"])
import pandas as pd,numpy as np,matplotlib.pyplot as plt
import sys,json,importlib.metadata as package_metadata
from functools import partial
import h5py
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
plt.rcParams.update(CONFIG["plot"])
root=Path(CONFIG["paths"]["root"]);base=root/CONFIG["paths"]["native"];out=root/CONFIG["paths"]["output"]
status_paths=list(base.glob("*/run_status.csv"))
statuses=pd.concat(map(pd.read_csv,status_paths),ignore_index=True)
statuses.to_csv(out/"TICCI_author_run_status.csv",index=False,encoding="utf-8")
done=statuses.loc[statuses.trajectory_exit.eq(0)&statuses.entropy_exit.eq(0)&statuses.prepare_exit.eq(0)].copy()
names=pd.Index(done.dataset)
paths=list(map(lambda d:base/d/"ticci_pseudotime.csv",names))
pred=pd.concat(map(partial(pd.read_csv,dtype={"cell_id":str}),paths),keys=names,names=["dataset","row"]).reset_index().drop(columns="row")
pred["source_file"]=pred.dataset.map(pd.Series(list(map(str,paths)),index=names))
pred["algorithm"]="TICCI_SLICE"
component_dataset=CONFIG["component_dataset"]
with h5py.File(base/component_dataset/"ticci_author_result.h5ad","r") as native:
    stored=native["obsp"]["connectivities"]
    adjacency=csr_matrix((stored["data"][:],stored["indices"][:],stored["indptr"][:]),shape=tuple(stored.attrs["shape"]))
    n_components,component_ids=connected_components(adjacency,directed=False)
    stored=native["obsp"]["distances"]
    original_knn=csr_matrix((np.ones_like(stored["data"][:]),stored["indices"][:],stored["indptr"][:]),shape=tuple(stored.attrs["shape"]))
    n_components_before,component_ids_before=connected_components(original_knn,directed=False)
    assert n_components==n_components_before and np.array_equal(component_ids,component_ids_before)
    cell_ids=native["obs"][native["obs"].attrs["_index"]].asstr()[:]
    root_index=int(native["uns"]["iroot"][()])
    population=native["obs"]["clusters"]
    population_names=population["categories"].asstr()[:][population["codes"][:]]
    native_time=native["obs"]["dpt_pseudotime"][:]
    component_cells=pd.DataFrame({"dataset":component_dataset,"cell_id":cell_ids,"component_id":component_ids,"population":population_names,
        "age_days":native["obs"]["age(days)"]["categories"][:][native["obs"]["age(days)"]["codes"][:]].astype(float),"native_pseudotime":native_time,"finite_prediction":np.isfinite(native_time)})
root_component=int(component_ids[root_index])
assert np.isfinite(native_time[component_ids==root_component]).all() and np.isinf(native_time[component_ids!=root_component]).all()
component_cells["root_connected"]=component_cells.component_id.eq(root_component)
component_cells["native_status"]=np.where(component_cells.root_connected,"finite_native_root_distance","native_infinite_distance_disconnected_from_root")
component_cells.to_csv(out/"TICCI_author_cell_components.csv",index=False,encoding="utf-8")
component_summary=component_cells.groupby(["dataset","component_id"]).agg(n_cells=("cell_id","size"),n_finite=("finite_prediction","sum"),root_connected=("root_connected","first"),native_status=("native_status","first"))
population_counts=component_cells.groupby(["dataset","component_id","population"]).size().rename("n").reset_index()
population_counts["description"]=population_counts.population+":"+population_counts.n.astype(str)
population_descriptions=population_counts.groupby(["dataset","component_id"]).description.agg(";".join).rename("populations")
age_counts=component_cells.groupby(["dataset","component_id","age_days"]).size().unstack(fill_value=0).add_prefix("n_age_")
component_summary=component_summary.join(population_descriptions).join(age_counts).reset_index()
component_summary["n_infinite"]=component_summary.n_cells-component_summary.n_finite
component_summary["native_root_cell_id"]=cell_ids[root_index]
component_summary["components_identical_before_after_CCI"]=True
component_summary["cross_component_time_alignment"]="undefined_by_native_single_root_distance;no_cross_component_edges_added"
component_summary["definition_source"]="https://academic.oup.com/bioinformatics/article/41/2/btaf027/7997266#2.5"
component_summary.to_csv(out/"TICCI_author_component_coverage.csv",index=False,encoding="utf-8")

old_binding=pd.read_csv(out/"native_TICCI_implementation_binding.csv")
old_paths=list(map(lambda p:Path(p).parent/"ticci_pseudotime.csv",old_binding.source_file))
old=pd.concat(map(partial(pd.read_csv,dtype={"cell_id":str}),old_paths),keys=old_binding.dataset,names=["dataset","row"]).reset_index().drop(columns="row")
old["source_file"]=old.dataset.map(pd.Series(list(map(str,old_paths)),index=old_binding.dataset))
old["algorithm"]="TICCI_historical_Shannon"
references={
 "dentategyrus":("ref_dentategyrus.csv","age(days)"),
 "gastrulation_erythroid":("ref_gastrulation_erythroid.csv","stage"),
 "zebrafish_neural_dev":("ref_zebrafish_neural_dev.csv","timepoint"),
 "lung":("ref_lung.csv","time_point"),
 "larry_invitro":("ref_larry_invitro.csv","time_info"),
 "human_fetal_liver_organoid":("ref_human_fetal_liver_organoid.csv","time_days"),
 "human_fetal_pancreas":("ref_human_fetal_pancreas.csv","pcw"),
 "zebrafish_tracerseq":("ref_zebrafish_tracerseq.csv","time_hours")}
synthetic_names=names[names.str.startswith("sim")]
references.update(dict(zip(synthetic_names,map(lambda d:("ref_"+d+".csv","true_pseudotime"),synthetic_names))))
refs=pd.concat(map(lambda k:pd.read_csv(root/CONFIG["paths"]["reference"]/references[k][0],dtype={"cell_id":str}).rename(columns={references[k][1]:"truth"})[["cell_id","truth"]],references),keys=references.keys(),names=["dataset","row"]).reset_index().drop(columns="row")
truth_numeric=pd.to_numeric(refs.truth,errors="coerce")
truth_label=refs.truth.astype("string").str.extract(r"([0-9]+(?:\.[0-9]+)?)",expand=False).astype(float)
refs["truth"]=truth_numeric.fillna(truth_label)
assert not refs.duplicated(["dataset","cell_id"]).any()
predictions=pd.concat([pred,old],ignore_index=True)
assert not predictions.duplicated(["dataset","algorithm","cell_id"]).any()
joined=predictions.merge(refs,on=["dataset","cell_id"],how="left",validate="many_to_one",indicator=True)
joined["matched"]=joined._merge.eq("both")
joined["finite_prediction"]=np.isfinite(joined.pseudotime)
keys=["dataset","algorithm"]
aligned=joined.groupby(keys).agg(n_prediction_cells=("cell_id","size"),n_aligned=("matched","sum"),n_finite_prediction=("finite_prediction","sum"),source_file=("source_file","first"))
aligned=aligned.join(refs.groupby("dataset").size().rename("n_reference"),on="dataset")
valid=joined.loc[joined.matched&np.isfinite(joined.pseudotime)&np.isfinite(joined.truth)].copy()
grouped=valid.groupby(keys)
valid["p01"]=(valid.pseudotime-grouped.pseudotime.transform("min"))/(grouped.pseudotime.transform("max")-grouped.pseudotime.transform("min"))
valid["t01"]=(valid.truth-grouped.truth.transform("min"))/(grouped.truth.transform("max")-grouped.truth.transform("min"))
valid["squared_error"]=(valid.p01-valid.t01)**2
spearman=valid.groupby(keys)[["p01","t01"]].corr(method="spearman").xs("p01",level=2).t01.rename("Spearman_signed")
pearson=valid.groupby(keys)[["p01","t01"]].corr(method="pearson").xs("p01",level=2).t01.rename("Pearson")
metrics=aligned.join(valid.groupby(keys).size().rename("n_used")).join(spearman).join(pearson).join(valid.groupby(keys).squared_error.mean().pow(.5).rename("RMSE")).reset_index()
metrics["cell_coverage"]=metrics.n_aligned/metrics.n_reference
metrics["output_status"]=np.select([metrics.n_reference.isna(),metrics.n_aligned.ne(metrics.n_reference)|metrics.n_prediction_cells.ne(metrics.n_reference),metrics.n_finite_prediction.ne(metrics.n_prediction_cells)],["no_independent_observed_time_reference","incomplete_input_scope","nonfinite_prediction"],default="scored")
native_disconnected=metrics.dataset.eq(component_dataset)&metrics.algorithm.eq("TICCI_SLICE")
metrics.loc[native_disconnected,"output_status"]="native_disconnected_infinite_distance"
metrics["nonfinite_reason"]=np.where(native_disconnected,"34_Cajal_Retzius_cells_disconnected_from_native_root;author_section2.5_defines_infinite_distance","")
metrics["root_input_mode"]="native_entropy_auto_root"
metrics["direction_policy"]="native_signed_no_truth_flip"
metrics["metric_source"]="native_author_TICCI_functional_entropy_full_input"
metrics["row_present"]=True
metrics["Spearman"]=metrics.Spearman_signed
metrics["Spearman_absolute"]=metrics.Spearman_signed.abs()
metrics["source_host"]="server"
original=pd.read_csv(root/"results/evaluation/final065/ordering_pseudotime_complete_Final065.csv").loc[lambda d:d.algorithm.eq("TICCI"),["dataset","Spearman_signed","Pearson","RMSE"]]
replay=metrics.loc[metrics.algorithm.eq("TICCI_historical_Shannon")].merge(original,on="dataset",suffixes=("_replay","_saved"))
replay["max_abs_metric_difference"]=np.nanmax(np.abs(replay[["Spearman_signed_replay","Pearson_replay","RMSE_replay"]].to_numpy()-replay[["Spearman_signed_saved","Pearson_saved","RMSE_saved"]].to_numpy()),axis=1)
assert replay.max_abs_metric_difference.max()<1e-10
replay.to_csv(out/"TICCI_metric_operator_readback.csv",index=False,encoding="utf-8")
new=metrics.loc[metrics.algorithm.eq("TICCI_SLICE")].copy()
new.to_csv(out/"classic_TICCI_author_ordering_metrics.csv",index=False,encoding="utf-8")
input_paths=list(map(lambda d:base/d/"input_binding.csv",names))
bindings=pd.concat(map(pd.read_csv,input_paths),ignore_index=True)
entropy_paths=list(map(lambda d:base/d/"ticci_entropy_summary.csv",names))
entropies=pd.concat(map(pd.read_csv,entropy_paths),ignore_index=True)
summaries=pd.concat(map(lambda d:pd.read_csv(base/d/"ticci_summary.csv").assign(dataset=d),names),ignore_index=True)
binding=bindings.merge(entropies,on=["dataset","n_cells","seed"],validate="one_to_one").merge(summaries.drop(columns=["method","seed","bootstrap_B","bootstrap_genes"]),on=["dataset","n_cells"],validate="one_to_one")
binding["entropy_mode"]=binding.entropy_mode_x.combine_first(binding.entropy_mode_y).fillna("go_kappa")
binding["implementation_label"]=np.where(binding.entropy_mode.eq("expression_gene_groups"),
 "TICCI SLICE supplied RNA gene groups B100;fixed RNA graph;zero matching CCI symbols",
 "TICCI SLICE GO-Kappa entropy B100;fixed CCI weight1 and resolution1")
binding["historical_graph_frozen_for_recompute"]=False
binding["python_runtime"]=sys.version.split()[0]
version_names=["numpy","pandas","scanpy","scipy","scikit-learn","igraph","networkx"]
binding["python_package_versions"]=json.dumps(dict(zip(version_names,map(package_metadata.version,version_names))),sort_keys=True)
binding["historical_comparison_limit"]="graph_reconstructed_in_current_environment;changes_not_attributed_to_entropy_alone"
binding.to_csv(out/"native_TICCI_author_implementation_binding.csv",index=False,encoding="utf-8")
cov=new.rename(columns={"algorithm":"method","output_status":"comparison_status"}).copy()
cov=cov.merge(bindings[["dataset","n_genes"]],on="dataset",validate="one_to_one")
cov["n_input_genes"]=cov.n_genes
cov["task"]="cell_ordering";cov["n_input_cells"]=cov.n_prediction_cells
cov["n_native_units"]=cov.n_prediction_cells;cov["n_reference_units"]=cov.n_reference
cov["n_scored_units"]=cov.n_used.fillna(0)
cov["classic_metrics"]="Spearman_signed;Pearson;normalized_RMSE"
cov["project_specific_metrics"]=""
cov["result_status"]="completed_author_functional_entropy_native_output"
cov["analysis_unit"]="cell_native_pseudotime"
cov["ground_truth_accuracy"]=cov.n_reference.notna()
cov["evidence_type"]="fixed_observed_collection_time;full_input_and_finite_metric_denominators_explicit"
cov=cov.merge(bindings[["dataset","entropy_mode"]],on="dataset",validate="one_to_one")
cov["entropy_mode"]=cov.entropy_mode.fillna("go_kappa")
cov["input_scope"]=np.where(cov.entropy_mode.eq("expression_gene_groups"),
 "full_original_cells_and_genes;SLICE_supplied_expression_gene_groups;B100_Bsize1000_seed12345;zero_matched_CCI",
 "full_original_native_input;SLICE_functional_Kappa;B100_Bsize1000_seed12345;CellChat_calls")

cov["denominator_note"]="all_original_cells_retained;metrics_only_finite_truth_and_finite_prediction;no_direction_flip"
cov.loc[cov.dataset.eq(component_dataset),"denominator_note"]="all2930_cells_retained;2896_finite_root_component;34_Cajal_Retzius_cells_native_infinite_distance;no_cross_component_order_imputed"
cov.loc[cov.dataset.eq(component_dataset),"result_status"]="completed_native_with_author_defined_disconnected_infinite_distances"
cov["implementation_label"]="TICCI SLICE entropy;CCI k1;Louvain resolution1;current environment"
cov["evaluation_source"]=(out/"classic_TICCI_author_ordering_metrics.csv").as_posix()
graph_cov=binding.assign(method="TICCI_SLICE",task="native_graph_topology",n_input_cells=binding.n_cells,n_native_units=binding.n_clusters,n_scored_units=0,analysis_unit="native_louvain_PAGA_and_author_directed_tree",classic_metrics="",project_specific_metrics="",ground_truth_accuracy=False,evidence_type="saved_native_graph_without_independent_GT_graph",comparison_status="descriptive_native_topology",result_status="saved_native_graph_read",evaluation_source=(out/"native_TICCI_author_implementation_binding.csv").as_posix(),denominator_note="all_native_input_cells;node_edge_counts_not_graph_accuracy")
pd.concat([cov,graph_cov],ignore_index=True,sort=False).to_csv(out/"native_TICCI_author_task_coverage.csv",index=False,encoding="utf-8")
print("TICCI_AUTHOR",len(done),"completed;",int(new.n_reference.notna().sum()),"fixed time references;historical_metric_max_difference",replay.max_abs_metric_difference.max(),flush=True)
print(new[["dataset","n_prediction_cells","n_reference","n_used","Spearman_signed","output_status"]].to_string(index=False),flush=True)
