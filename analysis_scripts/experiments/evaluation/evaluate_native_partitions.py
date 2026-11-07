import os
from pathlib import Path
CONFIG={"paths":{"root":str(Path(__file__).resolve().parents[2]),"native":"inputs/comparison/recovered_native","output":"results/evaluation/final065/saved_native","reference":"inputs/comparison/results/paper_tables/reference"},"threads":os.cpu_count(),"plot":{},"terminal_datasets":["simulated_branch","sim_2br_500c","sim_2br_2000c","sim_2br_5000c"],"terminal_min_time":0.4,"probability_tolerance":1e-6}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"]=str(CONFIG["threads"])
import operator as op
from functools import partial
import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse
from scipy.optimize import linear_sum_assignment
from sklearn.metrics import adjusted_rand_score,normalized_mutual_info_score
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG["plot"])
root=Path(CONFIG["paths"]["root"])
base=root/CONFIG["paths"]["native"]
out=root/CONFIG["paths"]["output"]
out.mkdir(parents=True,exist_ok=True)
paths=sorted((root/CONFIG["paths"]["reference"]).glob("ref_*.csv"))
names=pd.Index(map(op.attrgetter("stem"),paths)).str.removeprefix("ref_")
refs=pd.concat(map(op.methodcaller("rename",columns={"dataset":"reference_annotation_dataset","source_row":"reference_annotation_row"}),map(partial(pd.read_csv,dtype={"cell_id":"string","true_branch":"string"}),paths)),keys=names,names=["dataset","source_row"]).reset_index()
old_paths=sorted((base/"paper_tables/reference").glob("ref_*.csv"))
old_names=pd.Index(map(op.attrgetter("stem"),old_paths)).str.removeprefix("ref_")
old=pd.concat(map(op.methodcaller("rename",columns={"dataset":"reference_annotation_dataset","source_row":"reference_annotation_row"}),map(partial(pd.read_csv,dtype={"cell_id":"string","true_branch":"string"}),old_paths)),keys=old_names,names=["dataset","source_row"]).reset_index()
refs=pd.concat([old.loc[~old.dataset.isin(refs.dataset)],refs],ignore_index=True)
assert not refs.duplicated(["dataset","cell_id"]).any()
scope=refs.loc[refs.true_branch.notna()].copy()
scope=scope.loc[~scope.dataset.isin(CONFIG["terminal_datasets"])|scope.true_pseudotime.ge(CONFIG["terminal_min_time"])].copy()
scope["evaluation_scope"]=np.where(scope.dataset.isin(CONFIG["terminal_datasets"]),"submission_terminal_cells","all_cells_diagnostic")
planned=scope.groupby(["dataset","evaluation_scope"],observed=True).size().rename("n_scope_planned").reset_index()

paths=sorted((base/"by_algorithm/FateID").glob("*/fateid_bias.csv"))
names=pd.Index(map(lambda p:p.parent.name,paths))
headers=list(map(partial(pd.read_csv,nrows=0),paths))
columns=pd.DataFrame({"dataset":names,"native_target":list(map(op.attrgetter("columns"),headers))}).explode("native_target")
columns=columns.loc[columns.native_target.str.match(r"^t[0-9]+$")].copy()
columns["native_target_order"]=columns.groupby("dataset",sort=False).cumcount()
fi=pd.concat(map(partial(pd.read_csv,dtype={"cell_id":"string"}),paths),keys=names,names=["dataset","source_row"]).reset_index()
fi["source_file"]=fi.dataset.map(pd.Series(list(map(str,paths)),index=names))
fi_long=fi.melt(id_vars=["dataset","source_row","source_file","cell_id"],value_vars=columns.native_target.unique(),var_name="native_target",value_name="probability")
fi_long=fi_long.merge(columns,on=["dataset","native_target"],validate="many_to_one")
fi_long["method"]="FateID"
fi_long["input_condition"]="provided_input_clusters_used_for_native_targets"

paths=sorted((base/"misc").glob("cellrank2_*_result.h5ad"))
names=pd.Index(map(op.attrgetter("stem"),paths)).str.removeprefix("cellrank2_").str.removesuffix("_result")
objects=list(map(partial(ad.read_h5ad,backed="r"),paths))
obs=pd.concat(map(op.methodcaller("rename",columns={"dataset":"native_annotation_dataset","cell_id":"native_annotation_cell_id"}),map(op.attrgetter("obs"),objects)),keys=names,names=["dataset","cell_id"]).reset_index()
obs["source_row"]=obs.groupby("dataset",sort=False).cumcount()
obs["cell_id"]=obs.cell_id.astype("string")
matrices=list(map(op.itemgetter("lineages_fwd"),map(op.attrgetter("obsm"),objects)))
cr=pd.concat(map(pd.DataFrame,matrices),keys=names,names=["dataset","source_row"]).reset_index()
cr=cr.merge(obs[["dataset","source_row","cell_id"]],on=["dataset","source_row"],validate="one_to_one")
cr["source_file"]=cr.dataset.map(pd.Series(list(map(str,paths)),index=names))
counts=pd.Series(list(map(op.itemgetter(1),map(np.shape,matrices))),index=names)
cr_long=cr.melt(id_vars=["dataset","source_row","source_file","cell_id"],var_name="native_target",value_name="probability")
cr_long["native_target_order"]=pd.to_numeric(cr_long.native_target)
cr_long=cr_long.loc[cr_long.native_target_order.lt(cr_long.dataset.map(counts))].copy()
cr_long["native_target"]=cr_long.native_target.astype(str)
cr_long["method"]="CellRank2"
cr_long["input_condition"]="saved_native_ConnectivityKernel_GPCCA"
list(map(op.methodcaller("close"),map(op.attrgetter("file"),objects)))

paths=sorted((root/"inputs/comparison/results/by_algorithm/slingshot").glob("*/slingshot_weights.csv"))
names=pd.Index(map(lambda p:p.parent.name,paths))
weights=list(map(pd.read_csv,paths))
sh_columns=pd.DataFrame({"dataset":names,"native_target":list(map(op.attrgetter("columns"),weights))}).explode("native_target")
sh_columns["native_target_order"]=sh_columns.groupby("dataset",sort=False).cumcount()
sh=pd.concat(weights,keys=names,names=["dataset","source_row"]).reset_index()
sh_cells=pd.concat(map(partial(pd.read_csv,dtype={"cell_id":"string"}),map(lambda p:p.with_name("slingshot_pseudotime.csv"),paths)),keys=names,names=["dataset","source_row"]).reset_index()
sh=sh.merge(sh_cells[["dataset","source_row","cell_id"]],on=["dataset","source_row"],validate="one_to_one")
sh["source_file"]=sh.dataset.map(pd.Series(list(map(str,paths)),index=names))
sh_long=sh.melt(id_vars=["dataset","source_row","source_file","cell_id"],var_name="native_target",value_name="probability")
sh_long=sh_long.merge(sh_columns,on=["dataset","native_target"],validate="many_to_one")
sh_long["method"]="slingshot"
sh_long["input_condition"]="saved_slingCurveWeights_as_probs"

llr=pd.read_csv(out/"current_LineaLinkR_fates.csv.gz",dtype={"cell_id":"string","native_target":"string"})
llr["source_row"]=llr.groupby(["dataset","method","source_file","native_target"],sort=False).cumcount()
llr_targets=llr[["dataset","method","source_file","native_target"]].drop_duplicates()
llr_targets["native_target_order"]=llr_targets.groupby(["dataset","method","source_file"],sort=False).cumcount()
llr=llr.merge(llr_targets,on=["dataset","method","source_file","native_target"],validate="many_to_one")
llr["input_condition"]="saved_model_root_policy_recorded_in_R_object_catalog"
prob=pd.concat([fi_long,cr_long,sh_long,llr],ignore_index=True)
keys=["dataset","method","source_file"]
row_keys=keys+["cell_id"]
prob["finite"]=np.isfinite(prob.probability)
prob["outside_unit"]=prob.probability.lt(-CONFIG["probability_tolerance"])|prob.probability.gt(1+CONFIG["probability_tolerance"])
quality=prob.groupby(row_keys,observed=True).agg(n_fields=("probability","size"),n_finite=("finite","sum"),row_sum=("probability","sum"),outside_unit=("outside_unit","any")).reset_index()
quality["complete_finite"]=quality.n_fields.eq(quality.n_finite)
quality["normalized"]=quality.complete_finite&quality.row_sum.sub(1).abs().le(CONFIG["probability_tolerance"])
finite=prob.loc[prob.finite].sort_values(keys+["cell_id","native_target_order"])
index=finite.groupby(row_keys,observed=True).probability.idxmax()
hard=finite.loc[index,keys+["cell_id","native_target","input_condition"]].rename(columns={"native_target":"prediction"})
max_prob=finite.groupby(row_keys,observed=True).probability.transform("max")
ties=finite.loc[finite.probability.eq(max_prob)].groupby(row_keys,observed=True).size().rename("n_argmax_ties").reset_index()
hard=hard.merge(ties,on=row_keys,validate="one_to_one")
hard["prediction_unit"]="native_probability_argmax"

paths=sorted((base/"misc").glob("scfates_lineage_cell_segments_*.csv"))
names=pd.Index(map(op.attrgetter("stem"),paths)).str.removeprefix("scfates_lineage_cell_segments_")
sf=pd.concat(map(partial(pd.read_csv,dtype={"cell":"string"}),paths),keys=names,names=["dataset","source_row"]).reset_index().rename(columns={"cell":"cell_id","segment":"prediction"})
sf["source_file"]=sf.dataset.map(pd.Series(list(map(str,paths)),index=names))
sf["method"]="scFates_native_segments"
sf["input_condition"]="saved_native_segment_labels;root_policy_not_stored_in_export"
sf["prediction_unit"]="native_principal_tree_segment_partition"
sf["n_argmax_ties"]=1
sf["prediction"]=sf.prediction.astype("string")
hard=pd.concat([hard,sf[hard.columns]],ignore_index=True)
hard["cell_id"]=hard.cell_id.astype("string")
hard["prediction"]=hard.prediction.astype("string")
assert not hard.duplicated(row_keys).any()
reference_keys=refs[["dataset","cell_id"]].assign(reference_cell=True)
identity=hard.merge(reference_keys,on=["dataset","cell_id"],how="left",validate="many_to_one")
identity["reference_cell"]=identity.reference_cell.eq(True)
coverage=identity.groupby(keys,observed=True).agg(n_native_cells=("cell_id","size"),n_aligned_cells=("reference_cell","sum"),n_tied_cells=("n_argmax_ties",lambda x:x.gt(1).sum()),prediction_unit=("prediction_unit","first"),input_condition=("input_condition","first")).reset_index()
coverage["n_reference_cells"]=coverage.dataset.map(refs.groupby("dataset").size())
coverage["n_unexpected_cells"]=coverage.n_native_cells-coverage.n_aligned_cells
coverage["n_missing_cells"]=coverage.n_reference_cells-coverage.n_aligned_cells
coverage["exact_cell_set"]=coverage.n_unexpected_cells.eq(0)&coverage.n_missing_cells.eq(0)
pquality=quality.groupby(keys,observed=True).agg(n_probability_rows=("cell_id","size"),n_complete_finite=("complete_finite","sum"),n_normalized=("normalized","sum"),n_outside_unit=("outside_unit","sum")).reset_index()
coverage=coverage.merge(pquality,on=keys,how="left",validate="one_to_one")
matched=hard.merge(scope[["dataset","cell_id","true_branch","evaluation_scope"]],on=["dataset","cell_id"],validate="many_to_one")
groups=list(matched.groupby(keys+["evaluation_scope"],sort=True,observed=True))
labels=list(map(op.itemgetter(0),groups))
frames=list(map(op.itemgetter(1),groups))
truth_vectors=list(map(op.itemgetter("true_branch"),frames))
prediction_vectors=list(map(op.itemgetter("prediction"),frames))
classic=pd.DataFrame(labels,columns=keys+["evaluation_scope"])
classic["ARI"]=list(map(adjusted_rand_score,truth_vectors,prediction_vectors))
classic["NMI"]=list(map(normalized_mutual_info_score,truth_vectors,prediction_vectors))
sizes=matched.groupby(keys+["evaluation_scope"],observed=True).agg(n_scored=("cell_id","size"),n_truth_classes=("true_branch","nunique"),n_prediction_classes=("prediction","nunique")).reset_index()
classic=classic.merge(sizes,on=keys+["evaluation_scope"],validate="one_to_one").merge(planned,on=["dataset","evaluation_scope"],validate="many_to_one").merge(coverage,on=keys,validate="one_to_one")
classic["score_status"]=np.where(classic.exact_cell_set&classic.n_scored.eq(classic.n_scope_planned),"full_declared_cell_scope","partial_cell_scope_reported")
classic["evidence_type"]=np.where(classic.evaluation_scope.eq("submission_terminal_cells"),"declared_terminal_branch_reference","all_cell_segment_label_diagnostic")
classic["metric_family"]="classical_partition_agreement"
classic["ground_truth_accuracy"]=classic.evaluation_scope.eq("submission_terminal_cells")&classic.prediction_unit.eq("native_probability_argmax")
classic.to_csv(out/"classic_native_partition_metrics.csv",index=False,encoding="utf-8")
joint=matched.groupby(keys+["evaluation_scope","true_branch","prediction"],observed=True).size().rename("n").reset_index()
joint["truth_total"]=joint.groupby(keys+["evaluation_scope","true_branch"],observed=True).n.transform("sum")
joint["prediction_total"]=joint.groupby(keys+["evaluation_scope","prediction"],observed=True).n.transform("sum")
joint["pair_f1"]=2*joint.n/(joint.truth_total+joint.prediction_total)
truth_index=joint[keys+["evaluation_scope","true_branch"]].drop_duplicates()
truth_index["truth_index"]=truth_index.groupby(keys+["evaluation_scope"],observed=True).cumcount()
pred_index=joint[keys+["evaluation_scope","prediction"]].drop_duplicates()
pred_index["prediction_index"]=pred_index.groupby(keys+["evaluation_scope"],observed=True).cumcount()
joint=joint.merge(truth_index,on=keys+["evaluation_scope","true_branch"],validate="many_to_one").merge(pred_index,on=keys+["evaluation_scope","prediction"],validate="many_to_one")
fgroups=list(joint.groupby(keys+["evaluation_scope"],sort=True,observed=True))
fframes=list(map(op.itemgetter(1),fgroups))
coordinates=list(zip(map(op.itemgetter("truth_index"),fframes),map(op.itemgetter("prediction_index"),fframes)))
tables=list(map(sparse.coo_matrix,zip(map(op.itemgetter("pair_f1"),fframes),coordinates)))
dense=list(map(op.methodcaller("toarray"),tables))
assignments=list(map(partial(linear_sum_assignment,maximize=True),dense))
sums=np.fromiter(map(np.sum,map(op.getitem,dense,assignments)),dtype=float)
custom=pd.DataFrame(list(map(op.itemgetter(0),fgroups)),columns=keys+["evaluation_scope"])
custom["matched_f1_sum"]=sums
custom=custom.merge(classic.drop(columns=["ARI","NMI"]),on=keys+["evaluation_scope"],validate="one_to_one")
custom["Hungarian_matched_Branch_F1"]=custom.matched_f1_sum/np.maximum(custom.n_truth_classes,custom.n_prediction_classes)
custom["metric_family"]="project_specific_label_matching_protocol"
custom["denominator"]="max(number_of_truth_classes,number_of_native_prediction_classes)"
custom.to_csv(out/"project_native_matched_branch_F1.csv",index=False,encoding="utf-8")
coverage=coverage.merge(planned,on="dataset",how="left",validate="many_to_one")
coverage["reference_status"]=np.where(coverage.n_scope_planned.notna(),"explicit_branch_labels_available","no_explicit_branch_truth")
coverage.to_csv(out/"native_partition_input_coverage.csv",index=False,encoding="utf-8")
quality.to_csv(out/"native_probability_row_checks.csv.gz",index=False,encoding="utf-8")
matched.to_csv(out/"native_partition_matched_predictions.csv.gz",index=False,encoding="utf-8")
print(classic.groupby(["method","evaluation_scope"],observed=True).agg(n=("dataset","size"),ARI=("ARI","mean"),NMI=("NMI","mean")).to_string())

