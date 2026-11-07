import os
_VARIANT=os.environ.get("LLR_CMP_VARIANT","Final065")
CONFIG={"paths":{"root":os.environ.get("LLR_CMP_ROOT",os.environ.get("LLR_WORK_ROOT",os.getcwd())),
 "evaluation":os.environ.get("LLR_CMP_EVAL","results/evaluation/final065"),
 "ordering":"ordering_pseudotime_complete_"+_VARIANT+".csv",
 "author_ordering":"saved_native/classic_TICCI_author_ordering_metrics.csv",
 "partitions":"saved_native/classic_native_partition_metrics.csv",
 "matched_f1":"saved_native/project_native_matched_branch_F1.csv",
 "dynamic":"dynamic_methods.csv","lrsim":"lrsim_metrics.csv",
 "comparison":"main_task_comparison.csv","summary":"main_task_summary.csv"},
 "threads":os.cpu_count(),"plot":{},"main_algorithm":os.environ.get("LLR_CMP_MAIN_ALGORITHM","LLR_0.6.5_"+_VARIANT),
 "main_partition":"LineaLinkR_0.6.5_"+_VARIANT,"tolerance":1e-10,
 "ordering_metrics":["Spearman_signed","Spearman_absolute","Pearson","RMSE"],
 "partition_metrics":["ARI","NMI","Matched_F1"]}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"]=str(CONFIG["threads"])
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG["plot"])
CONFIG["paths"].update({"classic":"classic_metrics.csv","specific":"project_specific_metrics.csv"})
root=Path(CONFIG["paths"]["root"])
base=root/CONFIG["paths"]["evaluation"]

# Compare each actual method on its available native output; retain the fixed input grid.
detail=pd.read_csv(base/CONFIG["paths"]["ordering"])
own=detail.loc[detail.algorithm.eq(CONFIG["main_algorithm"])].copy()
author=pd.DataFrame({"dataset":own.dataset}).merge(
 pd.read_csv(base/CONFIG["paths"]["author_ordering"]),on="dataset",how="left",validate="one_to_one")
author["algorithm"]="TICCI_SLICE"
author["output_status"]=author.output_status.fillna("not_run_author_configuration")
external=pd.concat([detail.loc[~detail.algorithm.str.startswith("LLR_")&detail.algorithm.ne("TICCI")],
 author],ignore_index=True,sort=False)
own["Spearman_absolute"]=own.Spearman_signed.abs()
external["Spearman_absolute"]=external.Spearman_signed.abs()
main_order=own[["dataset","root_input_mode","output_status","n_aligned"]+CONFIG["ordering_metrics"]].rename(
 columns={"output_status":"main_status","n_aligned":"main_n_aligned"}).melt(
 id_vars=["dataset","root_input_mode","main_status","main_n_aligned"],
 value_vars=CONFIG["ordering_metrics"],var_name="metric",value_name="main_value")
other_order=external.reindex(columns=["dataset","algorithm","root_input_mode","output_status",
 "n_aligned","n_used","source_file"]+CONFIG["ordering_metrics"]).rename(
 columns={"algorithm":"comparator","root_input_mode":"comparator_root_input_mode",
 "output_status":"comparator_status","n_aligned":"comparator_n_aligned",
 "n_used":"comparator_n_used","source_file":"comparator_source_file"}).melt(
 id_vars=["dataset","comparator","comparator_root_input_mode","comparator_status",
 "comparator_n_aligned","comparator_n_used","comparator_source_file"],
 value_vars=CONFIG["ordering_metrics"],var_name="metric",value_name="comparator_value")
ordering=main_order.merge(other_order,on=["dataset","metric"],validate="one_to_many")
ordering["task"]="cell_ordering"
ordering["evaluation_scope"]="native_direction"
ordering["higher_is_better"]=ordering.metric.ne("RMSE")
ordering["same_candidate_denominator"]=ordering.main_n_aligned.eq(ordering.comparator_n_aligned)
ordering["status_eligible"]=ordering.main_status.eq("scored")&ordering.comparator_status.eq("scored")
ordering["comparison_policy"]="actual_method"
ordering["interpretation"]="native input conditions and full-input coverage accompany each method"
ordering["evaluation_source"]=CONFIG["paths"]["ordering"]

# State partitions use matched cells and reference labels, with scopes kept separate.
parts=pd.read_csv(base/CONFIG["paths"]["partitions"])
f1=pd.read_csv(base/CONFIG["paths"]["matched_f1"])
parts=parts.merge(f1[["dataset","method","evaluation_scope","Hungarian_matched_Branch_F1"]],
 on=["dataset","method","evaluation_scope"],validate="one_to_one").rename(
 columns={"Hungarian_matched_Branch_F1":"Matched_F1"})
# Keep the existing branch-table names synchronized with the canonical native partition evaluation.
branch_grid=parts.rename(columns={"method":"algorithm","evaluation_scope":"branch_scope","Matched_F1":"Branch_F1"}).copy()
branch_grid["algorithm"]=branch_grid.algorithm.replace({CONFIG["main_partition"]:CONFIG["main_algorithm"]})
branch_grid["branch_status"]="scored_native_partition"
branch_grid["branch_metric_source"]=CONFIG["paths"]["partitions"]
branch_grid["evaluation_object"]="state_partition"
branch_grid.to_csv(base/("ordering_branch_complete_"+_VARIANT+".csv"),index=False,encoding="utf-8")
part_own=parts.loc[parts.method.eq(CONFIG["main_partition"])]
part_other=parts.loc[~parts.method.str.startswith(("LineaLinkR","LLR_"))]
part_main=part_own[["dataset","evaluation_scope","n_scored","n_scope_planned","input_condition",
 "source_file"]+CONFIG["partition_metrics"]].rename(columns={"n_scored":"main_n_aligned",
 "n_scope_planned":"main_scope_planned","input_condition":"main_input_condition",
 "source_file":"main_source_file"}).melt(
 id_vars=["dataset","evaluation_scope","main_n_aligned","main_scope_planned",
 "main_input_condition","main_source_file"],value_vars=CONFIG["partition_metrics"],
 var_name="metric",value_name="main_value")
part_cmp=part_other[["dataset","method","evaluation_scope","n_scored","n_scope_planned",
 "input_condition","source_file","prediction_unit"]+CONFIG["partition_metrics"]].rename(
 columns={"method":"comparator","n_scored":"comparator_n_aligned",
 "n_scope_planned":"comparator_scope_planned","input_condition":"comparator_input_condition",
 "source_file":"comparator_source_file"}).melt(
 id_vars=["dataset","comparator","evaluation_scope","comparator_n_aligned",
 "comparator_scope_planned","comparator_input_condition","comparator_source_file","prediction_unit"],
 value_vars=CONFIG["partition_metrics"],var_name="metric",value_name="comparator_value")
branch=part_main.merge(part_cmp,on=["dataset","evaluation_scope","metric"],validate="one_to_many")
branch["task"]="state_partition"
branch["higher_is_better"]=True
branch["same_candidate_denominator"]=branch.main_n_aligned.eq(branch.comparator_n_aligned)&branch.main_scope_planned.eq(branch.comparator_scope_planned)
branch["status_eligible"]=branch.same_candidate_denominator
branch["main_status"]="saved_native_partition"
branch["comparator_status"]="saved_native_partition"
branch["comparison_policy"]="actual_method"
branch["interpretation"]="ARI/NMI and Hungarian-matched F1 on the declared identical cell scope and label reference"
branch["evaluation_source"]=CONFIG["paths"]["partitions"]

# Molecular nomination and controlled event projection have different main-method keys.
dynamic=pd.read_csv(base/CONFIG["paths"]["dynamic"])
is_main=(dynamic.task.eq("LR_pair_nomination")&dynamic.method.eq("LineaLinkR_event_max"))|(
 dynamic.task.eq("projected_fork_event")&dynamic.method.eq("LineaLinkR_event_score"))
dynamic_main=dynamic.loc[is_main].rename(columns={"method":"main_method",
 "n_candidates":"main_candidates","n_positive":"main_positive","n_targets":"main_targets",
 "precision_at_k":"main_precision_at_k","auprc":"main_auprc",
 "inference_status":"main_inference_status","evaluation_status":"main_status"})
dynamic_other=dynamic.loc[~is_main].rename(columns={"method":"comparator",
 "n_candidates":"comparator_candidates","n_positive":"comparator_positive",
 "precision_at_k":"comparator_precision_at_k","auprc":"comparator_auprc",
 "inference_status":"comparator_inference_status","evaluation_status":"comparator_status"})
paired=dynamic_main.merge(dynamic_other,on=["dataset","task","evaluation_scope"],validate="one_to_many")
paired=paired.melt(id_vars=["dataset","task","evaluation_scope","comparator","main_method",
 "main_candidates","comparator_candidates","main_positive","comparator_positive","main_targets",
 "main_inference_status","comparator_inference_status","main_status","comparator_status",
 "comparator_precision_at_k","comparator_auprc"],value_vars=["main_precision_at_k","main_auprc"],
 var_name="metric",value_name="main_value")
paired["comparator_value"]=np.where(paired.metric.eq("main_auprc"),paired.comparator_auprc,paired.comparator_precision_at_k)
paired["metric"]=paired.metric.replace({"main_auprc":"AUPRC_PRROC","main_precision_at_k":"precision_at_K"})
paired["higher_is_better"]=True
paired["same_candidate_denominator"]=paired.main_candidates.eq(paired.comparator_candidates)&paired.main_positive.eq(paired.comparator_positive)
paired["status_eligible"]=True
paired["comparison_policy"]=np.where(paired.comparator.str.startswith("LineaLinkR"),"matched_control","actual_method")
paired["interpretation"]=np.where(paired.task.eq("projected_fork_event"),
 "native score projected onto common LineaLinkR populations and windows","molecular pair nomination on the stated candidate universe")
paired["evaluation_source"]=CONFIG["paths"]["dynamic"]

# Availability and fixed-budget recovery have separate numerators and the same planned targets.
metrics=pd.read_csv(base/CONFIG["paths"]["lrsim"])
recovery=metrics[["dataset","n_targets","n_targets_available","n_targets_in_budget",
 "n_targets_correct_direction","target_availability","budget_recall","target_direction_accuracy",
 "n_total","n_evaluable","coverage","result_status"]].melt(
 id_vars=["dataset","n_targets","n_targets_available","n_targets_in_budget",
 "n_targets_correct_direction","n_total","n_evaluable","result_status"],
 value_vars=["target_availability","budget_recall","target_direction_accuracy","coverage"],
 var_name="metric",value_name="main_value")
recovery["main_count"]=np.select(
 [recovery.metric.eq("target_availability"),recovery.metric.eq("budget_recall"),
 recovery.metric.eq("target_direction_accuracy")],
 [recovery.n_targets_available,recovery.n_targets_in_budget,recovery.n_targets_correct_direction],
 default=recovery.n_evaluable)
recovery["main_denominator"]=np.where(recovery.metric.eq("coverage"),recovery.n_total,recovery.n_targets)
recovery["metric"]=recovery.metric.replace({"coverage":"candidate_coverage"})
recovery["task"]="native_event_target_recovery"
recovery["evaluation_scope"]="full_planned_targets"
recovery["comparator"]="descriptive_rate"
recovery["comparator_value"]=np.nan
recovery["higher_is_better"]=True
recovery["main_status"]=recovery.result_status
recovery["comparator_status"]="descriptive_rate"
recovery["same_candidate_denominator"]=True
recovery["status_eligible"]=True
recovery["comparison_policy"]="descriptive_rate"
recovery["interpretation"]="target availability; fixed-budget hits; correct local direction; native candidate coverage"
recovery["evaluation_source"]=CONFIG["paths"]["lrsim"]

# Biological readouts and peak localization retain their shared panels and units.
classic=pd.read_csv(base/CONFIG["paths"]["classic"])
specific=pd.read_csv(base/CONFIG["paths"]["specific"])
context_data=pd.concat([
 classic.loc[classic.task.isin(["clone_fate_composition","conditional_peak_localization"])],
 specific.loc[specific.task.isin(["fixed_depletion_response_enrichment",
 "known_IFN_axis_retrieval","known_IFN_ligand_retrieval"])]],ignore_index=True,sort=False)
main_by_task={"clone_fate_composition":"LLR_0.6.5_reference",
 "conditional_peak_localization":"LineaLinkR_native_pulse_fit",
 "fixed_depletion_response_enrichment":"llr_full",
 "known_IFN_axis_retrieval":"LineaLinkR_G_amplitude_Q_layout",
 "known_IFN_ligand_retrieval":"Assignment layout, geometric amplitude"}
context_data["value_upper"]=context_data.value_upper.fillna(context_data.value)
context_data["context_key"]=context_data.context_key.fillna("whole_receiver_scope")
context_data["evaluation_scope"]=context_data.evaluation_scope.fillna("declared_native_scope")
is_context_main=context_data.method.eq(context_data.task.map(main_by_task))
context_keys=["dataset","task","evaluation_scope","context_key","metric","n_reference_units"]
context_main=context_data.loc[is_context_main,context_keys+
 ["value","value_upper","n_scored_units","status","source_file","input_condition"]].rename(
 columns={"value":"main_value","value_upper":"main_value_upper","n_scored_units":"main_n_scored_units",
 "status":"main_status","source_file":"main_source_file","input_condition":"main_input_condition"})
context_other=context_data.loc[~is_context_main,context_keys+
 ["method","value","value_upper","n_scored_units","status","source_file","input_condition"]].rename(
 columns={"method":"comparator","value":"comparator_value","value_upper":"comparator_value_upper","n_scored_units":"comparator_n_scored_units",
 "status":"comparator_status","source_file":"comparator_source_file","input_condition":"comparator_input_condition"})
context_pairs=context_main.merge(context_other,on=context_keys,validate="one_to_many")
context_pairs["same_candidate_denominator"]=np.where(
 context_pairs.task.eq("clone_fate_composition")|context_pairs.metric.eq("peak_mean_distance_available"),
 context_pairs.main_n_scored_units.eq(context_pairs.comparator_n_scored_units),True)
context_pairs["rank_intervals_overlap"]=context_pairs.task.str.startswith("known_IFN")&(
 context_pairs.main_value_upper.gt(context_pairs.main_value)|context_pairs.comparator_value_upper.gt(context_pairs.comparator_value))&(
 context_pairs.main_value.le(context_pairs.comparator_value_upper)&context_pairs.comparator_value.le(context_pairs.main_value_upper))
context_pairs["status_eligible"]=~context_pairs.rank_intervals_overlap
context_pairs["higher_is_better"]=context_pairs.task.eq("fixed_depletion_response_enrichment")|(
 context_pairs.task.eq("conditional_peak_localization")&~context_pairs.metric.eq("peak_mean_distance_available"))
context_pairs["comparison_policy"]=np.where(
 context_pairs.comparator.str.lower().str.startswith(("llr_","linealinkr"))|
 context_pairs.comparator.isin(["Geometric population and layout","Geometric layout, assignment shape",
 "Assignment curve, mixed amplitude"]),"matched_control","actual_method")
context_pairs["interpretation"]=context_pairs.task.map({
 "clone_fate_composition":"shared terminal RNA decoder and matched clone-context-day units",
 "conditional_peak_localization":"available coverage and shared-target conditional localization",
 "fixed_depletion_response_enrichment":"same fixed ligand panel;receiver scope;K budget",
 "known_IFN_axis_retrieval":"rank of the experimentally supported IFN axis among 420 common pairs",
 "known_IFN_ligand_retrieval":"Ifnb1 rank among 138 common ligands"})
context_pairs["evaluation_source"]=np.where(
 context_pairs.task.isin(["clone_fate_composition","conditional_peak_localization"]),
 CONFIG["paths"]["classic"],CONFIG["paths"]["specific"])


combined=pd.concat([ordering,branch,paired,recovery,context_pairs],ignore_index=True,sort=False)
combined["dataset_scope"]=np.where(combined.dataset.str.startswith(("sim_","simulated_","lrsim_")),"synthetic","real")
combined["root_input_mode"]=combined.root_input_mode.fillna("not_applicable")
combined["comparable"]=np.isfinite(combined.main_value)&np.isfinite(combined.comparator_value)&combined.same_candidate_denominator&combined.status_eligible
# The per-dataset envelope is a separate diagnostic, chosen from actual eligible methods.
envelope=combined.loc[combined.task.isin(["cell_ordering","state_partition"])&combined.comparable].copy()
envelope["selection_value"]=np.where(envelope.higher_is_better,envelope.comparator_value,-envelope.comparator_value)
envelope=envelope.sort_values(["task","evaluation_scope","dataset","metric","selection_value","comparator"],
 ascending=[True,True,True,True,False,True]).drop_duplicates(["task","evaluation_scope","dataset","metric"])
envelope["comparison_policy"]="dataset_specific_envelope"
envelope["interpretation"]="diagnostic row-wise best external method within the declared scope"
branch_best=envelope.loc[envelope.task.eq("state_partition")].pivot(
 index=["dataset","evaluation_scope"],columns="metric",
 values=["main_value","comparator","comparator_value"])
branch_best.columns=branch_best.columns.map("_".join)
branch_best=branch_best.reset_index().rename(columns={
 "evaluation_scope":"branch_scope","main_value_ARI":"ARI","main_value_NMI":"NMI",
 "main_value_Matched_F1":"Branch_F1","comparator_ARI":"strongest_ari_method",
 "comparator_value_ARI":"strongest_ari","comparator_NMI":"strongest_nmi_method",
 "comparator_value_NMI":"strongest_nmi","comparator_Matched_F1":"strongest_f1_method",
 "comparator_value_Matched_F1":"strongest_branch_f1"})
branch_best["comparison_policy"]="dataset_specific_envelope"
branch_best.to_csv(base/("branch_vs_strongest_"+_VARIANT+".csv"),index=False,encoding="utf-8")
combined=pd.concat([combined,envelope],ignore_index=True,sort=False)
combined["main_method"]="LineaLinkR 0.6.5"
combined["raw_gap"]=combined.main_value-combined.comparator_value
combined["favorable_gap"]=np.where(combined.higher_is_better,combined.raw_gap,-combined.raw_gap)
combined["outcome"]=np.select([~combined.comparable,combined.favorable_gap.gt(CONFIG["tolerance"]),
 combined.favorable_gap.lt(-CONFIG["tolerance"])],["descriptive_or_unavailable","win","loss"],default="tie")
combined["win"]=combined.outcome.eq("win")
combined["loss"]=combined.outcome.eq("loss")
combined["tie"]=combined.outcome.eq("tie")
combined["metric_implementation"]=combined.metric.map({
 "AUPRC_PRROC":"PRROC::pr.curve auc.integral","Matched_F1":"Hungarian matching;divide by max class count",
 "RMSE":"normalized native coordinate RMSE","Spearman_signed":"native signed Spearman",
 "Spearman_absolute":"absolute native signed Spearman","Pearson":"native signed Pearson"})
combined["main_value_paired"]=combined.main_value.where(combined.comparable)
combined["comparator_value_paired"]=combined.comparator_value.where(combined.comparable)
combined["favorable_gap_paired"]=combined.favorable_gap.where(combined.comparable)
combined.sort_values(["task","evaluation_scope","metric","comparison_policy","comparator","dataset"]).to_csv(
 base/CONFIG["paths"]["comparison"],index=False,encoding="utf-8")

# Peak error is pooled by available target count, alongside dataset-level means.
peak_error=combined.task.eq("conditional_peak_localization")&combined.metric.eq("peak_mean_distance_available")
combined["main_metric_weight"]=combined.main_n_scored_units.where(peak_error)
combined["comparator_metric_weight"]=combined.comparator_n_scored_units.where(peak_error)
combined["main_metric_sum"]=combined.main_value*combined.main_metric_weight
combined["comparator_metric_sum"]=combined.comparator_value*combined.comparator_metric_weight
whole=combined.assign(summary_scope="whole_input_panel",root_input_mode="all_declared_inputs")
strata=combined.loc[combined.task.eq("cell_ordering")].assign(summary_scope="root_input_stratum")
summary_input=pd.concat([whole,strata],ignore_index=True,sort=False)
summary_input.loc[summary_input.comparison_policy.eq("dataset_specific_envelope"),"comparator"]="dataset_specific_envelope"
summary=summary_input.groupby(["task","metric","evaluation_scope","dataset_scope","summary_scope",
 "root_input_mode","comparison_policy","comparator"],dropna=False).agg(
 n_datasets=("dataset","nunique"),n_comparison_units=("dataset","size"),n_comparable=("comparable","sum"),wins=("win","sum"),ties=("tie","sum"),
 losses=("loss","sum"),mean_main_available=("main_value","mean"),mean_main_paired=("main_value_paired","mean"),
 mean_comparator_paired=("comparator_value_paired","mean"),median_favorable_gap=("favorable_gap_paired","median"),
 total_main_count=("main_count",lambda x:x.sum(min_count=1)),
 total_main_denominator=("main_denominator",lambda x:x.sum(min_count=1)),
 total_main_metric_sum=("main_metric_sum",lambda x:x.sum(min_count=1)),
 total_main_metric_weight=("main_metric_weight",lambda x:x.sum(min_count=1)),
 total_comparator_metric_sum=("comparator_metric_sum",lambda x:x.sum(min_count=1)),
 total_comparator_metric_weight=("comparator_metric_weight",lambda x:x.sum(min_count=1))).reset_index()
summary["pooled_main_rate"]=summary.total_main_count/summary.total_main_denominator
summary["pooled_main_metric"]=summary.total_main_metric_sum/summary.total_main_metric_weight
summary["pooled_comparator_metric"]=summary.total_comparator_metric_sum/summary.total_comparator_metric_weight
summary["n_not_comparable"]=summary.n_comparison_units-summary.n_comparable
summary.to_csv(base/CONFIG["paths"]["summary"],index=False,encoding="utf-8")
print(summary.loc[summary.task.eq("cell_ordering")&summary.metric.eq("Spearman_signed")&
 summary.summary_scope.eq("whole_input_panel")&summary.comparison_policy.eq("actual_method"),
 ["dataset_scope","comparator","n_datasets","n_comparable","mean_main_paired","mean_comparator_paired"]].to_string(index=False))
print(summary.loc[summary.task.eq("state_partition")&summary.summary_scope.eq("whole_input_panel")&
 summary.comparison_policy.eq("actual_method"),["metric","evaluation_scope","comparator","n_comparable"]].to_string(index=False))
print("TARGETS",int(metrics.n_targets.sum()),"SCORABLE",int(metrics.n_targets_available.sum()),
 "BUDGET_HITS",int(metrics.n_targets_in_budget.sum()))
