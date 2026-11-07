import os
CONFIG={"paths":{"root":os.environ.get("LLR_WORK_ROOT",os.getcwd()),"native":"results/by_algorithm","main_summary":"cache/main/linealinkr_summary.csv","evaluation":"results/evaluation/final065","output":"../feasible_clock/competitor_task_coverage.csv"},"threads":os.cpu_count(),"plot":{},"biological_datasets":["GSE274389_signal_recorder","GSE143437_muscle_regeneration"],"nonhuman_datasets":["celegans_embryo","monocle3_worm_l2","zebrafish_neural_dev","zebrafish_tracerseq","zebrafish_urd"],"stratified_diisco":["larry_invitro","lung","larry_cytokine_perturbation"],"excluded_methods":["FlowSig_NMF"],"own_prior_methods":["LLR_0.6.4","LineaLinkR_0.6.4","LineaLinkR_0.6.4_saved"],"exclusion_reason":"user_excluded_from_current_comparisons_and_all_main_supplementary_figures;native_files_preserved","core_tasks":{"ENTRAIN":"trajectory_associated_ligand_prioritization","NicheNet":"receiver_ligand_prioritization","scRICH_LR":"within_cell_type_communication_heterogeneity"}}
CONFIG["paths"].update({"rscript":r"D:\R\R-4.5.3\bin\Rscript.exe","native_manifest":"inputs/comparison/recovered_native/native_source_manifest.json"})
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"]=str(CONFIG["threads"])
from pathlib import Path
from itertools import chain
import numpy as np
import pandas as pd
import json
import subprocess
from io import StringIO
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG["plot"])
root=Path(CONFIG["paths"]["root"])
base=root/CONFIG["paths"]["native"]
out=root/CONFIG["paths"]["evaluation"]
diisco_whole=pd.Series(list(base.glob("DIISCO/*/coverage.csv")))
diisco_whole=diisco_whole.loc[~diisco_whole.map(lambda p:p.parent.name).isin(CONFIG["stratified_diisco"])]
diisco_units=pd.Series(list(base.glob("DIISCO/*/*/coverage.csv")))
diisco_units=diisco_units.loc[diisco_units.map(lambda p:p.parent.parent.name).isin(CONFIG["stratified_diisco"])&~diisco_units.map(lambda p:p.parent.name).str.startswith("failed_")]
paths=list(chain(base.glob("ENTRAIN/*/coverage.csv"),base.glob("NicheNet/*/coverage.csv"),base.glob("scRICH/*/coverage.csv"),diisco_whole,diisco_units,base.glob("FlowSig/*/coverage.csv"),base.glob("CellAgentChat/*/coverage.csv")))
index=pd.DataFrame({"source_id":np.arange(len(paths)),"source_file":pd.Series(paths,dtype="string")})
native=pd.concat(map(pd.read_csv,paths),keys=index.source_id,names=["source_id","native_row"]).reset_index().merge(index,on="source_id",validate="many_to_one")
native["scored_unit"]=native.result_status.str.startswith("scored")
native["clock_source"]=native.reindex(columns=["clock_source"]).iloc[:,0]
native["native_seconds"]=native.reindex(columns=["core_seconds"]).iloc[:,0].fillna(native.reindex(columns=["marker_seconds","activity_seconds"]).sum(axis=1,min_count=1)).fillna(native.reindex(columns=["preprocessing_seconds","network_seconds"]).sum(axis=1,min_count=1))
native["n_input_genes"]=native.reindex(columns=["n_input_genes"]).iloc[:,0]
status_strings=native[["dataset","method","result_status"]].drop_duplicates().sort_values(["dataset","method","result_status"]).groupby(["dataset","method"],observed=True).result_status.agg("|".join).rename("native_unit_statuses").reset_index()
completed=native.groupby(["dataset","method"],observed=True).agg(n_input_cells=("n_input_cells","first"),n_input_genes=("n_input_genes","first"),n_native_units=("result_status","size"),n_scored_units=("scored_unit","sum"),native_seconds=("native_seconds","first"),clock_source=("clock_source","first"),source_file=("source_file","first"),first_status=("result_status","first")).reset_index().merge(status_strings,on=["dataset","method"],validate="one_to_one")
completed["result_status"]=np.select([completed.n_scored_units.eq(completed.n_native_units),completed.n_scored_units.gt(0)],["completed_native","partial_native_coverage"],default=completed.first_status)
completed=completed.drop(columns="first_status")
repair_columns=["initial_result_status","initial_coverage_file","repair_detail","elbo_estimator","failure_detail"]
repair_details=native.reindex(columns=["dataset","method"]+repair_columns).groupby(["dataset","method"],observed=True).first().reset_index()
completed=completed.merge(repair_details,on=["dataset","method"],how="left",validate="one_to_one")
diisco=native.loc[native.method.eq("DIISCO")].copy()
diisco["n_scored_cells"]=diisco.n_input_cells.where(diisco.scored_unit,0)
diisco["unit_status"]=diisco.reindex(columns=["unit"]).iloc[:,0].fillna("whole_dataset").replace("","whole_dataset")+":"+diisco.result_status
diisco_counts=diisco.groupby("dataset",observed=True).agg(n_input_cells=("n_input_cells","sum"),n_scored_cells=("n_scored_cells","sum"),unit_statuses=("unit_status","|".join),native_seconds=("native_seconds","sum"),source_files=("source_file","|".join))
diisco_mask=completed.method.eq("DIISCO")
completed.loc[diisco_mask,"n_input_cells"]=completed.loc[diisco_mask,"dataset"].map(diisco_counts.n_input_cells)
completed.loc[diisco_mask,"native_seconds"]=completed.loc[diisco_mask,"dataset"].map(diisco_counts.native_seconds)
completed.loc[diisco_mask,"source_file"]=completed.loc[diisco_mask,"dataset"].map(diisco_counts.source_files)
completed["n_scored_cells"]=completed.dataset.where(diisco_mask).map(diisco_counts.n_scored_cells)
completed["unit_statuses"]=completed.dataset.where(diisco_mask).map(diisco_counts.unit_statuses)
main=pd.concat([pd.read_csv(root/CONFIG["paths"]["main_summary"]),pd.read_csv(root/"cache/final065/linealinkr_summary.csv")],ignore_index=True)[["dataset","n_cells","n_genes"]].drop_duplicates("dataset",keep="last")
bio=completed.loc[completed.dataset.isin(CONFIG["biological_datasets"]),["dataset","n_input_cells","n_input_genes"]].drop_duplicates().rename(columns={"n_input_cells":"n_cells","n_input_genes":"n_genes"})
native_datasets=completed[["dataset","n_input_cells","n_input_genes"]].rename(columns={"n_input_cells":"n_cells","n_input_genes":"n_genes"}).drop_duplicates("dataset")
datasets=pd.concat([main,bio,native_datasets],ignore_index=True).drop_duplicates("dataset")
datasets["dataset_scope"]=np.select([datasets.dataset.isin(CONFIG["biological_datasets"]),datasets.dataset.str.startswith(("sim_","simulated_","lrsim_"))],["biological_validation","synthetic_benchmark"],default="real_benchmark")
tasks=pd.Series(CONFIG["core_tasks"],name="task").rename_axis("method").reset_index()
core=datasets.merge(tasks,how="cross").merge(completed,on=["dataset","method"],how="left",validate="one_to_one")
unavailable=np.select([core.dataset_scope.eq("synthetic_benchmark"),core.dataset.isin(CONFIG["nonhuman_datasets"]),core.dataset.eq("human_fetal_pancreas"),core.dataset.eq("GSE210152_athero")&core.method.isin(["NicheNet","scRICH_LR"])],["not_applicable_synthetic_gene_ids_without_native_prior","not_applicable_native_prior_species","not_applicable_unmapped_gene_ids_and_missing_cell_type","missing_cell_type_annotations"],default="not_run_no_saved_native_output")
core["result_status"]=core.result_status.fillna(pd.Series(unavailable,index=core.index))
core["execution_scope"]="new_native_methods"
agreements=pd.read_csv(out/"signal_prioritization_agreement.csv")
agreement_status=agreements[["dataset","method","result_status"]].drop_duplicates().sort_values(["dataset","method","result_status"]).groupby(["dataset","method"],observed=True).result_status.agg("|".join).rename("comparison_status").reset_index()
core=core.merge(agreement_status,on=["dataset","method"],how="left",validate="one_to_one")
core.loc[core.method.eq("scRICH_LR")&core.result_status.eq("completed_native"),"comparison_status"]="native_heterogeneity_outputs_no_shared_label_truth"
core["score_interpretation"]=np.where(core.method.eq("scRICH_LR"),"native_communication_heterogeneity","descriptive_rank_agreement")
conditional_tasks=pd.DataFrame({"method":["FlowSig_NMF","DIISCO","CCCvelo"],"task":["condition_dependent_inflow_module_outflow","known_time_cell_type_interaction","spatial_dynamic_communication"]})
conditional=datasets.merge(conditional_tasks,how="cross").merge(completed,on=["dataset","method"],how="left",validate="one_to_one")
conditional["default_status"]=np.select([
conditional.method.eq("CCCvelo"),
conditional.method.eq("FlowSig_NMF")&conditional.dataset.isin(["GSE125527_IBD","GSE132465_CRC","GSE135893_IPF"]),
conditional.method.eq("FlowSig_NMF"),
conditional.method.eq("DIISCO")&conditional.dataset.eq("GSE274389_signal_recorder"),
conditional.method.eq("DIISCO")&conditional.dataset.eq("larry_invitro"),
conditional.method.eq("DIISCO")&conditional.dataset.eq("zebrafish_neural_dev"),
conditional.method.eq("DIISCO")&conditional.dataset.isin(["human_fetal_pancreas","zebrafish_tracerseq"])
],["missing_verified_spatial_measurements","not_run_no_saved_native_output","not_run_requires_explicit_condition_design_or_native_inputs","not_applicable_outcome_stratified_cell_sampling","requires_starting_population_stratification","not_applicable_single_cell_type","missing_cell_type_annotations"],default="missing_observed_sample_time_series")
conditional["result_status"]=conditional.result_status.fillna(conditional.default_status)
conditional["execution_scope"]="conditional_native_methods"
conditional=conditional.drop(columns="default_status")
flow_columns=["n_conditions","n_gems","n_input_flow_variables","n_retained_flow_variables","n_bootstraps","n_native_edges","n_retained_edges","condition_design"]
flow_details=native.loc[native.method.eq("FlowSig_NMF"),["dataset","method"]+flow_columns]
conditional=conditional.merge(flow_details,on=["dataset","method"],how="left",validate="one_to_one")
flow_mask=conditional.method.eq("FlowSig_NMF")
conditional.loc[flow_mask,"comparison_status"]=conditional.loc[flow_mask,"result_status"].replace({"completed_native":"native_network_summary_no_graph_ground_truth"})
conditional.loc[flow_mask,"score_interpretation"]="native_inflow_gem_outflow_network"
ticci=pd.read_csv(root/"results/evaluation/feasible_clock/ticci_reused_inputs.csv").rename(columns={"n_cells":"n_input_cells"})
ticci["method"]="TICCI"
ticci["task"]="cell_ordering"
ticci["execution_scope"]="reuse_existing_predictions"
ticci["result_status"]=np.where(ticci.n_finite.eq(ticci.n_input_cells),"reused_native_predictions","nonfinite_native_predictions")
ordering=pd.read_csv(out/"ordering_pseudotime_complete_Final065.csv")
ticci=ticci.merge(ordering.loc[ordering.algorithm.eq("TICCI"),["dataset","output_status","Spearman_signed","n_used"]],on="dataset",how="left",validate="one_to_one")
ticci["comparison_status"]=ticci.output_status.fillna("no_ordering_reference")
timetalk=pd.concat(map(pd.read_csv,(base/"reused").glob("*/TimeTalk_reused_coverage.csv")),ignore_index=True)
timetalk=timetalk.assign(task="dynamic_LR_TF_association",execution_scope="reuse_existing_predictions",comparison_status="native_output_current_reference_unavailable")
static=pd.concat(map(pd.read_csv,(base/"reused").glob("*/reused_coverage.csv")),ignore_index=True)
static=static.loc[static.method.isin(["CellChat","LIANA"])]
static=static.assign(task="static_ligand_receptor_ranking",execution_scope="reuse_existing_predictions",comparison_status="descriptive_rank_agreement")
static_scores=pd.read_csv(out/"static_ccc_agreement.csv")
liana_scores=static_scores.loc[static_scores.method.str.endswith("_via_LIANA")].copy()
liana_scores["scored_unit"]=liana_scores.result_status.eq("compared")
liana_methods=liana_scores.groupby(["dataset","method"],observed=True).agg(n_native_units=("score_type","size"),n_scored_units=("scored_unit","sum")).reset_index()
liana_methods=liana_methods.merge(static.loc[static.method.eq("LIANA"),["dataset","n_native_cells","n_native_genes","input_scope"]],on="dataset",how="left",validate="many_to_one")
liana_methods=liana_methods.assign(task="static_ligand_receptor_ranking",execution_scope="reuse_LIANA_method_scores",result_status="reused_LIANA_implementation_scores",comparison_status="descriptive_rank_agreement")
additional=completed.loc[completed.method.eq("CellAgentChat")].assign(task="static_ligand_receptor_ranking",execution_scope="new_native_method",score_interpretation="native_cellagentchat_CCI")
cellagent_comparisons=static_scores.loc[static_scores.method.eq("CellAgentChat"),["dataset","result_status"]].drop_duplicates().groupby("dataset",observed=True).result_status.agg("|".join).rename("comparison_status").reset_index()
additional=additional.merge(cellagent_comparisons,on="dataset",how="left",validate="one_to_one")
additional["comparison_status"]=additional.comparison_status.fillna("pending_native_scores")
trasig=pd.read_csv(out/"dynamic_methods.csv")
trasig=trasig.loc[trasig.task.eq("LR_pair_nomination")&trasig.method.eq("TraSig_native_max")&trasig.evaluation_scope.eq("common_candidates")].copy()
trasig=trasig.rename(columns={"evaluation_status":"comparison_status"})
trasig["method"]="TraSig"
trasig["dataset_scope"]="synthetic_benchmark"
trasig["n_input_cells"]=trasig.dataset.map(main.set_index("dataset").n_cells)
trasig["n_scored_cells"]=trasig.n_input_cells
trasig["n_native_units"]=1
trasig["n_scored_units"]=1
trasig["result_status"]=np.where(trasig.inference_status.eq("complete"),"completed_native","completed_observed_scores_only")
trasig["execution_scope"]=np.where(trasig.inference_status.eq("observed_new_fit"),"new_native_observed_scoring","reuse_existing_native_scores")
trasig["native_unit_statuses"]=trasig.inference_status
trasig["score_interpretation"]="native_LR_pair_raw_score_AP"
trasig["source_file"]=(out/"dynamic_methods").as_posix()+"/"+trasig.dataset+"/trasig_selected_task_pairs.csv"
combined=pd.concat([core,conditional,ticci,timetalk,static,liana_methods,additional,trasig],ignore_index=True,sort=False)
# Reuse the complete current ordering rows and saved author-module evaluations.
ordering_rows=pd.read_csv(out/"ordering_complete_comparison_Final065.csv")
ordering_rows=ordering_rows.loc[~ordering_rows.algorithm.str.startswith("LLR_")].copy()
ordering_rows=ordering_rows.rename(columns={"algorithm":"method","source_file":"original_prediction_source"})
ordering_rows=ordering_rows.assign(task="cell_ordering",execution_scope="reused_prediction_or_submission_metrics",active_comparison=True,scope_reason="",native_result_status="saved_ordering_metrics")
ordering_rows["native_result_status"]=np.select([ordering_rows.output_status.eq("missing_prediction"),ordering_rows.output_status.eq("nonfinite_prediction")],["not_available","saved_nonfinite_prediction"],default="saved_ordering_metrics")
ordering_rows["result_status"]=ordering_rows.output_status
ordering_rows["comparison_status"]=ordering_rows.output_status
ordering_rows["n_input_cells"]=ordering_rows.n_reference
ordering_rows["n_scored_cells"]=ordering_rows.n_used
ordering_rows["source_file"]=(out/"ordering_complete_comparison_Final065.csv").as_posix()
ordering_rows["score_interpretation"]="native_direction_correlation_and_normalized_RMSE"
branch_rows=pd.read_csv(out/"ordering_branch_complete_Final065.csv")
branch_rows=branch_rows.loc[branch_rows.algorithm.eq("slingshot")].rename(columns={"algorithm":"method"})
branch_rows=branch_rows.assign(task="branch_reconstruction",execution_scope="reused_native_lineage_weights",active_comparison=True,scope_reason="",native_result_status="saved_branch_metrics",score_interpretation="ARI_and_Hungarian_matched_branch_F1")
branch_rows["result_status"]=np.where(branch_rows.ARI.notna(),"scored","no_branch_metric_in_declared_panel")
branch_rows["comparison_status"]=branch_rows.result_status
branch_rows["source_file"]=(out/"ordering_branch_complete_Final065.csv").as_posix()
protocol_paths=list((out/"sccomm/no_TF_activity").glob("*/protocol.csv"))
sccomm=pd.concat(map(pd.read_csv,protocol_paths),ignore_index=True)
sccomm=sccomm.rename(columns={"n_cells":"n_input_cells","n_genes":"n_input_genes"})
sccomm=sccomm.assign(method="scComm",task="LR_pair_nomination",active_comparison=True,scope_reason="",execution_scope="reuse_saved_author_LR_module",native_result_status="completed_author_no_TF_activity_LR_module",result_status="completed_native")
sc_metrics=pd.read_csv(out/"dynamic_methods.csv")
sc_pair=sc_metrics.loc[sc_metrics.method.eq("scComm_no_TF_activity")&sc_metrics.task.eq("LR_pair_nomination")&sc_metrics.evaluation_scope.eq("common_candidates"),["dataset","evaluation_status"]]
sccomm=sccomm.merge(sc_pair.rename(columns={"evaluation_status":"comparison_status"}),on="dataset",how="left",validate="one_to_one")
sccomm["input_scope"]="same_cells_measured_genes_and_LR_resource; author_no_TF_activity_configuration"
sccomm["source_file"]=(out/"sccomm/no_TF_activity").as_posix()+"/"+sccomm.dataset+"/protocol.csv"
sc_projected=sccomm.copy().assign(task="projected_fork_event",execution_scope="saved_native_kernel_controlled_projection",result_status="completed_saved_projection")
sc_event=sc_metrics.loc[sc_metrics.method.eq("scComm_projected_native_kernel")&sc_metrics.evaluation_scope.eq("own_available"),["dataset","evaluation_status"]]
sc_projected=sc_projected.drop(columns="comparison_status").merge(sc_event.rename(columns={"evaluation_status":"comparison_status"}),on="dataset",how="left",validate="one_to_one")
sc_projected["source_file"]=(out/"dynamic_methods").as_posix()+"/"+sc_projected.dataset+"/projected_event_predictions.csv"
tra_projected=combined.loc[combined.method.eq("TraSig")&combined.task.eq("LR_pair_nomination")].copy()
tra_projection_status=sc_metrics.loc[sc_metrics.method.eq("TraSig_projected_contrast")&sc_metrics.evaluation_scope.eq("own_available"),["dataset","evaluation_status","n_candidates"]].rename(columns={"evaluation_status":"comparison_status","n_candidates":"n_projected_candidates"})
tra_projected=tra_projected.drop(columns="comparison_status").merge(tra_projection_status,on="dataset",how="left",validate="one_to_one")
tra_projected=tra_projected.assign(task="projected_fork_event",execution_scope="saved_native_curves_controlled_projection",native_result_status="saved_native_pair_fit_temporal_availability_per_projection")
tra_projected["result_status"]=np.where(tra_projected.n_projected_candidates.gt(0),"completed_saved_projection","no_evaluable_saved_temporal_curves")
tra_projected["n_scored_cells"]=np.nan
tra_projected["source_file"]=(out/"dynamic_methods").as_posix()+"/"+tra_projected.dataset+"/projected_event_predictions.csv"

combined=pd.concat([combined,ordering_rows,branch_rows,sccomm,sc_projected,tra_projected],ignore_index=True,sort=False).drop_duplicates(["dataset","method","task"],keep="last")

combined["n_input_cells"]=combined.n_input_cells.fillna(combined.reindex(columns=["n_cells"]).iloc[:,0]).fillna(combined.reindex(columns=["n_native_cells"]).iloc[:,0])
combined["n_input_genes"]=combined.n_input_genes.fillna(combined.reindex(columns=["n_genes"]).iloc[:,0]).fillna(combined.reindex(columns=["n_native_genes"]).iloc[:,0])
combined.loc[combined.method.eq("DIISCO"),"n_input_genes"]=np.nan
combined["dataset_scope"]=combined.dataset_scope.fillna(combined.dataset.map(datasets.set_index("dataset").dataset_scope)).fillna("historical_external")
combined["accuracy_scored"]=combined.task.isin(["cell_ordering","branch_reconstruction","LR_pair_nomination","projected_fork_event"])&combined.reindex(columns=["comparison_status"]).iloc[:,0].eq("scored")
combined["native_result_status"]=combined.native_result_status.fillna(combined.result_status)
combined["active_comparison"]=~combined.method.isin(CONFIG["excluded_methods"]+CONFIG["own_prior_methods"])
combined["scope_reason"]=np.where(combined.active_comparison,"",CONFIG["exclusion_reason"])
combined.loc[~combined.active_comparison,"execution_scope"]="historical_native_outputs_only"
combined["classic_metrics"]=combined.task.map({"cell_ordering":"Spearman_signed; Spearman_absolute; Pearson; normalized_RMSE","branch_reconstruction":"ARI; Hungarian_matched_branch_F1","LR_pair_nomination":"AUPRC; precision_at_original_K","projected_fork_event":"AUPRC; precision_at_original_K; recall_at_original_K","static_ligand_receptor_ranking":"descriptive_Spearman_rank_agreement; candidate_coverage","receiver_ligand_prioritization":"descriptive_Spearman_ligand_rank_agreement","trajectory_associated_ligand_prioritization":"descriptive_Spearman_ligand_rank_agreement"}).fillna("")
combined["context_metrics"]=np.where(combined.task.eq("projected_fork_event"),"fixed_42_target_recovery; local_direction; common_window_support","")
combined["analysis_unit"]=combined.task.map({"cell_ordering":"cell","branch_reconstruction":"cell_in_declared_all_or_terminal_subset","LR_pair_nomination":"LR_pair_per_simulation","projected_fork_event":"native_fork_context_x_LR_pair","static_ligand_receptor_ranking":"sender_receiver_LR_tuple","receiver_ligand_prioritization":"ligand_x_native_receiver","trajectory_associated_ligand_prioritization":"global_ligand","within_cell_type_communication_heterogeneity":"cell_pathway_mode_within_native_cell_type","known_time_cell_type_interaction":"sample_x_observed_cell_type","condition_dependent_inflow_module_outflow":"flow_variable_edge_by_observational_condition","spatial_dynamic_communication":"requires_verified_spatial_observations"}).fillna("")
combined["evaluation_source"]=combined.source_file
combined["unavailable_fields"]=""
combined["manuscript_use"]=np.select([combined.task.eq("cell_ordering"),combined.task.eq("branch_reconstruction"),combined.method.isin(["TraSig","scComm"]),combined.method.eq("TimeTalk"),combined.method.eq("ENTRAIN")],["ordering_Tables_1_3","branch_Table_4","simulation_Table_2_and_context_extension","fixed_receiver_response_Table_6","fixed_ligand_panels_Tables_6b_11"],default="saved_native_comparisons_pending_text_update")
combined.loc[combined.method.eq("FlowSig_NMF"),"unavailable_fields"]="native_nodes_are_inflow_GEM_outflow_variables; no_cell_population_sender_receiver_LR_fork_tuple"
combined.loc[combined.method.eq("DIISCO"),"unavailable_fields"]="native_W_f_is_sample_cell_type_interaction; no_molecular_LR_pair_or_per_cell_branch_membership"
combined.loc[combined.method.eq("ENTRAIN"),"unavailable_fields"]="global_feature_importance_has_no_native_receiver_fork_identifier_or_signed_response_window"
combined.loc[combined.method.eq("scRICH_LR"),"unavailable_fields"]="pathway_level_cell_modes_not_atomic_LR_events; canonical_WNT_and_CCL_absent_from_saved_biology_top10"

niche_response=combined.method.eq("NicheNet")&combined.dataset.eq("GSE143437_muscle_regeneration")
combined.loc[niche_response,"comparison_status"]="native_receiver_ranking_and_fixed_depletion_response_compared"
combined.loc[niche_response,"context_metrics"]="fixed_102_90_82_ligand_response_enrichment_K10_K20; rank_response_Spearman; six_balanced_sample_assignments; shared_regulatory_prior"
combined.loc[niche_response,"evaluation_source"]=(root/"results/validation/final065/tables/supp_muscle_depletion_validation.csv").as_posix()
scrich_contexts=pd.read_csv(out/"scrich/mode_context_agreement.csv")
scrich_mask=combined.method.eq("scRICH_LR")&combined.dataset.isin(scrich_contexts.dataset)
combined.loc[scrich_mask,"comparison_status"]="compared_native_mode_to_current_RNA_context"
combined.loc[scrich_mask,"classic_metrics"]="ARI; NMI_against_RNA_context_reference_descriptive"
combined.loc[scrich_mask,"context_metrics"]="within_native_type_G_Q_P_and_RNA_state_association; full_cell_matching"
combined.loc[scrich_mask,"evaluation_source"]=(out/"scrich/mode_context_agreement.csv").as_posix()
diisco_metrics=pd.read_csv(out/"diisco/native_reconstruction_metrics.csv")
diisco_mask=combined.method.eq("DIISCO")&combined.dataset.isin(diisco_metrics.dataset)
combined.loc[diisco_mask,"comparison_status"]="observed_sample_proportion_reconstruction_compared"
combined.loc[diisco_mask,"classic_metrics"]="in_sample_variational_Wf_RMSE; MAE; R2_relative_to_cell_type_mean"
combined.loc[diisco_mask,"context_metrics"]="known_observed_times; native_sample_cell_type_effects"
combined.loc[diisco_mask,"evaluation_source"]=(out/"diisco/native_reconstruction_metrics.csv").as_posix()
flow_mask=combined.method.eq("FlowSig_NMF")
combined.loc[flow_mask,"context_metrics"]="native_inflow_GEM_outflow_graph; bootstrap_edge_frequency; retained_edges; recorded_native_runtime"
combined.loc[combined.method.eq("scRICH_LR"),"unavailable_fields"]="pathway_level_cell_modes_not_atomic_LR_events"
combined.loc[combined.method.eq("scRICH_LR")&combined.dataset.isin(CONFIG["biological_datasets"]),"unavailable_fields"]="pathway_modes_not_atomic_LR_events; canonical_WNT_and_CCL_absent_from_saved_biology_top10"
new_time_inputs=combined.method.eq("DIISCO")&combined.dataset.isin(["larry_cytokine_perturbation","GSE253000_ulipstic_mLN_unsorted"])&combined.n_scored_units.isna()
combined.loc[new_time_inputs,"result_status"]="not_run_current_observed_time_input_not_adapted"
combined.loc[new_time_inputs,"native_result_status"]="no_saved_native_output"

reference_labels=pd.concat([pd.read_csv(out/"static_ccc_agreement.csv")[["dataset","reference_revision"]],pd.read_csv(out/"signal_prioritization_agreement.csv")[["dataset","reference_revision"]]],ignore_index=True).dropna().drop_duplicates("dataset").set_index("dataset").reference_revision
combined["reference_revision"]=combined.dataset.map(reference_labels)
combined.loc[combined.task.isin(["cell_ordering","branch_reconstruction"]),"reference_revision"]="fixed_time_or_branch_truth; current_Final065_comparison"
combined.loc[combined.task.eq("LR_pair_nomination"),"reference_revision"]="fixed_11_scenario_molecular_truth"
combined.loc[combined.task.eq("projected_fork_event"),"reference_revision"]="Final065_common_geometry_and_windows; fixed_42_targets"
combined.loc[scrich_mask,"reference_revision"]="Final065_G_Q_P_and_native_RNA_states"
combined.loc[diisco_mask,"reference_revision"]="native_observed_sample_proportions"
combined.loc[flow_mask,"reference_revision"]="native_observational_condition_graph"
combined.loc[combined.method.eq("scRICH_LR")&combined.dataset.eq("GSE143437_muscle_regeneration"),"unavailable_fields"]="no_atomic_LR_events; CCL_absent_from_saved_top10_pathways"
combined.loc[combined.method.eq("scRICH_LR")&combined.dataset.eq("GSE274389_signal_recorder"),"unavailable_fields"]="no_atomic_LR_events; canonical_WNT_not_selected; ncWNT_saved"
combined.loc[combined.method.eq("TimeTalk")&combined.dataset.eq("GSE143437_muscle_regeneration"),"comparison_status"]="fixed_depletion_response_compared"
combined.loc[combined.method.eq("TimeTalk")&combined.dataset.eq("GSE274389_signal_recorder"),"comparison_status"]="native_LR_TF_output_and_interpolation_verified"


# Supplied-hour BMDC native tasks use distinct LR-pair and ligand panels.
bmdc_tr_cov=pd.read_csv(base/"GSE48968_BMDC_WT/controlled_comparison_coverage.csv").iloc[0]
bmdc_en_cov=pd.read_csv(base/"ENTRAIN/GSE48968_BMDC_WT/coverage.csv").iloc[0]
bmdc_rows=pd.DataFrame([
 {"dataset":"GSE48968_BMDC_WT","dataset_scope":"controlled_paracrine_validation","method":"TraSig",
  "task":"known_IFN_axis_retrieval","active_comparison":True,"native_result_status":"completed_native",
  "execution_scope":"reuse_saved_native_scores","n_input_cells":bmdc_tr_cov.input_cells,
  "n_input_genes":bmdc_tr_cov.input_genes,"clock_source":"supplied_observed_hours_no_ordering_fit",
  "result_status":"completed_native","comparison_status":"compared_fixed_420_LR_panel",
  "score_interpretation":"rank_of_predefined_experimentally_established_LR_axis","accuracy_scored":False,
  "classic_metrics":"rank_in_fixed_420_LR_panel","context_metrics":"stimulus_specific_observed_hours; minimum_complex_subunit; author_nan2zero",
  "analysis_unit":"known_LR_axis_within_fixed_LR_candidate_panel",
  "source_file":(base/"GSE48968_BMDC_WT/controlled_comparison_coverage.csv").as_posix(),
  "evaluation_source":(root/"results/paper_tables/current/gse48968_axis_retrieval.csv").as_posix(),
  "unavailable_fields":"trajectory_ordering_not_fitted; source_or_target_type_contrast_not_available",
  "manuscript_use":"Results_2_6_Table_11_and_Figure_6","reference_revision":"Final065_known_IFN_axis_fixed_420_LR_panel"},
 {"dataset":"GSE48968_BMDC_WT","dataset_scope":"controlled_paracrine_validation","method":"ENTRAIN",
  "task":"known_IFN_ligand_retrieval","active_comparison":True,"native_result_status":"completed_native",
  "execution_scope":"reuse_saved_native_fit","n_input_cells":bmdc_en_cov.n_cells,
  "n_input_genes":bmdc_en_cov.n_genes,"native_seconds":bmdc_en_cov.seconds,
  "clock_source":bmdc_en_cov.clock_source,"input_scope":"log1p_author_TPM_over_100;399_native_ligands;500_trees;seed_12345",
  "result_status":"completed_native","comparison_status":"compared_fixed_138_ligand_panel",
  "score_interpretation":"rank_of_predefined_experimentally_established_ligand","accuracy_scored":False,
  "classic_metrics":"rank_in_fixed_138_common_ligand_panel",
  "context_metrics":"supplied_observed_hours; native_trajectory_associated_ligand_importance",
  "analysis_unit":"ligand_nomination_within_fixed_ligand_panel",
  "source_file":(base/"ENTRAIN/GSE48968_BMDC_WT/coverage.csv").as_posix(),
  "evaluation_source":(root/"results/paper_tables/current/table11.csv").as_posix(),
  "unavailable_fields":"receptor_specific_rank_not_returned; secretion_onset_not_estimated",
  "manuscript_use":"Results_2_6_and_Table_11_ligand_panel","reference_revision":"native_ENTRAIN_fixed_138_ligand_panel"}])
tt_cov=pd.read_csv(base/"TimeTalk/GSE48968_BMDC_WT/coverage.csv").iloc[0]
assert tt_cov.result_status=="completed_BMDC_MRA_autocrine"
tt_bmdc=pd.DataFrame({"task":["known_IFN_axis_retrieval","known_IFN_ligand_retrieval"],
 "n_native_units":[tt_cov.n_common_pairs,tt_cov.n_common_ligands],
 "analysis_unit":["known_LR_axis_within_fixed_LR_candidate_panel","ligand_nomination_within_fixed_ligand_panel"]})
tt_bmdc=tt_bmdc.assign(dataset="GSE48968_BMDC_WT",dataset_scope="controlled_paracrine_validation",
 method="TimeTalk",active_comparison=True,native_result_status=tt_cov.result_status,
 execution_scope="native_BMDC_RNA_TRN_MRA_and_TimeTalk",n_input_cells=tt_cov.n_cells,
 n_input_genes=tt_cov.n_genes,clock_source=tt_cov.time_input,result_status="completed_native",
 comparison_status="common_panel_target_rank_with_tie_bounds",accuracy_scored=False,
 score_interpretation="predefined_IFN_circuit_retrieval;native_SCC_and_Granger_filters",
 context_metrics="BMDC_wide_MRA;RNA_stimulus_time_phenotype;shared_zero_hour_baseline;condition_maximum",
 source_file=(base/"TimeTalk/GSE48968_BMDC_WT/coverage.csv").as_posix(),
 evaluation_source=(root/"results/paper_tables/current/table11.csv").as_posix(),
 manuscript_use="Results_2_5;Supplementary_Table_S18;Supplementary_Methods_S11",
 reference_revision="BMDC_native_TRN_MRA_fixed_420_pair_and_138_ligand_panels")
bmdc_rows=pd.concat([bmdc_rows,tt_bmdc],ignore_index=True,sort=False)
combined=pd.concat([combined,bmdc_rows],ignore_index=True).drop_duplicates(["dataset","method","task"],keep="last")
combined["manuscript_use"]=combined.manuscript_use.fillna("").str.replace("pending_text_update","reported_current_main_text_and_supplementary_data",regex=False)

# Integrate task-matched readers without refitting any saved method.
saved=out/"saved_native"
binding=pd.read_csv(saved/"native_TICCI_implementation_binding.csv")
larry_summary_path=root/"inputs/comparison/recovered_native/by_algorithm/TICCI/full_from_revision_20260919/larry_invitro/ticci_summary.csv"
larry_summary=pd.read_csv(larry_summary_path)
larry_summary["dataset"]="larry_invitro"
larry_summary["row"]=0
larry_summary["source_file"]=larry_summary_path.as_posix()
binding=pd.concat([binding.loc[~binding.dataset.eq("larry_invitro")],larry_summary],ignore_index=True,sort=False)
binding[["implementation_label","author_entropy_protocol","entropy_equivalent_to_author_protocol","saved_formula_probe_dataset","saved_formula_probe_max_abs_error"]]=binding[["implementation_label","author_entropy_protocol","entropy_equivalent_to_author_protocol","saved_formula_probe_dataset","saved_formula_probe_max_abs_error"]].ffill()
binding.to_csv(saved/"native_TICCI_implementation_binding.csv",index=False,encoding="utf-8")
part_input=pd.read_csv(saved/"native_partition_input_coverage.csv")
part_classic=pd.read_csv(saved/"classic_native_partition_metrics.csv")
part_custom=pd.read_csv(saved/"project_native_matched_branch_F1.csv")
gene_metrics=pd.read_csv(saved/"classic_native_gene_metrics.csv")
ccc_metrics=pd.read_csv(saved/"native_CCC_comparison_metrics.csv")
graph_metrics=pd.read_csv(saved/"native_graph_structure_readback.csv")
remaining=pd.read_csv(saved/"native_remaining_artifact_coverage.csv")
additional_objects=pd.read_csv(saved/"native_additional_object_readback.csv")
combined["source_series"]="existing_server_saved"
combined["project_specific_metrics"]=combined.task.map({
 "branch_reconstruction":"Hungarian_matched_Branch_F1",
 "LR_pair_nomination":"precision_at_original_K",
 "projected_fork_event":"precision_at_original_K;recall_at_original_K;fixed_42_target_recovery",
 "known_IFN_axis_retrieval":"rank_in_fixed_420_LR_panel",
 "known_IFN_ligand_retrieval":"rank_in_fixed_138_ligand_panel"}).fillna("")
combined.loc[combined.task.eq("branch_reconstruction"),"classic_metrics"]="ARI"
combined.loc[combined.task.eq("LR_pair_nomination"),"classic_metrics"]="AUPRC_PRROC_auc_integral"
combined.loc[combined.task.eq("projected_fork_event"),"classic_metrics"]="AUPRC_PRROC_auc_integral"
combined.loc[combined.task.isin(["known_IFN_axis_retrieval","known_IFN_ligand_retrieval"]),"classic_metrics"]=""
combined.loc[combined.method.eq("NicheNet")&combined.dataset.eq("GSE143437_muscle_regeneration"),"project_specific_metrics"]="fixed_depletion_response_topK_enrichment;exact_six_sample_assignments"
combined["ground_truth_accuracy"]=combined.task.isin(["LR_pair_nomination","projected_fork_event"])&combined.accuracy_scored
combined["evidence_type"]=np.select([combined.task.eq("cell_ordering"),combined.task.eq("branch_reconstruction"),combined.ground_truth_accuracy],
 ["fixed_simulated_time_or_observed_collection_time","declared_branch_labels_with_scope_in_source","simulation_truth_with_declared_candidate_scope"],default="descriptive_or_project_specific_biological_readout")
combined["denominator_note"]="see_original_evaluation_source_and_task_specific_metric_tables"
combined.loc[combined.method.eq("FlowSig_NMF")&combined.dataset.eq("GSE135893_IPF"),["result_status","comparison_status","native_result_status"]]=[
 "saved_partial_bootstrap_checkpoint","checkpoint_only_no_final_graph","partial_native_output"]
combined.loc[combined.method.eq("FlowSig_NMF"),"execution_scope"]="reuse_saved_graph_or_checkpoint;no_bootstrap_or_NMF_rerun"

part_cov=part_input.copy()
part_cov["branch_scope"]=part_cov.evaluation_scope
part_cov["task"]=np.where(part_cov.prediction_unit.eq("native_principal_tree_segment_partition"),"native_segment_partition","branch_reconstruction")
part_cov["ground_truth_accuracy"]=part_cov.evaluation_scope.eq("submission_terminal_cells")&part_cov.prediction_unit.eq("native_probability_argmax")&part_cov.exact_cell_set
part_cov["evidence_type"]=np.where(part_cov.evaluation_scope.eq("submission_terminal_cells"),"declared_terminal_branch_reference","all_cell_segment_diagnostic_or_unlabelled_real_input")
part_cov["n_input_cells"]=part_cov.n_native_cells
part_cov["n_reference_units"]=part_cov.n_scope_planned
part_cov["n_scored_cells"]=part_cov.n_scope_planned.where(part_cov.exact_cell_set)
part_cov["n_native_units"]=1
part_cov["n_scored_units"]=part_cov.n_scope_planned.notna().astype(int)
part_cov["result_status"]=np.where(part_cov.n_scope_planned.notna(),"evaluated_saved_native_partition","saved_partition_without_independent_branch_truth")
part_cov["comparison_status"]=np.where(part_cov.n_scope_planned.notna(),"full_declared_scope_classic_and_project_metrics","no_explicit_branch_truth")
part_cov["classic_metrics"]="ARI;NMI"
part_cov["project_specific_metrics"]="Hungarian_matched_Branch_F1"
part_cov["analysis_unit"]=part_cov.prediction_unit
part_cov["reference_revision"]="explicit_true_branch;terminal_time_ge_0.4_for_four_submission_datasets"
part_cov["evaluation_source"]=(saved/"classic_native_partition_metrics.csv").as_posix()
part_cov["denominator_note"]="planned_terminal_cells_or_declared_all_cell_diagnostic;missing_cells_never_silently_dropped"

gene_cov=gene_metrics.copy()
gene_cov["task"]="trajectory_gene_association"
gene_cov["n_input_cells"]=gene_cov.n_cells
gene_cov["n_input_genes"]=gene_cov.n_snapshot_genes
gene_cov["n_native_genes"]=gene_cov.n_genes
gene_cov["n_native_units"]=gene_cov.n_native_test_rows
gene_cov["n_scored_units"]=gene_cov.n_finite_p
gene_cov["n_reference_units"]=gene_cov.n_truth_genes
gene_cov["result_status"]="saved_native_gene_test_read"
gene_cov["comparison_status"]=gene_cov.truth_accuracy_status
gene_cov["classic_metrics"]="AUROC;average_precision_step;BH_precision;BH_recall_fixed_truth;BH_F1_fixed_truth;empirical_FDP_BH"
gene_cov["project_specific_metrics"]=""
gene_cov["analysis_unit"]="gene"
gene_cov["evaluation_source"]=(saved/"classic_native_gene_metrics.csv").as_posix()
gene_cov["denominator_note"]="AUC_finite_saved_p_with_truth;BH_recall_all_original_dynamic_genes;BH_over_finite_p"
gene_cov["ground_truth_accuracy"]=gene_cov.n_truth_genes.notna()
gene_cov["evidence_type"]="actual_generator_gene_class_and_present_branch_set"

ccc_cov=ccc_metrics.assign(compared=ccc_metrics.comparison_status.str.startswith("compared_descriptive",na=False)).groupby(["dataset","method"],observed=True).agg(
 n_native_score_sets=("score_field","size"),n_compared_score_sets=("compared","sum"),
 n_native_units=("n_distinct_native_returned_keys","max"),n_scored_units=("n_common_atomic_keys","max"),
 comparison_status=("comparison_status",lambda s:"|".join(sorted(s.unique()))),
 n_common_atomic_keys=("n_common_atomic_keys","max"),n_reference_units=("n_reference_atomic_keys","first")).reset_index()
native_catalog=pd.read_csv(saved/"native_R_object_catalog.csv").groupby(["dataset","method"],observed=True)[["n_cells","n_genes"]].first().reset_index()
ccc_cov=ccc_cov.merge(native_catalog,on=["dataset","method"],how="left",validate="one_to_one")
snapshot_paths=sorted((root/"inputs/comparison/recovered_native/input_metadata").glob("*/cells.txt"))
snapshot_keys=pd.Index(map(lambda p:p.parent.name,snapshot_paths))
snapshot_cells=pd.Series(list(map(lambda p:len(p.read_text(encoding="utf-8-sig").splitlines()),snapshot_paths)),index=snapshot_keys)
snapshot_genes=pd.Series(list(map(lambda p:len(p.with_name("genes.txt").read_text(encoding="utf-8-sig").splitlines()),snapshot_paths)),index=snapshot_keys)
ccc_cov["n_input_cells"]=ccc_cov.n_cells.fillna(ccc_cov.dataset.map(snapshot_cells))
ccc_cov["n_native_genes"]=ccc_cov.n_genes
ccc_cov["n_input_genes"]=ccc_cov.dataset.map(snapshot_genes)
ccc_cov["input_scope"]="native_object_cell_ids_where_stored;otherwise_declared_original_snapshot;matrix_equality_not_assumed"
ccc_binding=ccc_metrics.groupby(["dataset","method"],observed=True).agg(reference_dataset=("reference_dataset","first"),
 planned_native_population_pairs=("planned_native_population_pairs","max"),completed_native_population_pairs=("completed_native_population_pairs","max"),
 n_full_native_input_cells=("n_full_native_input_cells","first"),n_full_native_input_gene_rows=("n_full_native_input_gene_rows","first"),
 n_distinct_native_gene_names=("n_distinct_native_gene_names","first"),n_ambiguous_native_gene_names=("n_ambiguous_native_gene_names","first"),
 n_unresolved_native_gene_name_rows=("n_unresolved_native_gene_name_rows","sum"),n_invalid_native_model_score_rows=("n_invalid_native_model_score_rows","sum")).reset_index()
ccc_cov=ccc_cov.merge(ccc_binding,on=["dataset","method"],how="left",validate="one_to_one")
ccc_cov["n_input_cells"]=ccc_cov.n_input_cells.fillna(ccc_cov.n_full_native_input_cells)
ccc_cov["n_input_genes"]=ccc_cov.n_input_genes.fillna(ccc_cov.n_full_native_input_gene_rows)
ccc_cov["n_native_genes"]=ccc_cov.n_native_genes.fillna(ccc_cov.n_full_native_input_gene_rows)
sc_ccc_mask=ccc_cov.method.eq("scTenifoldXct_native")
ccc_cov.loc[sc_ccc_mask,"input_scope"]="full_original_2730_cell_3451_feature_input;native_pair_cell_IDs_and_complete_population_membership_verified;native_uppercase_gene_order_exact;original_alias_collision_retained;"+ccc_cov.loc[sc_ccc_mask,"completed_native_population_pairs"].astype(int).astype(str)+"_of_"+ccc_cov.loc[sc_ccc_mask,"planned_native_population_pairs"].astype(int).astype(str)+"_original_ordered_pairs"
model_revisions=pd.read_csv(saved/"static_reference_model_binding.csv").set_index("dataset").result_revision
ccc_cov["reference_revision"]=ccc_cov.reference_dataset.map(model_revisions)
ccc_cov["primary_dataset_key"]=ccc_cov.dataset.eq(ccc_cov.reference_dataset)
ccc_cov["task"]="static_ligand_receptor_ranking"
ccc_cov["source_file"]=(saved/"native_CCC_comparison_metrics.csv").as_posix()
ccc_cov["result_status"]="read_saved_native_scores"
ccc_cov["classic_metrics"]="descriptive_Spearman_rank_agreement"
ccc_cov["project_specific_metrics"]=""
ccc_cov["analysis_unit"]="sender_receiver_ligand_receptor_atomic_key"
ccc_cov["score_interpretation"]="agreement_with_saved_static_reference_without_GT"
ccc_cov["evaluation_source"]=ccc_cov.source_file
ccc_cov["denominator_note"]="atomic_keys_are_analysis_units;coverage_counts_are_maximum_available_key_counts_across_native_runs_and_score_fields;each_classic_metric_has_exact_native_run_and_score_field_denominators;score_set_counts_separate;partial_native_pair_plan_retained"
ccc_cov["ground_truth_accuracy"]=False
ccc_cov["evidence_type"]="descriptive_rank_agreement"

graph_cov=graph_metrics.groupby(["dataset","method"],observed=True).agg(n_native_units=("source_file","size"),n_native_rows=("n_saved_edges","sum"),source_file=("source_file","first")).reset_index()
graph_cov["method"]=graph_cov.method.replace({"CytoTalk_native_PCST":"CytoTalk_native","TICCI_native_graph":"TICCI","scFates_principal_graph":"scFates"})
graph_cov["task"]="native_graph_topology"
graph_cov["n_scored_units"]=0
graph_cov["result_status"]="saved_native_graph_read"
graph_cov["comparison_status"]="descriptive_topology_only;no_explicit_GT_graph_and_node_map"
graph_cov["classic_metrics"]=""
graph_cov["project_specific_metrics"]=""
graph_cov["analysis_unit"]="native_graph_node_and_edge"
graph_cov["evaluation_source"]=(saved/"native_graph_structure_readback.csv").as_posix()
graph_cov["denominator_note"]="all_saved_graphs;observed_nodes_only_for_edge_list_exports"
graph_cov["ground_truth_accuracy"]=False
graph_cov["evidence_type"]="native_graph_structure_without_accuracy_truth"

remaining["task"]=remaining.artifact.map({"saved_transport_map":"cell_transport_between_timepoints","max_fate_summary":"fate_confidence_summary","gene_pseudoorder":"gene_ordering","cell_pseudotime":"cell_ordering"})
other_cov=remaining.groupby(["dataset","method","task"],observed=True).agg(n_native_units=("source_file","size"),source_file=("source_file","first"),n_saved_cells=("n_saved_cells","first"),n_finite_native_scores=("n_finite_native_scores","first"),n_saved_genes=("n_saved_genes","first"),comparison_status=("evaluation_status",lambda s:"|".join(sorted(s.unique())))).reset_index()
other_cov["n_input_cells"]=other_cov.n_saved_cells
other_cov["n_input_cells"]=other_cov["n_input_cells"].where(~other_cov.method.eq("WOT"),other_cov.dataset.map(remaining.loc[remaining.method.eq("WOT")].groupby("dataset").n_input_cohort_cells.first()))
other_cov["n_native_genes"]=other_cov.n_saved_genes
other_cov["n_scored_units"]=0
other_cov["result_status"]="saved_native_smoke_artifact_read"
other_cov["classic_metrics"]=""
other_cov["project_specific_metrics"]=""
other_cov["analysis_unit"]=other_cov.task
other_cov["evaluation_source"]=(saved/"native_remaining_artifact_coverage.csv").as_posix()
other_cov["ground_truth_accuracy"]=False
other_cov["evidence_type"]="native_smoke_artifact_with_explicit_missing_truth_or_prediction_fields"
other_cov["denominator_note"]="smoke_scope_preserved;WOT_max_fate_is_not_a_full_target_probability_vector"
other_cov.loc[other_cov.method.eq("SLICER")&other_cov.dataset.eq("sim_2br_500c"),["n_scored_units","comparison_status","classic_metrics","ground_truth_accuracy"]]=[
 1,"scored_original_simulation_all_500_cells","Spearman_signed;Pearson;normalized_RMSE",True]
object_cov=additional_objects.copy()
object_cov["task"]=np.where(object_cov.method.eq("Tempora"),"cluster_pathway_trajectory","static_ligand_receptor_ranking")
object_cov["n_input_cells"]=object_cov.n_cells
object_cov["n_input_genes"]=object_cov.n_genes
object_cov["n_native_units"]=1
object_cov["n_scored_units"]=0
object_cov["result_status"]=object_cov.object_status
object_cov["comparison_status"]=object_cov.accuracy_status
object_cov["classic_metrics"]=""
object_cov["project_specific_metrics"]=""
object_cov["analysis_unit"]="saved_native_object"
object_cov["ground_truth_accuracy"]=False
object_cov["evaluation_source"]=(saved/"native_additional_object_readback.csv").as_posix()
clone_metrics=pd.read_csv(saved/"classic_native_clonal_metrics.csv")
clone_cov=clone_metrics.groupby(["dataset","method"],observed=True).agg(n_saved_cells=("n_saved_cells","first"),n_native_units=("n_planned_clones","sum"),n_context_day_groups=("context","size"),n_planned_clones=("n_planned_clones","sum"),n_truth_available=("n_truth_available","sum"),n_matched_clone_units=("n_scored","sum"),reference_protocol=("reference_protocol","first"),native_probability_policy=("native_probability_policy","first")).reset_index()
clone_cov["task"]="clone_fate_composition"
clone_cov["n_input_cells"]=clone_cov.n_saved_cells
clone_cov["n_reference_units"]=clone_cov.n_planned_clones
clone_cov["n_scored_units"]=clone_cov.n_matched_clone_units
clone_cov["result_status"]="saved_fate_probabilities_evaluated"
clone_cov["comparison_status"]="scored_on_observed_future_clonal_compositions"
clone_cov["classic_metrics"]="composition_Brier;expected_descendant_Brier;log_loss"
clone_cov["project_specific_metrics"]=""
clone_cov["analysis_unit"]="clone_within_context_and_current_day;repeated_clones_not_independent_studies"
clone_cov["input_scope"]=clone_cov.reference_protocol+";"+clone_cov.native_probability_policy+";LLR_root=observed_timepoint;CellRank2_native_connectivity_kernel"
clone_cov["denominator_note"]="planned_clone_context_day_units;future_truth_available_units;all_49302_native_input_cells"
clone_cov["ground_truth_accuracy"]=True
clone_cov["evidence_type"]="observed_future_clonal_composition;terminal_RNA_label_reference_shared"
clone_cov["evaluation_source"]=(saved/"classic_native_clonal_metrics.csv").as_posix()
phlower_log=root/"inputs/comparison/recovered_native/native_logs/misc/phlower_cellrank_stemnet.log"
phlower_report=pd.Series(phlower_log.read_text(encoding="utf-8-sig").splitlines()).str.extract(r"PHLOWER OK: ([0-9]+) cells, ([0-9]+) assigned").dropna().iloc[0].astype(int)
phlower_cov=pd.DataFrame([{"dataset":"cellrank_stemnet","method":"PHLOWER_ddhodge_component","task":"cell_ordering","n_input_cells":phlower_report.iloc[0],"n_native_units":0,"n_scored_units":0,"result_status":"native_component_success_log_only","comparison_status":"saved_prediction_file_unavailable;cell_identity_binding_unverified","classic_metrics":"","project_specific_metrics":"","analysis_unit":"expected_cell_pseudotime;log_only","ground_truth_accuracy":False,"evidence_type":"saved_stdout_of_native_ddhodge_component;not_complete_PHLOWER_workflow","input_scope":"2531_reported_cells;provided_or_fallback_groups;low_PC1_group_root;native_prediction_csv_not_recovered","denominator_note":"assigned_count_from_log;no_per_cell_values_or_ids_available_for_evaluation","source_file":phlower_log.as_posix(),"evaluation_source":phlower_log.as_posix()}])
ticci_author_cov=pd.read_csv(saved/"native_TICCI_author_task_coverage.csv")
new_cov=pd.concat([part_cov,gene_cov,ccc_cov,graph_cov,other_cov,object_cov,clone_cov,ticci_author_cov,phlower_cov],ignore_index=True,sort=False)
new_cov["source_series"]=np.where(new_cov.method.eq("TICCI_SLICE"),"authorized_native_TICCI_recompute","verified_saved_native_readback")
new_cov.loc[new_cov.branch_scope.notna(),"source_series"]="saved_native_partition_"+new_cov.loc[new_cov.branch_scope.notna(),"branch_scope"]
new_cov["execution_scope"]=np.where(new_cov.method.eq("TICCI_SLICE"),"authorized_native_TICCI_recompute_then_read_saved_outputs","read_saved_native_output_only")
new_cov["active_comparison"]=True
new_cov["scope_reason"]=""
new_cov["native_result_status"]=new_cov.result_status
new_cov["source_file"]=new_cov.source_file.fillna(new_cov.evaluation_source)
new_cov["manuscript_use"]="evaluation_only;manuscript_unchanged"
new_cov["accuracy_scored"]=new_cov.ground_truth_accuracy.fillna(False)
new_cov["dataset_scope"]=np.where(new_cov.dataset.str.startswith(("sim_","simulated_")),"synthetic_benchmark","real_or_native_smoke")
combined=pd.concat([combined,new_cov],ignore_index=True,sort=False).drop_duplicates(["dataset","method","task","source_series"],keep="last")
reference_counts=part_input.groupby("dataset",observed=True).n_reference_cells.first()
ordering_mask=combined.task.eq("cell_ordering")
combined.loc[ordering_mask,"n_input_cells"]=combined.loc[ordering_mask,"n_input_cells"].fillna(combined.loc[ordering_mask,"dataset"].map(reference_counts))
combined.loc[ordering_mask,"n_scored_cells"]=combined.loc[ordering_mask,"n_scored_cells"].fillna(combined.loc[ordering_mask,"n_aligned"])


# Long metric tables keep standard statistics separate from project-designed protocols.
order=pd.read_csv(out/"ordering_pseudotime_complete_Final065.csv").rename(columns={"algorithm":"method"})
order["task"]="cell_ordering"
order["evaluation_scope"]="fixed_ordering24_reference"
order["analysis_unit"]="cell"
fixed_reference_sizes=part_input.groupby("dataset",observed=True).n_reference_cells.first()
order["n_reference_units"]=order.n_reference.fillna(order.dataset.map(fixed_reference_sizes))
order["n_native_units"]=order.n_prediction_cells
order["n_scored_units"]=order.n_used.fillna(order.n_aligned)
order["identity_basis"]=np.where(order.n_prediction_cells.notna(),"saved_full_prediction_cell_ids_read_back","saved_metric_n_aligned;native_object_counts_preserved_when_available")
order["status"]=order.output_status
order["evidence_type"]=np.where(order.dataset.str.startswith(("sim_","simulated_")),"simulated_true_time","observed_collection_time")
order["ground_truth_accuracy"]=True
order["metric_implementation"]="native_signed_direction;separate_absolute_correlation;minmax_normalized_RMSE"
order["denominator"]="n_used_finite_aligned_cells;n_reference_retained_as_planned_denominator"
order["source_series"]="existing_server_saved"
order["source_file"]=(out/"ordering_pseudotime_complete_Final065.csv").as_posix()
order_long=order.melt(id_vars=list(order.columns.difference(["Spearman_signed","Spearman_absolute","Pearson","RMSE"])),value_vars=["Spearman_signed","Spearman_absolute","Pearson","RMSE"],var_name="metric",value_name="value")
part_classic["task"]=np.where(part_classic.prediction_unit.eq("native_principal_tree_segment_partition"),"native_segment_partition","branch_reconstruction")
part_classic["n_reference_units"]=part_classic.n_scope_planned
part_classic["n_native_units"]=part_classic.n_native_cells
part_classic["n_scored_units"]=part_classic.n_scored
part_classic["status"]=part_classic.score_status
part_classic["metric_implementation"]="sklearn_adjusted_rand_score_and_normalized_mutual_info_score"
part_classic["denominator"]="full_planned_terminal_cells_or_all_cell_diagnostic"
part_classic["source_series"]="verified_saved_native_readback"
part_long=part_classic.melt(id_vars=list(part_classic.columns.difference(["ARI","NMI"])),value_vars=["ARI","NMI"],var_name="metric",value_name="value")
gene_metrics["task"]="trajectory_gene_association"
gene_metrics["evaluation_scope"]="saved_associationTest;finite_score_AUC_and_full_truth_recall"
gene_metrics["n_reference_units"]=gene_metrics.n_truth_genes
gene_metrics["n_native_units"]=gene_metrics.n_native_test_rows
gene_metrics["n_scored_units"]=gene_metrics.n_finite_p
gene_metrics["n_full_positive_targets"]=gene_metrics.n_dynamic_truth
gene_metrics["n_scored_positive"]=gene_metrics.n_positive_finite_scores
gene_metrics["n_scored_null"]=gene_metrics.n_null_finite_scores
gene_metrics["analysis_unit"]="gene"
gene_metrics["status"]=gene_metrics.truth_accuracy_status
gene_metrics["ground_truth_accuracy"]=gene_metrics.n_truth_genes.notna()
gene_metrics["evidence_type"]="generator_dynamic_gene_truth;unobserved_branches_are_static_null"
gene_metrics["metric_implementation"]="sklearn_ROC_AUC_and_step_average_precision;BH_on_finite_saved_p;alpha_0.05"
gene_metrics["denominator"]="AUC_finite_scores_with_truth;recall_all_original_dynamic_genes;precision_all_BH_calls"
gene_metrics["source_series"]="verified_saved_native_readback"
gene_fields=["AUROC","average_precision_step","BH_precision","BH_recall_fixed_truth","BH_F1_fixed_truth","empirical_FDP_BH"]
gene_long=gene_metrics.melt(id_vars=list(gene_metrics.columns.difference(gene_fields)),value_vars=gene_fields,var_name="metric",value_name="value")
ccc_metrics["task"]="static_ligand_receptor_ranking"
ccc_metrics["reference_revision"]=ccc_metrics.reference_dataset.map(model_revisions)
ccc_metrics["context_key"]=ccc_metrics.native_run+"|"+ccc_metrics.score_field
ccc_metrics["evaluation_scope"]="exact_complete_cohort_and_population_binding;common_atomic_keys"
ccc_metrics["n_reference_units"]=ccc_metrics.n_reference_atomic_keys
ccc_metrics["n_native_units"]=ccc_metrics.n_distinct_native_returned_keys
ccc_metrics["n_scored_units"]=ccc_metrics.n_common_atomic_keys
ccc_metrics["analysis_unit"]="sender_receiver_ligand_receptor_atomic_key"
ccc_metrics["status"]=ccc_metrics.comparison_status
ccc_metrics["evidence_type"]="descriptive_agreement_with_saved_LLR_static_scores"
ccc_metrics["metric_implementation"]="Pearson_correlation_of_average_tie_ranks"
ccc_metrics["denominator"]="common_atomic_keys;both_full_native_and_reference_key_denominators_saved"
ccc_metrics["source_series"]="verified_saved_native_readback"
ccc_metrics["source_file"]=(saved/"native_CCC_comparison_metrics.csv").as_posix()
ccc_long=ccc_metrics.rename(columns={"Spearman_rank_agreement":"value"}).assign(metric="Spearman_rank_agreement")
signal=pd.read_csv(out/"signal_prioritization_agreement.csv")
signal["context_key"]=signal.receiver.astype(str)+"|"+signal.score_field
signal["analysis_unit"]="ligand_within_native_receiver_or_global_ligand"
signal["n_reference_units"]=signal.n_reference_ligands
signal["n_native_units"]=signal.n_native_ligands
signal["n_scored_units"]=signal.n_common
signal["status"]=signal.result_status
signal["ground_truth_accuracy"]=False
signal["evidence_type"]="descriptive_ligand_rank_agreement"
signal["metric_implementation"]="Spearman_on_common_native_ligands"
signal["denominator"]="common_native_ligands;full_native_and_reference_ligand_counts_retained"
signal["source_file"]=(out/"signal_prioritization_agreement.csv").as_posix()
signal["source_series"]="existing_server_saved"
signal_long=signal.rename(columns={"spearman_agreement":"value"}).assign(metric="Spearman_ligand_rank_agreement")
scr=scrich_contexts.copy()
scr["context_key"]=scr.pathway+"|"+scr.scRICH_group+"|"+scr.reference
scr["task"]="within_cell_type_communication_heterogeneity"
scr["n_reference_units"]=scr.n_cells
scr["n_native_units"]=scr.n_cells
scr["n_scored_units"]=scr.n_cells
scr["analysis_unit"]="cell_within_native_population_and_pathway"
scr["status"]=scr.result_status
scr["evidence_type"]="descriptive_mode_to_RNA_context_association"
scr["metric_implementation"]="ARI_and_NMI"
scr["denominator"]="all_matched_cells_in_named_native_population_and_pathway"
scr["source_series"]="existing_server_saved"
scr["source_file"]=(out/"scrich/mode_context_agreement.csv").as_posix()
scr_long=scr.melt(id_vars=list(scr.columns.difference(["ARI","NMI"])),value_vars=["ARI","NMI"],var_name="metric",value_name="value")
di=diisco_metrics.copy()
di["method"]="DIISCO"
di["context_key"]=di.unit_dir
di["n_reference_units"]=di.n_rows
di["n_native_units"]=di.n_rows
di["n_scored_units"]=di.n_rows
di["analysis_unit"]="sample_by_cell_type"
di["status"]="in_sample_reconstruction"
di["ground_truth_accuracy"]=False
di["evidence_type"]="in_sample_factorized_variational_W_times_f"
di["metric_implementation"]="RMSE_MAE_and_R2_relative_to_cell_type_mean"
di["denominator"]="all_observed_sample_cell_type_proportions;no_heldout_prediction"
di["source_series"]="existing_server_saved"
di["source_file"]=(out/"diisco/native_reconstruction_metrics.csv").as_posix()
di_fields=["fitted_rmse","fitted_mae","r2_over_cell_type_mean"]
di_long=di.melt(id_vars=list(di.columns.difference(di_fields)),value_vars=di_fields,var_name="metric",value_name="value")
dyn=pd.read_csv(out/"dynamic_methods.csv")
dyn["n_reference_units"]=dyn.n_candidates
dyn["n_native_units"]=dyn.n_candidates
dyn["n_scored_units"]=dyn.n_candidates
dyn["n_full_positive_targets"]=dyn.n_targets
dyn["n_scored_positive"]=dyn.n_positive
dyn["n_scored_null"]=dyn.n_null
dyn["analysis_unit"]=np.where(dyn.task.eq("LR_pair_nomination"),"LR_pair","native_fork_context_by_LR_pair")
dyn["status"]=dyn.evaluation_status
dyn["ground_truth_accuracy"]=True
dyn["evidence_type"]="simulation_truth_on_declared_own_or_common_candidate_scope"
dyn["metric_implementation"]="PRROC_pr_curve_auc_integral;distinct_from_step_average_precision"
dyn["denominator"]="AUPRC_on_scored_candidate_scope;positive_and_null_counts_explicit"
dyn["source_series"]="existing_server_saved"
dyn["source_file"]=(out/"dynamic_methods.csv").as_posix()
dyn_long=dyn.rename(columns={"auprc":"value"}).assign(metric="AUPRC_PRROC_auc_integral")
sl=pd.read_csv(saved/"classic_native_SLICER_ordering.csv")
sl["task"]="cell_ordering"
sl["n_reference_units"]=sl.n_reference
sl["n_native_units"]=sl.n_reference
sl["n_scored_units"]=sl.n_used
sl["status"]="full_original_500_cell_simulation"
sl["analysis_unit"]="cell"
sl["ground_truth_accuracy"]=True
sl["evidence_type"]="simulated_true_time"
sl["denominator"]="all_500_original_simulated_cells;native_signed_no_truth_flip"
sl["source_series"]="verified_saved_native_readback"
sl["source_file"]=(saved/"classic_native_SLICER_ordering.csv").as_posix()
sl_long=sl.melt(id_vars=list(sl.columns.difference(["Spearman_signed","Pearson","RMSE"])),value_vars=["Spearman_signed","Pearson","RMSE"],var_name="metric",value_name="value")
old_static=static_scores.copy()
old_static["context_key"]=old_static.score_type+"|"+old_static.score_field
old_static["n_reference_units"]=old_static.n_reference_keys
old_static["n_native_units"]=old_static.n_native_keys
old_static["n_scored_units"]=old_static.n_common
old_static["analysis_unit"]="sender_receiver_ligand_receptor_atomic_key"
old_static["status"]=old_static.result_status
old_static["ground_truth_accuracy"]=False
old_static["evidence_type"]="descriptive_saved_static_rank_agreement"
old_static["denominator"]="common_atomic_keys;full_native_and_reference_key_counts_retained"
old_static["source_series"]="existing_server_saved"
old_static["source_file"]=(out/"static_ccc_agreement.csv").as_posix()
old_static_long=old_static.rename(columns={"spearman_agreement":"value"}).assign(metric="Spearman_rank_agreement")
clone_metrics["context_key"]=clone_metrics.context+"|day="+clone_metrics.day.astype(str)
clone_metrics["analysis_unit"]=clone_metrics.unit
clone_metrics["n_reference_units"]=clone_metrics.n_planned_clones
clone_metrics["n_native_units"]=clone_metrics.n_planned_clones
clone_metrics["n_scored_units"]=clone_metrics.n_scored
clone_metrics["ground_truth_accuracy"]=True
clone_metrics["status"]="scored_on_matched_observed_future_clones"
clone_metrics["evidence_type"]=clone_metrics.ground_truth_type
clone_metrics["metric_implementation"]="multiclass_composition_Brier;expected_descendant_multiclass_Brier;cross_entropy"
clone_metrics["denominator"]=clone_metrics.reference_protocol+";"+clone_metrics.native_probability_policy+";planned_and_truth_available_counts_retained"
clone_metrics["source_series"]="verified_saved_native_readback"
clone_metrics["source_file"]=(saved/"classic_native_clonal_metrics.csv").as_posix()
clone_long=clone_metrics.melt(id_vars=list(clone_metrics.columns.difference(["soft_brier","descendant_brier","log_loss"])),value_vars=["soft_brier","descendant_brier","log_loss"],var_name="metric",value_name="value")
clone_long["metric"]=clone_long.metric.replace({"soft_brier":"composition_Brier","descendant_brier":"expected_descendant_Brier"})
ticci_author=pd.read_csv(saved/"classic_TICCI_author_ordering_metrics.csv").rename(columns={"algorithm":"method"})
ticci_author["task"]="cell_ordering"
ticci_author["analysis_unit"]="cell_native_pseudotime"
ticci_author["n_reference_units"]=ticci_author.n_reference
ticci_author["n_native_units"]=ticci_author.n_prediction_cells
ticci_author["n_scored_units"]=ticci_author.n_used
ticci_author["status"]=ticci_author.output_status
ticci_author["ground_truth_accuracy"]=ticci_author.n_reference.notna()
ticci_author["metric_implementation"]="signed_Spearman;Pearson;RMSE_after_separate_minmax_on_finite_matched_cells"
ticci_author["denominator"]="original_native_cells;reference_cells;finite_matched_cells_reported_separately;no_truth_flip"
ticci_author["evidence_type"]="fixed_observed_collection_time"
ticci_author["input_condition"]="provided_input_population_labels_for_saved_CellChat;functional_entropy_auto_root"
ticci_author["source_series"]="authorized_native_TICCI_recompute"
ticci_author_long=ticci_author.melt(id_vars=list(ticci_author.columns.difference(["Spearman_signed","Pearson","RMSE"])),value_vars=["Spearman_signed","Pearson","RMSE"],var_name="metric",value_name="value")
full_order=pd.read_csv(saved/"classic_full_native_ordering_metrics.csv").rename(columns={"algorithm":"method"})
full_time_reference=pd.read_csv(root/"inputs/comparison/results/paper_tables/completed_real_ordering.csv").groupby("dataset")[["n_reference","n_reference_known_time"]].first()
full_order["task"]="cell_ordering"
full_order["analysis_unit"]="cell_native_full_workflow_pseudotime"
full_order["n_reference_units"]=full_order.n_reference
full_order["n_native_units"]=full_order.n_prediction_cells
full_order["n_scored_units"]=full_order.n_used
full_order["status"]=full_order.output_status
full_order["ground_truth_accuracy"]=full_order.n_reference_known_time.fillna(0).gt(0)
full_order["metric_implementation"]="same_saved_prediction_scoring_operator;signed_Spearman;absolute_correlation_separate;Pearson;separate_minmax_RMSE;no_truth_flip"
full_order["denominator"]="all_original_prediction_cells;all_reference_cells;known_time_reference_cells;finite_matched_cells_reported_separately"
full_order["evidence_type"]=np.where(full_order.n_reference_known_time.fillna(0).gt(0),"observed_collection_time","no_independent_observed_time_reference")
full_order["source_series"]="current_full_native_task_workflow"
full_order_long=full_order.melt(id_vars=list(full_order.columns.difference(["Spearman_signed","Spearman_absolute","Pearson","RMSE"])),value_vars=["Spearman_signed","Spearman_absolute","Pearson","RMSE"],var_name="metric",value_name="value")
classic_long=pd.concat([order_long,part_long,gene_long,ccc_long,signal_long,scr_long,di_long,dyn_long,sl_long,old_static_long,clone_long,ticci_author_long,full_order_long],ignore_index=True,sort=False)
time_files=list((out/"dynamic_methods").glob("*/temporal_peak_window_metrics.csv"))
time_native=pd.concat(map(pd.read_csv,time_files),ignore_index=True)
time_native["n_full_positive_targets"]=time_native.groupby(["dataset","method"]).n_targets.transform("max")
time_native["peak_available_fraction"]=time_native.n_available/time_native.n_targets
time_native["peak_plateau_hit_fraction"]=time_native.n_hit/time_native.n_targets
time_native=time_native.rename(columns={"scope":"evaluation_scope","n_targets":"n_reference_units",
 "n_available":"n_scored_units","mean_distance":"peak_mean_distance_available"})
time_long=time_native.melt(id_vars=["dataset","method","evaluation_scope","n_reference_units",
 "n_scored_units","n_hit","n_full_positive_targets"],
 value_vars=["peak_available_fraction","peak_plateau_hit_fraction","peak_mean_distance_available"],
 var_name="metric",value_name="value")
time_long["task"]="conditional_peak_localization"
time_long["context_key"]="generating_plateau_0.62_0.68"
time_long["analysis_unit"]="planted_transient_target_in_truth_mapped_native_branch"
time_long["status"]="available_native_peak_summary"
time_long["evidence_type"]="conditional_localization_with_declared_branch_and_time_reference"
time_long["ground_truth_accuracy"]=True
time_long["metric_implementation"]="coverage;plateau hits;distance to generating plateau"
time_long["denominator"]="scope-specific targets;full positive-target count and available estimates reported separately"
time_long["input_condition"]="same_G_populations_Q_bins_and_evaluation_time_calibration"
time_long["source_file"]=(out/"dynamic_methods").as_posix()+"/"+time_long.dataset+"/temporal_peak_window_recovery.csv"
time_long["source_series"]="native_temporal_peak_readout"
classic_long=pd.concat([classic_long,time_long],ignore_index=True,sort=False)

metric_columns=["dataset","method","task","evaluation_scope","context_key","metric","value","value_upper","analysis_unit","n_reference_units","n_native_units","n_scored_units","n_full_positive_targets","n_scored_positive","n_scored_null","status","evidence_type","ground_truth_accuracy","metric_implementation","denominator","source_series","source_file","reference_revision","input_condition","root_input_mode","direction_policy","identity_basis","reference_dataset","primary_dataset_key","native_run","score_field","n_truth_classes","n_prediction_classes","matched_f1_sum","metric_denominator","implementation_label"]
classic_long=classic_long.loc[~classic_long.method.isin(CONFIG["own_prior_methods"])].copy()
classic_long["implementation_label"]=classic_long.method.map({"TICCI":"TICCI expression-Shannon adapter + CellChat + author getLineage + Scanpy DPT","TICCI_SLICE":"TICCI SLICE functional entropy B100;CCI k1,resolution1;current environment"}).fillna(classic_long.method)
metric_columns+=["implementation_repair","n_empty_suffix_sentinels","runtime_scope","method_family","implementation_id","comparison_stage","family_result_status","authoritative_implementation"]
repaired_classic=classic_long.reindex(columns=["implementation_repair"]).iloc[:,0].fillna("none").ne("none")
classic_long.loc[repaired_classic,"implementation_label"]=classic_long.loc[repaired_classic,"implementation_label"]+";explicit_empty_suffix_sentinel_repair"
classic_long.reindex(columns=metric_columns).to_csv(out/"classic_metrics.csv",index=False,encoding="utf-8")
part_custom["task"]=np.where(part_custom.prediction_unit.eq("native_principal_tree_segment_partition"),"native_segment_partition","branch_reconstruction")
part_custom["analysis_unit"]="cell_partition_with_optimal_label_matching"
part_custom["n_reference_units"]=part_custom.n_scope_planned
part_custom["n_native_units"]=part_custom.n_native_cells
part_custom["n_scored_units"]=part_custom.n_scored
part_custom["status"]=part_custom.score_status
part_custom["metric_denominator"]=np.maximum(part_custom.n_truth_classes,part_custom.n_prediction_classes)
part_custom["source_series"]="verified_saved_native_readback"
part_custom_long=part_custom.rename(columns={"Hungarian_matched_Branch_F1":"value"}).assign(metric="Hungarian_matched_Branch_F1")
dyn["metric_implementation"]="original_budget_K_equals_full_positive_target_count;unavailable_candidates_not_zero_filled"
dyn["denominator"]="original_K_and_all_original_targets;scored_candidate_scope_preserved"
dyn_custom=dyn.melt(id_vars=list(dyn.columns.difference(["precision_at_k","recall_at_k"])),value_vars=["precision_at_k","recall_at_k"],var_name="metric",value_name="value")
muscle_reference_file=root/"results/validation/final065/tables/supp_muscle_depletion_validation.csv"
entrain_response_file=muscle_reference_file.with_name("supp_muscle_depletion_ligand_level.csv")
entrain_response=pd.read_csv(entrain_response_file).loc[lambda panel: panel.method.eq("entrain")].assign(
 method="ENTRAIN",context="common_ligands",source_file=entrain_response_file.as_posix())
mus=pd.concat([pd.read_csv(muscle_reference_file).assign(source_file=muscle_reference_file.as_posix()),
 entrain_response],ignore_index=True,sort=False)
mus["dataset"]="GSE143437_muscle_regeneration"
mus["task"]="fixed_depletion_response_enrichment"
mus["context_key"]=mus.context+"|K="+mus.k.astype(str)
mus["analysis_unit"]="ligand_in_fixed_receiver_panel"
mus["n_reference_units"]=mus.n_candidates
mus["n_native_units"]=mus.n_available
mus["n_scored_units"]=mus.k
mus["ground_truth_accuracy"]=False
mus["status"]="fixed_response_panel;four_biological_samples;six_balanced_assignments"
mus["evidence_type"]="experimental_response_association_with_shared_regulatory_prior"
mus["denominator"]="fixed_panel_candidates_and_original_K;available_ligands_and_biological_sample_counts_in_source"
mus["source_series"]="existing_server_saved"
mus_custom=mus.melt(id_vars=list(mus.columns.difference(["enrichment","exact_sample_p"])),value_vars=["enrichment","exact_sample_p"],var_name="metric",value_name="value")
entrain_native_dir=root/"results/by_algorithm/ENTRAIN/GSE143437_muscle_regeneration"
entrain_native=combined.loc[combined.dataset.eq("GSE143437_muscle_regeneration")
 &combined.method.eq("ENTRAIN")&combined.task.eq("trajectory_associated_ligand_prioritization")]
entrain_panel=entrain_native.assign(
 task="fixed_depletion_response_enrichment",analysis_unit="ligand_in_fixed_receiver_panel",
 source_file=entrain_response_file.as_posix(),evaluation_source=entrain_response_file.as_posix(),
 original_prediction_source=(entrain_native_dir/"ligand_importance.csv").as_posix(),
 n_native_units=pd.read_csv(entrain_native_dir/"coverage.csv").n_ligands.iloc[0],
 n_scored_units=entrain_response.n_available.iloc[0],n_reference_units=entrain_response.n_candidates.iloc[0],
 input_scope="GSE143437_native_ligand_ranks;GSE291091_day3_FAP_depletion",
 reference_revision="GSE291091_day3_FAP_depletion",clock_source="TimeTalk",native_seconds=np.nan,
 classic_metrics="",context_metrics="",project_specific_metrics="enrichment;exact_sample_p",
 score_interpretation="top_ligand_enrichment_for_independent_depletion_response",
 comparison_status="fixed_response_panel;four_biological_samples;six_balanced_assignments",
 evidence_type="experimental_response_association_with_shared_regulatory_prior",
 accuracy_scored=False,ground_truth_accuracy=False,
 denominator_note="75 common ligands; K=10,20; four biological samples; six balanced assignments",
 manuscript_use="Table 2;Table S11")
combined=pd.concat([combined,entrain_panel],ignore_index=True,sort=False)
axis=pd.read_csv(root/"results/paper_tables/current/gse48968_axis_retrieval.csv")
axis=axis.loc[axis.ligand_genes.str.upper().eq("IFNB1")&axis.receptor_genes.str.upper().str.contains("IFNAR1")&axis.receptor_genes.str.upper().str.contains("IFNAR2")].copy()
axis["dataset"]="GSE48968_BMDC_WT"
axis["task"]="known_IFN_axis_retrieval"
axis["context_key"]=axis.stimulation.astype("string").fillna("not_stored")+"|"+axis.ligand_genes+"|"+axis.receptor_genes
axis["analysis_unit"]="known_LR_axis"
axis["n_reference_units"]=axis.n_common_candidates
axis["n_native_units"]=axis.n_candidates
axis["n_scored_units"]=axis.n_common_candidates.where(axis.rank_common.notna())
axis["ground_truth_accuracy"]=False
axis["status"]="fixed_known_axis_retrieval;not_global_causal_accuracy"
axis["evidence_type"]="predefined_experimentally_established_IFN_axis"
axis["denominator"]="fixed_420_common_LR_panel;full_native_candidate_count_retained"
axis["source_file"]=(root/"results/paper_tables/current/gse48968_axis_retrieval.csv").as_posix()
axis["source_series"]="existing_server_saved"
axis_custom=axis.rename(columns={"rank_common":"value"}).assign(metric="known_axis_rank_in_common_420_panel")
entrain_axis=pd.read_csv(root/"results/paper_tables/current/table11.csv").loc[lambda d:d["Ligand rank / 138"].notna()].rename(columns={"Method":"method","Ligand rank / 138":"value"})
entrain_axis=entrain_axis.assign(dataset="GSE48968_BMDC_WT",task="known_IFN_ligand_retrieval",metric="known_ligand_rank_in_common_138_panel",analysis_unit="known_ligand",n_reference_units=138,ground_truth_accuracy=False,status="fixed_known_IFNB1_ligand_panel",source_series="existing_server_saved",source_file=(root/"results/paper_tables/current/table11.csv").as_posix(),denominator="fixed_138_common_ligands",evidence_type="predefined_experimentally_established_IFN_axis")
tt_axis=pd.read_csv(root/"results/paper_tables/current/table11.csv").loc[lambda d:d.Method.eq("TimeTalk")].rename(
 columns={"Method":"method","Pair rank / 420":"value","Pair rank upper":"value_upper"})
tt_axis=tt_axis.assign(dataset="GSE48968_BMDC_WT",task="known_IFN_axis_retrieval",
 metric="known_axis_rank_in_common_420_panel",context_key="all_primary_conditions|IFNB1|IFNAR1_IFNAR2",analysis_unit="known_LR_axis",n_reference_units=420,
 n_native_units=tt_cov.n_input_molecular_pairs,n_scored_units=420,ground_truth_accuracy=False,
 status="fixed_IFNB1_IFNAR1_IFNAR2_rank_bounds",source_series="native_BMDC_TRN_MRA",
 source_file=(base/"TimeTalk/GSE48968_BMDC_WT/TimeTalk_common_pair_ranking.csv").as_posix(),
 denominator="fixed_420_common_pairs;value_and_value_upper_bound_tied_ranks",
 evidence_type="predefined_experimentally_established_IFN_axis")
entrain_axis["value_upper"]=entrain_axis["Ligand rank upper"]
custom_long=pd.concat([part_custom_long,dyn_custom,mus_custom,axis_custom,entrain_axis,tt_axis],ignore_index=True,sort=False)
custom_long=custom_long.loc[~custom_long.method.isin(CONFIG["own_prior_methods"])].copy()
custom_long["implementation_label"]=custom_long.method.replace({"llr_previous":"current-fit matched scoring control"})
custom_long.reindex(columns=metric_columns).to_csv(out/"project_specific_metrics.csv",index=False,encoding="utf-8")
print("classic_metric_rows",len(classic_long),"project_specific_metric_rows",len(custom_long),flush=True)

combined.loc[combined.task.eq("cell_ordering"),"ground_truth_accuracy"]=combined.loc[combined.task.eq("cell_ordering"),"comparison_status"].str.startswith("scored",na=False)
combined.loc[combined.source_series.eq("existing_server_saved")&combined.task.eq("branch_reconstruction"),"ground_truth_accuracy"]=combined.loc[combined.source_series.eq("existing_server_saved")&combined.task.eq("branch_reconstruction"),"branch_scope"].eq("submission_terminal_cells")

# Pending scopes are registered without launching or scoring a new method run.
pending_ticci=pd.DataFrame([dict(dataset="human_fetal_pancreas",method="TICCI_SLICE",
 task="cell_ordering",source_series="current_author_protocol_pending_full_input",execution_scope="ongoing_author_protocol_full_input",
 n_input_cells=28430,n_native_units=0,n_scored_units=0,n_reference_units=26290,
 result_status="running_author_protocol_full_input",
 native_result_status="running_author_protocol_full_input",
 comparison_status="full_author_trajectory_pending;not_mouse_pancreas_3696",
 input_scope="full_28430_cell_human_fetal_pancreas;32955_original_gene_rows_preserved;unique_resource_ID_binding;native_Gap_K10_B50",
 ground_truth_accuracy=False,evaluation_source=(saved/"native_TICCI_author_task_coverage.csv").as_posix(),
 source_file=(root/"inputs/comparison/recovered_native/input_metadata/human_fetal_pancreas/meta.csv").as_posix(),
 next_action="finish_active_author_Gap_and_CellChat_then_native_trajectory_and_evaluation",
 classic_metrics="",project_specific_metrics="")])
pending_input_binding=pd.read_csv(base/"TICCI_author_functional_entropy/human_fetal_pancreas/input_binding.csv")
pending_ticci["n_input_genes"]=int(pending_input_binding.n_genes.iloc[0])
pending_ticci["n_reference_units"]=pending_ticci.dataset.map(full_time_reference.n_reference)
pending_ticci["n_reference_known_time"]=pending_ticci.dataset.map(full_time_reference.n_reference_known_time)
completed_ticci_datasets=pd.Index(pd.Series(list((base/"TICCI_author_functional_entropy").glob("*/ticci_author_result.h5ad"))).map(lambda p:p.parent.name))
pending_ticci=pending_ticci.loc[~pending_ticci.dataset.isin(completed_ticci_datasets)]
cyto_pairs=pd.read_csv(saved/"native_CytoTalk_pair_coverage.csv")
cyto_pairs["planned_pair"]=cyto_pairs.planned_pair.fillna(False).astype(bool)
cyto_pairs["saved_network"]=cyto_pairs.saved_network.fillna(False).astype(bool)
cyto_pairs["valid_network"]=cyto_pairs.valid_native_network.eq(True)
cyto_pairs["failed_model"]=cyto_pairs.native_model_status.str.startswith("failed_",na=False)
cyto_pairs["unverified_model"]=cyto_pairs.native_model_status.eq("saved_native_network_finiteness_unverified")
cyto_pairs["missing_model"]=cyto_pairs.native_model_status.eq("missing_native_network")
cyto_scope=cyto_pairs.groupby("dataset",observed=True).agg(
 planned_units=("planned_pair","sum"),saved_units=("saved_network","sum"),
 valid_units=("valid_network","sum"),failed_units=("failed_model","sum"),
 unverified_units=("unverified_model","sum"),missing_units=("missing_model","sum"),
 n_input_cells=("input_cells","first"),input_scope=("scope_status","first")).reset_index()
missing_cyto=cyto_pairs.loc[cyto_pairs.planned_pair&cyto_pairs.missing_model].groupby("dataset").pair_key.agg(";".join)
failed_cyto=cyto_pairs.loc[cyto_pairs.planned_pair&cyto_pairs.failed_model].groupby("dataset").pair_key.agg(";".join)
cyto_scope["missing_unit_keys"]=cyto_scope.dataset.map(missing_cyto).fillna("")
cyto_scope["failed_unit_keys"]=cyto_scope.dataset.map(failed_cyto).fillna("")
cyto_scope["planned_units"]=cyto_scope.planned_units.mask(cyto_scope.dataset.eq("paul15_smoke"))
cyto_scope["missing_units"]=cyto_scope.missing_units.mask(cyto_scope.planned_units.isna())
repaired_cyto=cyto_pairs.loc[cyto_pairs.valid_network&cyto_pairs.shared_failure_status.notna()].groupby("dataset",observed=True).size()
cyto_scope["numerically_repaired_units"]=cyto_scope.dataset.map(repaired_cyto).fillna(0).astype(int)
cyto_scope["implementation_repair"]=np.where(cyto_scope.numerically_repaired_units.gt(0),"zero_degree_generalized_inverse;all_gene_nodes_retained;positive_degree_formula_unchanged","none")
cyto_scope["method"]="CytoTalk_native"
cyto_scope["task"]="native_cell_type_pair_network"
cyto_scope["source_series"]="original_native_population_pair_plan_readback"
cyto_scope["n_native_units"]=cyto_scope.valid_units
cyto_scope["n_scored_units"]=0
cyto_scope["n_reference_units"]=cyto_scope.planned_units
cyto_scope["result_status"]=np.select([
 cyto_scope.planned_units.isna(),cyto_scope.unverified_units.gt(0),cyto_scope.missing_units.gt(0),cyto_scope.failed_units.gt(0)],
 ["smoke_original_plan_unverified","saved_native_finiteness_unverified","partial_native_pair_plan","complete_original_plan_accounted_with_native_failures"],
 default="complete_valid_native_pair_plan")
cyto_scope["native_result_status"]=cyto_scope.result_status
cyto_scope["comparison_status"]=cyto_scope.valid_units.astype(str)+"_valid_native_networks;"+cyto_scope.failed_units.astype(str)+"_native_failures;"+cyto_scope.saved_units.astype(str)+"_files_saved;no_graph_accuracy_truth"
cyto_scope["ground_truth_accuracy"]=False
cyto_scope["denominator_note"]="original unordered population-pair plan;gene_population_graph_vertices_are_not_cells;file_presence_and_valid_model_outputs_counted_separately"
cyto_scope["evaluation_source"]=(saved/"native_CytoTalk_pair_coverage.csv").as_posix()
cyto_scope["source_file"]=cyto_scope.evaluation_source
cyto_scope["classic_metrics"]=""
cyto_scope["project_specific_metrics"]=""
cyto_scope["next_action"]=np.select([cyto_scope.missing_units.gt(0),cyto_scope.failed_units.gt(0),cyto_scope.planned_units.isna()],
 ["complete_only_missing_original_pairs","retain_failed_native_normalization_status_in_full_plan_denominator","recover_original_smoke_input_plan"],
 default="")
resource_check=pd.read_csv(saved/"ticci_nonhuman_input_resource_check.csv")
resource_check["method"]="TICCI_SLICE"
resource_check["task"]="cell_ordering"
resource_check["source_series"]="current_author_resource_applicability"
resource_check["n_native_units"]=0
resource_check["n_reference_units"]=resource_check.n_reference_cells
resource_check["n_scored_units"]=0
resource_check["ground_truth_accuracy"]=False
resource_check["result_status"]="not_applicable_with_distributed_author_functional_resources"
resource_check["native_result_status"]=resource_check.result_status
resource_check["comparison_status"]=resource_check.reason
resource_check["input_scope"]="all_original_cells_and_gene_IDs_checked;no_human_mouse_Kappa_substitution_or_Shannon_replacement"
resource_check["classic_metrics"]=""
resource_check["project_specific_metrics"]=""
resource_check["analysis_unit"]="author_functional_entropy_and_species_matched_LR_resource_contract"
resource_check["evaluation_source"]=(saved/"ticci_nonhuman_input_resource_check.csv").as_posix()
resource_check["denominator_note"]="native_full_input_size_retained;not_run_or_scored_when_functional_gene_annotation_contract_unavailable"
resource_check["next_action"]=""
combined=pd.concat([combined,pending_ticci,cyto_scope,resource_check],ignore_index=True,sort=False)
wot_files=list((base/"WOT").glob("*/coverage.csv"))
wot_index=pd.DataFrame({"source_id":np.arange(len(wot_files)),"evaluation_source":pd.Series(wot_files,dtype="string")})
wot_native=pd.concat(map(pd.read_csv,wot_files),keys=wot_index.source_id,names=["source_id","native_row"]).reset_index().merge(wot_index,on="source_id",validate="many_to_one")
wot_native["n_native_genes"]=wot_native.reindex(columns=["n_model_genes"]).iloc[:,0]
wot_native["source_series"]="current_author_WOT_protocol"
wot_native["execution_scope"]="new_native_author_protocol_after_adapter_repair"
wot_native["native_result_status"]=wot_native.result_status
wot_native["input_scope"]="full_input_cell_ID_binding;native_gene_filter;unit_growth_multiplier_1;time_origin="+wot_native.time_origin
wot_native["ground_truth_accuracy"]=False
wot_native["n_scored_units"]=0
wot_native["classic_metrics"]=""
wot_native["project_specific_metrics"]=""
wot_transport=wot_native.assign(task="cell_transport_between_timepoints",
 n_native_units=wot_native.n_transport_maps,analysis_unit="native_map_between_two_observed_times",
 comparison_status="all_saved_transport_entries_finite;no_cross_time_cell_identity_GT",
 denominator_note="native_maps_on_full_time_observed_cells;mass_entries_are_not_independent_biological_replicates")
wot_fates=wot_native.assign(task="target_identified_fate_probability",
 n_native_units=wot_native.n_native_fates,analysis_unit="input_cell_by_terminal_RNA_target",
 comparison_status="full_named_target_probability_vector_saved;independent_fate_GT_not_scored",
 denominator_note="all_time_observed_cells_by_all_named_terminal_targets;unknown_times_preserved_in_full_input_binding")
combined=pd.concat([combined,wot_transport,wot_fates],ignore_index=True,sort=False)
full_files=list((base/"PHLOWER").glob("*/coverage.csv"))+list((base/"Tempora").glob("*/coverage.csv"))+list((base/"SLICER").glob("*/coverage.csv"))+list((base/"GeneTrajectory").glob("*/coverage.csv"))
full_index=pd.DataFrame({"source_id":np.arange(len(full_files)),"evaluation_source":pd.Series(full_files,dtype="string")})
full_native=pd.concat(map(pd.read_csv,full_files),keys=full_index.source_id,names=["source_id","native_row"]).reset_index().merge(full_index,on="source_id",validate="many_to_one")
full_native["source_series"]="current_full_native_task_workflow"
full_native["execution_scope"]="new_native_full_input_workflow"
full_native["native_result_status"]=full_native.result_status
full_native["ground_truth_accuracy"]=False
full_native["n_scored_units"]=0
full_native["classic_metrics"]=""
full_native["project_specific_metrics"]=""
full_native["source_file"]=full_native.evaluation_source
ph_full=full_native.loc[full_native.method.eq("PHLOWER")].copy()
ph_full["input_scope"]="all_original_cells_and_genes_bound;author_HVG2000;10000_native_trajectories;root_lowest_median_PC1_provided_group;no_reference_orientation"
ph_complete=ph_full.result_status.eq("completed_full_native_workflow")
ph_repaired=ph_full.reindex(columns=["n_empty_suffix_sentinels"]).iloc[:,0].fillna(0).gt(0)
ph_full["initial_result_status"]=np.where(ph_repaired,"failed_native_harmonic_tree_zero_weight_node","")
ph_full["initial_coverage_file"]=np.where(ph_repaired,(base/"PHLOWER/dentategyrus/phlower_original_failure_coverage.csv").as_posix(),"")
ph_full["repair_detail"]=ph_full.reindex(columns=["implementation_repair"]).iloc[:,0].fillna("none")
ph_full["input_scope"]=ph_full.input_scope+";implementation_repair="+ph_full.repair_detail
ph_full["comparison_status"]=np.where(ph_complete,"complete_native_tree_and_projection;tree_reference_unavailable;clock_reference_checked_separately","native_harmonic_tree_failure;ddhodge_checkpoint_not_complete_PHLOWER_result")
ph_full["next_action"]=np.where(ph_complete,"","retain_native_zero_weight_tree_node_failure;do_not_impute_tree_coordinates_or_score_partial_component")
ph_order=ph_full.assign(task="cell_ordering",n_native_units=ph_full.n_finite_pseudotime.fillna(0),
 analysis_unit="all_original_cell_IDs_with_native_S0_pseudotime",denominator_note="full_input_cells;provided_RNA_groups_are_root_inputs_not_independent_time_truth")
ph_tree=ph_full.assign(task="native_lineage_tree",n_native_units=ph_full.n_native_tree_edges.fillna(0),
 analysis_unit="native_harmonic_tree_edge",denominator_note="all_original_cells_and_native_tree;tree_segments_are_not_terminal_fates_or_independent_biological_replicates")
te_full=full_native.loc[full_native.method.eq("Tempora")].copy()
te_full["task"]="cluster_pathway_trajectory"
te_full["n_native_units"]=te_full.n_native_edges
te_full["analysis_unit"]="provided_cluster_by_pathway_trajectory_with_observed_time_input"
te_full["clock_source"]="observed_collection_stage_used_as_input"
te_full["input_scope"]="all_9815_cells_and_53801_genes;5_provided_celltypes;7_observed_stages;50_Mouse_Hallmark_pathways;native_GSVA_and_3PC_trajectory"
te_full["comparison_status"]="complete_native_cluster_pathway_graph;known_time_input;no_independent_graph_accuracy_truth"
te_full["denominator_note"]="full_input_cells;5_RNA_groups_are_analysis_nodes;observed_time_used_by_method_not_blind_cell_ordering_truth"
te_full["next_action"]=""
sc_native=pd.read_csv(base/"scTenifoldXct/paul15/coverage.csv")
sc_scope=sc_native.groupby(["dataset","method"],observed=True).agg(n_input_cells=("n_input_cells","first"),
 n_input_genes=("n_input_genes","first"),valid_units=("source_celltype","size"),
 n_native_units=("n_native_candidates","sum"),native_seconds=("native_seconds","sum")).reset_index()
sc_scope["planned_units"]=342
sc_scope["missing_units"]=sc_scope.planned_units-sc_scope.valid_units
sc_scope["task"]="ordered_population_pair_LR_nomination"
sc_scope["source_series"]="current_full_native_ordered_pair_workflow"
sc_scope["execution_scope"]="new_native_author_workflow_after_adapter_repair"
sc_scope["result_status"]=np.where(sc_scope.missing_units.eq(0),"complete_native_ordered_pair_plan","partial_native_ordered_pair_plan_running")
sc_scope["native_result_status"]=sc_scope.result_status
sc_scope["n_scored_units"]=0
sc_scope["ground_truth_accuracy"]=False
sc_scope["classic_metrics"]=""
sc_scope["project_specific_metrics"]=""
sc_scope["input_scope"]="full_2730_cells_and_3451_genes;19_provided_populations;342_ordered_pairs;native_nComp5;1000_steps;seed12345;native_seed0_preserved"
sc_scope["analysis_unit"]="ordered_population_pair_by_native_LR_candidate"
sc_scope["comparison_status"]="saved_native_aligned_distance_and_significance;no_independent_molecular_truth"
sc_scope["denominator_note"]="342_ordered_population_pairs;original_all_cells_and_genes_retained;candidate_scores_and_significant_calls_saved_separately"
sc_scope["evaluation_source"]=(base/"scTenifoldXct/paul15/coverage.csv").as_posix()
sc_scope["source_file"]=sc_scope.evaluation_source
sc_scope["next_action"]=np.where(sc_scope.missing_units.gt(0),"continue_active_remaining_ordered_pairs_reusing_saved_GRNs_and_completed_pair_bundles","evaluate_only_verified_common_molecular_candidate_keys")
sl_full=full_native.loc[full_native.method.eq("SLICER")].copy()
sl_full["task"]="cell_ordering"
sl_full["n_native_units"]=sl_full.n_finite_pseudotime
sl_full["input_scope"]="all_original_cells_and_genes_bound;same_saved_native_lognorm;native_gene_selection;complete_K5_to_500_by5_search"
sl_full["analysis_unit"]="all_original_cell_IDs_with_native_cell_order"
sl_full["comparison_status"]="complete_native_cell_order;reference_availability_checked_separately"
sl_full["denominator_note"]="full_input_cells;native_extreme_cell_root;no_reference_orientation"
sl_full["next_action"]=""
native_order=pd.concat([ph_order,sl_full],ignore_index=True,sort=False)
native_order=native_order.merge(full_order[["dataset","method","n_reference","n_reference_known_time","n_used","output_status"]],on=["dataset","method"],how="left",validate="one_to_one")
native_order["n_reference_units"]=native_order.n_reference.fillna(native_order.dataset.map(full_time_reference.n_reference))
native_order["n_reference_known_time"]=native_order.n_reference_known_time.fillna(native_order.dataset.map(full_time_reference.n_reference_known_time))
native_order["n_scored_units"]=native_order.n_used.fillna(0)
native_order["ground_truth_accuracy"]=native_order.n_reference_known_time.fillna(0).gt(0)&native_order.n_used.fillna(0).gt(0)
native_order["comparison_status"]=native_order.output_status.fillna(native_order.comparison_status)
native_order["classic_metrics"]=np.where(native_order.n_used.fillna(0).gt(0),"Spearman_signed;Spearman_absolute;Pearson;normalized_RMSE","")
gene_full=full_native.loc[full_native.method.eq("GeneTrajectory")].copy()
gene_full["task"]="gene_ordering"
gene_full["n_native_units"]=gene_full.reindex(columns=["n_returned_genes"]).iloc[:,0]
gene_full["input_scope"]="all_original_cells_and_genes_bound;author_2000_gene_selection;all_1999000_native_pair_transports;unchanged_EMD50000_budget"
gene_full["analysis_unit"]="gene_with_native_trajectory_pseudoorder"
gene_full["comparison_status"]="full_native_gene_trajectory;no_independent_gene_onset_truth;not_cell_ordering_or_gene_DE_pvalues"
gene_full["denominator_note"]="native_2000_gene_panel_and_full_original_input_genes_retained;pair_transports_not_independent_biological_samples"
gene_full["next_action"]=""
pending_full=pd.DataFrame([
 dict(dataset="paul15",method="GeneTrajectory",task="gene_ordering",n_input_cells=2730,n_input_genes=3451,n_native_genes=2000,
  planned_units=2000,n_native_gene_pairs=1999000,analysis_unit="gene_native_pseudoorder",
  input_scope="full_2730_cell_input;3451_original_features;author_2000_gene_panel;all_1999000_OT_pairs;50000_iteration_budget;full_result_not_saved"),
 dict(dataset="paul15",method="SLICER",task="cell_ordering",n_input_cells=2730,n_input_genes=3451,n_native_genes=370,
  planned_units=2730,analysis_unit="cell_native_order",
  input_scope="full_saved_native_input;native_370_gene_selection;original_K5_to_500_by5_search;full_result_not_saved"),
 dict(dataset="dentategyrus",method="SLICER",task="cell_ordering",n_input_cells=2930,n_input_genes=13913,
  planned_units=2930,analysis_unit="cell_native_order",
  input_scope="full_saved_native_input;queued_after_paul15;native_gene_selection;original_K5_to_500_by5_search;full_result_not_saved")])
pending_full=pending_full.merge(full_native[["dataset","method"]].drop_duplicates().assign(native_result_saved=True),on=["dataset","method"],how="left",validate="one_to_one")
pending_full=pending_full.loc[pending_full.native_result_saved.isna()].copy()
pending_full["source_series"]="current_full_native_scope_pending"
pending_full["execution_scope"]="authorized_full_native_run_in_progress_or_queued"
pending_full["n_native_units"]=0
pending_full["n_scored_units"]=0
pending_full["n_reference_units"]=pending_full.dataset.map(full_time_reference.n_reference).where(pending_full.task.eq("cell_ordering"))
pending_full["n_reference_known_time"]=pending_full.dataset.map(full_time_reference.n_reference_known_time).where(pending_full.task.eq("cell_ordering"))
pending_full["ground_truth_accuracy"]=False
pending_full["result_status"]="pending_complete_native_output"
pending_full["native_result_status"]=pending_full.result_status
pending_full["comparison_status"]="no_complete_native_prediction_saved;not_scored"
pending_full["classic_metrics"]=""
pending_full["project_specific_metrics"]=""
pending_full["denominator_note"]="original_full_input_and_declared_native_model_panel_retained;intermediate_computations_not_counted_as_completed_predictions"
pending_full["evaluation_source"]=pending_full.method.map({"GeneTrajectory":(root/"experiments/evaluation/run_genetrajectory.py").as_posix(),"SLICER":(root/"experiments/evaluation/run_slicer_full.R").as_posix()})
pending_full["source_file"]=pending_full.evaluation_source
pending_full["next_action"]="await_existing_full_native_job_then_read_saved_result;do_not_restart_completed_work"
combined=pd.concat([combined,native_order,ph_tree,te_full,gene_full,sc_scope,pending_full],ignore_index=True,sort=False)
combined["next_action"]=combined.reindex(columns=["next_action"]).iloc[:,0].fillna("")
combined.loc[combined.method.eq("PHLOWER_ddhodge_component"),"next_action"]="recover_full_native_prediction_and_cell_ID_binding;component_log_is_not_full_PHLOWER_workflow"
combined.loc[combined.method.eq("WOT")&~combined.source_series.eq("current_author_WOT_protocol"),"next_action"]="historical_nonfinite_maps_preserved;adapter_repair_completed_on_same_bound_500_cell_probe_and_full_observed_time_input"
combined.loc[combined.method.isin(["SLICER","GeneTrajectory","Tempora"])&~combined.source_series.isin(["current_full_native_task_workflow","current_full_native_scope_pending"]),"next_action"]="saved_smoke_scope_only;retain_declared_cells_genes_and_task;not_full_main_panel_completion"
combined.loc[combined.method.eq("DIISCO")&combined.dataset.isin(["larry_cytokine_perturbation","GSE253000_ulipstic_mLN_unsorted"])&combined.n_scored_units.isna(),"next_action"]="applicable_observed_time_input_not_adapted_or_run;sample_by_cell_type_output_only"
combined.loc[combined.method.eq("DIISCO")&combined.result_status.isin(["missing_cell_type_annotations","missing_observed_time"]),"next_action"]="required_native_input_information_missing;not_replace_with_pseudotime_or_inferred_labels"

# Match competitors to the object actually returned, including implementation variants.
objects=pd.DataFrame([
 ["clone_fate_composition","clone_context_day_fate_distribution","same_future_clonal_composition_with_shared_terminal_RNA_decoder","Do saved native fates predict the same observed future clone compositions?"],
 ["cell_ordering","cell_time","same_output","Is the native ordering consistent with observed collection or simulated time?"],
 ["branch_reconstruction","cell_fate_membership","same_output_with_declared_terminal_scope","Does the method recover the same terminal branch labels?"],
 ["native_segment_partition","principal_tree_segment","different_from_terminal_fate","How does the native segment partition relate to declared branch labels?"],
 ["trajectory_gene_association","gene_association_pvalue","auxiliary_gene_task_LinkR_has_no_gene_pvalue_output","Which genes vary along the saved trajectory, with original gene-truth recall?"],
 ["native_cell_type_pair_network","native_population_pair_PCST_network","original_unordered_population_pair_plan;descriptive_without_graph_GT","Which networks exist for the original complete cell-type pair plan?"],
 ["native_lineage_tree","native_harmonic_principal_tree","full_native_tree_requires_independent_matched_tree_GT","What complete native tree and cell projections were saved?"],
 ["ordered_population_pair_LR_nomination","ordered_population_pair_by_native_LR_candidate","native_molecular_candidate_distance_and_significance","Which native LR candidates and calls exist for the complete ordered population-pair plan?"],
 ["native_graph_topology","native_graph_nodes_and_edges","descriptive_without_matched_GT_graph","What native topology was saved and which cells or nodes can be bound?"],
 ["LR_pair_nomination","molecular_LR_pair","same_output_on_common_candidate_panel","Which true molecular pairs are ranked within the original candidate and budget scope?"],
 ["projected_fork_event","fork_context_by_LR_pair","controlled_projection_not_native_competitor_event_output","Can saved native pair curves support the same fixed fork-event targets after declared projection?"],
 ["static_ligand_receptor_ranking","sender_receiver_LR_tuple","descriptive_common_keys_without_accuracy_truth","How do native static scores order the same verified atomic LR keys?"],
 ["fixed_depletion_response_enrichment","ligand_within_fixed_receiver_response_panel","project_specific_experimental_response_with_four_biological_samples","Do fixed ligand priorities enrich measured depletion response under all six balanced sample assignments?"],
 ["receiver_ligand_prioritization","ligand_within_receiver","fixed_receiver_response_panel","Do ligand priorities agree with measured receiver response in the same fixed panel?"],
 ["trajectory_associated_ligand_prioritization","global_ligand","global_priority_control_requires_fixed_response_panel","Does the native global ligand ranking recover the fixed biological response or known-axis panel?"],
 ["within_cell_type_communication_heterogeneity","cell_pathway_mode","auxiliary_mode_not_LR_event_or_fate_prediction","What communication roles coexist within the same annotated cell population?"],
 ["known_time_cell_type_interaction","sample_by_cell_type_proportion","in_sample_reconstruction_not_LR_accuracy","How well does the saved variational W-times-f reconstruct observed sample compositions?"],
 ["spatial_dynamic_communication","spatial_dynamic_LR_field","not_applicable_without_verified_spatial_measurements","Are the required measured spatial coordinates present?"],
 ["dynamic_LR_TF_association","native_LR_TF_association","native_response_readout_not_fork_event_output","Does the saved native association support the fixed response or known signaling axis?"],
 ["known_IFN_axis_retrieval","predefined_IFNB1_IFNAR_axis","project_specific_known_axis_rank","Where is the established IFN axis ranked in the fixed 420-pair panel?"],
 ["known_IFN_ligand_retrieval","predefined_IFNB1_ligand","project_specific_known_ligand_rank","Where is the established IFNB1 ligand ranked in the fixed 138-ligand panel?"],
 ["cell_transport_between_timepoints","cell_to_cell_transport_mass","saved_smoke_nonfinite_transport","Is a finite native transport matrix available on the original input cells?"],
 ["fate_confidence_summary","maximum_fate_only","insufficient_output_for_fate_accuracy","Was a full target-identified fate probability vector saved?"],
 ["target_identified_fate_probability","cell_by_named_terminal_target_probability","full_native_probability_requires_matched_future_outcome_GT","Are full target-identified native fate vectors available on the declared input cells?"],
 ["gene_ordering","gene_pseudoorder","smoke_subset_without_gene_onset_GT","Is an independent gene-onset reference available for the saved 16-gene subset?"],
 ["cluster_pathway_trajectory","cluster_pathway_graph","smoke_without_independent_trajectory_GT","What cluster-level pathway trajectory was saved?"]
],columns=["task","predicted_object","comparison_role","question_answered"])
combined=combined.merge(objects,on="task",how="left",validate="many_to_one")
combined["implementation_label"]=combined.method
combined.loc[combined.source_series.eq("current_author_WOT_protocol"),"implementation_label"]="WOT author OTModel and full native fates;growth multiplier 1;unchanged native convergence budget"
combined.loc[combined.source_series.eq("current_author_WOT_protocol")&combined.task.eq("cell_transport_between_timepoints"),"comparison_role"]="finite_native_transport;no_matched_cross_time_cell_identity_GT"
combined.loc[combined.method.eq("TICCI"),"implementation_label"]="TICCI expression-Shannon adapter + CellChat + author getLineage + Scanpy DPT"
combined.loc[combined.method.eq("TICCI_SLICE"),"implementation_label"]="TICCI SLICE functional entropy B100;CCI k1,resolution1;current environment"
combined.loc[combined.method.eq("TICCI_CellChat"),"implementation_label"]="CellChat stage saved by TICCI adapter"
combined.loc[combined.method.eq("TICCI_topology_only_historical"),"implementation_label"]="TICCI historical graph-only stage;no saved cell time"
combined.loc[combined.method.eq("TICCI"),"input_scope"]=combined.loc[combined.method.eq("TICCI"),"input_scope"].fillna("")+";expression_Shannon_entropy_replaces_author_SLICE_functional_entropy"
history_datasets=pd.read_csv(saved/"native_graph_structure_readback.csv").loc[lambda d:d.method.eq("TICCI_topology_only_historical"),"dataset"]
history_time=combined.method.eq("TICCI")&combined.task.eq("cell_ordering")&combined.dataset.isin(history_datasets)&combined.comparison_status.eq("missing_prediction")
combined.loc[history_time,"native_result_status"]="historical_graph_only_output"
combined.loc[history_time,"comparison_status"]="missing_full_workflow_pseudotime;historical_topology_saved"
combined.loc[combined.method.eq("TICCI")&combined.dataset.eq("dentategyrus")&combined.task.eq("cell_ordering"),"unavailable_fields"]="34_infinite_DPT_distances_preserved;not_a_missing_time_parser"
combined.loc[combined.method.eq("PHLOWER_ddhodge_component"),"implementation_label"]="PHLOWER native ddhodge component;complete workflow not established"
combined.loc[combined.method.eq("PHLOWER_ddhodge_component"),"comparison_role"]="success_log_only;no_saved_prediction_to_score"
current_ph_datasets=pd.Index(ph_full.loc[ph_complete,"dataset"])
historical_ph_component=combined.method.eq("PHLOWER_ddhodge_component")&combined.dataset.isin(current_ph_datasets)
combined.loc[historical_ph_component,"comparison_role"]="historical_ddhodge_log;current_full_author_workflow_available"
combined.loc[historical_ph_component,"next_action"]="current_full_native_result_available;old_component_CSV_not_recovered"
combined.loc[combined.method.eq("CytoTalk_native")&combined.dataset.isin(repaired_cyto.index),"implementation_label"]="CytoTalk author PCST with explicit zero-degree generalized inverse on repaired pairs;all original nodes and null budgets retained"
combined.loc[combined.method.eq("TICCI_SLICE")&combined.dataset.eq("dentategyrus")&combined.task.eq("cell_ordering"),"unavailable_fields"]="34_Cajal_Retzius_cells_have_author_defined_infinite_DPT_distance;2896_native_root_component_finite"

combined.loc[combined.method.eq("GeneTrajectory")&combined.source_series.isin(["current_full_native_task_workflow","current_full_native_scope_pending"]),"comparison_role"]="native_gene_ordering;independent_gene_onset_reference_unavailable"
combined.loc[combined.method.eq("GeneTrajectory")&combined.source_series.isin(["current_full_native_task_workflow","current_full_native_scope_pending"]),"question_answered"]="What native gene order is available for the full original input and declared 2000-gene panel?"
combined.loc[combined.source_series.eq("current_full_native_task_workflow")&combined.method.eq("Tempora"),"comparison_role"]="full_known_time_cluster_pathway_graph;not_blind_cell_clock_recovery"
combined.loc[combined.source_series.eq("current_full_native_task_workflow")&combined.method.eq("PHLOWER"),"implementation_label"]="PHLOWER full author Hodge trajectories, harmonic tree and STREAM cell projection"
ph_repaired_scope=combined.method.eq("PHLOWER")&combined.reindex(columns=["n_empty_suffix_sentinels"]).iloc[:,0].fillna(0).gt(0)
combined.loc[ph_repaired_scope,"implementation_label"]=combined.loc[ph_repaired_scope,"implementation_label"]+";explicit repair of unobserved suffix sentinel leaves"
combined["active_comparison"]=~combined.method.isin(CONFIG["excluded_methods"]+CONFIG["own_prior_methods"])
combined.loc[historical_ph_component,"active_comparison"]=False
repaired_wot_datasets=pd.Index(wot_native.loc[wot_native.result_status.eq("completed_native"),"dataset"])
historical_wot=combined.method.eq("WOT")&~combined.source_series.eq("current_author_WOT_protocol")&combined.dataset.isin(repaired_wot_datasets)
combined.loc[historical_wot,"active_comparison"]=False
combined.loc[historical_wot,"comparison_role"]="historical_invalid_smoke_artifact;current_finite_transport_available"
combined.loc[historical_wot,"next_action"]="use_current_author_WOT_protocol_record;historical_transport_files_preserved"
diagnostic_wot=combined.source_series.eq("current_author_WOT_protocol")&combined.reindex(columns=["time_origin"]).iloc[:,0].str.startswith("historical_",na=False)
combined.loc[diagnostic_wot,"active_comparison"]=False
combined.loc[diagnostic_wot,"scope_reason"]="adapter_repair_probe_only;historical_time_bins_not_verified_as_observed_experimental_time"
combined.loc[combined.method.isin(CONFIG["excluded_methods"]),"scope_reason"]=CONFIG["exclusion_reason"]
combined.loc[combined.method.isin(CONFIG["excluded_methods"]),"comparison_role"]="excluded_by_user;native_files_preserved"
combined["biological_line"]=combined.dataset.map({
 "GSE143437_muscle_regeneration":"muscle_regeneration_and_independent_FAP_depletion",
 "GSE274389_signal_recorder":"WNT_reporter_and_signal_history",
 "GSE48968_BMDC_WT":"BMDC_IFN_axis_and_perturbation",
 "lung":"lung_repair_with_protein_and_aging_followups",
 "larry_invitro":"LARRY_clonal_fate",
 "larry_cytokine_perturbation":"LARRY_clonal_fate",
 "GSE253000_ulipstic_mLN_unsorted":"uLIPSTIC_independent_contact"})
# Two saved lung exports are supplementary contexts, separate from the fixed RNA benchmark.
lung_joint=pd.read_csv(root/"results/by_algorithm/TraSig/GSE141259_joint_lung_epithelial/coverage.csv").iloc[0]
lung_whole=pd.read_csv(root/"results/by_algorithm/TraSig/GSE141259_lung_whole/coverage.csv").iloc[0]
lung_auxiliary=pd.DataFrame([
 dict(dataset="GSE141259_joint_lung_epithelial",n_input_cells=lung_joint.all_input_cells,n_native_units=lung_joint.all_contexts,
      source_file="results/by_algorithm/TraSig/GSE141259_joint_lung_epithelial/coverage.csv",
      original_prediction_source="results/by_algorithm/TraSig/GSE141259_joint_lung_epithelial/all_native_dot_scores.npz",
      input_scope="all_61856_cells;900_role_state_contexts;30_observed_time_paths;1160_compatible_LR_pairs"),
 dict(dataset="GSE141259_lung_whole",n_input_cells=lung_whole.input_cells,n_input_genes=lung_whole.input_genes,n_native_units=lung_whole.native_contexts,
      source_file="results/by_algorithm/TraSig/GSE141259_lung_whole/coverage.csv",
      original_prediction_source="results/by_algorithm/TraSig/GSE141259_lung_whole/all_native_dot_scores.npz",
      input_scope="all_29297_cells;6400_condition_state_contexts;80_observed_time_paths;1177_monomer_LR_pairs;window20")])
lung_auxiliary=lung_auxiliary.assign(dataset_scope="supplementary_lung_context",method="TraSig_lung_observed_dot",
 task="auxiliary_lung_observed_dot",source_series="saved_lung_observed_dot_exports",
 result_status="saved_observed_dot_export;permutations_run0",native_result_status="saved_observed_dot_export;permutations_run0",
 active_comparison=True,ground_truth_accuracy=False,comparison_role="supplementary_native_context",
 comparison_status="descriptive_native_context_coverage",score_interpretation="observed_TraSig_paired_dot",
 analysis_unit="author_state_condition_path_context",clock_source="observed_collection_dates",
 predicted_object="observed_dot_scores_on_author_state_paths",implementation_label="TraSig observed-dot; calendar paths; zero permutations",
 manuscript_use="supplementary_native_coverage",scope_reason="auxiliary_GSE141259_inputs_separate_from_lung_and_primary39",
 biological_line="lung_RNA_context",next_action="retain_saved_native_export_and_record_scope")
combined=pd.concat([combined,lung_auxiliary],ignore_index=True,sort=False)
# Method families identify authors' algorithms; implementation IDs retain native stages and projections.
family_alias={"TICCI_SLICE":"TICCI","TICCI_topology_only_historical":"TICCI",
 "TICCI_CellChat":"CellChat","PHLOWER_ddhodge_component":"PHLOWER",
 "CellPhoneDB_native":"CellPhoneDB","CellPhoneDB_via_LIANA":"CellPhoneDB",
 "Connectome_via_LIANA":"Connectome","NATMI_native":"NATMI","NATMI_via_LIANA":"NATMI",
 "SingleCellSignalR_native":"SingleCellSignalR","SingleCellSignalR_via_LIANA":"SingleCellSignalR",
 "CytoTalk_native":"CytoTalk","LIANA_rank_aggregate":"LIANA","scFates_native_segments":"scFates",
 "scRICH_LR":"scRICH","scTenifoldXct_native":"scTenifoldXct","FlowSig_NMF":"FlowSig"}
family_alias.update({"TraSig_native_max":"TraSig","TraSig_projected_contrast":"TraSig",
 "TraSig_native_dot_same_resource":"TraSig","TraSig paired dot":"TraSig","TraSig-style paired dot":"TraSig",
 "TraSig_native_curve":"TraSig","TraSig_lung_observed_dot":"TraSig","scComm_native_kernel":"scComm",
 "scComm_no_TF_activity":"scComm","scComm_projected_native_kernel":"scComm",
 "timetalk":"TimeTalk","timetalk_scc":"TimeTalk","nichenet":"NicheNet",
 "Assignment layout, geometric amplitude":"LineaLinkR",
 "Geometric layout, assignment shape":"LineaLinkR",
 "Geometric population and layout":"LineaLinkR","Assignment curve, mixed amplitude":"LineaLinkR",
 "pooled_LR":"LR score reference","maximum_observed_mean_LR_product":"LR score reference",
 "pooled_mean_LR_product":"LR score reference","Peak mean product":"LR score reference",
 "Pooled mean product":"LR score reference"})
combined["method_family"]=combined.method.replace(family_alias)
combined.loc[combined.method.str.lower().str.startswith(("llr_","linealinkr"),na=False),"method_family"]="LineaLinkR"
combined["implementation_id"]=combined.method
historical_ticci=combined.method.isin(["TICCI","TICCI_topology_only_historical"])
combined.loc[historical_ticci,"active_comparison"]=False
combined.loc[historical_ticci,"scope_reason"]="retained_expression_Shannon_or_graph_stage;author_SLICE_records_define_current_TICCI_coverage"
combined.loc[historical_ticci,"next_action"]="use_author_SLICE_implementation_record"
historical_partition=combined.source_series.eq("existing_server_saved")&combined.task.eq("branch_reconstruction")
combined.loc[historical_partition,"active_comparison"]=False
combined.loc[historical_partition,"next_action"]="use_scope_specific_saved_native_partition_record"
combined["comparison_stage"]=np.select(
 [combined.method.isin(CONFIG["excluded_methods"]),~combined.active_comparison,
  combined.task.eq("projected_fork_event")],
 ["excluded","retained_history","controlled_projection"],default="current_native_or_task_adapter")
current_ticci=combined.loc[combined.method.eq("TICCI_SLICE")&combined.task.eq("cell_ordering"),
 ["dataset","result_status","source_file"]].drop_duplicates("dataset").set_index("dataset")
combined["family_result_status"]=combined.result_status
combined["authoritative_implementation"]=combined.method
ticci_order=combined.method_family.eq("TICCI")&combined.task.eq("cell_ordering")&combined.dataset.isin(current_ticci.index)
combined.loc[ticci_order,"family_result_status"]=combined.loc[ticci_order,"dataset"].map(current_ticci.result_status)
combined.loc[ticci_order,"authoritative_implementation"]="TICCI_SLICE"

combined["status_record"]=combined.dataset+":"+combined.comparison_status.fillna(combined.result_status)
task_table=combined.loc[combined.active_comparison].groupby(
 ["task","method","method_family","implementation_id","comparison_stage","branch_scope","source_series","predicted_object","comparison_role","implementation_label","question_answered"],
 observed=True,dropna=False).agg(
 n_dataset_keys=("dataset","nunique"),dataset_keys=("dataset",lambda s:";".join(sorted(s.unique()))),
 status_by_dataset=("status_record",lambda s:";".join(sorted(s.unique()))),
 classic_metrics=("classic_metrics",lambda s:";".join(sorted(s.dropna().unique()))),
 project_specific_metrics=("project_specific_metrics",lambda s:";".join(sorted(s.dropna().unique()))),
 input_conditions=("input_scope",lambda s:";".join(sorted(s.dropna().unique()))),
 source_files=("evaluation_source",lambda s:";".join(sorted(s.dropna().unique()))),
 next_actions=("next_action",lambda s:";".join(sorted(set(s.dropna())-{""})))).reset_index()
# Include the scored LinkR outputs and fixed-panel adapters omitted by the native-coverage registry.
score_map=pd.concat([
 classic_long.assign(classic_metric_name=classic_long.metric,custom_metric_name=""),
 custom_long.assign(classic_metric_name="",custom_metric_name=custom_long.metric)
],ignore_index=True,sort=False)
score_map=score_map.loc[
 score_map.task.isin(["LR_pair_nomination","projected_fork_event","fixed_depletion_response_enrichment","known_IFN_axis_retrieval","known_IFN_ligand_retrieval","conditional_peak_localization"])
 |(score_map.task.eq("cell_ordering")&score_map.method.eq("LLR_0.6.5_Final065"))
].copy()
score_map["status_record"]=score_map.dataset.astype("string")+":"+score_map.status.astype("string").fillna("saved_metric")
input_parts=score_map.reindex(columns=["input_condition","root_input_mode","direction_policy","denominator"]).fillna("").astype("string")
score_map["input_record"]=input_parts.input_condition+";"+input_parts.root_input_mode+";"+input_parts.direction_policy+";"+input_parts.denominator
scored_task_map=score_map.groupby(["task","method"],observed=True,dropna=False).agg(
 n_dataset_keys=("dataset","nunique"),dataset_keys=("dataset",lambda s:";".join(sorted(s.dropna().unique()))),
 status_by_dataset=("status_record",lambda s:";".join(sorted(s.dropna().unique()))),
 classic_metrics=("classic_metric_name",lambda s:";".join(sorted(set(s.dropna())-{""}))),
 project_specific_metrics=("custom_metric_name",lambda s:";".join(sorted(set(s.dropna())-{""}))),
 input_conditions=("input_record",lambda s:";".join(sorted(s.dropna().unique()))),
 source_files=("source_file",lambda s:";".join(sorted(s.dropna().unique())))).reset_index()
scored_task_map=scored_task_map.merge(task_table[["task","method"]].drop_duplicates().assign(already_listed=True),on=["task","method"],how="left",validate="one_to_one")
scored_task_map=scored_task_map.loc[scored_task_map.already_listed.isna()].drop(columns="already_listed").merge(objects,on="task",how="left",validate="many_to_one")
scored_task_map["source_series"]="existing_server_saved_metric_table"
scored_task_map["implementation_label"]=scored_task_map.method
scored_task_map.loc[scored_task_map.method.eq("LineaLinkR_event_score"),"comparison_role"]="native_fork_event_output_against_fixed_simulation_truth"
scored_task_map.loc[scored_task_map.task.eq("known_IFN_ligand_retrieval")&~scored_task_map.method.eq("ENTRAIN"),"input_conditions"]=scored_task_map.loc[scored_task_map.task.eq("known_IFN_ligand_retrieval")&~scored_task_map.method.eq("ENTRAIN"),"input_conditions"]+";saved_LR_aggregated_ligand_rank_from_table11_not_native_global_ligand_output"
task_table=pd.concat([task_table,scored_task_map],ignore_index=True,sort=False)
historical_failed_task=pd.DataFrame([dict(task="static_ligand_receptor_ranking",method="scTenifoldXct",
 source_series="historical_adapter_error_reported_not_replayed",
 predicted_object="ordered_population_pair_LR_significance",
 comparison_role="historical_unbound_adapter_failure;not_current_native_coverage",
 implementation_label="scTenifoldXct historical adapter error;unbound input;current full native result reported separately",
 question_answered="What unbound historical adapter error was reported before the verified current native run?",
 n_dataset_keys=0,dataset_keys="",status_by_dataset="unbound_historical_input:duplicate_nComp_argument_error_reported;error_log_not_recovered",
 classic_metrics="",project_specific_metrics="",
 input_conditions="historical_input_identity_and_error_log_not_recovered;not_scored_as_method_accuracy",
 source_files=(root/"inputs/comparison/recovered_native/source_adapters/run_sctenifoldxct.py").as_posix(),
 next_actions="historical_unbound_failure_preserved;current_full_paul15_native_run_has_repaired_nComp_and_explicit_ID_binding")])
task_table=pd.concat([task_table,historical_failed_task],ignore_index=True,sort=False)
task_table=task_table.loc[~task_table.method.isin(CONFIG["own_prior_methods"])].copy()
task_table["method_family"]=task_table.method.replace(family_alias)
task_table.loc[task_table.method.str.lower().str.startswith(("llr_","linealinkr"),na=False),"method_family"]="LineaLinkR"
task_table["implementation_id"]=task_table.implementation_id.fillna(task_table.method)
# One task navigator links author code, original runners, inputs, native outputs and paper locations.
source_rows=[
 ("LineaLinkR","experiments/run/run_linealinkr.R","software/LineaLinkR","LineaLinkR"),
 ("LR score reference","experiments/evaluation/gse48968_comparators.py","",""),
 ("CCCvelo","","","cccvelo"),
 ("CellAgentChat","experiments/evaluation/run_cellagentchat.py","tools/competitors/original/CellAgentChat","CellAgentChat"),
 ("CellChat","run_cellchat.R","","CellChat"),
 ("CellPhoneDB","run_cellphonedb.py","","cellphonedb"),
 ("CellRank2","01_run_other_algos.R#RunCellRank2","","cellrank"),
 ("Connectome","run_liana.py","","liana"),
 ("CytoTalk","experiments/evaluation/run_cytotalk_resume.R","tools/competitors/dependencies/CytoTalk","CytoTalk"),
 ("DIISCO","experiments/evaluation/run_diisco.py","tools/competitors/original/DIISCO","diisco"),
 ("DPT","run_dpt.py","","scanpy"),
 ("ENTRAIN","experiments/evaluation/run_niche_entrain.R","tools/competitors/original/ENTRAIN","ENTRAIN"),
 ("FateID","run_fateid.R","","FateID"),
 ("GeneTrajectory","experiments/evaluation/run_genetrajectory.py","tools/competitors/original/GeneTrajectory","GeneTrajectory"),
 ("LIANA","run_liana.py","","liana"),
 ("MultistageOT","run_multistageot.py","","multistageot"),
 ("NATMI","run_natmi.py","","NATMI"),
 ("NicheNet","experiments/evaluation/run_niche_entrain.R","tools/competitors/original/NicheNet","nichenetr"),
 ("PAGA","run_paga.py","","scanpy"),
 ("PHLOWER","experiments/evaluation/run_phlower_full.py","tools/competitors/original/PHLOWER","phlower"),
 ("Palantir","run_palantir.py","","palantir"),
 ("SCORPIUS","run_scorpius.R","","SCORPIUS"),
 ("SLICER","experiments/evaluation/run_slicer_full.R","tools/competitors/original/SLICER","SLICER"),
 ("SingleCellSignalR","run_scsignalr.R","","SingleCellSignalR"),
 ("StaVia","run_stavia.py","","pyVIA"),
 ("TICCI","experiments/evaluation/run_ticci_author_queue.R","inputs/comparison/external_methods/TICCI","TICCI"),
 ("TSCAN","run_tscan_official_batch.R","","TSCAN"),
 ("Tempora","experiments/evaluation/run_tempora_full.R","tools/competitors/original/Tempora","Tempora"),
 ("TimeTalk","run_timetalk.R","","TimeTalk"),
 ("TraSig","run_trasig.R","","TraSig"),
 ("WOT","experiments/evaluation/run_wot.py","tools/competitors/original/wot","wot"),
 ("monocle3","run_monocle3.R","","monocle3"),
 ("pCreode","run_pcreode.py","","pcreode"),
 ("scComm","experiments/evaluation/score_sccomm_native.R","tools/competitors/original/scComm","scComm"),
 ("scFates","run_scfates.py","","scFates"),
 ("scRICH","experiments/evaluation/run_scrich.py","tools/competitors/original/scRICH","scRICH"),
 ("scTenifoldXct","experiments/evaluation/run_sctenifoldxct.py","tools/competitors/original/scTenifoldXct","scTenifoldXct"),
 ("scTour","run_sctour.py","","sctour"),
 ("slingshot","run_slingshot.R","","slingshot"),
 ("tradeSeq","run_tradeseq.R","","tradeSeq")]
source_index=pd.DataFrame(source_rows,columns=["method_family","native_run_entry","author_source_snapshot","package_name"])
bare_runner=source_index.native_run_entry.ne("")&~source_index.native_run_entry.str.contains("/",regex=False)
source_index.loc[bare_runner,"native_run_entry"]="inputs/comparison/recovered_native/source_adapters/"+source_index.loc[bare_runner,"native_run_entry"]
registry_code='root<-Sys.getenv("LLR_WORK_ROOT");Sys.setenv(LLR_WORK_ROOT=root,LLR_ASSETS=root);.libPaths(c(file.path(root,"library"),.libPaths()));invisible(capture.output(source(file.path(root,"experiments/config/00_config.R"))));p<-c(CONFIG$dataset_registry,CONFIG$dataset_narrative_registry);d<-data.frame(kind="input",key=names(p),path=unname(p),version="");ip<-as.data.frame(installed.packages()[,c("Package","LibPath","Version")],stringsAsFactors=FALSE);v<-data.frame(kind="package",key=ip$Package,path=paste0(ip$LibPath,"/",ip$Package),version=ip$Version);s<-paste(capture.output(write.table(rbind(d,v),file="",sep="\t",row.names=FALSE,quote=FALSE)),collapse="\n");cat(paste(sprintf("%02x",as.integer(charToRaw(enc2utf8(s)))),collapse=""))'
registry_env=dict(os.environ,LLR_WORK_ROOT=str(root),LLR_ASSETS=str(root))
registry_run=subprocess.run([CONFIG["paths"]["rscript"],"--encoding=UTF-8","-"],input=registry_code,
 cwd=root,env=registry_env,capture_output=True,text=True,encoding="utf-8",check=True)
registry=pd.read_csv(StringIO(bytes.fromhex(registry_run.stdout).decode("utf-8")),sep="\t").fillna("")
input_registry=registry.loc[registry.kind.eq("input")].drop_duplicates("key",keep="last").set_index("key")
package_registry=registry.loc[registry.kind.eq("package")].drop_duplicates("key").set_index("key")
source_index["available_implementation_version"]=source_index.package_name.map(package_registry.version).fillna("source_snapshot_version_unrecorded")
source_index["available_package_code"]=source_index.package_name.map(package_registry.path).fillna("")
source_index["available_author_source"]=source_index.author_source_snapshot.where(source_index.author_source_snapshot.ne(""),source_index.available_package_code)
source_index["entry_code_file"]=source_index.native_run_entry.str.split("#").str[0]
source_index["entry_code_available"]=source_index.entry_code_file.map(lambda p:(root/p).is_file() if p else False)
source_index["author_source_available"]=source_index.available_author_source.map(lambda p:Path(p).is_dir() if Path(p).is_absolute() else (root/p).is_dir() if p else False)
source_index["adapted_source"]="tools/competitors/optimized/"+source_index.method_family
source_index["adapted_source"]=source_index.adapted_source.where(source_index.adapted_source.map(lambda p:(root/p).is_dir()),"")
task_table=task_table.merge(source_index,on="method_family",how="left",validate="many_to_one")
liana_adapter=task_table.method.str.endswith("_via_LIANA",na=False)
task_table.loc[liana_adapter,"native_run_entry"]="inputs/comparison/recovered_native/source_adapters/run_liana.py"
task_table.loc[task_table.method.eq("TICCI_CellChat"),"native_run_entry"]="experiments/evaluation/run_ticci_author_cellchat.R"
task_table["comparison_stage"]=task_table.comparison_stage.fillna("current_native_or_task_adapter")
task_table.loc[task_table.source_series.str.contains("historical",case=False,na=False),"comparison_stage"]="retained_history"
task_table["active_comparison"]=~task_table.comparison_stage.isin(["retained_history","excluded"])
task_table["method_role"]=np.select([task_table.method_family.eq("LineaLinkR"),task_table.method_family.eq("LR score reference")],["current_method_or_component","simple_reference"],default="external_method")
task_table["method_name"]=task_table.method_family.replace({"slingshot":"Slingshot","monocle3":"Monocle3","LineaLinkR":"LineaLinkR 0.6.5"})
task_table["paper_location"]=task_table.task.map({
 "cell_ordering":"Results 2.3;Figure 3a,c,d;Supplementary Tables S1,S5;Figures S7-S9",
 "branch_reconstruction":"Results 2.3;Figure 3b;Supplementary Table S6;Figure S10d",
 "native_segment_partition":"Results 2.3;Supplementary Table S6;Figure S10d",
 "LR_pair_nomination":"Results 2.2;Figure 2;Table 1;Supplementary Table S2;Figures S2-S3",
 "projected_fork_event":"Results 2.2;Figure 2b;Table 1;Supplementary Table S3;Figures S2-S3",
 "conditional_peak_localization":"Results 2.2;Figure 2c;Supplementary Table S4",
 "fixed_depletion_response_enrichment":"Results 2.4;Figure 4e;Table 2;Supplementary Tables S10-S11;Figure S12c",
 "clone_fate_composition":"Results 2.3;Figure 3e;Supplementary Tables S13,S15-S17;Figures S5,S10",
 "known_IFN_axis_retrieval":"Results 2.5;Figure 5;Supplementary Tables S18,S22",
 "known_IFN_ligand_retrieval":"Results 2.5;Figure 5;Supplementary Tables S18,S22"}).fillna("Supplementary Data; task-specific native object and metric tables")
task_table["evaluation_code"]=task_table.task.map({
 "cell_ordering":"experiments/evaluation/evaluate_ordering.R;experiments/evaluation/evaluate_ticci_author_results.py",
 "branch_reconstruction":"experiments/evaluation/evaluate_native_partitions.py;experiments/evaluation/compare_current_main.py",
 "native_segment_partition":"experiments/evaluation/evaluate_native_partitions.py",
 "LR_pair_nomination":"experiments/evaluation/compare_lrsim_dynamic_methods.R;experiments/evaluation/compare_sccomm_pairs.R",
 "projected_fork_event":"experiments/evaluation/project_sccomm_events.R",
 "conditional_peak_localization":"experiments/evaluation/export_common_temporal_curves.R",
 "fixed_depletion_response_enrichment":"experiments/evaluation/validate_muscle_depletion.R;experiments/evaluation/compare_receiver_nomination.R",
 "clone_fate_composition":"experiments/evaluation/evaluate_native_clones.R;experiments/evaluation/evaluate_clone_composition.R",
 "known_IFN_axis_retrieval":"experiments/evaluation/gse48968_comparators.py",
 "known_IFN_ligand_retrieval":"experiments/evaluation/gse48968_comparators.py;experiments/evaluation/gse48968_entrain.R",
 "static_ligand_receptor_ranking":"experiments/evaluation/compare_native_ccc.py",
 "receiver_ligand_prioritization":"experiments/evaluation/compare_signal_prioritization.py",
 "trajectory_associated_ligand_prioritization":"experiments/evaluation/compare_signal_prioritization.py",
 "within_cell_type_communication_heterogeneity":"experiments/evaluation/compare_scrich_contexts.R",
 "observed_cell_type_proportion_reconstruction":"experiments/evaluation/evaluate_diisco_saved.R",
 "trajectory_gene_association":"experiments/evaluation/evaluate_native_genes_graphs.py",
 "native_graph_topology":"experiments/evaluation/evaluate_native_genes_graphs.py"
}).fillna("experiments/evaluation/evaluate_native_genes_graphs.py")
input_long=task_table[["dataset_keys"]].copy()
input_long["row_id"]=np.arange(len(input_long))
input_long=input_long.assign(dataset=input_long.dataset_keys.fillna("").str.split(";")).explode("dataset")
input_long["input_file"]=input_long.dataset.map(input_registry.path)
input_files=input_long.groupby("row_id",sort=False).input_file.agg(lambda x:";".join(sorted(x.dropna().unique())))
task_table["input_files"]=pd.Series(np.arange(len(task_table))).map(input_files)
manifest=pd.DataFrame(json.loads((root/CONFIG["paths"]["native_manifest"]).read_text(encoding="utf-8-sig")))
manifest["relative_path"]=manifest.relative_path.str.replace("\\","/",regex=False)
manifest["source_norm"]=manifest.source_file.str.replace("\\\\","\\",regex=False).str.replace("\\","/",regex=False)
source_lookup=manifest.drop_duplicates("source_norm",keep="last").set_index("source_norm").server_file
relative_lookup=manifest.drop_duplicates("relative_path",keep="last").set_index("relative_path").server_file
native_links=combined.loc[combined.active_comparison,["dataset","method","task","source_series","original_prediction_source","source_file"]].copy()
native_links["native_path"]=native_links.original_prediction_source.fillna(native_links.source_file).astype("string").str.replace("\\\\","\\",regex=False).str.replace("\\","/",regex=False)
native_links=native_links.assign(native_path=native_links.native_path.str.split(r"[;|]")).explode("native_path")
native_links["relative_path"]=native_links.native_path.str.extract(r"/(?:results|recovered_native)/(by_algorithm/[^;|]+|misc/[^;|]+)$",expand=False)
native_links["resolved_native_path"]=native_links.native_path.map(source_lookup).fillna(native_links.relative_path.map(relative_lookup)).fillna(native_links.native_path)
native_links["resolved_native_path"]=native_links.resolved_native_path.str.replace("D:/biotrade/LineaLinkR_0.6.5",root.as_posix(),regex=False)
native_group=native_links.groupby(["method","task","source_series"],dropna=False).resolved_native_path.agg(
 lambda x:";".join(sorted(x.dropna().unique()))).rename("native_output_files").reset_index()
task_table=task_table.merge(native_group,on=["method","task","source_series"],how="left",validate="many_to_one")
task_table["native_asset_manifest"]=CONFIG["paths"]["native_manifest"]
task_table["input_metadata_manifest"]="inputs/comparison/recovered_native/metadata_source_manifest.json"
task_table["metric_source_files"]=task_table.source_files
task_table.loc[task_table.task.eq("LR_pair_nomination")&task_table.method_family.eq("LR score reference"),"native_run_entry"]="experiments/evaluation/compare_lrsim_dynamic_methods.R"
task_table.loc[task_table.task.str.startswith("known_IFN")&task_table.method_family.eq("TraSig"),"native_run_entry"]="experiments/evaluation/gse48968_comparators.py"
task_table.loc[task_table.task.str.startswith("known_IFN")&task_table.method_family.eq("TimeTalk"),"native_run_entry"]="experiments/evaluation/gse48968_timetalk_applicability.R"
task_table.loc[task_table.task.str.startswith("known_IFN")&task_table.method_family.eq("TimeTalk"),"evaluation_code"]="experiments/evaluation/gse48968_timetalk_applicability.R;writing/update_current_manuscript.py"


# Display identities follow the evaluated object and actual implementation.
task_table.loc[task_table.method_family.eq("TraSig")&task_table.task.eq("LR_pair_nomination"),"method_name"]="TraSig native trajectories"
task_table.loc[task_table.method_family.eq("TraSig")&task_table.task.eq("projected_fork_event"),"method_name"]="TraSig projected contrast"
task_table.loc[task_table.method_family.eq("TraSig")&task_table.task.str.startswith("known_IFN"),"method_name"]="TraSig-style paired dot (one measurement)"
task_table.loc[task_table.method_family.eq("scComm"),"method_name"]="scComm LR (no TF activity)"
task_table.loc[task_table.method_family.eq("scRICH"),"method_name"]="scRICH LR (10 pathways; target=False)"
auxiliary_lung=task_table.task.eq("auxiliary_lung_observed_dot")
task_table.loc[auxiliary_lung,"method_name"]="TraSig observed dot (calendar paths)"
task_table.loc[auxiliary_lung,"paper_location"]="Supplementary Data; supplementary lung native coverage"
task_table.loc[auxiliary_lung,"native_run_entry"]=""
task_table.loc[auxiliary_lung,"evaluation_code"]="experiments/evaluation/summarize_competitor_tasks.py"
task_table.loc[task_table.method_family.eq("tradeSeq")&task_table.task.eq("trajectory_gene_association"),"paper_location"]="Supplementary Results S1; Supplementary Methods S13; Supplementary Table S19"
task_table["entry_code_file"]=task_table.native_run_entry.fillna("").str.split("#").str[0]
task_table["entry_function_anchor"]=task_table.native_run_entry.fillna("").str.partition("#")[2]
task_table["entry_code_available"]=task_table.entry_code_file.map(lambda p:(root/p).is_file() if p else False)
entry_source_text=task_table.entry_code_file.map(lambda p:(root/p).read_text(encoding="utf-8-sig") if p.endswith(".R") else "")
snapshot_bootstrap=entry_source_text.str.contains("../config/llr_bootstrap.R",regex=False)
task_table["native_entry_kind"]=np.select([task_table.entry_function_anchor.ne(""),snapshot_bootstrap,task_table.entry_code_available],
 ["archived_function_snapshot","archived_adapter_snapshot","project_runner_or_task_adapter"],default="saved_output_or_input_requirement_record")
task_table["native_entry_dependencies"]=np.where(snapshot_bootstrap,
 "original_pipeline/config/llr_bootstrap.R;original_pipeline/config/00_config.R;original_pipeline/adapters/_mem_guard.R;original_pipeline/config/panel_skip_rules.R;recorded_author_package_runtime_and_dataset_registry",
 "declared_input_files;available_package_code;available_author_source")
task_table["native_entry_completeness"]=np.select([snapshot_bootstrap,task_table.entry_code_available],
 ["snapshot_requires_original_bootstrap_context","project_entry_and_dependencies_indexed"],default="saved_output_or_applicability_record")
task_table["evaluation_replay_entry"]="experiments/launchers/run_final_evaluation.ps1 -SummariesOnly"
task_table["native_adoption"]=np.where(auxiliary_lung,"auxiliary_context_export",np.where(task_table.active_comparison,"task_specific_current_output","retained_source_record"))

task_table.to_csv(out/"task_matched_comparison.csv",index=False,encoding="utf-8")

combined["n_scored_units"]=pd.to_numeric(combined.n_scored_units.replace({True:1,False:0,"True":1,"False":0}),errors="raise")
combined["native_execution_origin"]=combined.execution_scope
combined["execution_scope"]="read_saved_results_no_refit"
combined["accuracy_scored"]=combined.ground_truth_accuracy.fillna(False)&combined.reindex(columns=["comparison_status"]).iloc[:,0].str.contains("scored|full_declared_scope|conditional_AUC",na=False)
columns=["dataset","dataset_scope","method","task","active_comparison","scope_reason","native_result_status","execution_scope","n_input_cells","n_input_genes","n_native_genes","input_scope","n_native_units","n_scored_units","n_scored_cells","n_targets","n_timepoints","n_time_observed","n_model_genes","unit_statuses","native_seconds","clock_source","n_native_clusters","n_observed_timepoints","n_native_pathways","n_native_edges","n_native_tree_nodes","n_native_tree_edges","n_native_tree_leaves","n_trajectories","n_finite_pseudotime","n_reference_known_time","n_native_gene_pairs","n_returned_genes","n_nonoptimal_OT_pairs","author_species_functional_kappa_available","CellChat_species_database_available","direct_mouse_kappa_ID_matches","direct_human_kappa_ID_matches","direct_zebrafish_LR_endpoint_ID_matches",*flow_columns,"result_status","native_unit_statuses","comparison_status","score_interpretation","accuracy_scored",*repair_columns,"source_file","classic_metrics","context_metrics","analysis_unit","evaluation_source","unavailable_fields","manuscript_use","branch_scope","root_input_mode","direction_policy","cell_coverage","original_prediction_source","reference_revision","source_series","project_specific_metrics","evidence_type","ground_truth_accuracy","denominator_note","n_native_rows","n_reference_units","n_native_score_sets","n_compared_score_sets","planned_native_population_pairs","completed_native_population_pairs","n_distinct_native_gene_names","n_ambiguous_native_gene_names","n_unresolved_native_gene_name_rows","n_invalid_native_model_score_rows","n_common_atomic_keys","n_finite_native_scores","n_saved_genes","n_saved_cells","native_execution_origin","reference_dataset","primary_dataset_key","predicted_object","comparison_role","question_answered","implementation_label","biological_line","planned_units","saved_units","valid_units","failed_units","unverified_units","missing_units","missing_unit_keys","failed_unit_keys","next_action"]
columns+=["implementation_repair","n_empty_suffix_sentinels","runtime_scope"]
combined=combined.loc[~combined.method.isin(CONFIG["own_prior_methods"])].copy()
columns=pd.Index(columns+["method_family","implementation_id","comparison_stage","family_result_status","authoritative_implementation"]).drop_duplicates().tolist()
combined.reindex(columns=columns).sort_values(["task","method","dataset"]).to_csv(out/CONFIG["paths"]["output"],index=False,encoding="utf-8")
assert not combined.duplicated(["dataset","method","task","source_series"]).any()
print(combined.groupby(["method","result_status"],dropna=False).size().to_string(),flush=True)
print("datasets",combined.dataset.nunique(),"rows",len(combined),flush=True)



