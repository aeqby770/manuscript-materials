import os
from pathlib import Path
CONFIG={"paths":{"root":str(Path(__file__).resolve().parents[2]),"native":"inputs/comparison/recovered_native","output":"results/evaluation/final065/saved_native","canonical":"inputs/benchmark/canonical"},"threads":os.cpu_count(),"plot":{},"alpha":0.05,"bifurcation_time":0.35,"dynamic_classes":["early","late","late_shared","branch_specific","bifurcation"]}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"]=str(CONFIG["threads"])
import operator as op
from functools import partial
import numpy as np
import pandas as pd
import anndata as ad
import networkx as nx
from scipy.stats import spearmanr,pearsonr
from sklearn.metrics import roc_auc_score,average_precision_score
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG["plot"])
root=Path(CONFIG["paths"]["root"])
base=root/CONFIG["paths"]["native"]
out=root/CONFIG["paths"]["output"]
out.mkdir(parents=True,exist_ok=True)

paths=sorted((root/CONFIG["paths"]["canonical"]).glob("sim_*_canonical.h5ad"))+[root/CONFIG["paths"]["canonical"]/"simulated_complex_3branch_canonical.h5ad",root/"inputs/benchmark/literature_benchmarks/simulated_branch/simulated_branch.h5ad"]
names=pd.Index(map(op.attrgetter("stem"),paths)).str.removesuffix("_canonical")
objects=list(map(partial(ad.read_h5ad,backed="r"),paths))
genes=pd.concat(map(op.methodcaller("rename",columns={"dataset":"native_annotation_dataset","gene":"native_annotation_gene"}),map(op.attrgetter("var"),objects)),keys=names,names=["dataset","gene"]).reset_index()
cells=pd.concat(map(op.methodcaller("rename",columns={"dataset":"native_annotation_dataset","cell_id":"native_annotation_cell_id"}),map(op.attrgetter("obs"),objects)),keys=names,names=["dataset","cell_id"]).reset_index()
branches=cells[["dataset","true_branch"]].drop_duplicates().astype({"true_branch":str})
branches["supported_class"]="branch"+branches.true_branch
support=branches[["dataset","supported_class"]].assign(class_present=True)
genes["branch_class"]=genes.gene_class.str.match(r"^branch[0-9]+$")
genes=genes.merge(support,left_on=["dataset","gene_class"],right_on=["dataset","supported_class"],how="left",validate="many_to_one")
genes["truth_status"]=np.select([genes.gene_class.isin(CONFIG["dynamic_classes"]),genes.branch_class&genes.class_present.eq(True),genes.gene_class.eq("noise")|genes.branch_class&~genes.class_present.eq(True)],["dynamic","dynamic","static_null"],default="unlabelled")
genes["truth_reason"]=np.select([genes.gene_class.eq("noise"),genes.branch_class&~genes.class_present.eq(True),genes.branch_class&genes.class_present.eq(True)],["generator_constant_noise","generator_branch_absent_in_this_input","generator_branch_present"],default="generator_time_dependent_expression")
genes["truth_dynamic"]=genes.truth_status.eq("dynamic")
genes["truth_source"]="saved_var.gene_class_plus_actual_obs.true_branch; verified_source_generators"
genes[["dataset","gene","gene_class","truth_status","truth_dynamic","truth_reason","truth_source"]].to_csv(out/"simulation_gene_truth_readback.csv",index=False,encoding="utf-8")
list(map(op.methodcaller("close"),map(op.attrgetter("file"),objects)))

paths=sorted((base/"by_algorithm/tradeSeq").glob("*/tradeseq_association_test.csv"))
names=pd.Index(map(lambda p:p.parent.name,paths))
tests=pd.concat(map(pd.read_csv,paths),keys=names,names=["dataset","source_row"]).reset_index()
tests["source_file"]=tests.dataset.map(pd.Series(list(map(str,paths)),index=names))
assert not tests.duplicated(["dataset","gene"]).any()
tests["finite_p"]=np.isfinite(tests.pvalue)
assert tests.loc[tests.finite_p,"pvalue"].between(0,1).all()
keys=["dataset","source_file"]
valid=tests.loc[tests.finite_p].sort_values(keys+["pvalue","gene"]).copy()
valid["m"]=valid.groupby(keys,observed=True).pvalue.transform("size")
valid["p_rank"]=valid.groupby(keys,observed=True).cumcount()+1
valid["bh_raw"]=valid.pvalue*valid.m/valid.p_rank
reverse=valid.sort_values(keys+["p_rank"],ascending=[True,True,False]).copy()
valid["padj_BH"]=reverse.groupby(keys,observed=True).bh_raw.cummin().clip(upper=1)
tests=tests.merge(valid[keys+["gene","padj_BH"]],on=keys+["gene"],how="left",validate="one_to_one")
tests=tests.merge(genes[["dataset","gene","gene_class","truth_status","truth_dynamic","truth_reason"]],on=["dataset","gene"],how="left",validate="one_to_one")
tests["called_BH"]=tests.padj_BH.le(CONFIG["alpha"])
tests["called_raw_p"]=tests.pvalue.lt(CONFIG["alpha"])
tests["true_call"]=tests.called_BH&tests.truth_status.eq("dynamic")
tests["false_call"]=tests.called_BH&tests.truth_status.eq("static_null")
tests["scored_positive"]=tests.finite_p&tests.truth_status.eq("dynamic")
tests["scored_null"]=tests.finite_p&tests.truth_status.eq("static_null")
summary=tests.groupby(keys,observed=True).agg(n_native_test_rows=("gene","size"),n_finite_p=("finite_p","sum"),n_raw_p05=("called_raw_p","sum"),n_BH05=("called_BH","sum"),n_true_BH_calls=("true_call","sum"),n_false_BH_calls=("false_call","sum"),n_positive_finite_scores=("scored_positive","sum"),n_null_finite_scores=("scored_null","sum")).reset_index()
truth_sizes=genes.groupby("dataset",observed=True).agg(n_truth_genes=("gene","size"),n_dynamic_truth=("truth_dynamic","sum")).reset_index()
summary=summary.merge(truth_sizes,on="dataset",how="left",validate="many_to_one")
summary["n_unreported_genes"]=summary.n_truth_genes-summary.n_native_test_rows
summary["positive_score_coverage"]=summary.n_positive_finite_scores/summary.n_dynamic_truth
summary["BH_precision"]=summary.n_true_BH_calls/summary.n_BH05
summary["BH_recall_fixed_truth"]=summary.n_true_BH_calls/summary.n_dynamic_truth
summary["BH_F1_fixed_truth"]=2*summary.n_true_BH_calls/(summary.n_BH05+summary.n_dynamic_truth)
summary["empirical_FDP_BH"]=np.where(summary.n_BH05.gt(0),summary.n_false_BH_calls/summary.n_BH05,0)
eligible=summary.loc[summary.n_positive_finite_scores.gt(0)&summary.n_null_finite_scores.gt(0),keys]
ranked=tests.loc[tests.finite_p&tests.truth_status.isin(["dynamic","static_null"])].merge(eligible,on=keys,validate="many_to_one")
ranked["truth_dynamic"]=ranked.truth_status.eq("dynamic")
groups=list(ranked.groupby(keys,sort=True,observed=True))
frames=list(map(op.itemgetter(1),groups))
areas=pd.DataFrame(list(map(op.itemgetter(0),groups)),columns=keys)
areas["AUROC"]=list(map(roc_auc_score,map(op.itemgetter("truth_dynamic"),frames),map(op.neg,map(op.itemgetter("pvalue"),frames))))
areas["average_precision_step"]=list(map(average_precision_score,map(op.itemgetter("truth_dynamic"),frames),map(op.neg,map(op.itemgetter("pvalue"),frames))))
summary=summary.merge(areas,on=keys,how="left",validate="one_to_one")
summary["method"]="tradeSeq"
summary["metric_family"]="classical_gene_association"
summary["evaluation_unit"]="gene_in_saved_associationTest"
summary["truth_accuracy_status"]=np.select([summary.n_truth_genes.isna(),summary.n_null_finite_scores.eq(0),summary.n_positive_finite_scores.eq(0)],["no_explicit_dynamic_gene_truth","no_finite_null_scores_for_conditional_AUC","no_finite_positive_scores_for_conditional_AUC"],default="conditional_AUC_and_fixed_truth_recall")
summary["AUC_denominator"]="finite_saved_pvalues_with_explicit_generator_truth"
summary["recall_denominator"]="all_generator_dynamic_genes_in_original_input"
summary["p_adjustment"]="BH_over_finite_saved_pvalues; original_fit_unchanged"
summary.loc[summary.n_truth_genes.isna(),["n_true_BH_calls","n_false_BH_calls","n_positive_finite_scores","n_null_finite_scores","BH_precision","BH_recall_fixed_truth","BH_F1_fixed_truth","empirical_FDP_BH"]]=np.nan
binding=pd.read_csv(out/"native_gene_input_binding.csv")
summary=summary.merge(binding.drop(columns=["source_file","method"]),on="dataset",how="left",validate="one_to_one")
summary.to_csv(out/"classic_native_gene_metrics.csv",index=False,encoding="utf-8")
tests.to_csv(out/"native_gene_test_readback.csv.gz",index=False,encoding="utf-8")

paga_paths=sorted((base/"misc").glob("paga_*_connectivities.csv"))
paga_names=pd.Index(map(op.attrgetter("stem"),paga_paths)).str.removeprefix("paga_").str.removesuffix("_connectivities")
adj=list(map(partial(pd.read_csv,index_col=0),paga_paths))
adj=list(map(lambda a:a.set_axis(a.index.astype(str),axis=0).set_axis(a.columns.astype(str),axis=1),adj))
paga_graphs=list(map(nx.from_pandas_adjacency,adj))
sf_paths=sorted((base/"misc").glob("scfates_lineage_edges_*.csv"))
sf_names=pd.Index(map(op.attrgetter("stem"),sf_paths)).str.removeprefix("scfates_lineage_edges_")
edges=list(map(pd.read_csv,sf_paths))
sf_graphs=list(map(partial(nx.from_pandas_edgelist,source="from_node",target="to_node"),edges))
ticci_paths=sorted(list((root/"inputs/comparison/results/by_algorithm/TICCI").glob("*/ticci_lineage_edges.csv"))+list((base/"by_algorithm/TICCI/full_from_revision_20260919").glob("*/ticci_lineage_edges.csv")))
ticci_names=pd.Index(map(lambda p:p.parent.name,ticci_paths))
ticci_frames=list(map(pd.read_csv,ticci_paths))
ticci_source=pd.DataFrame({"dataset":ticci_names,"columns":list(map(op.attrgetter("columns"),ticci_frames))})
ticci_source["columns"]=ticci_source["columns"].astype(str)
ticci_source.to_csv(out/"native_TICCI_graph_schema.csv",index=False,encoding="utf-8")
graph_manifest=pd.DataFrame({"dataset":np.r_[paga_names,sf_names],"method":np.r_[np.repeat("PAGA",len(paga_paths)),np.repeat("scFates_principal_graph",len(sf_paths))],"source_file":list(map(str,paga_paths+sf_paths))})
graphs=paga_graphs+sf_graphs
graph_manifest["n_observed_nodes"]=list(map(nx.number_of_nodes,graphs))
graph_manifest["n_saved_edges"]=list(map(nx.number_of_edges,graphs))
graph_manifest["n_components"]=list(map(nx.number_connected_components,graphs))
graph_manifest["n_self_loops"]=list(map(nx.number_of_selfloops,graphs))
degrees=list(map(lambda g:pd.Series(dict(g.degree())),graphs))
graph_manifest["n_degree_one_nodes"]=list(map(lambda x:x.eq(1).sum(),degrees))
graph_manifest["n_degree_at_least_three"]=list(map(lambda x:x.ge(3).sum(),degrees))
graph_manifest["metric_family"]="descriptive_native_graph_structure"
graph_manifest["node_scope"]=np.r_[np.repeat("all_rows_of_saved_adjacency",len(paga_paths)),np.repeat("nodes_present_in_saved_edge_list",len(sf_paths))]
graph_manifest["input_annotation_condition"]=np.where(graph_manifest.method.eq("PAGA"),"native_PAGA_uses_provided_metadata_groups","saved_principal_points")
graph_manifest["classic_graph_accuracy_status"]="no_saved_explicit_GT_graph_and_node_correspondence"
graph_manifest.to_csv(out/"native_graph_structure_readback.csv",index=False,encoding="utf-8")

paths=sorted((base/"misc").glob("scfates_lineage_segments_*.csv"))
names=pd.Index(map(op.attrgetter("stem"),paths)).str.removeprefix("scfates_lineage_segments_")
segments=pd.concat(map(pd.read_csv,paths),keys=names,names=["dataset","source_row"]).reset_index()
paths=sorted((base/"misc").glob("scfates_lineage_cell_segments_*.csv"))
names=pd.Index(map(op.attrgetter("stem"),paths)).str.removeprefix("scfates_lineage_cell_segments_")
cell_seg=pd.concat(map(pd.read_csv,paths),keys=names,names=["dataset","source_row"]).reset_index()
labels=cell_seg[["dataset","segment"]].drop_duplicates()
named=segments[["dataset","segment_id"]].drop_duplicates().assign(named_column_match=True)
third=segments[["dataset","length"]].drop_duplicates().assign(third_column_match=True)
labels=labels.merge(named,left_on=["dataset","segment"],right_on=["dataset","segment_id"],how="left",validate="many_to_one").merge(third,left_on=["dataset","segment"],right_on=["dataset","length"],how="left",validate="many_to_one")
labels["named_column_match"]=labels.named_column_match.eq(True)
labels["third_column_match"]=labels.third_column_match.eq(True)
segment_check=labels.groupby("dataset",observed=True).agg(n_cell_segment_labels=("segment","size"),n_match_named_segment_id=("named_column_match","sum"),n_match_third_export_column=("third_column_match","sum")).reset_index()
segment_check["evaluation_reader"]="cell_segment_labels_read_directly; no_role_inference_from_misnamed_segment_table"
segment_check.to_csv(out/"native_scfates_segment_schema_check.csv",index=False,encoding="utf-8")

paths=sorted((base/"smoke/by_algorithm/SLICER").glob("*/slicer_pseudotime.csv"))
names=pd.Index(map(lambda p:p.parent.name,paths))
slicer=pd.concat(map(partial(pd.read_csv,dtype={"cell_id":"string"}),paths),keys=names,names=["dataset","source_row"]).reset_index()
ref=cells[["dataset","cell_id","true_pseudotime"]].astype({"cell_id":"string"})
sl=slicer.merge(ref,on=["dataset","cell_id"],how="inner",validate="one_to_one")
sl=sl.loc[np.isfinite(sl.pseudotime)&np.isfinite(sl.true_pseudotime)].copy()
sl["p01"]=(sl.pseudotime-sl.groupby("dataset").pseudotime.transform("min"))/(sl.groupby("dataset").pseudotime.transform("max")-sl.groupby("dataset").pseudotime.transform("min"))
sl["t01"]=(sl.true_pseudotime-sl.groupby("dataset").true_pseudotime.transform("min"))/(sl.groupby("dataset").true_pseudotime.transform("max")-sl.groupby("dataset").true_pseudotime.transform("min"))
sl["sq_error"]=(sl.p01-sl.t01)**2
groups=list(sl.groupby("dataset",sort=True))
frames=list(map(op.itemgetter(1),groups))
sl_metrics=pd.DataFrame({"dataset":list(map(op.itemgetter(0),groups))})
sl_metrics["Spearman_signed"]=list(map(op.itemgetter(0),map(spearmanr,map(op.itemgetter("p01"),frames),map(op.itemgetter("t01"),frames))))
sl_metrics["Pearson"]=list(map(op.itemgetter(0),map(pearsonr,map(op.itemgetter("p01"),frames),map(op.itemgetter("t01"),frames))))
sl_sizes=sl.groupby("dataset").agg(n_used=("cell_id","size"),MSE=("sq_error","mean")).reset_index()
sl_metrics=sl_metrics.merge(sl_sizes,on="dataset",validate="one_to_one")
sl_metrics["RMSE"]=np.sqrt(sl_metrics.MSE)
sl_metrics["n_reference"]=sl_metrics.dataset.map(ref.groupby("dataset").size())
sl_metrics["cell_coverage"]=sl_metrics.n_used/sl_metrics.n_reference
sl_metrics["method"]="SLICER"
sl_metrics["metric_family"]="classical_ordering"
sl_metrics["execution_scope"]="saved_smoke_artifact_on_original_500_cell_simulation"
sl_metrics.to_csv(out/"classic_native_SLICER_ordering.csv",index=False,encoding="utf-8")
print(summary[["dataset","n_native_test_rows","n_finite_p","n_dynamic_truth","n_BH05","BH_recall_fixed_truth","AUROC","truth_accuracy_status"]].to_string(index=False))
print(graph_manifest.groupby("method").size().to_string())

# Saved topology objects remain descriptive without a matched GT graph.
from scipy import sparse
cy_paths=sorted((base/"by_algorithm/CytoTalk").rglob("FinalNetwork.txt"))
cy_names=pd.Index(map(lambda p:p.parent.parent.name,cy_paths))
cy_frames=list(map(partial(pd.read_csv,sep="\t"),cy_paths))
cy_frames=list(map(lambda d:d.assign(node_a=d.node1_type.astype(str)+"|"+d.node1.astype(str),node_b=d.node2_type.astype(str)+"|"+d.node2.astype(str)),cy_frames))
cy_graphs=list(map(partial(nx.from_pandas_edgelist,source="node_a",target="node_b"),cy_frames))
ticci_graphs=list(map(partial(nx.from_pandas_edgelist,source="from_cluster",target="to_cluster"),ticci_frames))
all_graphs=cy_graphs+ticci_graphs
extra_graph=pd.DataFrame({"dataset":np.r_[cy_names,ticci_names],"method":np.r_[np.repeat("CytoTalk_native_PCST",len(cy_paths)),np.repeat("TICCI_native_graph",len(ticci_paths))],"source_file":list(map(str,cy_paths+ticci_paths))})
extra_graph["n_observed_nodes"]=list(map(nx.number_of_nodes,all_graphs))
extra_graph["n_saved_edges"]=list(map(nx.number_of_edges,all_graphs))
extra_graph["n_components"]=list(map(nx.number_connected_components,all_graphs))
extra_graph["n_self_loops"]=list(map(nx.number_of_selfloops,all_graphs))
extra_degrees=list(map(lambda g:pd.Series(dict(g.degree()),dtype=float),all_graphs))
extra_graph["n_degree_one_nodes"]=list(map(lambda d:d.eq(1).sum(),extra_degrees))
extra_graph["n_degree_at_least_three"]=list(map(lambda d:d.ge(3).sum(),extra_degrees))
extra_graph["n_native_cross_population_LR_edges"]=np.r_[list(map(lambda d:d.is_ct_edge.eq(True).sum(),cy_frames)),np.repeat(np.nan,len(ticci_paths))]
extra_graph["metric_family"]="descriptive_native_graph_structure"
extra_graph["node_scope"]="nodes_present_in_saved_edges; topology_summary_uses_undirected_projection"
extra_graph["input_annotation_condition"]=np.where(extra_graph.method.eq("CytoTalk_native_PCST"),"native_population_pair_and_intracellular_genes","native_cluster_graph")
extra_graph["classic_graph_accuracy_status"]="no_saved_explicit_GT_graph_and_node_correspondence"
graph_manifest=pd.concat([graph_manifest,extra_graph],ignore_index=True)
history_paths=sorted((base/"by_algorithm/TICCI/topology_only_20260731").glob("*/ticci_lineage_edges.csv"))
history_edges=list(map(pd.read_csv,history_paths))
history_graphs=list(map(partial(nx.from_pandas_edgelist,source="from_cluster",target="to_cluster"),history_edges))
history_summary=pd.concat(map(pd.read_csv,map(lambda p:p.with_name("ticci_summary.csv"),history_paths)),keys=list(map(lambda p:p.parent.name,history_paths)),names=["dataset","row"]).reset_index()
history=pd.DataFrame({"dataset":list(map(lambda p:p.parent.name,history_paths)),"method":"TICCI_topology_only_historical","source_file":list(map(str,history_paths))})
history["n_observed_nodes"]=list(map(nx.number_of_nodes,history_graphs))
history["n_saved_edges"]=list(map(nx.number_of_edges,history_graphs))
history["n_components"]=list(map(nx.number_connected_components,history_graphs))
history["n_self_loops"]=list(map(nx.number_of_selfloops,history_graphs))
history["n_declared_input_cells"]=history.dataset.map(history_summary.set_index("dataset").n_cells)
history["n_declared_input_clusters"]=history.dataset.map(history_summary.set_index("dataset").n_clusters)
history["metric_family"]="descriptive_saved_graph_structure"
history["node_scope"]="only_saved_tree_nodes;no_per_cell_cluster_membership_saved"
history["input_annotation_condition"]="historical_ticci_lineage_topology_only;not_current_CellChat_CCI_workflow"
history["classic_graph_accuracy_status"]="no_matched_GT_node_correspondence;no_saved_pseudotime"
graph_manifest=pd.concat([graph_manifest,history],ignore_index=True)
graph_manifest.to_csv(out/"native_graph_structure_readback.csv",index=False,encoding="utf-8")
wot_paths=sorted((base/"smoke/by_algorithm/wot").rglob("*.h5ad"))
wot_objects=list(map(ad.read_h5ad,wot_paths))
wot_values=list(map(lambda x:sparse.csr_matrix(x.X).toarray(),wot_objects))
wot=pd.DataFrame({"dataset":list(map(lambda p:p.parent.name,wot_paths)),"method":"WOT","source_file":list(map(str,wot_paths)),"artifact":"saved_transport_map"})
wot["n_source_cells"]=list(map(op.attrgetter("n_obs"),wot_objects))
wot["n_target_cells"]=list(map(op.attrgetter("n_vars"),wot_objects))
wot["n_transport_entries"]=wot.n_source_cells*wot.n_target_cells
wot["n_finite_entries"]=list(map(lambda x:np.isfinite(x).sum(),wot_values))
wot["n_negative_entries"]=list(map(lambda x:(x<0).sum(),wot_values))
wot["native_score_min"]=list(map(np.nanmin,wot_values))
wot["native_score_max"]=list(map(np.nanmax,wot_values))
wot["evaluation_status"]=np.where(wot.n_finite_entries.eq(wot.n_transport_entries),"saved_finite_transport_matrix;no_explicit_matched_transport_GT","nonfinite_saved_transport_entries;no_valid_transport_accuracy")
wot_cells=pd.DataFrame({"dataset":wot.dataset,"cell_id":list(map(lambda x:list(pd.Index(x.obs_names).union(x.var_names)),wot_objects))}).explode("cell_id")
wot["n_input_cohort_cells"]=wot.dataset.map(wot_cells.groupby("dataset").cell_id.nunique())
wot_fates_paths=sorted((base/"smoke/by_algorithm/wot").rglob("wot_fates_summary.csv"))
wot_fates_frames=list(map(partial(pd.read_csv,dtype={"cell_id":"string"}),wot_fates_paths))
wot_fates=pd.DataFrame({"dataset":list(map(lambda p:p.parent.name,wot_fates_paths)),"method":"WOT","source_file":list(map(str,wot_fates_paths)),"artifact":"max_fate_summary"})
wot_fates["n_saved_cells"]=list(map(len,wot_fates_frames))
wot_fates["n_finite_native_scores"]=list(map(lambda d:np.isfinite(d.max_fate).sum(),wot_fates_frames))
wot_fates["evaluation_status"]=np.where(wot_fates.n_finite_native_scores.eq(wot_fates.n_saved_cells),"all_finite_max_fate_only;no_target_identity_or_full_probability_vector","partial_nonfinite_max_fate_summary;no_target_identity_or_full_probability_vector")
wot_fates["n_input_cohort_cells"]=wot_fates.dataset.map(wot.groupby("dataset").n_input_cohort_cells.first())
gt_paths=sorted((base/"smoke/misc").glob("genetrajectory_*_gene_pt.csv"))
gt_frames=list(map(pd.read_csv,gt_paths))
gt=pd.DataFrame({"dataset":pd.Index(map(op.attrgetter("stem"),gt_paths)).str.removeprefix("genetrajectory_").str.removesuffix("_gene_pt"),"method":"GeneTrajectory","source_file":list(map(str,gt_paths)),"artifact":"gene_pseudoorder"})
gt["n_saved_genes"]=list(map(len,gt_frames))
gt["n_finite_native_scores"]=list(map(lambda d:np.isfinite(d["Pseudoorder-1"]).sum(),gt_frames))
gt["evaluation_status"]="native_gene_order_on_declared_16_gene_smoke_subset;no_independent_gene_onset_rank_GT"
slicer_all=pd.DataFrame({"dataset":slicer.dataset.unique(),"method":"SLICER","artifact":"cell_pseudotime"})
slicer_all["source_file"]=slicer_all.dataset.map(pd.Series(list(map(str,paths)),index=names))
slicer_all["n_saved_cells"]=slicer_all.dataset.map(slicer.groupby("dataset").size())
slicer_all["evaluation_status"]=np.where(slicer_all.dataset.isin(sl_metrics.dataset),"scored_original_simulation_all_500_cells","saved_smoke_cell_pseudotime;no_independent_ordering_reference")
pd.concat([wot,wot_fates,gt,slicer_all],ignore_index=True).to_csv(out/"native_remaining_artifact_coverage.csv",index=False,encoding="utf-8")
print("native_graph_objects",len(graph_manifest),"WOT_maps",len(wot),"WOT_fate_rows",int(wot_fates.n_saved_cells.sum()),"GeneTrajectory_inputs",len(gt))

flow_paths=sorted((root/"results/by_algorithm/FlowSig").glob("*/flow_edges.csv.gz"))
flow_frames=list(map(pd.read_csv,flow_paths))
flow_summary=pd.DataFrame({"dataset":list(map(lambda p:p.parent.name,flow_paths)),"method":"FlowSig_NMF","source_file":list(map(str,flow_paths))})
flow_variables=list(map(lambda p:pd.read_csv(p.with_name("retained_flow_variables.csv")),flow_paths))
flow_summary["n_native_variables"]=list(map(len,flow_variables))
flow_summary["n_distinct_saved_node_labels"]=list(map(lambda d:pd.Index(d.source).union(d.target).nunique(),flow_frames))
flow_summary["n_distinct_typed_node_labels"]=list(map(lambda d:d[["Unnamed: 0","Type"]].drop_duplicates().shape[0],flow_variables))
flow_summary["node_identity_basis"]="row_position_and_Type_in_saved_retained_flow_variables;bare_name_not_unique"
flow_summary["n_saved_possible_edge_entries"]=list(map(len,flow_frames))
flow_summary["n_positive_bootstrap_edges"]=list(map(lambda d:d.bootstrap_frequency.gt(0).sum(),flow_frames))
flow_summary["n_biological_flow_edges"]=list(map(lambda d:d.biological_flow_weight.ne(0).sum(),flow_frames))
flow_summary["n_retained_edges"]=list(map(lambda d:d.retained_weight.ne(0).sum(),flow_frames))
flow_summary["n_nonfinite_retained_entries"]=list(map(lambda d:(~np.isfinite(d.retained_weight)).sum(),flow_frames))
assert np.array_equal(flow_summary.n_native_variables.to_numpy()**2,flow_summary.n_saved_possible_edge_entries.to_numpy())
flow_summary["comparison_status"]="saved_native_inflow_module_outflow_graph;no_explicit_GT_graph"
flow_summary.to_csv(out/"native_FlowSig_graph_readback.csv",index=False,encoding="utf-8")


