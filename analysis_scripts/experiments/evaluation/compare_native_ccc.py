import os
import pickle
from pathlib import Path
CONFIG={"paths":{"root":str(Path(__file__).resolve().parents[2]),"native":"inputs/comparison/recovered_native","output":"results/evaluation/final065/saved_native"},"threads":os.cpu_count(),"plot":{},"population_columns":["cell_type","clusters","cluster","celltype","seurat_clusters","auto_clusters","state_info","cell_ontology_class","ct","ct_subtype","paul15_clusters","louvain_cluster","meta_celltype"],"LIANA_fields":[["LIANA_rank_aggregate","magnitude_rank",-1],["LIANA_rank_aggregate","specificity_rank",-1],["CellPhoneDB_via_LIANA","lr_means",1],["CellPhoneDB_via_LIANA","cellphone_pvals",-1],["NATMI_via_LIANA","expr_prod",1],["NATMI_via_LIANA","spec_weight",1],["SingleCellSignalR_via_LIANA","lrscore",1],["Connectome_via_LIANA","expr_prod",1],["Connectome_via_LIANA","scaled_weight",1]],"NATMI_fields":["Edge average expression weight","Edge average expression derived specificity","Edge total expression weight","Edge total expression derived specificity"]}
os.environ["OMP_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["OPENBLAS_NUM_THREADS"]=str(CONFIG["threads"])
os.environ["MKL_NUM_THREADS"]=str(CONFIG["threads"])
import operator as op
from functools import partial
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update(CONFIG["plot"])
root=Path(CONFIG["paths"]["root"])
base=root/CONFIG["paths"]["native"]
out=root/CONFIG["paths"]["output"]
ref=pd.read_csv(out/"static_reference_scores.csv.gz")
refcells=pd.read_csv(out/"static_reference_cells.csv.gz",dtype={"cell_id":"string","reference_label":"string"})
assert not refcells.duplicated(["dataset","cell_id"]).any()
refsign=refcells.sort_values(["dataset","cell_id"]).groupby("dataset",observed=True).cell_id.agg(tuple).rename("cell_signature").reset_index().rename(columns={"dataset":"reference_dataset"})
reference_sizes=refcells.groupby("dataset",observed=True).size().rename("n_reference_cells")
object_cells=pd.read_csv(out/"native_object_cells.csv.gz",dtype={"cell_id":"string","native_label":"string"})
ticci_root=root/"inputs/comparison/results/by_algorithm/TICCI"
ticci_cell_paths=sorted(list(ticci_root.glob("*/ticci_pseudotime.csv"))+list((base/"by_algorithm/TICCI/full_from_revision_20260919").glob("*/ticci_pseudotime.csv")))
prior_ticci_datasets=pd.Index(map(lambda p:p.parent.name,ticci_cell_paths))
author_ticci_cells=list((root/"results/by_algorithm/TICCI_author_functional_entropy").glob("*/ticci_pseudotime.csv"))
ticci_cell_paths+=list(filter(lambda p:p.parent.name not in prior_ticci_datasets,author_ticci_cells))
ticci_datasets=pd.Index(map(lambda p:p.parent.name,ticci_cell_paths))
ticci_cells=pd.concat(map(partial(pd.read_csv,dtype={"cell_id":"string"}),ticci_cell_paths),keys=ticci_datasets,names=["dataset","row"]).reset_index()
ticci_cells=ticci_cells.rename(columns={"classification":"native_label"})
ticci_cells["method"]="TICCI_CellChat"
ticci_cells["source_file"]=ticci_cells.dataset.map(pd.Series(list(map(lambda p:str(p.with_name("ticci_cellchat_interactions.csv")).replace("\\","/"),ticci_cell_paths)),index=ticci_datasets))
object_cells=pd.concat([object_cells,ticci_cells[object_cells.columns]],ignore_index=True)
object_ids=object_cells.sort_values(["dataset","method","source_file","cell_id"]).groupby(["dataset","method","source_file"],observed=True).cell_id.agg(tuple).rename("cell_signature").reset_index()
object_cohort=object_ids.merge(refsign,on="cell_signature",how="left",validate="many_to_one").drop(columns="cell_signature")
object_cohort["identity_basis"]=np.where(object_cohort.method.eq("TICCI_CellChat"),"exact_companion_TICCI_cell_ids_and_classification;summary_counts_not_used","exact_cell_set_in_saved_native_object")
object_join=object_cells.merge(object_cohort,on=["dataset","method","source_file"],validate="many_to_one")
object_join=object_join.merge(refcells[["dataset","cell_id","reference_label"]].rename(columns={"dataset":"reference_dataset"}),on=["reference_dataset","cell_id"],how="left",validate="many_to_one")
object_join["reference_missing"]=object_join.reference_label.isna()
maps=object_join.groupby(["dataset","method","source_file","native_label"],observed=True).agg(reference_dataset=("reference_dataset","first"),reference_label=("reference_label","first"),n_native_cells=("cell_id","size"),n_reference_labels=("reference_label","nunique"),n_missing_reference=("reference_missing","sum")).reset_index()
pop_sizes=refcells.groupby(["dataset","reference_label"],observed=True).size().rename("n_reference_population").reset_index().rename(columns={"dataset":"reference_dataset"})
maps=maps.merge(pop_sizes,on=["reference_dataset","reference_label"],how="left",validate="many_to_one")
maps["complete_reference_population"]=maps.n_reference_labels.eq(1)&maps.n_missing_reference.eq(0)&maps.n_native_cells.eq(maps.n_reference_population)
reverse=maps.loc[maps.complete_reference_population].groupby(["dataset","method","source_file","reference_label"],observed=True).native_label.nunique().rename("n_native_labels_for_reference").reset_index()
maps=maps.merge(reverse,on=["dataset","method","source_file","reference_label"],how="left",validate="many_to_one")
maps["mapping_accepted"]=maps.complete_reference_population&maps.n_native_labels_for_reference.eq(1)
object_map=maps.loc[maps.mapping_accepted,["dataset","method","source_file","native_label","reference_dataset","reference_label"]].copy()
object_map["identity_basis"]="exact_native_cells_and_complete_one_to_one_population"
maps.to_csv(out/"native_CCC_object_population_mapping.csv",index=False,encoding="utf-8")

meta_paths=sorted((base/"input_metadata").glob("*/meta.csv"))
meta_names=pd.Index(map(lambda p:p.parent.name,meta_paths))
meta=pd.concat(map(op.methodcaller("rename",columns={"dataset":"metadata_dataset_annotation","source_row":"metadata_row_annotation"}),map(partial(pd.read_csv,dtype=str),meta_paths)),keys=meta_names,names=["dataset","source_row"]).reset_index()
old_paths=sorted((base/"paper_tables/reference").glob("ref_*.csv"))
old_names=pd.Index(map(op.attrgetter("stem"),old_paths)).str.removeprefix("ref_")
old_meta=pd.concat(map(op.methodcaller("rename",columns={"dataset":"metadata_dataset_annotation","source_row":"metadata_row_annotation"}),map(partial(pd.read_csv,dtype=str),old_paths)),keys=old_names,names=["dataset","source_row"]).reset_index()
meta=pd.concat([meta,old_meta.loc[~old_meta.dataset.isin(meta.dataset)]],ignore_index=True)
meta["cell_id"]=meta.cell_id.astype("string")
assert not meta.duplicated(["dataset","cell_id"]).any()
meta_ids=meta.sort_values(["dataset","cell_id"]).groupby("dataset",observed=True).cell_id.agg(tuple).rename("cell_signature").reset_index()
snapshot_cohort=meta_ids.merge(refsign,on="cell_signature",how="left",validate="many_to_one").drop(columns="cell_signature")
snapshot_cohort["identity_basis"]="exact_declared_snapshot_cell_set; cell_ids_not_in_score_file"
meta=meta.merge(snapshot_cohort,on="dataset",validate="many_to_one").merge(refcells[["dataset","cell_id","reference_label"]].rename(columns={"dataset":"reference_dataset"}),on=["reference_dataset","cell_id"],how="left",validate="many_to_one")
meta_long=meta.melt(id_vars=["dataset","cell_id","reference_dataset","reference_label"],value_vars=meta.columns.intersection(CONFIG["population_columns"]),var_name="metadata_column",value_name="native_label").dropna(subset=["native_label"])
meta_long["reference_missing"]=meta_long.reference_label.isna()
candidate=meta_long.groupby(["dataset","metadata_column","native_label"],observed=True).agg(reference_dataset=("reference_dataset","first"),reference_label=("reference_label","first"),n_native_cells=("cell_id","size"),n_reference_labels=("reference_label","nunique"),n_missing_reference=("reference_missing","sum")).reset_index()
candidate=candidate.merge(pop_sizes,on=["reference_dataset","reference_label"],how="left",validate="many_to_one")
candidate["complete_population"]=candidate.n_reference_labels.eq(1)&candidate.n_missing_reference.eq(0)&candidate.n_native_cells.eq(candidate.n_reference_population)
usable=candidate.loc[candidate.complete_population,["dataset","native_label","reference_dataset","reference_label"]].drop_duplicates()
ambiguity=usable.groupby(["dataset","native_label"],observed=True).reference_label.nunique().rename("n_possible_reference_labels").reset_index()
usable=usable.merge(ambiguity,on=["dataset","native_label"],validate="many_to_one")
snapshot_map=usable.loc[usable.n_possible_reference_labels.eq(1),["dataset","native_label","reference_dataset","reference_label"]].drop_duplicates()
snapshot_map["identity_basis"]="declared_input_snapshot_and_complete_population_correspondence"
candidate.to_csv(out/"native_CCC_snapshot_population_candidates.csv",index=False,encoding="utf-8")

r_scores=pd.read_csv(out/"native_object_CCC_scores.csv.gz",dtype={"source_a":"string","source_b":"string","entity_a":"string","entity_b":"string"})
r_scores["native_run"]=r_scores.source_file
r_scores.loc[r_scores.method.eq("CytoTalk_native"),"native_run"]=r_scores.loc[r_scores.method.eq("CytoTalk_native"),"source_file"].str.replace(r"/[^/]+/cytotalk_result[.]rds$","",regex=True)
r_scores["returned_scope"]=np.where(r_scores.method.eq("CellChat"),"all_saved_net_probability_entries","selected_native_network_cross_edges")
scsr=pd.read_csv(out/"native_SingleCellSignalR_returned_scores.csv.gz")
parts=scsr.native_pair.str.partition("_vs_")
scsr["source_a"]=parts[0]
scsr["source_b"]=np.where(scsr.kind.eq("autocrine"),parts[0],parts[2])
assert (scsr.kind.eq("autocrine")|parts[1].eq("_vs_")).all()
scsr=scsr.rename(columns={"L":"entity_a","R":"entity_b","LR.score":"native_score","pval":"pvalue"})
scsr["score_field"]="LR.score"
scsr["direction_scope"]="native_LR"
scsr["native_run"]=scsr.source_file
scsr["returned_scope"]="only_native_returned_significant_interactions; absent_not_zero"

text_paths=sorted((base/"by_algorithm/CytoTalk").rglob("FinalNetwork.txt"))
rds_paths=set(pd.read_csv(out/"native_R_object_catalog.csv").loc[lambda d:d.method.eq("CytoTalk_native"),"source_file"])
only_text=list(filter(lambda p:str(p.with_name("cytotalk_result.rds")).replace("\\","/") not in rds_paths,text_paths))
empty_network=pd.DataFrame({"is_ct_edge":pd.Series(dtype=bool),"node1_type":pd.Series(dtype="string"),"node2_type":pd.Series(dtype="string"),"node1":pd.Series(dtype="string"),"node2":pd.Series(dtype="string"),"cost":pd.Series(dtype=float)})
cy=pd.concat([empty_network,*map(partial(pd.read_csv,sep="\t"),only_text)],keys=["",*map(str,only_text)],names=["source_file","source_row"]).reset_index()
cy["source_file"]=cy.source_file.str.replace("\\","/",regex=False)
cy["dataset"]=cy.source_file.str.split("/").str[-3]
cy["method"]="CytoTalk_native"
cy=cy.loc[cy.is_ct_edge.eq(True)].copy().rename(columns={"node1_type":"source_a","node2_type":"source_b","node1":"entity_a","node2":"entity_b"})
cy["native_score"]=-cy.cost
cy["pvalue"]=np.nan
cy["score_field"]="negative_native_PCST_edge_cost"
cy["direction_scope"]="native_cross_population_unoriented"
cy["native_pair"]=cy.entity_a+"|"+cy.entity_b+"|"+cy.source_a+"|"+cy.source_b
cy["native_run"]=cy.source_file.str.replace(r"/[^/]+/FinalNetwork[.]txt$","",regex=True)
cy["returned_scope"]="selected_network_text_without_saved_RDS"

paths=sorted((base/"misc").glob("natmi_*_edges.csv"))
names=pd.Index(map(op.attrgetter("stem"),paths)).str.removeprefix("natmi_").str.removesuffix("_edges")
nat=pd.concat(map(pd.read_csv,paths),keys=names,names=["dataset","source_row"]).reset_index()
nat["source_file"]=nat.dataset.map(pd.Series(list(map(str,paths)),index=names))
nat=nat.rename(columns={"Sending cluster":"source_a","Target cluster":"source_b","Ligand symbol":"entity_a","Receptor symbol":"entity_b"})
nat=nat.melt(id_vars=["dataset","source_file","source_row","source_a","source_b","entity_a","entity_b"],value_vars=CONFIG["NATMI_fields"],var_name="score_field",value_name="native_score")
nat["method"]="NATMI_native"
nat["pvalue"]=np.nan
nat["direction_scope"]="native_LR"
nat["native_pair"]=nat.entity_a+"|"+nat.entity_b
nat["native_run"]=nat.source_file
nat["returned_scope"]="all_native_edge_table_entries"

paths=sorted((base/"misc").glob("liana_*_result.csv"))
names=pd.Index(map(op.attrgetter("stem"),paths)).str.removeprefix("liana_").str.removesuffix("_result")
li=pd.concat(map(pd.read_csv,paths),keys=names,names=["dataset","source_row"]).reset_index()
li["source_file"]=li.dataset.map(pd.Series(list(map(str,paths)),index=names))
li=li.rename(columns={"source":"source_a","target":"source_b","ligand_complex":"entity_a","receptor_complex":"entity_b"})
definitions=pd.DataFrame(CONFIG["LIANA_fields"],columns=["method","score_field","direction"])
li=li.melt(id_vars=["dataset","source_file","source_row","source_a","source_b","entity_a","entity_b"],value_vars=definitions.score_field.unique(),var_name="score_field",value_name="native_score").merge(definitions,on="score_field",validate="many_to_many")
li["native_score"]=li.native_score*li.direction
li["pvalue"]=np.nan
li["direction_scope"]="native_LR"
li["native_pair"]=li.entity_a+"|"+li.entity_b
li["native_run"]=li.source_file
li["returned_scope"]="all_saved_LIANA_rows_and_named_score_implementations"

mean_paths=sorted((base/"by_algorithm/CellPhoneDB").glob("*/statistical_analysis_means_*.txt"))
pval_paths=sorted((base/"by_algorithm/CellPhoneDB").glob("*/statistical_analysis_pvalues_*.txt"))
paths=mean_paths+pval_paths
artifact_names=list(map(lambda p:str(p).replace("\\","/"),paths))
headers=list(map(partial(pd.read_csv,sep="\t",nrows=0),paths))
schema=pd.DataFrame({"native_artifact":artifact_names,"celltype_pair":list(map(op.attrgetter("columns"),headers))}).explode("celltype_pair")
schema=schema.loc[schema.celltype_pair.str.contains("|",regex=False)]
cp=pd.concat(map(partial(pd.read_csv,sep="\t"),paths),keys=artifact_names,names=["native_artifact","source_row"]).reset_index()
cp["dataset"]=cp.native_artifact.str.split("/").str[-2]
cp["native_run"]=cp.native_artifact.str.replace(r"statistical_analysis_(means|pvalues)_","statistical_analysis_run_",regex=True)
cp["score_field"]=np.where(cp.native_artifact.str.contains("statistical_analysis_means_"),"native_mean","negative_native_permutation_pvalue")
cp["direction"]=np.where(cp.score_field.eq("native_mean"),1,-1)
deconv_paths=sorted((base/"by_algorithm/CellPhoneDB").glob("*/statistical_analysis_deconvoluted_*.txt"))
deconv=pd.concat(map(partial(pd.read_csv,sep="\t"),deconv_paths),keys=list(map(lambda p:str(p).replace("\\","/"),deconv_paths)),names=["native_artifact","source_row"]).reset_index()
deconv["native_run"]=deconv.native_artifact.str.replace("statistical_analysis_deconvoluted_","statistical_analysis_run_",regex=False)
components=deconv.loc[deconv.is_complex.eq(True)&deconv.gene_name.notna(),["native_run","complex_name","gene_name"]].drop_duplicates().sort_values(["native_run","complex_name","gene_name"])
complexes=components.groupby(["native_run","complex_name"],observed=True).gene_name.agg("_".join).rename("components").reset_index()
cp["partner_a_name"]=cp.partner_a.str.removeprefix("complex:")
cp["partner_b_name"]=cp.partner_b.str.removeprefix("complex:")
cp=cp.merge(complexes.rename(columns={"complex_name":"partner_a_name","components":"complex_a"}),on=["native_run","partner_a_name"],how="left",validate="many_to_one")
cp=cp.merge(complexes.rename(columns={"complex_name":"partner_b_name","components":"complex_b"}),on=["native_run","partner_b_name"],how="left",validate="many_to_one")
cp["entity_a"]=cp.gene_a.fillna(cp.complex_a)
cp["entity_b"]=cp.gene_b.fillna(cp.complex_b)
cp=cp.melt(id_vars=["dataset","native_artifact","source_row","native_run","score_field","direction","id_cp_interaction","entity_a","entity_b"],value_vars=schema.celltype_pair.unique(),var_name="celltype_pair",value_name="native_score")
cp=cp.merge(schema,on=["native_artifact","celltype_pair"],validate="many_to_one")
pair=cp.celltype_pair.str.split("|",n=1,expand=True,regex=False)
cp["source_a"]=pair[0]
cp["source_b"]=pair[1]
cp["native_score"]=pd.to_numeric(cp.native_score)*cp.direction
cp["method"]="CellPhoneDB_native"
cp["pvalue"]=np.where(cp.score_field.eq("negative_native_permutation_pvalue"),-cp.native_score,np.nan)
cp["direction_scope"]="native_partner_A_B_unoriented"
cp["source_file"]=cp.native_artifact
cp["native_pair"]=cp.id_cp_interaction.astype(str)
cp["returned_scope"]="all_native_partner_pair_tests; missing_scores_preserved"

ticci_score_paths=list(map(lambda p:p.with_name("ticci_cellchat_interactions.csv"),ticci_cell_paths))
ticci_scores=pd.concat(map(pd.read_csv,ticci_score_paths),keys=ticci_datasets,names=["dataset","row"]).reset_index()
ticci_scores=ticci_scores.rename(columns={"source":"source_a","target":"source_b","prob":"native_score","pval":"pvalue","interaction_name":"native_pair"})
molecule_dictionary=r_scores.loc[r_scores.method.eq("CellChat"),["dataset","native_pair","entity_a","entity_b"]].drop_duplicates()
assert not molecule_dictionary.duplicated(["dataset","native_pair"]).any()
ticci_scores=ticci_scores.merge(molecule_dictionary,on=["dataset","native_pair"],how="left",validate="many_to_one")
ticci_display=ticci_scores.interaction_name_2.str.replace(r"[()]","",regex=True).str.split(r"\s+-\s+",n=1,expand=True)
ticci_scores["entity_a"]=ticci_scores.entity_a.fillna(ticci_display[0].str.replace("+","_",regex=False))
ticci_scores["entity_b"]=ticci_scores.entity_b.fillna(ticci_display[1].str.replace("+","_",regex=False))
ticci_scores["source_file"]=ticci_scores.dataset.map(pd.Series(list(map(lambda p:str(p).replace("\\","/"),ticci_score_paths)),index=ticci_datasets))
ticci_scores["native_run"]=ticci_scores.source_file
ticci_scores["method"]="TICCI_CellChat"
ticci_scores["score_field"]="native_CellChat_selected_probability"
ticci_scores["direction_scope"]="native_LR"
ticci_scores["returned_scope"]="native_TICCI_CellChat_stage_calls;complexes_from_saved_CellChat_DB_or_explicit_native_interaction_name_2"
sc_dir=root/"results/by_algorithm/scTenifoldXct/paul15"
sc_paths=sorted((sc_dir/"sctenifoldxct_pairs").glob("*.pkl"))
sc_bundles=list(map(lambda p:pickle.loads(p.read_bytes()),sc_paths))
sc_pair_meta=pd.DataFrame(map(op.itemgetter("coverage"),sc_bundles))
assert len(sc_bundles)==len(sc_pair_meta[["source_celltype","target_celltype"]].drop_duplicates())
assert sc_pair_meta.result_status.eq("completed_native_ordered_pair").all()
sc_cells=pd.DataFrame({"source_file":pd.Series(sc_paths,dtype="string"),"source_label":sc_pair_meta.source_celltype,
 "target_label":sc_pair_meta.target_celltype,"source_cells":list(map(op.itemgetter("source_cell_ids"),sc_bundles)),
 "target_cells":list(map(op.itemgetter("target_cell_ids"),sc_bundles))})
sc_cells=sc_cells.melt(id_vars=["source_file","source_label","target_label"],value_vars=["source_cells","target_cells"],var_name="role",value_name="cell_id").explode("cell_id")
sc_cells["native_label"]=np.where(sc_cells.role.eq("source_cells"),sc_cells.source_label,sc_cells.target_label)
sc_cells["cell_id"]=sc_cells.cell_id.astype("string")
sc_full_meta=meta.loc[meta.dataset.eq("paul15"),["cell_id","paul15_clusters"]]
assert not sc_full_meta.cell_id.duplicated().any()
sc_bound=sc_cells.merge(sc_full_meta,on="cell_id",how="left",validate="many_to_one")
assert sc_bound.native_label.eq(sc_bound.paul15_clusters).all()
sc_signatures=sc_cells.sort_values("cell_id").groupby(["source_file","native_label"],observed=True).cell_id.agg(tuple).reset_index()
sc_expected=sc_full_meta.sort_values("cell_id").groupby("paul15_clusters",observed=True).cell_id.agg(tuple)
assert sc_signatures.cell_id.eq(sc_signatures.native_label.map(sc_expected)).all()
sc_genes=pd.read_csv(sc_dir/"native_gene_id_binding.csv")
sc_ambiguous_gene_names=pd.Index(sc_genes.loc[sc_genes.native_gene_id.duplicated(keep=False),"native_gene_id"].unique())
sc_expected_gene_order=sc_genes.native_gene_id.to_numpy()
assert np.array_equal(np.vstack(list(map(op.itemgetter("source_gene_names"),sc_bundles))),np.broadcast_to(sc_expected_gene_order,(len(sc_bundles),len(sc_expected_gene_order))))
assert np.array_equal(np.vstack(list(map(op.itemgetter("target_gene_names"),sc_bundles))),np.broadcast_to(sc_expected_gene_order,(len(sc_bundles),len(sc_expected_gene_order))))
sc_gene_values=np.unique(np.concatenate(list(map(op.itemgetter("source_gene_names"),sc_bundles))+list(map(op.itemgetter("target_gene_names"),sc_bundles))))
assert pd.Index(sc_gene_values).difference(sc_genes.native_gene_id).empty
sc_scores=pd.concat(map(op.itemgetter("aligned_candidate_distances"),sc_bundles),keys=np.arange(len(sc_paths)),names=["source_id","native_row"]).reset_index()
sc_scores["dataset"]="paul15"
sc_scores["method"]="scTenifoldXct_native"
sc_scores["native_run"]=(sc_dir/"coverage.csv").as_posix()
sc_scores["source_file"]=sc_scores.source_id.map(pd.Series(list(map(str,sc_paths))))
sc_scores=sc_scores.rename(columns={"source_celltype":"source_a","target_celltype":"source_b","ligand":"entity_a","receptor":"entity_b","index":"native_pair"})
sc_scores["native_score"]=-sc_scores.dist
sc_scores["pvalue"]=np.nan
sc_scores["score_field"]="negative_native_aligned_distance"
sc_scores["direction_scope"]="native_LR"
sc_scores["returned_scope"]="all_native_eligible_cartesian_LR_candidates_of_completed_ordered_pairs;distance_not_pvalue;"+str(len(sc_paths))+"_of_342_original_ordered_pairs"
columns=["dataset","method","native_run","source_file","source_a","source_b","entity_a","entity_b","native_score","pvalue","score_field","direction_scope","native_pair","returned_scope"]
scores=pd.concat([r_scores[columns],scsr[columns],cy[columns],nat[columns],li[columns],cp[columns],ticci_scores[columns],sc_scores[columns]],ignore_index=True)
scores[["source_a","source_b","entity_a","entity_b"]]=scores[["source_a","source_b","entity_a","entity_b"]].astype("string")
scores["finite_score"]=np.isfinite(scores.native_score)
scores["native_molecule_identity_resolved"]=~(scores.method.eq("scTenifoldXct_native")&(scores.entity_a.isin(sc_ambiguous_gene_names)|scores.entity_b.isin(sc_ambiguous_gene_names)))
cyto_pair_validity=pd.read_csv(out/"native_CytoTalk_pair_coverage.csv").loc[lambda d:d.saved_network.eq(True)]
assert not cyto_pair_validity.network_file.duplicated().any()
cyto_model_valid=pd.Series(cyto_pair_validity.valid_native_network.eq(True).to_numpy(),index=cyto_pair_validity.network_file)
cyto_model_key=scores.source_file.str.replace(r"/cytotalk_result[.]rds$","/FinalNetwork.txt",regex=True)
assert cyto_model_key.loc[scores.method.eq("CytoTalk_native")].isin(cyto_model_valid.index).all()
scores["native_model_valid"]=~scores.method.eq("CytoTalk_native")|cyto_model_key.map(cyto_model_valid).fillna(False)
groups=["dataset","method","native_run","score_field"]
native_counts=scores.groupby(groups,observed=True).agg(n_native_rows=("native_score","size"),n_finite_native_rows=("finite_score","sum"),returned_scope=("returned_scope","first")).reset_index()
direct=scores.method.isin(["CellChat","SingleCellSignalR_native","TICCI_CellChat"])
direct_a=object_map.rename(columns={"native_label":"source_a","reference_label":"mapped_source"})
direct_b=object_map.rename(columns={"native_label":"source_b","reference_label":"mapped_target"}).drop(columns=["reference_dataset","identity_basis"])
direct_a=direct_a.drop(columns=["reference_dataset","identity_basis"])
direct_scores=scores.loc[direct].merge(object_cohort,on=["dataset","method","source_file"],how="left",validate="many_to_one").merge(direct_a,on=["dataset","method","source_file","source_a"],how="left",validate="many_to_one").merge(direct_b,on=["dataset","method","source_file","source_b"],how="left",validate="many_to_one")
snapshot_a=snapshot_map.rename(columns={"native_label":"source_a","reference_label":"mapped_source"}).drop(columns=["reference_dataset","identity_basis"])
snapshot_b=snapshot_map.rename(columns={"native_label":"source_b","reference_label":"mapped_target"}).drop(columns=["reference_dataset","identity_basis"])
prefixed_a=snapshot_a.loc[snapshot_a.source_a.str.fullmatch(r"[0-9]+")].copy()
prefixed_a["source_a"]="C"+prefixed_a.source_a
prefixed_b=snapshot_b.loc[snapshot_b.source_b.str.fullmatch(r"[0-9]+")].copy()
prefixed_b["source_b"]="C"+prefixed_b.source_b
nat_a=pd.concat([snapshot_a,prefixed_a],ignore_index=True).drop_duplicates()
nat_b=pd.concat([snapshot_b,prefixed_b],ignore_index=True).drop_duplicates()
snap_regular=scores.loc[~direct&~scores.method.eq("NATMI_native")].merge(snapshot_cohort,on="dataset",how="left",validate="many_to_one").merge(snapshot_a,on=["dataset","source_a"],how="left",validate="many_to_one").merge(snapshot_b,on=["dataset","source_b"],how="left",validate="many_to_one")
snap_nat=scores.loc[scores.method.eq("NATMI_native")].merge(snapshot_cohort,on="dataset",how="left",validate="many_to_one").merge(nat_a,on=["dataset","source_a"],how="left",validate="many_to_one").merge(nat_b,on=["dataset","source_b"],how="left",validate="many_to_one")
mapped=pd.concat([direct_scores,snap_regular,snap_nat],ignore_index=True)
entities=pd.concat([mapped.entity_a,mapped.entity_b,ref.ligand,ref.receptor],ignore_index=True).dropna().drop_duplicates().reset_index(drop=True).astype(str).rename("entity")
parts=entities.str.split("_").explode().rename("gene").to_frame().join(entities)
parts["gene"]=parts.gene.str.upper()
canonical=parts.drop_duplicates().sort_values(["entity","gene"]).groupby("entity",sort=False).gene.agg("_".join)
mapped["entity_a"]=mapped.entity_a.map(canonical)
mapped["entity_b"]=mapped.entity_b.map(canonical)
ref["ligand"]=ref.ligand.map(canonical)
ref["receptor"]=ref.receptor.map(canonical)
ref=ref.rename(columns={"dataset":"reference_dataset"})
ref_keys=["reference_dataset","sender","receiver","ligand","receptor"]
ref_audit=ref.groupby(ref_keys,observed=True).reference_score.agg(["size","min","max"]).reset_index()
ref_audit["unique_score"]=ref_audit["min"].eq(ref_audit["max"])
ref_unique=ref_audit.loc[ref_audit.unique_score,ref_keys+["min"]].rename(columns={"min":"reference_score"})
direction_pairs=ref_unique[["reference_dataset","ligand","receptor"]].drop_duplicates()
forward=direction_pairs.rename(columns={"ligand":"entity_a","receptor":"entity_b"}).assign(orientation=1)
backward=direction_pairs.loc[direction_pairs.ligand.ne(direction_pairs.receptor)].rename(columns={"ligand":"entity_b","receptor":"entity_a"}).assign(orientation=-1)
orientations=pd.concat([forward,backward],ignore_index=True).drop_duplicates()
orientation_counts=orientations.groupby(["reference_dataset","entity_a","entity_b"],observed=True).orientation.agg(["size","first"]).reset_index()
mapped=mapped.merge(orientation_counts,on=["reference_dataset","entity_a","entity_b"],how="left",validate="many_to_one").rename(columns={"size":"n_resource_directions","first":"resource_orientation"})
native_lr=mapped.direction_scope.eq("native_LR")
mapped["orientation_resolved"]=native_lr|mapped.n_resource_directions.eq(1)
flip=~native_lr&mapped.resource_orientation.eq(-1)
mapped["sender"]=np.where(flip,mapped.mapped_target,mapped.mapped_source)
mapped["receiver"]=np.where(flip,mapped.mapped_source,mapped.mapped_target)
mapped["ligand"]=np.where(flip,mapped.entity_b,mapped.entity_a)
mapped["receptor"]=np.where(flip,mapped.entity_a,mapped.entity_b)
mapped["population_mapped"]=mapped.mapped_source.notna()&mapped.mapped_target.notna()
mapped["molecule_mapped"]=mapped.entity_a.notna()&mapped.entity_b.notna()&mapped.native_molecule_identity_resolved
mapped["eligible_mapping"]=mapped.finite_score&mapped.population_mapped&mapped.molecule_mapped&mapped.orientation_resolved&mapped.native_model_valid
mapping_counts=mapped.groupby(groups,observed=True).agg(n_population_mapped=("population_mapped","sum"),n_molecule_resolved=("molecule_mapped","sum"),n_direction_resolved=("orientation_resolved","sum"),n_eligible_mapping_rows=("eligible_mapping","sum"),reference_dataset=("reference_dataset","first"),identity_basis=("identity_basis","first")).reset_index()
valid=mapped.loc[mapped.eligible_mapping].copy()
atomic=valid.groupby(groups+ref_keys,observed=True).native_score.agg(["size","min","max"]).reset_index()
atomic["score_unambiguous"]=atomic["min"].eq(atomic["max"])
ambiguities=atomic.loc[~atomic.score_unambiguous].copy()
atomic=atomic.loc[atomic.score_unambiguous].rename(columns={"min":"native_score"})
aligned=atomic.merge(ref_unique,on=ref_keys,how="inner",validate="many_to_one")
aligned["native_rank"]=aligned.groupby(groups,observed=True).native_score.rank(method="average")
aligned["reference_rank"]=aligned.groupby(groups,observed=True).reference_score.rank(method="average")
aligned["native_centered"]=aligned.native_rank-aligned.groupby(groups,observed=True).native_rank.transform("mean")
aligned["reference_centered"]=aligned.reference_rank-aligned.groupby(groups,observed=True).reference_rank.transform("mean")
aligned["cross"]=aligned.native_centered*aligned.reference_centered
aligned["native_ss"]=aligned.native_centered**2
aligned["reference_ss"]=aligned.reference_centered**2
agreement=aligned.groupby(groups,observed=True).agg(n_common_atomic_keys=("native_score","size"),rank_cross=("cross","sum"),native_ss=("native_ss","sum"),reference_ss=("reference_ss","sum")).reset_index()
agreement["Spearman_rank_agreement"]=agreement.rank_cross/np.sqrt(agreement.native_ss*agreement.reference_ss)
summary=native_counts.merge(mapping_counts,on=groups,how="left",validate="one_to_one").merge(agreement,on=groups,how="left",validate="one_to_one")
summary["n_common_atomic_keys"]=summary.n_common_atomic_keys.fillna(0).astype(int)
summary["n_reference_atomic_keys"]=summary.reference_dataset.map(ref_unique.groupby("reference_dataset").size())
summary["reference_key_coverage"]=summary.n_common_atomic_keys/summary.n_reference_atomic_keys
summary["n_ambiguous_native_atomic_keys"]=summary.set_index(groups).index.map(ambiguities.groupby(groups,observed=True).size()).fillna(0).astype(int)
summary["comparison_status"]=np.select([summary.reference_dataset.isna(),summary.n_eligible_mapping_rows.eq(0),summary.n_common_atomic_keys.lt(2),summary.native_ss.eq(0)|summary.reference_ss.eq(0)],["no_verified_complete_cohort_binding","population_molecule_or_direction_unresolved","insufficient_common_atomic_keys","constant_score_on_common_keys"],default="compared_descriptive")
unresolved_native_gene_rows=scores.assign(unresolved_native_gene_name=~scores.native_molecule_identity_resolved).groupby(groups,observed=True).unresolved_native_gene_name.sum()
summary["n_unresolved_native_gene_name_rows"]=summary.set_index(groups).index.map(unresolved_native_gene_rows)
invalid_native_model_rows=scores.assign(invalid_native_model=~scores.native_model_valid).groupby(groups,observed=True).invalid_native_model.sum()
summary["n_invalid_native_model_score_rows"]=summary.set_index(groups).index.map(invalid_native_model_rows)
summary.loc[summary.n_invalid_native_model_score_rows.gt(0),"comparison_status"]+=";invalid_native_models_excluded_from_metric"
summary["native_model_validation_scope"]=np.where(summary.method.eq("CytoTalk_native"),"original_population_pair_plan_and_native_core_input_finiteness","not_this_CytoTalk_failure_check")
summary["n_full_native_input_cells"]=np.where(summary.method.eq("scTenifoldXct_native"),int(sc_pair_meta.n_input_cells.iloc[0]),np.nan)
summary["n_full_native_input_gene_rows"]=np.where(summary.method.eq("scTenifoldXct_native"),len(sc_genes),np.nan)
summary["n_distinct_native_gene_names"]=np.where(summary.method.eq("scTenifoldXct_native"),sc_genes.native_gene_id.nunique(),np.nan)
summary["n_ambiguous_native_gene_names"]=np.where(summary.method.eq("scTenifoldXct_native"),len(sc_ambiguous_gene_names),np.nan)
summary["planned_native_population_pairs"]=np.where(summary.method.eq("scTenifoldXct_native"),342,np.nan)
summary["completed_native_population_pairs"]=np.where(summary.method.eq("scTenifoldXct_native"),len(sc_paths),np.nan)
summary["native_population_pair_scope"]=np.where(summary.method.eq("scTenifoldXct_native"),
 str(len(sc_paths))+"_of_342_original_ordered_pairs","original_saved_native_scope")
summary.loc[summary.method.eq("scTenifoldXct_native")&summary.completed_native_population_pairs.lt(summary.planned_native_population_pairs),"comparison_status"]+=";partial_native_ordered_pair_scope"
summary["metric_family"]="classical_rank_statistic_descriptive_without_GT"
summary["ground_truth_accuracy"]=False
summary["reference"]="unchanged_saved_native_static_LR_scores; source_model_binding_saved"
summary["primary_dataset_key"]=summary.dataset.eq(summary.reference_dataset)
summary["mapping_protocol"]="exact_full_cohort; complete_population_correspondence; explicit_complex_subunits; unique_LR_direction_for_unoriented_native_edges"

returned_keys=mapped.drop_duplicates(groups+["source_a","source_b","entity_a","entity_b"]).groupby(groups,observed=True).size()
summary["n_distinct_native_returned_keys"]=summary.set_index(groups).index.map(returned_keys)
summary["native_returned_key_coverage"]=summary.n_common_atomic_keys/summary.n_distinct_native_returned_keys
summary["n_nonfinite_native_scores"]=summary.n_native_rows-summary.n_finite_native_rows
summary["n_population_unmapped_rows"]=summary.n_native_rows-summary.n_population_mapped
summary["mapping_failure_detail"]=np.select([
summary.comparison_status.eq("no_verified_complete_cohort_binding"),
summary.method.eq("NATMI_native")&summary.dataset.isin(["gastrulation_erythroid","larry_invitro"]),
summary.n_population_mapped.eq(0)&summary.method.str.contains("LIANA"),
summary.n_direction_resolved.eq(0),
summary.comparison_status.eq("insufficient_common_atomic_keys")],
["no_saved_reference_model_with_identical_complete_cell_cohort",
"native_grouping_uses_sample_or_Library_field;not_current_cell_type_partition",
"native_Leiden_population_membership_not_saved;cannot_map_group_ids_to_reference_cell_types",
"no_unique_molecular_LR_direction_in_saved_reference_resource",
"native_selected_interactions_have_less_than_two_common_atomic_keys"],
default="")

summary.drop(columns=["rank_cross","native_ss","reference_ss"]).to_csv(out/"native_CCC_comparison_metrics.csv",index=False,encoding="utf-8")
aligned[groups+ref_keys+["native_score","reference_score","size"]].to_csv(out/"native_CCC_aligned_scores.csv.gz",index=False,encoding="utf-8")
ambiguities.to_csv(out/"native_CCC_ambiguous_atomic_scores.csv",index=False,encoding="utf-8")
ref_audit.loc[~ref_audit.unique_score].to_csv(out/"static_reference_ambiguous_atomic_scores.csv",index=False,encoding="utf-8")
object_cohort.to_csv(out/"native_CCC_object_cohort_binding.csv",index=False,encoding="utf-8")
snapshot_cohort.to_csv(out/"native_CCC_snapshot_cohort_binding.csv",index=False,encoding="utf-8")
print(summary.groupby(["method","comparison_status"],observed=True).agg(n_runs_scores=("dataset","size"),n_datasets=("dataset","nunique"),n_common_keys=("n_common_atomic_keys","sum")).to_string())
