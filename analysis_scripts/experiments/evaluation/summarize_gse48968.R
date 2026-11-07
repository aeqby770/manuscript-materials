CONFIG <- list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),
 input="results/validation/gse48968",comparison="results/by_algorithm/GSE48968_BMDC_WT/same_resource_pair_rankings.csv",
 output="results/paper_tables/current",timetalk="results/by_algorithm/TimeTalk/GSE48968_BMDC_WT",threads=parallel::detectCores())
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads)
library(data.table)
setDTthreads(CONFIG$threads)
IN <- file.path(CONFIG$root,CONFIG$input)
OUT <- file.path(CONFIG$root,CONFIG$output)
pair_source <- fread(file.path(CONFIG$root,CONFIG$comparison))
pmap <- unique(pair_source[,.(pair_idx,ligand_genes,receptor_genes,ligand_id,receptor_id,interaction_name,pathway_name,interaction_aliases)])
base <- pair_source[!startsWith(method,"LineaLinkR"),.(pair_idx,stimulation,score,method,time_window_measurements)]
gg <- fread(file.path(IN,"external_validation_time_root_legacy/native_pair_ranking.csv"))
gq <- fread(file.path(IN,"external_validation_time_root/native_pair_ranking.csv"))
shape <- fread(file.path(IN,"external_validation_time_root_shape/native_pair_ranking.csv"))
layout <- fread(file.path(IN,"external_validation_time_root_layout/native_pair_ranking.csv"))
gg[,`:=`(method="LineaLinkR_G_point_G_curve",stimulation="all_primary_conditions")]
gq[,`:=`(method="LineaLinkR_G_point_Q_curve",stimulation="all_primary_conditions")]
shape[,`:=`(method="LineaLinkR_G_amplitude_Q_shape",stimulation="all_primary_conditions")]
layout[,`:=`(method="LineaLinkR_G_amplitude_Q_layout",stimulation="all_primary_conditions")]
own <- rbindlist(list(gg,gq,shape,layout),use.names=TRUE,fill=TRUE)
shared <- own[is.finite(score),.(n_methods=uniqueN(method)),by=pair_idx][n_methods==4L,pair_idx]
ranks <- merge(rbindlist(list(base,own[,.(pair_idx,stimulation,score,method)]),fill=TRUE),pmap,by="pair_idx",all.x=TRUE)
ranks[,`:=`(in_common_evaluable_panel=pair_idx%in%shared,n_common_candidates=length(shared))]
ranks[,score_common:=fifelse(in_common_evaluable_panel,score,NA_real_)]
ranks[,rank_common:=frank(-score_common,ties.method="min",na.last="keep"),by=.(method,stimulation)]
ranks[,rank:=frank(-score,ties.method="min",na.last="keep"),by=.(method,stimulation)]
ranks[,n_candidates:=.N,by=.(method,stimulation)]
ranks[,n_scored:=sum(is.finite(score)),by=.(method,stimulation)]
fwrite(ranks,file.path(CONFIG$root,CONFIG$comparison))
fwrite(ranks[stimulation=="all_primary_conditions"&ligand_genes=="IFNB1"&receptor_genes=="IFNAR1_IFNAR2"],file.path(OUT,"gse48968_axis_retrieval.csv"))
pair_axis <- ranks[stimulation=="all_primary_conditions"&ligand_genes=="IFNB1"&receptor_genes=="IFNAR1_IFNAR2",
 .(method,unit="molecular_pair",rank_all=rank,n_all=n_scored,rank_common,n_common=n_common_candidates,
 score,output_status="scored")]
entrain <- fread(file.path(CONFIG$root,"results/by_algorithm/ENTRAIN/GSE48968_BMDC_WT/ligand_importance.csv"))
entrain[,method:="ENTRAIN_native_regression"]
ligands <- ranks[stimulation=="all_primary_conditions",
 .(score=sort(score,decreasing=TRUE,na.last=TRUE)[1L]),by=.(method,ligand=ligand_genes)]
ligands <- rbindlist(list(ligands,entrain[,.(method,ligand,score)]))
common_ligands <- ligands[is.finite(score),.(n_methods=uniqueN(method)),by=ligand][n_methods==uniqueN(ligands$method),ligand]
ligands[,score_common:=fifelse(ligand%in%common_ligands,score,NA_real_)]
ligands[,rank_common:=frank(-score_common,ties.method="min",na.last="keep"),by=method]
ligands[,rank_all:=frank(-score,ties.method="min",na.last="keep"),by=method]
ligands[,n_all:=sum(is.finite(score)),by=method]
ligands[,n_common:=length(common_ligands)]
fwrite(ligands,file.path(IN,"all_method_ligand_rankings.csv"))
ligand_axis <- ligands[ligand=="IFNB1",
 .(method,unit="ligand",rank_all,n_all,rank_common,n_common,score,output_status="scored")]
tt_dir <- file.path(CONFIG$root,CONFIG$timetalk)
timetalk <- fread(file.path(tt_dir,"coverage.csv"))
tt_strength <- fread(file.path(tt_dir,"TimeTalk_condition_pair_strength.csv"))[
 ,.(score=max(score)),by=.(ligand_genes=L,receptor_genes=R)]
tt_pairs <- merge(pmap,tt_strength,by=c("ligand_genes","receptor_genes"),all.x=TRUE)
tt_pairs[,score:=fcoalesce(score,0)]
tt_pairs[,rank_all:=frank(-score,ties.method="min")]
tt_common <- fread(file.path(tt_dir,"TimeTalk_common_pair_ranking.csv"))
tt_ligand_common <- fread(file.path(tt_dir,"TimeTalk_common_ligand_ranking.csv"))
tt_pair_target <- merge(tt_pairs[ligand_genes=="IFNB1"&receptor_genes=="IFNAR1_IFNAR2"],
 tt_common[,.(pair_idx,rank_common=rank_lower,common_score=score)],by="pair_idx")
stopifnot(nrow(tt_pair_target)==1L,tt_pair_target$score==tt_pair_target$common_score)
tt_ligands <- tt_pairs[,.(score=max(score)),by=.(ligand=ligand_genes)]
tt_ligands[,rank_all:=frank(-score,ties.method="min")]
tt_ligand_target <- merge(tt_ligands[ligand=="IFNB1"],
 tt_ligand_common[,.(ligand,rank_common=rank_lower,common_score=score)],by="ligand")
stopifnot(nrow(tt_ligand_target)==1L,tt_ligand_target$score==tt_ligand_target$common_score)
tt_axis <- rbindlist(list(
 tt_pair_target[,.(method="TimeTalk",unit="molecular_pair",rank_all,n_all=nrow(tt_pairs),
 rank_common,n_common=nrow(tt_common),score,output_status=timetalk$result_status)],
 tt_ligand_target[,.(method="TimeTalk",unit="ligand",rank_all,n_all=nrow(tt_ligands),
 rank_common,n_common=nrow(tt_ligand_common),score,output_status=timetalk$result_status)]))
table <- rbindlist(list(pair_axis,ligand_axis,tt_axis),use.names=TRUE)
fwrite(table,file.path(OUT,"table11_independent_paracrine_axis.csv"))
profiles <- fread(file.path(IN,"all_experiment_module_profiles.csv"))[
 gene_universe=="same_genes_measured_in_both_matrices"&CLUSTER%chin%c("Id","IIIc","IIId")]
fwrite(profiles,file.path(OUT,"supp_independent_paracrine_modules.csv"))
meta <- fread(file.path(CONFIG$root,"inputs/external_validation/GSE48968/GSE48968_all_library_metadata.csv"))
design <- meta[,.(n_libraries=.N,n_author_quality_passed=sum(author_quality_passed),n_single_cells=sum(single_cell),
 n_used_in_graph=sum(model_input)),by=.(matrix_source,experimental_role)]
fwrite(design,file.path(OUT,"supp_independent_paracrine_design.csv"))
print(table)
