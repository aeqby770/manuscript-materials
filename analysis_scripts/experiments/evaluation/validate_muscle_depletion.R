# Sample-level perturbation readout for fixed FAP-origin predictions.
# A ligand-target program measures transcriptional response in the receiver;
# neither receptor downregulation nor ligand downregulation is required.
CONFIG <- list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE), threads=parallel::detectCores(),
  k=c(10L,20L), tolerance=1e-12,
  nichenet_native="results/by_algorithm/NicheNet/GSE143437_muscle_regeneration/nichenet_native.rds",
  workspace=Sys.getenv("LLR_BIO_WS","cache/final065"),
  frozen_panel="results/validation/sample_dynamic/predictions/native_pair_union.csv",
  candidate_dir=Sys.getenv("LLR_BIO_EVAL_DIR","results/evaluation/biology_final065"),
  baseline_dir=Sys.getenv("LLR_BIO_BASELINE_DIR","results/evaluation/biology_final065_ablation/legacy_matched"),
  output=Sys.getenv("LLR_VALIDATION_OUT_DIR","results/validation/final065/tables"),
  data_output=Sys.getenv("LLR_VALIDATION_DATA_OUT","results/validation/final065/predictions"),
  main_table=Sys.getenv("LLR_VALIDATION_MAIN_TABLE","results/validation/final065/table6_real_communication.csv"))
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads, OPENBLAS_NUM_THREADS=CONFIG$threads, MKL_NUM_THREADS=CONFIG$threads)
library(data.table)
.libPaths(c(Sys.getenv("LLR_PACKAGE_LIBRARY",file.path(CONFIG$root,"library")),file.path(CONFIG$root,"library"),.libPaths()))
library(LineaLinkR)
methods <- asNamespace("LineaLinkR")
setDTthreads(CONFIG$threads)
RAW <- file.path(CONFIG$root,"results/validation/gse291091")
OUT <- file.path(CONFIG$root,CONFIG$output)
RESULT <- file.path(CONFIG$root,CONFIG$data_output)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
dir.create(RESULT,recursive=TRUE,showWarnings=FALSE)
cand <- fread(file.path(CONFIG$root,CONFIG$candidate_dir,"muscle_native_fap_candidates.csv"))
previous <- fread(file.path(CONFIG$root,CONFIG$baseline_dir,"muscle_native_fap_candidates.csv"))
unadjusted <- copy(cand)
unadjusted[,event_ranking_score := sample_unadjusted_ranking_score]
tt <- fread(file.path(CONFIG$root,"results/by_algorithm/reused/GSE143437_muscle_regeneration/TimeTalk_native_lr_summary.csv.gz"))
receiver_map <- data.table(target_label=c("Anti-inflammatory macrophages","Resident Macrophages/APCs","Monocytes/Macrophages/Platelets","Pro-inflammatory macrophages"),
  context=c("C1qa macrophages","APC","Gpnmb macrophages","Inflammatory MNP"))
llr <- merge(cand[evaluable %in% TRUE],receiver_map,by="target_label")
llr <- rbindlist(list(llr,copy(llr)[,context:="All MNP"]))
llr <- llr[,.(llr_full=sort(event_ranking_score,decreasing=TRUE,na.last=TRUE)[1L],llr_nofork=sort(no_fork_ranking_score,decreasing=TRUE,na.last=TRUE)[1L],n_events=.N),
  by=.(context,ligand=ligand_genes,receptor=receptor_genes)]

# Compare scores on the same molecular candidates, without selecting the
# receiver or the evaluation budget from its observed perturbation response.
reference_events <- rbindlist(list(copy(previous)[,clock:="llr_previous"],
  copy(unadjusted)[,clock:="llr_unadjusted"]),fill=TRUE)
reference_events <- merge(reference_events[evaluable %in% TRUE],receiver_map,by="target_label")
reference_events <- rbindlist(list(reference_events,copy(reference_events)[,context:="All MNP"]))
reference_scores <- reference_events[,.(score=sort(event_ranking_score,decreasing=TRUE,na.last=TRUE)[1L]),
  by=.(clock,context,ligand=ligand_genes,receptor=receptor_genes)]
reference_scores[,clock:=factor(clock,levels=c("llr_previous","llr_unadjusted"))]
reference_scores <- dcast(reference_scores,context+ligand+receptor~clock,value.var="score",drop=c(TRUE,FALSE))

tc <- merge(tt[sender=="FAPs"],receiver_map,by.x="receiver",by.y="target_label")
tc[, native_score:=native_timetalk_score]
tc <- rbindlist(list(tc,copy(tc)[,context:="All MNP"]))
tc <- tc[,.(timetalk=max(native_score),timetalk_scc=max(mean_signed_scc),n_native_eLR=sum(n_native_eLR)),
  by=.(context,ligand,receptor)]
pair_union <- merge(llr,tc,by=c("context","ligand","receptor"),all=TRUE)
pair_union <- merge(pair_union,reference_scores,by=c("context","ligand","receptor"),all=TRUE)
pair_union[,common_pair:=is.finite(llr_full)&is.finite(llr_nofork)&is.finite(timetalk)&
  is.finite(llr_previous)&is.finite(llr_unadjusted)]
fwrite(pair_union,file.path(RESULT,"native_pair_union.csv"))
coverage <- pair_union[,.(n_union_pairs=.N,n_llr_pairs=sum(is.finite(llr_full)),
  n_previous_pairs=sum(is.finite(llr_previous)),n_unadjusted_pairs=sum(is.finite(llr_unadjusted)),
  n_new_timetalk_pairs=sum(is.finite(llr_full)&is.finite(timetalk)),
  n_timetalk_pairs=sum(is.finite(timetalk)),n_timetalk_eLR_pairs=sum(timetalk>0,na.rm=TRUE),
  n_common_pairs=sum(common_pair)),by=context]
fwrite(coverage,file.path(OUT,"supp_muscle_depletion_coverage.csv"))

# The unsigned NicheNet prior is used only as an independent readout of response
# magnitude, never as a causal sign or a receptor-specific ground truth.
program <- fread(file.path(RAW,"ligand_target_response_assignments.csv"))
response <- melt(program,id.vars="ligand",variable.name="context_assignment",value.name="response")
response[, c("context","assignment"):=tstrsplit(context_assignment,"|",fixed=TRUE)]
response[,assignment:=as.integer(assignment)]
background <- unique(response[assignment==1L & is.finite(response),.(context,ligand)])
frozen <- fread(file.path(CONFIG$root,CONFIG$frozen_panel))[common_pair %in% TRUE,.(context,ligand,receptor)]
common <- merge(frozen,pair_union,by=c("context","ligand","receptor"),all.x=TRUE)
score_fields <- c("llr_full","llr_nofork","llr_previous","llr_unadjusted","timetalk","timetalk_scc")
ligand_scores <- melt(common,id.vars=c("context","ligand","receptor"),measure.vars=score_fields,
  variable.name="method",value.name="score",variable.factor=FALSE)
ligand_scores <- ligand_scores[,.(score=fifelse(any(is.finite(score)),max(score,na.rm=TRUE),NA_real_)),
  by=.(context,ligand,method)]
ligands <- dcast(ligand_scores,context+ligand~method,value.var="score")
ligands <- merge(ligands,common[,.(n_common_pairs=.N),by=.(context,ligand)],by=c("context","ligand"))
ligands <- merge(ligands,background,by=c("context","ligand"))
# Reuse independent atlas fits. Pooled receivers use the same maximum rule as
# the native TimeTalk adapter; no depletion response enters this ranking.
niche_fit <- readRDS(file.path(CONFIG$root,CONFIG$nichenet_native))
niche <- merge(as.data.table(niche_fit$activities),receiver_map,
  by.x="receiver",by.y="target_label")
niche <- rbindlist(list(niche,copy(niche)[,context:="All MNP"]))
niche_scores <- niche[,.(nichenet=max(pearson)),by=.(context,ligand=test_ligand)]
ligands <- merge(ligands,niche_scores,by=c("context","ligand"),all.x=TRUE)
# Query the actual native nominator on the pre-existing molecular/background panel.
# Receiver weights average observed type fractions over biological samples.
fwrite(ligands,file.path(RESULT,"peak_ligand_predictions.csv"))
native <- readRDS(file.path(CONFIG$root,CONFIG$workspace,"GSE143437_muscle_regeneration","linealinkr_complete.rds"))
meta <- fread(file.path(native$input_extract_dir,"GSE143437_muscle_regeneration_meta.csv"))
stopifnot(identical(meta$barcode_id,native$cells))
niche_cells <- readLines(file.path(niche_fit$config$input,"cells.txt"),encoding="UTF-8")
stopifnot(setequal(native$cells,niche_cells))
types <- unique(data.table(target=native$celltypes,target_label=native$celltype_labels))
query_groups <- merge(types,receiver_map,by="target_label")[,.(context,target)]
query_groups <- rbindlist(list(query_groups,copy(query_groups)[,context:="All MNP"]))
sample_cells <- data.table(target=native$celltypes,sample=as.character(meta[[native$fork_sample_column]]))
sample_counts <- merge(sample_cells,query_groups,by="target",allow.cartesian=TRUE)[,
 .(n_cells=.N),by=.(context,sample,target)]
totals <- sample_counts[,.(total=sum(n_cells)),by=.(context,sample)]
planned_counts <- merge(totals,query_groups,by="context",allow.cartesian=TRUE)
planned_counts <- merge(planned_counts,sample_counts,by=c("context","sample","target"),all.x=TRUE)
nomination_groups <- planned_counts[,.(weight=mean(fcoalesce(n_cells,0L)/total)),by=.(context,target)]
nomination_panel <- merge(frozen,background,by=c("context","ligand"))
current_evidence <- methods$ComposeForkEvidence(cand,nomination_groups=nomination_groups,
 pair_gene_map=unique(cand[,.(pair_idx,ligand_genes,receptor_genes)]),
 nomination_candidates=nomination_panel)
legacy_evidence <- methods$ComposeForkEvidence(previous,nomination_groups=nomination_groups,
 pair_gene_map=unique(previous[,.(pair_idx,ligand_genes,receptor_genes)]),
 nomination_candidates=nomination_panel)
current_nomination <- as.data.table(attr(current_evidence,"ligand_nominations"))
legacy_nomination <- as.data.table(attr(legacy_evidence,"ligand_nominations"))
fwrite(current_nomination,file.path(RESULT,"population_ligand_nominations.csv"))
fwrite(legacy_nomination,file.path(RESULT,"population_legacy_nominations.csv"))
fwrite(as.data.table(attr(current_evidence,"nomination_receiver_evidence")),
 file.path(RESULT,"population_receiver_evidence.csv"))
fwrite(nomination_groups,file.path(RESULT,"population_receiver_weights.csv"))
current_nomination[,method:=c(event_ranking_score="llr_full",no_fork_ranking_score="llr_nofork",
 sample_unadjusted_ranking_score="llr_unadjusted")[method]]
current_scores <- dcast(current_nomination,context+ligand~method,value.var="nomination_order_score")
legacy_scores <- legacy_nomination[method=="event_ranking_score",
 .(context,ligand,llr_previous=nomination_order_score)]
ligands[,c("llr_full","llr_nofork","llr_unadjusted","llr_previous"):=NULL]
ligands <- merge(merge(ligands,current_scores,by=c("context","ligand"),all.x=TRUE),
 legacy_scores,by=c("context","ligand"),all.x=TRUE)
fwrite(ligands,file.path(RESULT,"common_ligand_predictions.csv"))
ranks <- melt(ligands,id.vars=c("context","ligand","n_common_pairs"),
  measure.vars=c("llr_previous","llr_unadjusted","llr_full","llr_nofork","timetalk","timetalk_scc","nichenet"),
  variable.name="method",value.name="score",variable.factor=FALSE)
ranks[,rank_score:=fcoalesce(score,-Inf)]
ranks[,n_available:=sum(is.finite(score)),by=.(context,method)]
ranks[, rank_lo:=frank(-rank_score,ties.method="min"),by=.(context,method)]
ranks[, n_ties:=.N,by=.(context,method,rank_score)]
ranks[, n_candidates:=.N,by=.(context,method)]
ranks[, join_key:=1L]
budgets <- data.table(k=CONFIG$k,join_key=1L)
ranks <- merge(ranks,budgets,by="join_key",allow.cartesian=TRUE)
ranks[, weight:=fifelse(is.finite(score),pmin(1,pmax(0,(k-rank_lo+1)/n_ties)),0)]
stopifnot(all(abs(ranks[,sum(weight)-pmin(k[1L],n_available[1L]),by=.(context,method,k)]$V1)<CONFIG$tolerance))
readout <- merge(ranks,response[,.(context,ligand,assignment,response)],by=c("context","ligand"),allow.cartesian=TRUE)
stats <- readout[,.(n_candidates=unique(n_candidates),mean_response_top=sum(weight*response)/unique(k),
  mean_response_rest=sum((1-weight)*response)/(unique(n_candidates)-unique(k))),
  by=.(context,method,k,assignment)]
availability <- unique(ranks[,.(context,method,k,n_available)])
stats <- merge(stats,availability,by=c("context","method","k"))
stats[,enrichment:=fifelse(n_available>=k,mean_response_top-mean_response_rest,NA_real_)]
stats[n_available<k,c("mean_response_top","mean_response_rest"):=.(NA_real_,NA_real_)]
summ <- stats[,.(n_candidates=unique(n_candidates),n_available=unique(n_available),mean_response_top=mean_response_top[assignment==1L],
  mean_response_rest=mean_response_rest[assignment==1L],enrichment=enrichment[assignment==1L],
  exact_sample_p=sum(enrichment>=enrichment[assignment==1L]-CONFIG$tolerance)/.N,
  n_biological_samples=4L,n_label_assignments=.N),by=.(context,method,k)]
fwrite(summ,file.path(OUT,"supp_muscle_depletion_validation.csv"))
fwrite(stats,file.path(RESULT,"program_enrichment_assignments.csv"))
rank_response <- merge(unique(ranks[,.(context,ligand,method,rank_score)]),
  response[assignment==1L,.(context,ligand,response)],by=c("context","ligand"))
rank_agreement <- rank_response[,.(n_candidates=.N,n_available=sum(is.finite(rank_score)),
  spearman_response=cor(rank_score,response,method="spearman")),
  by=.(context,method)]
fwrite(rank_agreement,file.path(OUT,"supp_muscle_depletion_rank_agreement.csv"))
fwrite(ranks[weight>0],file.path(OUT,"supp_muscle_depletion_top_ligands.csv"))

wide <- dcast(stats,context+k+assignment~method,value.var="enrichment")
diffs <- rbindlist(list(
  wide[,.(context,k,assignment,comparison="Full minus legacy event baseline",difference=llr_full-llr_previous)],
  wide[,.(context,k,assignment,comparison="Full minus no sample adjustment",difference=llr_full-llr_unadjusted)],
  wide[,.(context,k,assignment,comparison="No sample adjustment minus legacy event baseline",difference=llr_unadjusted-llr_previous)],
  wide[,.(context,k,assignment,comparison="Full minus no fork",difference=llr_full-llr_nofork)],
  wide[,.(context,k,assignment,comparison="Full minus TimeTalk",difference=llr_full-timetalk)],
  wide[,.(context,k,assignment,comparison="Full minus NicheNet",difference=llr_full-nichenet)],
  wide[,.(context,k,assignment,comparison="No fork minus TimeTalk",difference=llr_nofork-timetalk)]))
direct <- diffs[,.(difference=difference[assignment==1L],
  exact_two_sided_p=sum(abs(difference)>=abs(difference[assignment==1L])-CONFIG$tolerance)/.N,
  n_biological_samples=4L,n_label_assignments=.N),by=.(context,k,comparison)]
fwrite(direct,file.path(OUT,"supp_muscle_depletion_direct_differences.csv"))
fwrite(diffs,file.path(RESULT,"method_difference_assignments.csv"))

# ENTRAIN produces a global ligand ranking. Compare it only on ligands, on a
# matched ligand panel for the pooled receiver, retaining its original clock.
ent <- fread(file.path(CONFIG$root,"results/by_algorithm/ENTRAIN/GSE143437_muscle_regeneration/ligand_importance.csv"))
four <- merge(ligands[context=="All MNP"],ent[,.(ligand,entrain=importance)],by="ligand")
# Recompute the same nomination fronts on the fixed ENTRAIN shared panel.
panel75 <- nomination_panel[context=="All MNP" & ligand %chin% four$ligand]
current75 <- methods$ComposeForkEvidence(cand,
 nomination_groups=nomination_groups[context=="All MNP"],
 pair_gene_map=unique(cand[,.(pair_idx,ligand_genes,receptor_genes)]),
 nomination_candidates=panel75)
legacy75 <- methods$ComposeForkEvidence(previous,
 nomination_groups=nomination_groups[context=="All MNP"],
 pair_gene_map=unique(previous[,.(pair_idx,ligand_genes,receptor_genes)]),
 nomination_candidates=panel75)
nom75 <- as.data.table(attr(current75,"ligand_nominations"))
nom75[,method:=c(event_ranking_score="llr_full",no_fork_ranking_score="llr_nofork",
 sample_unadjusted_ranking_score="llr_unadjusted")[method]]
ref75 <- as.data.table(attr(legacy75,"ligand_nominations"))[method=="event_ranking_score",
 .(context,ligand,llr_previous=nomination_order_score)]
four[,c("llr_full","llr_nofork","llr_unadjusted","llr_previous"):=NULL]
four <- merge(merge(four,dcast(nom75,context+ligand~method,value.var="nomination_order_score"),
 by=c("context","ligand"),all.x=TRUE),ref75,by=c("context","ligand"),all.x=TRUE)
fwrite(four,file.path(RESULT,"entrain_common_ligand_predictions.csv"))
four <- melt(four,id.vars=c("ligand","context"),measure.vars=c("llr_previous","llr_unadjusted","llr_full","llr_nofork","timetalk","entrain","nichenet"),
  variable.name="method",value.name="score",variable.factor=FALSE)
four[,rank_score:=fcoalesce(score,-Inf)]
four[,n_available:=sum(is.finite(score)),by=method]
four[,rank_lo:=frank(-rank_score,ties.method="min"),by=method]
four[,n_ties:=.N,by=.(method,rank_score)]
four[,n_candidates:=.N,by=method]
four[,join_key:=1L]
four <- merge(four,budgets,by="join_key",allow.cartesian=TRUE)
four[,weight:=fifelse(is.finite(score),pmin(1,pmax(0,(k-rank_lo+1)/n_ties)),0)]
four <- merge(four,response[context=="All MNP",.(ligand,assignment,response)],by="ligand",allow.cartesian=TRUE)
four_stats <- four[,.(n_candidates=unique(n_candidates),mean_response_top=sum(weight*response)/unique(k),
  mean_response_rest=sum((1-weight)*response)/(unique(n_candidates)-unique(k))),by=.(method,k,assignment)]
four_availability <- unique(four[,.(method,k,n_available)])
four_stats <- merge(four_stats,four_availability,by=c("method","k"))
four_stats[,enrichment:=fifelse(n_available>=k,mean_response_top-mean_response_rest,NA_real_)]
four_stats[n_available<k,c("mean_response_top","mean_response_rest"):=.(NA_real_,NA_real_)]
four_summ <- four_stats[,.(n_candidates=unique(n_candidates),n_available=unique(n_available),mean_response_top=mean_response_top[assignment==1L],
  enrichment=enrichment[assignment==1L],exact_sample_p=sum(enrichment>=enrichment[assignment==1L]-CONFIG$tolerance)/.N,
  n_biological_samples=4L,n_label_assignments=.N),by=.(method,k)]
fwrite(four_summ,file.path(OUT,"supp_muscle_depletion_ligand_level.csv"))
fwrite(four_stats,file.path(RESULT,"entrain_common_assignments.csv"))
four_wide <- dcast(four_stats,k+assignment~method,value.var="enrichment")
four_diffs <- rbindlist(list(
 four_wide[,.(k,assignment,comparison="Full minus legacy event baseline",difference=llr_full-llr_previous)],
 four_wide[,.(k,assignment,comparison="Full minus no sample adjustment",difference=llr_full-llr_unadjusted)],
 four_wide[,.(k,assignment,comparison="Full minus ENTRAIN",difference=llr_full-entrain)],
 four_wide[,.(k,assignment,comparison="Full minus NicheNet",difference=llr_full-nichenet)],
 four_wide[,.(k,assignment,comparison="Full minus TimeTalk",difference=llr_full-timetalk)],
 four_wide[,.(k,assignment,comparison="Full minus no fork",difference=llr_full-llr_nofork)]))
four_direct <- four_diffs[,.(difference=difference[assignment==1L],
 exact_two_sided_p=sum(abs(difference)>=abs(difference[assignment==1L])-CONFIG$tolerance)/.N,
 n_biological_samples=4L,n_label_assignments=.N),by=.(k,comparison)]
fwrite(four_direct,file.path(OUT,"supp_muscle_depletion_matched_method_differences.csv"))
main_table <- dcast(summ[method %chin% c("llr_previous","llr_unadjusted","llr_full","llr_nofork","timetalk","nichenet")],
  context+k+n_candidates~method,value.var="enrichment")
fwrite(main_table,file.path(CONFIG$root,CONFIG$main_table))


print(coverage)
print(summ[method!="timetalk_scc"])
print(direct)
print(four_summ)
