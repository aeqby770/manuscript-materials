CONFIG<-list(root=normalizePath(Sys.getenv("LLR_WORK_ROOT",getwd()),winslash="/",mustWork=TRUE),threads=parallel::detectCores(),
 output="results/evaluation/final065/key_ablation",k=c(10L,20L),tolerance=1e-12,
 models=c("full","legacy_matched"),workspaces=c("final065","final065_ablation/legacy_matched"),
 query="Same current point scores, fixed molecular candidates and sample-balanced receiver weights; compare peak versus population-front nomination",
 units="Four biological samples and all six balanced assignments",
 response="results/validation/gse291091/ligand_target_response_assignments.csv")
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads)
library(data.table)
setDTthreads(CONFIG$threads)
VALIDATION<-file.path(CONFIG$root,"results/validation")
CONFIG$refit_root<-Sys.getenv("LLR_ABLATION_REFIT_ROOT","")
switch(as.character(nzchar(CONFIG$refit_root)),"TRUE"={
 CONFIG$output<-file.path("results/evaluation",CONFIG$refit_root,"key_ablation")
 VALIDATION<-file.path(CONFIG$root,"results/evaluation",CONFIG$refit_root,"validation")
 CONFIG$workspaces<-CONFIG$models
})
OUT<-file.path(CONFIG$root,CONFIG$output)
dir.create(OUT,recursive=TRUE,showWarnings=FALSE)
paths<-file.path(VALIDATION,CONFIG$workspaces,"predictions")
peak<-rbindlist(setNames(lapply(file.path(paths,"peak_ligand_predictions.csv"),fread),CONFIG$models),idcol="model")
front<-rbindlist(setNames(lapply(file.path(paths,"common_ligand_predictions.csv"),fread),CONFIG$models),idcol="model")
front75<-rbindlist(setNames(lapply(file.path(paths,"entrain_common_ligand_predictions.csv"),fread),CONFIG$models),idcol="model")
peak75<-merge(peak[context=="All MNP"],unique(front75[,.(model,context,ligand)]),by=c("model","context","ligand"))
peak[,rule:="peak"]
peak75[,rule:="peak"]
front[,rule:="front"]
front75[,rule:="front"]
scores<-rbindlist(list(copy(peak)[,panel:="Full fixed panel"],
 copy(front)[,panel:="Full fixed panel"],copy(peak75)[,panel:="ENTRAIN75"],
 copy(front75)[,panel:="ENTRAIN75"]),use.names=TRUE,fill=TRUE)
ranks<-melt(scores,id.vars=c("model","rule","panel","context","ligand"),
 measure.vars=c("llr_full","llr_nofork"),variable.name="method",value.name="score",variable.factor=FALSE)
ranks[,rank_score:=fcoalesce(score,-Inf)]
ranks[,rank:=frank(-rank_score,ties.method="min"),by=.(model,rule,panel,context,method)]
ranks[,n_ties:=.N,by=.(model,rule,panel,context,method,rank_score)]
ranks[,n_candidates:=.N,by=.(model,rule,panel,context,method)]
ranks[,n_available:=sum(is.finite(score)),by=.(model,rule,panel,context,method)]
ranks<-merge(ranks[,join_key:=1L],data.table(k=CONFIG$k,join_key=1L),by="join_key",allow.cartesian=TRUE)[,join_key:=NULL]
ranks[,weight:=fifelse(is.finite(score),pmin(1,pmax(0,(k-rank+1)/n_ties)),0)]
stopifnot(all(abs(ranks[,sum(weight)-pmin(k[1L],n_available[1L]),
 by=.(model,rule,panel,context,method,k)]$V1)<CONFIG$tolerance))
response<-melt(fread(file.path(CONFIG$root,CONFIG$response)),id.vars="ligand",
 variable.name="context_assignment",value.name="response")
response[,c("context","assignment"):=tstrsplit(context_assignment,"|",fixed=TRUE)]
response[,assignment:=as.integer(assignment)]
readout<-merge(ranks,response,by=c("context","ligand"),allow.cartesian=TRUE)
stats<-readout[,.(n_candidates=unique(n_candidates),n_available=unique(n_available),
 enrichment=fifelse(unique(n_available)>=unique(k),
 sum(weight*response)/unique(k)-sum((1-weight)*response)/(unique(n_candidates)-unique(k)),NA_real_)),
 by=.(model,rule,panel,context,method,k,assignment)]
summary<-stats[,.(n_candidates=unique(n_candidates),n_available=unique(n_available),
 enrichment=enrichment[assignment==1L],n_label_assignments=.N,n_biological_samples=4L),
 by=.(model,rule,panel,context,method,k)]
wide<-dcast(stats,panel+context+method+k+assignment~model+rule,value.var="enrichment")
differences<-rbindlist(list(
 wide[,.(panel,context,method,k,assignment,comparison="Full: front minus peak",difference=full_front-full_peak)],
 wide[,.(panel,context,method,k,assignment,comparison="Matched old: front minus peak",difference=legacy_matched_front-legacy_matched_peak)],
 wide[,.(panel,context,method,k,assignment,comparison="Full front minus matched old front",difference=full_front-legacy_matched_front)]))
direct<-differences[,.(difference=difference[assignment==1L],
 exact_two_sided_p=mean(abs(difference)>=abs(difference[assignment==1L])-CONFIG$tolerance),
 n_biological_samples=4L,n_label_assignments=.N),by=.(panel,context,method,k,comparison)]
expected<-rbindlist(setNames(lapply(file.path(VALIDATION,CONFIG$workspaces,
 "tables/supp_muscle_depletion_validation.csv"),fread),CONFIG$models),idcol="model")
check<-merge(summary[rule=="front"&panel=="Full fixed panel"],expected,
 by=c("model","context","method","k"),suffixes=c("_replay","_saved"))
supported<-is.finite(check$enrichment_replay)
stopifnot(nrow(check)==24L,identical(supported,is.finite(check$enrichment_saved)),
  all(check$n_available_replay==check$n_available_saved),
  all(abs(check$enrichment_replay[supported]-check$enrichment_saved[supported])<CONFIG$tolerance))
fwrite(summary,file.path(OUT,"receiver_nomination_ablation_summary.csv"))
fwrite(ranks,file.path(OUT,"receiver_nomination_all_scores_ranks.csv"))
fwrite(stats,file.path(OUT,"receiver_nomination_all_assignments.csv"))
fwrite(direct,file.path(OUT,"receiver_nomination_direct_differences.csv"))
print(summary[context=="All MNP"&method=="llr_full"])
print(direct[context=="All MNP"&method=="llr_full"])
