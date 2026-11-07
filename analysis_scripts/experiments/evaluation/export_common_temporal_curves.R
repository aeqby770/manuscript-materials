# Conditional pulse localisation from saved native activities on one common layout.
# True times and the generating plateau enter evaluation only.
TIME_CONFIG <- list(min_pct=0.05,min_bins=3L,plateau_lo=0.62,plateau_hi=0.68,
  reader="fixed_truth_mapped_context; saved_G_memberships_and_Q_bin_layout")
time_detail <- fread(file.path(base,"lrsim_contrast_detail.csv"))[
  dataset==CONFIG$dataset & sender_direction%in%TRUE & truth_row%in%TRUE &
    role=="transient" & duplicate_emission%in%c(FALSE,NA)]
time_targets <- unique(time_detail[,.(target_key=hit_target_key,unit_id,child,pair_idx)])
setorder(time_targets,target_key,unit_id,child,pair_idx)
time_targets <- unique(time_targets,by="target_key")
time_units <- attr(x$fork_contrasts,"event_units")$units
time_unit_table <- data.table(unit_id=seq_along(time_units),unit=time_units)
time_children <- time_unit_table[,.(child=unit[[1L]]$children,
  descendants=unit[[1L]]$desc),by=unit_id]
time_leaves <- time_children[,.(leaf=unlist(descendants)),by=.(unit_id,child)]
time_bounds <- x$temporal_lr$intermediate$bin_bounds[[1L]]
time_n_bins <- length(time_bounds)-1L
time_windows <- time_children[,.(bin=seq_len(time_n_bins)),by=.(unit_id,child)]
time_windows[,c("window_id","mid"):=list(.I,
  (time_bounds[bin]+time_bounds[bin+1L])/2)]
time_meta <- fread(file.path(CONFIG$root,"cache/final065",CONFIG$dataset,"extract",paste0(CONFIG$dataset,"_meta.csv")))[
  match(x$cells,barcode_id)]
stopifnot(identical(time_meta$barcode_id,x$cells))
time_inputs <- fread(file.path(CONFIG$root,"results/by_algorithm/TraSig",
  CONFIG$dataset,"trasig_cell_inputs.csv"))[match(x$cells,cell_id)]
stopifnot(identical(time_inputs$cell_id,x$cells))
time_cells <- data.table(cell_index=seq_along(x$cells),cell_type=x$celltypes,
  bin=findInterval(x$pseudotime,time_bounds,all.inside=TRUE),
  true_time=time_meta$true_pseudotime,native_path=time_inputs$cell_path,
  native_time=round(time_inputs$cell_time,2))
time_weights <- as.data.table(x$lineage_reference_weights/rowSums(x$lineage_reference_weights))
setnames(time_weights,as.character(seq_len(ncol(time_weights))))
time_weights[,cell_index:=.I]
time_weights <- melt(time_weights,id.vars="cell_index",variable.name="leaf",value.name="weight")
time_weights[,leaf:=as.integer(as.character(leaf))]
time_weights <- merge(time_weights[weight>0],time_leaves,by="leaf",allow.cartesian=TRUE)
time_weights <- time_weights[,.(weight=sum(weight)),by=.(unit_id,child,cell_index)][
  weight>=CONFIG$min_weight]
time_weights <- merge(merge(time_weights,time_cells,by="cell_index"),time_windows,
  by=c("unit_id","child","bin"))
time_support <- time_weights[,.(mass=sum(weight),n_eff=sum(weight)^2/sum(weight^2)),
  by=.(window_id,cell_type)]
time_calibration <- time_weights[cell_type==target_type,
  .(true_center=sum(true_time*weight)/sum(weight)),by=window_id]
time_weights[,normalized_weight:=weight/sum(weight),by=.(window_id,cell_type)]
time_sender <- time_weights[cell_type==source_type]
time_receiver <- time_weights[cell_type==target_type]
time_source_matrix <- sparseMatrix(i=time_sender$cell_index,j=time_sender$window_id,
  x=time_sender$normalized_weight,dims=c(n_cells,nrow(time_windows)))
time_target_matrix <- sparseMatrix(i=time_receiver$cell_index,j=time_receiver$window_id,
  x=time_receiver$normalized_weight,dims=c(n_cells,nrow(time_windows)))
time_source_matrix <- as.matrix(time_source_matrix)
time_target_matrix <- as.matrix(time_target_matrix)
time_llr <- sqrt((ligand%*%time_source_matrix)*(receptor%*%time_target_matrix))
time_llr[(ligand>0)%*%time_source_matrix<TIME_CONFIG$min_pct |
  (receptor>0)%*%time_target_matrix<TIME_CONFIG$min_pct] <- 0
# The pre-Hill geometric expression has exactly the native Hill peak location.
time_kernel_weights <- time_source_matrix[source_index,,drop=FALSE]*
  time_target_matrix[target_index,,drop=FALSE]
time_sccomm <- (kernel%*%time_kernel_weights)*native_weight
time_score_info <- data.table(method=rep(c("LineaLinkR_G_activity","scComm_native_kernel"),
  each=nrow(pairs)),pair_key=rep(pair_keys,2L))
time_matrix <- rbind(time_llr,time_sccomm)
colnames(time_matrix) <- as.character(time_windows$window_id)
time_scores <- melt(cbind(time_score_info,as.data.table(time_matrix)),
  id.vars=c("method","pair_key"),variable.name="window_id",value.name="score")
time_scores[,window_id:=as.integer(as.character(window_id))]
time_scores <- merge(time_scores,pairs[,.(pair_key,pair_idx)],by="pair_key")
time_scores[,pair_key:=NULL]
time_sender_mix <- time_sender[,.(source_mass=sum(weight)),
  by=.(window_id,source_path=native_path)]
time_sender_mix[,source_fraction:=source_mass/sum(source_mass),by=window_id]
time_receiver_mix <- time_receiver[,.(receiver_mass=sum(weight)),
  by=.(window_id,target_path=native_path,time=native_time)]
time_mixture <- merge(time_sender_mix,time_receiver_mix,by="window_id",allow.cartesian=TRUE)
time_mixture[,path_pair:=paste(source_path,target_path,sep="_")]
time_curves <- fread(file.path(CONFIG$root,"results/by_algorithm/TraSig",CONFIG$dataset,"trasig_temporal.csv"))
time_curves <- merge(time_curves,pairs,by=c("ligand","receptor"))
time_values <- time_curves[,.(path_pair,time,pair_idx,score,
  source_cells_in_window,target_cells_in_window)]
time_projection <- merge(time_mixture,time_values,by=c("path_pair","time"),
  all.x=TRUE,allow.cartesian=TRUE)
time_projection[,available:=is.finite(score)&source_cells_in_window>0&target_cells_in_window>0]
time_trasig <- time_projection[,.(weighted_score=sum(fifelse(available,score,0)*
  source_fraction*receiver_mass),supported_mass=sum(fifelse(available,
  source_fraction*receiver_mass,0))),by=.(window_id,pair_idx)]
time_receiver_mass <- time_receiver[,.(receiver_mass=sum(weight)),by=window_id]
time_trasig <- merge(time_trasig,time_receiver_mass,by="window_id")
time_trasig[,score:=weighted_score/receiver_mass]
time_trasig[supported_mass/receiver_mass<1-CONFIG$tolerance,score:=NA_real_]
time_scores <- rbindlist(list(time_scores,time_trasig[,.(window_id,pair_idx,
  method="TraSig_native_curve",score)]),use.names=TRUE)
time_window_info <- merge(time_windows,time_calibration,by="window_id",all.x=TRUE)
time_window_info <- merge(time_window_info,time_support[cell_type==source_type,
  .(window_id,source_n_eff=n_eff)],by="window_id",all.x=TRUE)
time_window_info <- merge(time_window_info,time_support[cell_type==target_type,
  .(window_id,target_n_eff=n_eff)],by="window_id",all.x=TRUE)
time_window_info[,supported:=is.finite(source_n_eff)&is.finite(target_n_eff)&
  pmin(source_n_eff,target_n_eff)>=CONFIG$min_n_eff]
time_scores <- merge(time_scores,time_window_info,by="window_id")
time_scores[supported%in%FALSE,score:=NA_real_]
time_scores <- merge(time_scores,pairs[,.(pair_idx,ligand,receptor)],by="pair_idx")
time_scores[,dataset:=CONFIG$dataset]
fwrite(time_scores,file.path(out,"common_temporal_activity_profiles.csv"))
time_peaks <- time_scores[is.finite(score)&is.finite(true_center),.(
  n_bins=.N,span=max(score)-min(score),peak_bin=bin[which.max(score)],
  peak_clock=mid[which.max(score)],peak_true= true_center[which.max(score)]),
  by=.(method,unit_id,child,pair_idx)]
time_peaks[n_bins<TIME_CONFIG$min_bins | span<=CONFIG$tolerance,
  c("peak_bin","peak_clock","peak_true"):=list(NA_integer_,NA_real_,NA_real_)]
time_native <- unique(time_detail[,.(target_key=hit_target_key,unit_id,child,pair_idx,
  native_peak=focal_peak_time)])
time_native <- merge(time_targets,time_native,by=c("target_key","unit_id","child","pair_idx"))
time_native_cal <- merge(time_native,time_window_info[supported%in%TRUE & is.finite(true_center)],
  by=c("unit_id","child"),allow.cartesian=TRUE)
time_native_cal <- time_native_cal[,.SD[.N>=TIME_CONFIG$min_bins],by=target_key]
time_native_peak <- time_native_cal[,.(method="LineaLinkR_native_pulse_fit",
  n_bins=.N,span=NA_real_,peak_bin=NA_integer_,peak_clock=native_peak[1L],
  peak_true=approx(mid,true_center,xout=native_peak[1L])$y),
  by=.(unit_id,child,pair_idx)]
time_peaks <- rbindlist(list(time_peaks,time_native_peak),use.names=TRUE)
time_methods <- data.table(method=c("LineaLinkR_G_activity","scComm_native_kernel",
  "TraSig_native_curve","LineaLinkR_native_pulse_fit"),join_key=1L)
time_planned <- merge(copy(time_targets)[,join_key:=1L],time_methods,
  by="join_key",allow.cartesian=TRUE)[,join_key:=NULL]
time_recovery <- merge(time_planned,time_peaks,
  by=c("method","unit_id","child","pair_idx"),all.x=TRUE)
time_recovery[,c("dataset","plateau_lo","plateau_hi","reader"):=list(CONFIG$dataset,
  TIME_CONFIG$plateau_lo,TIME_CONFIG$plateau_hi,TIME_CONFIG$reader)]
time_recovery[,available:=is.finite(peak_true)]
time_recovery[,distance_to_plateau:=pmax(plateau_lo-peak_true,peak_true-plateau_hi,0)]
time_recovery[,hit_plateau:=available&distance_to_plateau==0]
time_joint <- time_recovery[method!="LineaLinkR_native_pulse_fit",
  .(common_profile_support=all(available)),by=target_key]
time_recovery <- merge(time_recovery,time_joint,by="target_key")
time_all_readers <- time_recovery[,.(common_all_reader_support=all(available)),by=target_key]
time_recovery <- merge(time_recovery,time_all_readers,by="target_key")
fwrite(time_recovery,file.path(out,"temporal_peak_window_recovery.csv"))
time_scopes <- rbindlist(list(copy(time_recovery)[,scope:="all_fixed_transient_targets"],
  copy(time_recovery[common_profile_support%in%TRUE])[,scope:="common_profile_support"],
  copy(time_recovery[common_all_reader_support%in%TRUE])[,scope:="common_all_reader_support"]))
time_metrics <- time_scopes[,.(n_targets=.N,n_available=sum(available),n_hit=sum(hit_plateau),
  mean_distance=mean(distance_to_plateau,na.rm=TRUE)),by=.(dataset,method,scope)]
fwrite(time_metrics,file.path(out,"temporal_peak_window_metrics.csv"))
print(time_metrics)
