CONFIG <- list(paths=list(root=normalizePath(commandArgs(TRUE)[1L],winslash="/",mustWork=TRUE),
  native="inputs/comparison/recovered_native",output="results/evaluation/final065/saved_native"),
  threads=parallel::detectCores(),plot=list())
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,MKL_NUM_THREADS=CONFIG$threads)
.libPaths(c(file.path(CONFIG$paths$root,"library"),.libPaths()))
library(data.table)
library(Matrix)
library(SingleCellExperiment)
setDTthreads(CONFIG$threads)
base <- file.path(CONFIG$paths$root,CONFIG$paths$native)
out <- file.path(CONFIG$paths$root,CONFIG$paths$output)
dir.create(out,recursive=TRUE,showWarnings=FALSE)

p <- list.files(file.path(base,"by_algorithm/CellChat"),pattern="^cellchat[.]rds$",recursive=TRUE,full.names=TRUE)
manifest <- data.table(source_file=p,dataset=basename(dirname(p)),method="CellChat")
cc <- manifest[,{
  x <- readRDS(source_file)
  ids <- colnames(x@data)
  obs <- data.table(cell_id=ids,native_label=as.character(x@meta[ids,"labels"]))
  map <- as.data.table(x@LR$LRsig,keep.rownames="interaction")
  parts <- melt(as.data.table(x@DB$complex,keep.rownames="partner"),id.vars="partner",value.name="gene")
  parts <- parts[!is.na(gene)&nzchar(gene),.(components=paste(sort(unique(as.character(gene))),collapse="_")),by=partner]
  map[parts,on=.(ligand=partner),gene_a:=i.components]
  map[parts,on=.(receptor=partner),gene_b:=i.components]
  map[,gene_a:=fcoalesce(gene_a,as.character(ligand))]
  map[,gene_b:=fcoalesce(gene_b,as.character(receptor))]
  scores <- switch(as.character(length(x@net$prob)>0L),"TRUE"={
    d <- as.data.table(as.table(x@net$prob))
    setnames(d,c("sender","receiver","interaction","native_score"))
    d[,pvalue:=as.numeric(x@net$pval)]
    d <- merge(d,map[,.(interaction,gene_a,gene_b)],by="interaction",all.x=TRUE)
    d[,.(source_a=as.character(sender),source_b=as.character(receiver),
      entity_a=gene_a,entity_b=gene_b,native_score,pvalue,
      score_field="probability",direction_scope="native_LR",native_pair=interaction)]
  },"FALSE"=data.table())
  .(cells=list(obs),scores=list(scores),n_cells=length(ids),n_genes=nrow(x@data),
    n_native_rows=nrow(scores),object_status=fifelse(nrow(scores)>0,"saved_scores","saved_object_without_scores"),
    identity_basis="saved_expression_cell_ids",native_parameters=paste(names(x@options),collapse=";"))
},by=.(dataset,method,source_file)]
cc_cells <- cc[,cells[[1L]],by=.(dataset,method,source_file)]
cc_scores <- cc[n_native_rows>0,scores[[1L]],by=.(dataset,method,source_file)]

p <- list.files(file.path(base,"by_algorithm/SingleCellSignalR"),pattern="^scsignalr_scsr[.]rds$",recursive=TRUE,full.names=TRUE)
manifest <- data.table(source_file=p,dataset=basename(dirname(p)),method="SingleCellSignalR_native")
scsr <- manifest[,{
  x <- attributes(readRDS(source_file))
  raw <- attributes(x$bsrdm.comp)
  counts <- raw$ncounts
  obs <- data.table(cell_id=colnames(counts),native_label=as.character(x$populations))
  para <- rbindlist(x$paracrines,idcol="native_pair",fill=TRUE)
  auto <- rbindlist(x$autocrines,idcol="native_pair",fill=TRUE)
  para[,kind:=rep("paracrine",.N)]
  auto[,kind:=rep("autocrine",.N)]
  d <- rbindlist(list(para,auto),fill=TRUE)
  .(cells=list(obs),scores=list(d),n_cells=ncol(counts),n_genes=nrow(counts),
    n_native_rows=nrow(d),object_status=fifelse(nrow(d)>0,"saved_returned_interactions","saved_object_without_returned_interactions"),
    identity_basis="saved_native_counts_cell_ids",
    native_parameters=paste0("n.resample=",x$inf.param$n.resample,";subsample.size=",x$inf.param$subsample.size))
},by=.(dataset,method,source_file)]
scsr_cells <- scsr[,cells[[1L]],by=.(dataset,method,source_file)]
scsr_scores <- rbindlist(setNames(scsr$scores,scsr$source_file),idcol="source_file",use.names=TRUE,fill=TRUE)
scsr_scores <- merge(scsr_scores,scsr[,.(dataset,method,source_file)],by="source_file",all.x=TRUE)

p <- list.files(file.path(base,"by_algorithm/CytoTalk"),pattern="^cytotalk_result[.]rds$",recursive=TRUE,full.names=TRUE)
manifest <- data.table(source_file=p,dataset=basename(dirname(dirname(p))),method="CytoTalk_native")
cyto <- manifest[,{
  x <- readRDS(source_file)
  d <- as.data.table(x$pcst$final_network)
  cross <- d[is_ct_edge%in%TRUE,.(source_a=as.character(node1_type),source_b=as.character(node2_type),
    entity_a=as.character(node1),entity_b=as.character(node2),native_score=-cost,pvalue=NA_real_,
    score_field="negative_native_PCST_edge_cost",direction_scope="native_cross_population_unoriented",
    native_pair=paste(node1,node2,node1_type,node2_type,sep="|"))]
  .(scores=list(cross),n_cells=NA_integer_,n_genes=NA_integer_,n_native_rows=nrow(cross),
    n_graph_edges=nrow(d),object_status="saved_selected_network",
    identity_basis="cell_ids_not_stored_in_network_object",
    native_parameters=paste0("population_a=",x$params$cell_type_a,";population_b=",x$params$cell_type_b,
      ";ntrial=",x$params$ntrial))
},by=.(dataset,method,source_file)]
cyto_scores <- cyto[n_native_rows>0,scores[[1L]],by=.(dataset,method,source_file)]

p <- list.files(file.path(base,"by_algorithm/tradeSeq"),pattern="^tradeseq_sce[.]rds$",recursive=TRUE,full.names=TRUE)
manifest <- data.table(source_file=p,dataset=basename(dirname(p)),method="tradeSeq")
ts <- manifest[,{
  x <- readRDS(source_file)
  d <- fread(file.path(dirname(source_file),"tradeseq_association_test.csv"))
  stopifnot(setequal(rownames(x),as.character(d$gene)),!anyDuplicated(colnames(x)))
  .(cells=list(data.table(cell_id=colnames(x),native_label=NA_character_)),n_cells=ncol(x),n_genes=nrow(x),
    n_native_rows=nrow(d),n_finite_p=sum(is.finite(d$pvalue)),object_status="saved_GAM_and_association_tests",
    identity_basis="saved_SCE_cell_ids",native_parameters="saved_associationTest;no_fit_or_test_recomputed")
},by=.(dataset,method,source_file)]
ts_cells <- ts[,cells[[1L]],by=.(dataset,method,source_file)]

p <- list.files(file.path(base,"by_algorithm/FateID"),pattern="^fateid_bias[.]rds$",recursive=TRUE,full.names=TRUE)
manifest <- data.table(source_file=p,dataset=basename(dirname(p)),method="FateID")
fi <- manifest[,{
  x <- readRDS(source_file)
  d <- fread(file.path(dirname(source_file),"fateid_bias.csv"),colClasses=c(cell_id="character"))
  target_names <- colnames(x$probs)
  difference <- max(abs(as.matrix(d[,target_names,with=FALSE])-as.matrix(x$probs)))
  stopifnot(nrow(d)==nrow(x$probs),difference==0,!anyDuplicated(d$cell_id))
  .(n_cells=nrow(d),n_genes=NA_integer_,n_native_rows=nrow(d),probability_max_difference=difference,
    object_status="RDS_and_exported_probabilities_identical",identity_basis="saved_cell_id_export_and_native_zero_based_rows",
    native_parameters=paste0("n_targets=",ncol(x$probs),";native_row_min=",min(as.integer(rownames(x$probs)))))
},by=.(dataset,method,source_file)]

datasets <- readLines(file.path(CONFIG$paths$root,"experiments/config/ordering24_datasets.txt"),encoding="UTF-8")
manifest <- data.table(dataset=datasets,method="LineaLinkR_0.6.5_Final065",
  source_file=file.path(CONFIG$paths$root,"cache/final065",datasets,"linealinkr_complete.rds"))
llr <- manifest[,{
  x <- readRDS(source_file)
  w <- as.data.table(x$lineage_weights)
  w[,cell_id:=as.character(x$cells)]
  d <- melt(w,id.vars="cell_id",variable.name="native_target",value.name="probability")
  .(cells=list(data.table(cell_id=as.character(x$cells),native_label=as.character(x$celltype_labels))),
    probabilities=list(d),n_cells=x$n_cells,n_genes=x$n_genes,n_native_rows=nrow(d),
    object_status="saved_current_lineage_weights",identity_basis="saved_current_model_cell_ids",
    native_parameters=paste0("revision=",x$result_revision,";root=",x$root_evidence_mode))
},by=.(dataset,method,source_file)]
llr_cells <- llr[,cells[[1L]],by=.(dataset,method,source_file)]
llr_probabilities <- llr[,probabilities[[1L]],by=.(dataset,method,source_file)]
fwrite(llr_probabilities,file.path(out,"current_LineaLinkR_fates.csv.gz"))

cells <- rbindlist(list(cc_cells,scsr_cells,ts_cells,llr_cells),fill=TRUE)
stopifnot(!anyDuplicated(cells[,.(dataset,method,source_file,cell_id)]))
fwrite(cells,file.path(out,"native_object_cells.csv.gz"))
fwrite(rbindlist(list(cc_scores,cyto_scores),fill=TRUE),file.path(out,"native_object_CCC_scores.csv.gz"))
fwrite(scsr_scores,file.path(out,"native_SingleCellSignalR_returned_scores.csv.gz"))
catalog <- rbindlist(list(cc[,!c("cells","scores")],scsr[,!c("cells","scores")],
  cyto[,!"scores"],ts[,!"cells"],fi,llr[,!c("cells","probabilities")]),fill=TRUE)
catalog[,evaluation_scope:="read_saved_native_objects_only"]
fwrite(catalog,file.path(out,"native_R_object_catalog.csv"))
print(catalog[,.(n_objects=.N,n_saved_score_rows=sum(n_native_rows)),by=method])

old_paths <- list.files(file.path(CONFIG$paths$root,"results/by_algorithm/reused"),
  pattern="^LineaLinkR_static_native_scores[.]csv[.]gz$",recursive=TRUE,full.names=TRUE)
new_paths <- list.files(file.path(CONFIG$paths$root,"results/evaluation/final065/static_ccc"),
  pattern="^LineaLinkR_static_native_scores[.]csv[.]gz$",recursive=TRUE,full.names=TRUE)
manifest <- rbindlist(list(
  data.table(dataset=basename(dirname(old_paths)),workspace="feasible_clock")[
    !dataset%in%basename(dirname(new_paths))],
  data.table(dataset=basename(dirname(new_paths)),workspace="final065")))
manifest[,source_file:=file.path(CONFIG$paths$root,"cache",workspace,dataset,"linealinkr_complete.rds")]
reference <- manifest[,{
  x <- readRDS(source_file)
  type_map <- unique(data.table(type_id=as.character(x$celltypes),label=as.character(x$celltype_labels)))
  stopifnot(!anyDuplicated(type_map$type_id),!anyDuplicated(x$cells))
  d <- as.data.table(x$lr_scores)
  d[,c("source","target"):=list(as.character(source),as.character(target))]
  d[type_map,on=.(source=type_id),sender:=i.label]
  d[type_map,on=.(target=type_id),receiver:=i.label]
  d[as.data.table(x$pair_gene_map),on=.(pair_idx),
    c("ligand","receptor"):=list(i.ligand_genes,i.receptor_genes)]
  stopifnot(!anyNA(d$sender),!anyNA(d$receiver),!anyNA(d$ligand),!anyNA(d$receptor))
  .(cells=list(data.table(cell_id=as.character(x$cells),reference_label=as.character(x$celltype_labels))),
    scores=list(d[,.(sender,receiver,ligand,receptor,reference_score=score)]),
    n_cells=x$n_cells,n_genes=x$n_genes,n_saved_score_rows=nrow(d),
    result_revision=x$result_revision)
},by=.(dataset,workspace,source_file)]
reference_cells <- reference[,cells[[1L]],by=.(dataset,workspace,source_file)]
reference_scores <- reference[,scores[[1L]],by=.(dataset,workspace,source_file)]
fwrite(reference_cells,file.path(out,"static_reference_cells.csv.gz"))
fwrite(reference_scores,file.path(out,"static_reference_scores.csv.gz"))
fwrite(reference[,!c("cells","scores")],file.path(out,"static_reference_model_binding.csv"))
print(reference[,.(dataset,workspace,n_cells,n_saved_score_rows)])

# Read saved input identities for gene tests and the remaining native smoke objects.
gene_paths <- list.files(file.path(base,"by_algorithm/tradeSeq"),
  pattern="^tradeseq_sce[.]rds$",recursive=TRUE,full.names=TRUE)
gene_binding <- data.table(dataset=basename(dirname(gene_paths)),source_file=gene_paths)[
  dataset%in%basename(list.dirs(file.path(base,"input_metadata"),recursive=FALSE)),{
  x <- readRDS(source_file)
  input_dir <- file.path(base,"input_metadata",dataset)
  input_cells <- readLines(file.path(input_dir,"cells.txt"),encoding="UTF-8")
  input_genes <- readLines(file.path(input_dir,"genes.txt"),encoding="UTF-8")
  .(method="tradeSeq",n_saved_cells=ncol(x),n_snapshot_cells=length(input_cells),
    n_missing_input_cells=length(setdiff(input_cells,colnames(x))),
    n_unexpected_cells=length(setdiff(colnames(x),input_cells)),
    exact_input_cell_set=setequal(input_cells,colnames(x)),
    n_saved_genes=nrow(x),n_snapshot_genes=length(input_genes),
    n_saved_genes_outside_input=length(setdiff(rownames(x),input_genes)),
    n_input_genes_without_saved_test=length(setdiff(input_genes,rownames(x))),
    identity_basis="saved_SCE_ids_and_original_exported_input_ids",
    matrix_equality_status="native_preprocessed_expression;not_assumed_equal_to_input_matrix")
},by=.(dataset,source_file)]
gene_binding <- merge(
  fread(file.path(out,"native_R_object_catalog.csv"))[method=="tradeSeq",
    .(dataset,source_file,n_cells,n_genes)],gene_binding,
  by=c("dataset","source_file"),all.x=TRUE)
gene_binding[,input_binding_status:=fifelse(is.na(exact_input_cell_set),
  "no_saved_original_input_snapshot","saved_SCE_and_original_input_snapshot_read")]
fwrite(gene_binding,file.path(out,"native_gene_input_binding.csv"))
tempora_source <- file.path(base,"smoke/by_algorithm/Tempora/smoke_human_cd34_500c/tempora_result.rds")
tx <- attributes(readRDS(tempora_source))
stopifnot(identical(colnames(tx$data),rownames(tx$meta.data)))
fwrite(data.table(cell_id=colnames(tx$data),as.data.table(tx$meta.data)),
  file.path(out,"native_Tempora_cells.csv.gz"))
sx_source <- file.path(base,"own_saved_references/cellrank_stemnet/linealinkr_complete.rds")
sx <- readRDS(sx_source)
additional <- rbindlist(list(
  data.table(dataset="smoke_human_cd34_500c",method="Tempora",source_file=tempora_source,
    n_cells=ncol(tx$data),n_genes=nrow(tx$data),n_native_rows=nrow(tx$trajectory),
    n_populations=uniqueN(tx$meta.data$Clusters),n_timepoints=uniqueN(tx$meta.data$Timepoints),
    n_pathways=nrow(tx$cluster.pathways),object_status="saved_cluster_pathway_trajectory",
    accuracy_status="no_independent_GT_graph_or_cell_pseudotime;smoke_timepoint_provenance_unverified"),
  data.table(dataset="cellrank_stemnet",method="LineaLinkR_0.6.4_saved",
    source_file=sx_source,n_cells=sx$n_cells,n_genes=sx$n_genes,
    n_native_rows=nrow(sx$lr_scores),n_populations=sx$n_types,
    object_status="saved_static_LR_scores_and_trajectories",
    accuracy_status="historical_model_has_no_full_saved_celltype_id_label_map",
    result_revision=sx$result_revision)),fill=TRUE)
additional<-additional[method!="LineaLinkR_0.6.4_saved"]
fwrite(additional,file.path(out,"native_additional_object_readback.csv"))
print(gene_binding[,.(n_datasets=.N,n_exact_cell_sets=sum(exact_input_cell_set),
  n_saved_genes_outside_input=sum(n_saved_genes_outside_input))])
print(additional[,.(dataset,method,n_cells,n_native_rows,accuracy_status)])



# Read only CytoTalk original population-pair plans and saved networks.
cyto_lung_meta <- fread(file.path(base,"input_metadata/lung/meta.csv"))
cyto_paul_meta <- fread(file.path(base,"input_metadata/paul15/meta.csv"))
cyto_labels <- rbindlist(list(
 cyto_lung_meta[,.(dataset="lung",population=fifelse(is.na(cell_type),"unknown",as.character(cell_type)))],
 cyto_paul_meta[,.(dataset="paul15",population=fifelse(is.na(paul15_clusters),"unknown",as.character(paul15_clusters)))]))
cyto_counts <- cyto_labels[,.N,by=.(dataset,population)]
cyto_counts[,native_label:=make.names(population)]
cyto_pairs <- merge(cyto_counts,cyto_counts,by="dataset",suffixes=c("_a","_b"),allow.cartesian=TRUE)
cyto_pairs <- cyto_pairs[population_a<population_b]
cyto_pairs[,pair_key:=paste(pmin(native_label_a,native_label_b),pmax(native_label_a,native_label_b),sep="__")]
cyto_saved <- data.table(network_file=list.files(file.path(base,"by_algorithm/CytoTalk"),
 pattern="^FinalNetwork[.]txt$",recursive=TRUE,full.names=TRUE))
cyto_saved[,c("dataset","directory"):=list(basename(dirname(dirname(network_file))),basename(dirname(network_file)))]
cyto_saved[,c("native_label_a","native_label_b"):=tstrsplit(directory,"__",fixed=TRUE)]
cyto_saved[,pair_key:=paste(pmin(native_label_a,native_label_b),pmax(native_label_a,native_label_b),sep="__")]
cyto_saved[,saved_pair_rds:=file.exists(file.path(dirname(network_file),"cytotalk_result.rds"))]
cyto_pairs <- merge(cyto_pairs,cyto_saved[,.(dataset,pair_key,network_file,saved_pair_rds)],
 by=c("dataset","pair_key"),all.x=TRUE,sort=FALSE)
cyto_pairs[,c("planned_pair","saved_network","scope_status"):=list(TRUE,!is.na(network_file),
 "original_metadata_all_unordered_population_pairs")]
cyto_pairs[,input_cells:=cyto_labels[,.N,by=dataset]$N[match(dataset,cyto_labels[,.N,by=dataset]$dataset)]]
cyto_smoke <- cyto_saved[dataset=="paul15_smoke"]
cyto_smoke[,c("planned_pair","saved_network","input_cells","scope_status"):=list(
 NA,TRUE,NA_integer_,"saved_smoke_networks;original_input_population_plan_unverified")]
cyto_scope <- rbindlist(list(cyto_pairs,cyto_smoke),use.names=TRUE,fill=TRUE)
cyto_core <- cyto_saved[,{
 folder<-dirname(network_file)
 rds_path<-file.path(folder,"cytotalk_result.rds")
 node_path<-file.path(folder,"IntegratedNodes.txt")
 edge_path<-file.path(folder,"IntegratedEdges.txt")
 kind<-fcase(all(file.exists(c(node_path,edge_path))),"native_tables",
  file.exists(rds_path),"native_rds",default="unavailable")
 net<-switch(kind,
  native_tables=list(nodes=fread(node_path),edges=fread(edge_path)),
  native_rds=readRDS(rds_path)$integrated_net,
  unavailable=NULL)
 switch(as.character(!is.null(net)),
  "TRUE"=list(core_format=kind,n_core_gene_population_vertices=nrow(net$nodes),
   n_core_edges=nrow(net$edges),n_nonfinite_prizes=sum(!is.finite(net$nodes$prize)),
   n_nonfinite_costs=sum(!is.finite(net$edges$cost)),
   model_inputs_finite=all(is.finite(net$nodes$prize))&&all(is.finite(net$edges$cost))),
  "FALSE"=list(core_format=kind,n_core_gene_population_vertices=NA_integer_,n_core_edges=NA_integer_,
   n_nonfinite_prizes=NA_integer_,n_nonfinite_costs=NA_integer_,model_inputs_finite=NA))
},by=.(dataset,pair_key,network_file)]
cyto_scope<-merge(cyto_scope,cyto_core,by=c("dataset","pair_key","network_file"),all.x=TRUE,sort=FALSE)
cyto_failures<-fread(file.path(CONFIG$paths$root,"results/by_algorithm/CytoTalk/paul15/cytotalk_pair_failures.csv"))
cyto_scope<-merge(cyto_scope,cyto_failures[,.(dataset,pair_key,shared_failure_status=result_status,shared_failure_reason=failure_reason)],
 by=c("dataset","pair_key"),all.x=TRUE,sort=FALSE)
cyto_scope[,valid_native_network:=saved_network&model_inputs_finite]
cyto_scope[,native_model_status:=fcase(
 valid_native_network%in%TRUE,"valid_native_network",
 !is.na(shared_failure_status),shared_failure_status,
 saved_network&model_inputs_finite%in%FALSE,"failed_saved_nonfinite_native_model_inputs",
 saved_network,"saved_native_network_finiteness_unverified",default="missing_native_network")]
setorder(cyto_scope,dataset,pair_key)
fwrite(cyto_scope,file.path(out,"native_CytoTalk_pair_coverage.csv"))
print(cyto_scope[,.(n_planned=sum(planned_pair%in%TRUE),n_saved_networks=sum(saved_network),
 n_missing_networks=sum(planned_pair%in%TRUE&!saved_network),
 input_population_count=uniqueN(c(population_a,population_b),na.rm=TRUE)),by=dataset])
