CONFIG <- list(root=Sys.getenv("LLR_WORK_ROOT"),dataset=commandArgs(TRUE)[2L],threads=parallel::detectCores())
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads)
library(data.table)
library(AnnotationDbi)
library(org.Hs.eg.db)
setDTthreads(CONFIG$threads)
out <- file.path(CONFIG$root,"results/by_algorithm/TICCI_author_functional_entropy",CONFIG$dataset)
genes <- readLines(file.path(out,"entropy_genes.txt"),encoding="UTF-8")
primary <- as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=genes,keytype="ENSEMBL",columns="SYMBOL"))
primary <- unique(primary[!is.na(SYMBOL),.(original_gene_id=ENSEMBL,resource_symbol=SYMBOL)])
primary[,n_symbols_per_input:=uniqueN(resource_symbol),by=original_gene_id]
primary[,n_inputs_per_symbol:=uniqueN(original_gene_id),by=resource_symbol]
unique_primary <- primary[n_symbols_per_input==1L&n_inputs_per_symbol==1L,.(original_gene_id,resource_symbol)]
gene_binding <- merge(data.table(original_gene_id=genes,input_column=seq_along(genes)),unique_primary,
 by="original_gene_id",all.x=TRUE,sort=FALSE)
setorder(gene_binding,input_column)
gene_binding[,mapping_status:=fifelse(is.na(resource_symbol),"unmapped_or_ambiguous_retained_original_ID","unique_Ensembl_to_primary_symbol")]
gene_binding[,resource_symbol:=fcoalesce(resource_symbol,original_gene_id)]
gene_binding[,annotation_package:=paste0("org.Hs.eg.db_",packageVersion("org.Hs.eg.db"))]
stopifnot(nrow(gene_binding)==length(genes),!anyDuplicated(gene_binding$resource_symbol))
fwrite(gene_binding,file.path(out,"ticci_gene_id_binding.csv"))
load(file.path(CONFIG$root,"inputs/comparison/external_methods/TICCI/data/hs_km.Rda"))
kappa_symbols <- data.table(kappa_symbol=rownames(km),kappa_row=seq_len(nrow(km)))
direct <- merge(kappa_symbols,unique_primary,by.x="kappa_symbol",by.y="resource_symbol")
aliases <- as.data.table(AnnotationDbi::select(org.Hs.eg.db,keys=kappa_symbols$kappa_symbol,
 keytype="ALIAS",columns="ENSEMBL"))
aliases <- unique(aliases[ENSEMBL%in%genes,.(kappa_symbol=ALIAS,original_gene_id=ENSEMBL)])
aliases <- aliases[!kappa_symbol%in%direct$kappa_symbol]
candidates <- unique(rbindlist(list(direct[,.(kappa_symbol,original_gene_id)],aliases)))
candidates[,n_inputs:=uniqueN(original_gene_id),by=kappa_symbol]
candidates[,n_kappa_symbols:=uniqueN(kappa_symbol),by=original_gene_id]
matched <- candidates[n_inputs==1L&n_kappa_symbols==1L,.(kappa_symbol,original_gene_id)]
kappa_binding <- merge(kappa_symbols,matched,by="kappa_symbol",all.x=TRUE,sort=FALSE)
setorder(kappa_binding,kappa_row)
kappa_binding[,prior_gene_label:=fcoalesce(original_gene_id,paste0("unmatched_resource::",kappa_symbol))]
kappa_binding[,numeric_prior_modified:=FALSE]
stopifnot(!anyDuplicated(kappa_binding$prior_gene_label),sum(!is.na(kappa_binding$original_gene_id))>1000L)
fwrite(kappa_binding,file.path(out,"ticci_kappa_id_binding.csv"))
fwrite(data.table(dataset=CONFIG$dataset,n_input_genes=length(genes),
 n_input_cells=nrow(fread(file.path(out,"entropy_cells.txt"),header=FALSE)),
 n_primary_symbol_mapped=sum(gene_binding$mapping_status=="unique_Ensembl_to_primary_symbol"),
 n_kappa_rows_bound=sum(!is.na(kappa_binding$original_gene_id)),
 input_gene_rows_dropped=0L,numeric_prior_modified=FALSE,input_expression_modified=FALSE,
 bootstrap_gene_ID_order_modified=FALSE),
 file.path(out,"ticci_resource_binding_summary.csv"))
print(fread(file.path(out,"ticci_resource_binding_summary.csv")))
