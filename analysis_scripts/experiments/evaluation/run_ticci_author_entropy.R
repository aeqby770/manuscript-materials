CONFIG <- list(paths=list(root=normalizePath(commandArgs(TRUE)[1L],winslash="/",mustWork=TRUE),
 input="results/by_algorithm/TICCI_author_functional_entropy",
 vendor="inputs/comparison/external_methods/TICCI"),threads=as.integer(Sys.getenv("LLR_NATIVE_THREADS",parallel::detectCores())),plot=list(),
 seed=12345L,B=100L,B_size=1000L,exp_cutoff=1)
Sys.setenv(OMP_NUM_THREADS=CONFIG$threads,OPENBLAS_NUM_THREADS=CONFIG$threads,MKL_NUM_THREADS=CONFIG$threads)
.libPaths(c(file.path(CONFIG$paths$root,"library"),.libPaths()))
library(data.table)
library(Matrix)
library(entropy)
setDTthreads(CONFIG$threads)
dataset <- commandArgs(TRUE)[2L]
out <- file.path(CONFIG$paths$root,CONFIG$paths$input,dataset)
vendor <- file.path(CONFIG$paths$root,CONFIG$paths$vendor)
binding <- fread(file.path(out,"input_binding.csv"))
entropy_mode <- fcoalesce(as.character(binding$entropy_mode)[1L],"go_kappa")
switch(entropy_mode,
 go_kappa=load(file.path(vendor,"data",paste0(c(mouse="mm",human="hs")[[binding$species]],"_km.Rda"))),
 expression_gene_groups={km<-NULL})
x <- switch(as.character(file.exists(file.path(out,"entropy_input.h5"))),
 "TRUE"={
  h5file <- file.path(out,"entropy_input.h5")
  dimensions <- rhdf5::h5read(h5file,"shape")
  new("dgCMatrix",i=as.integer(rhdf5::h5read(h5file,"indices")),
   p=as.integer(rhdf5::h5read(h5file,"indptr")),
   x=as.numeric(rhdf5::h5read(h5file,"data")),Dim=as.integer(rev(dimensions)))
 },
 "FALSE"=as(readMM(file.path(out,"entropy_input.mtx")),"CsparseMatrix"))
rownames(x) <- readLines(file.path(out,"entropy_genes.txt"),encoding="UTF-8")
switch(as.character(file.exists(file.path(out,"ticci_kappa_id_binding.csv"))),"TRUE"={
 resource_ids <- fread(file.path(out,"ticci_kappa_id_binding.csv"))
 stopifnot(identical(rownames(km),resource_ids$kappa_symbol),!anyDuplicated(resource_ids$prior_gene_label))
 dimnames(km) <- list(resource_ids$prior_gene_label,resource_ids$prior_gene_label)
})

colnames(x) <- readLines(file.path(out,"entropy_cells.txt"),encoding="UTF-8")
gene_groups <- switch(entropy_mode,go_kappa=NULL,expression_gene_groups={
 normalized_expression<-as.matrix(x%*%Diagonal(x=10000/pmax(Matrix::colSums(x),.Machine$double.eps)))
 normalized_expression<-log1p(normalized_expression)
 centered_expression<-normalized_expression-rowMeans(normalized_expression)
 centered_expression<-centered_expression/pmax(sqrt(rowSums(centered_expression^2)),.Machine$double.eps)
 similarity<-tcrossprod(centered_expression)
 similarity[similarity>1]<-1
 similarity[similarity< -1]<- -1
 diag(similarity)<-1
 set.seed(CONFIG$seed)
 groups<-stats::kmeans(as.dist(1-similarity),centers=floor(sqrt(min(CONFIG$B_size,nrow(x))/2)),
  iter.max=10,nstart=1)$cluster
 groups<-data.table(GENE=rownames(x),CLUSTER=as.integer(groups))
 fwrite(groups,file.path(out,"ticci_expression_gene_groups.csv"))
 rm(normalized_expression,centered_expression,similarity)
 groups
})
source_lines <- readLines(file.path(vendor,"R_tools/slice.R"),encoding="UTF-8")
# Materialize only the sampled genes; retain the author's calculation and random draws.
source_lines <- sub("g.exp.m <- exp.m[as.character(gclus$GENE), ]",
 "g.exp.m <- as.matrix(exp.m[as.character(gclus$GENE), , drop=FALSE])",
 source_lines,fixed=TRUE)
source_lines <- sub("y <- aggregate(g.exp.m.b, by=list(gclus$CLUSTER), sum)",
 "y <- cbind(Group.1=sort(unique(gclus$CLUSTER)),rowsum(g.exp.m.b,gclus$CLUSTER,reorder=TRUE))",
 source_lines,fixed=TRUE)
source_text <- paste(source_lines,collapse="\n")
source_text <- sub(
 "for (i in 1:ncells) {\n        entropy[i] <- entropy::entropy.Dirichlet(y[,i], a=a)\n      }",
 "posterior_counts <- sweep(as.matrix(y),1L,a,\"+\")\n      posterior_probability <- sweep(posterior_counts,2L,colSums(posterior_counts),\"/\")\n      entropy <- -colSums(posterior_probability*log(posterior_probability))",
 source_text,fixed=TRUE)
stopifnot(grepl("posterior_probability <-",source_text,fixed=TRUE))
author_expressions <- parse(text=source_text)
author_declarations <- vapply(author_expressions,deparse1,character(1L))
eval(author_expressions[[which(startsWith(author_declarations,"scEntropy <- function"))]])
eval(author_expressions[[which(startsWith(author_declarations,"expressed.all <- function"))]])
eval(author_expressions[[which(startsWith(author_declarations,"expressed.random <- function"))]])
# The author assigns one shared bootstrap gene draw to every cell; only its union is used.
cached_allgenes <- sort(rownames(x)[Matrix::rowSums(x>CONFIG$exp_cutoff)>0])
random_declaration <- deparse(expressed.random)
random_declaration <- sub("genelists <- expressed.all(exp.m, exp.cutoff = exp.cutoff)",
 "genelists <- NULL",random_declaration,fixed=TRUE)
random_declaration <- sub("allgenes <- sort(unique(as.character(unlist(genelists))))",
 "allgenes <- cached_allgenes",random_declaration,fixed=TRUE)
random_text <- paste(random_declaration,collapse="\n")
random_text <- sub(
 "for (i in 1:n) {\n        genelists[[i]] <- rndgenes\n        names(genelists)[i] <- cells[i]\n    }",
 "genelists[[1L]] <- rndgenes",random_text,fixed=TRUE)
expressed.random <- eval(parse(text=random_text))
started <- proc.time()[["elapsed"]]
ret <- scEntropy(x,clusters=gene_groups,km=km,calculation="bootstrap",r=CONFIG$B,
 exp.cutoff=CONFIG$exp_cutoff,n=CONFIG$B_size,p=0,
 k=floor(sqrt(CONFIG$B_size/2)),cmethod="kmeans",prior="pr",
 random.seed=CONFIG$seed,context_str=paste0(dataset,"_"),verbose=FALSE)
stopifnot(length(ret$entropy)==ncol(x),all(is.finite(ret$entropy)))
fwrite(data.table(cell_id=colnames(x),entropy=as.numeric(ret$entropy)),
 file.path(out,"ticci_functional_entropy.csv"))
saveRDS(ret,file.path(out,"ticci_SLICE_entropy_bootstrap.rds"))
fwrite(data.table(dataset,n_cells=ncol(x),n_input_genes=nrow(x),
 n_functional_kappa_genes=sum(rownames(x)%in%rownames(km)),entropy_mode=entropy_mode,
 n_expression_gene_groups=uniqueN(gene_groups$CLUSTER),
 B=CONFIG$B,B_size=CONFIG$B_size,exp_cutoff=CONFIG$exp_cutoff,
 clustering_k=floor(sqrt(CONFIG$B_size/2)),seed=CONFIG$seed,
 n_finite=sum(is.finite(ret$entropy)),entropy_min=min(ret$entropy),
 entropy_max=max(ret$entropy),seconds=proc.time()[["elapsed"]]-started,
 engine="author_SLICE_bootstrap_draws;shared_gene_union;matrix_counts;Dirichlet_plugin_entropy"),
 file.path(out,"ticci_entropy_summary.csv"))
cat("TICCI_FUNCTIONAL_ENTROPY",dataset,ncol(x),"cells;",CONFIG$B,"bootstrap draws\n")
