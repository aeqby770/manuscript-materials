# Formal API contracts for null calibration and adaptive methods.

test_that("ComputeFakeLRNull rejects n_fake_per_real < 100", {
  skip_if_not_installed("Matrix")
  set.seed(1)
  n_cells <- 50L
  x <- Matrix::rsparsematrix(15L, n_cells, density = 0.3)
  x@x <- abs(x@x) + 0.05
  ct <- rep(1:2, each = 25L)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(1, n_cells, 1L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  expect_error(
    ComputeFakeLRNull(
      x, ct, lr_db, pt, lw,
      candidate_fake_genes = 3:8,
      n_fake_per_real = 10L,
      n_bins = 4L, min_pct = 0, min_n_eff = 0.1,
      min_valid_bins = 2L, min_valid_fraction = 0,
      seed = 1L
    ),
    "n_fake_per_real >= 100"
  )
})

test_that("ComputeTrajectoryFakeLRNull rejects n_fake_per_real < 100", {
  skip_if_not_installed("Matrix")
  set.seed(2)
  n_cells <- 50L
  x <- Matrix::rsparsematrix(15L, n_cells, density = 0.3)
  x@x <- abs(x@x) + 0.05
  ct <- rep(1:2, each = 25L)
  pt <- seq(0, 1, length.out = n_cells)
  lw <- matrix(1, n_cells, 1L)
  lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
  expect_error(
    ComputeTrajectoryFakeLRNull(
      x, ct, lr_db, 3:8, pt, lw,
      n_fake_per_real = 50L,
      n_bins = 4L, min_pct = 0, min_n_eff = 0.1,
      min_valid_bins = 2L, min_valid_fraction = 0,
      min_fake_observed = 0L, min_fake_coverage = 0,
      seed = 1L
    ),
    "n_fake_per_real >= 100"
  )
})

test_that("adaptive APIs run without experimental warning", {
  f <- local({
    skip_if_not_installed("Matrix")
    set.seed(3)
    n_cells <- 120L
    x <- Matrix::rsparsematrix(25L, n_cells, density = 0.3)
    x@x <- abs(x@x)
    ct <- sample.int(2L, n_cells, replace = TRUE)
    pt <- seq(0, 1, length.out = n_cells)
    lw <- matrix(1, n_cells, 1L)
    lr_db <- list(ligand_rows = list(1L), receptor_rows = list(2L))
    list(x = x, ct = ct, pt = pt, lw = lw, lr_db = lr_db)
  })
  expect_no_warning(
    ComputeAdaptiveTemporalPermutation(
      f$x, f$ct, f$lr_db, f$pt, f$lw,
      n_bins = 4L, min_pct = 0, min_n_eff = 0.1,
      min_valid_bins = 2L, min_valid_fraction = 0,
      n_perms_screen = 19L, n_perms_confirm = 29L,
      confirm_top_k = 2L, q_confirm_threshold = 1,
      min_perm_observed = 1L, min_perm_coverage = 0, seed = 4L
    )
  )
})

test_that("disconnected forest and zero parent spread give valid local memberships", {
  emb <- rbind(c(0,0),c(0,0),c(1,1),c(1,1),c(1,-1),c(1,-1),
    c(10,0),c(10,0),c(11,0),c(11,0))
  labels <- rep(1:5,each=2L)
  paths <- list(c(1L,2L),c(1L,3L),c(4L,5L))
  q <- ComputeLineageWeightsLocalBranchResidual(emb,labels,paths)
  expect_true(all(is.finite(q)))
  expect_equal(rowSums(q),rep(1,10),tolerance=1e-6)
  expect_equal(q[1:6,3],rep(0,6))
  expect_equal(q[7:10,1:2],matrix(0,4,2))
  expect_equal(q[1:2,1:2],matrix(.5,2,2))
})

test_that("state refinement retains weak observed bridges without joining separate components", {
  emb <- cbind(c(seq(0,.5,length.out=6),seq(10,10.5,length.out=6)),
    rep(c(0,.01),6))
  edge <- seq_len(11L)
  support <- rep(1,11)
  support[6] <- 1e-4
  conn <- Matrix::sparseMatrix(i=c(edge,edge+1L),j=c(edge+1L,edge),
    x=rep(support,2),dims=c(12,12))
  connected <- RefineStateGraph(emb,conn,1L,k=4L,n_neighbors=2L,project_root=FALSE)
  expect_equal(connected$tree$n_components,1L)
  expect_gt(connected$diagnostics$native_bridge_edges,0)
  conn[6,7] <- conn[7,6] <- 0
  forest <- RefineStateGraph(emb,conn,1L,k=4L,n_neighbors=2L,project_root=FALSE)
  expect_equal(forest$tree$n_components,2L)
  expect_true(all(is.finite(forest$lineage_weights)))
  expect_equal(rowSums(forest$lineage_weights),rep(1,12),tolerance=1e-6)
  expect_setequal(unique(unlist(forest$lineage_paths)),unique(forest$labels))
})

test_that("singleton components retain unit mass without artificial connections", {
  emb <- cbind(c(0,1,50),c(0,0,0))
  conn <- Matrix::sparseMatrix(i=c(1,2),j=c(2,1),x=1,dims=c(3,3))
  fit <- RefineStateGraph(emb,conn,1L,k=2L,n_neighbors=1L,project_root=FALSE)
  expect_equal(fit$tree$n_components,2L)
  expect_equal(length(fit$lineage_paths),2L)
  expect_equal(rowSums(fit$lineage_weights),rep(1,3))
  q <- ComputeLineageWeightsLocalBranchResidual(fit$embedding,fit$labels,fit$lineage_paths)
  expect_true(all(is.finite(q)))
  expect_equal(rowSums(q),rep(1,3))
})

test_that("requested cell-type labels survive gaps in within-sample numeric IDs", {
  expr <- Matrix::Matrix(matrix(1,2,8),sparse=TRUE)
  celltypes <- rep(c(2L,4L),4)
  pt <- rep(c(.1,.3,.7,.9),each=2)
  db <- list(ligand_rows=list(1L),receptor_rows=list(2L))
  request <- data.frame(pair_idx=1L,source=2L,target=4L)
  result <- ComputeTemporalLRFixedBounds(expr,celltypes,db,pt,matrix(1,8,1),
    list(c(0,.5,1)),1,min_n_eff=1,requested_contexts=request)
  expect_equal(nrow(result),2L)
  expect_equal(result$source,rep(2L,2))
  expect_equal(result$target,rep(4L,2))
  expect_equal(result$score,rep(.5,2))
})

test_that("observed sibling levels preserve paired observations and missing contexts", {
  profiles <- list(scores=data.frame(pair_idx=1L,source=1L,target=2L,
    biological_sample=rep(c("s1","s2"),each=2),observed_time=rep(1:2,each=2),
    condition="all",lineage=rep(c(2L,3L),2),
    score=c(.5,.2,.8,.6),structural_valid=TRUE,expression_supported=TRUE),
    observations=data.frame(biological_sample=c("s1","s2"),observed_time=1:2,condition="all"),
    requested_contexts=data.frame(pair_idx=1:2,source=1L,target=2L),
    profile_nodes=1:3,sampling_unit="biological_sample")
  units <- list(units=list(list(event_node=1L,children=c(2L,3L))))
  result <- SummarizeObservedSiblingLevels(profiles,units)
  expected <- mean(c(log(1.6),log(2)))
  expect_equal(result$effects[pair_idx==1L & child==2L,level_effect],expected)
  expect_equal(result$effects[pair_idx==1L & child==3L,level_effect],-expected)
  expect_true(all(result$effects[pair_idx==2L,n_samples]==0L))
  expect_true(all(is.na(result$effects[pair_idx==2L,level_effect])))
  zero <- profiles
  zero$scores$score <- rep(c(.5,0),2)
  zero$scores$expression_supported <- rep(c(TRUE,FALSE),2)
  observed_zero <- SummarizeObservedSiblingLevels(zero,units)
  expect_equal(observed_zero$effects[pair_idx==1L & child==2L,level_effect],log(2))
  profiles$scores <- profiles$scores[-4L,]
  partial <- SummarizeObservedSiblingLevels(profiles,units)
  expect_equal(partial$effects[pair_idx==1L & child==2L,level_effect],log(1.6))
  expect_equal(partial$effects[pair_idx==1L & child==2L,paired_sample_fraction],.5)
})
