# Shared by the full runner and LR-only recomputation. CellChat partner names
# can be complex identifiers; only the complex table defines their subunits.
PrepareLRExpression <- function(expr, scale_factor = 10000, input_scale = "counts") {
  switch(input_scale,
    tpm = {
      # TPM is already measured per one million transcripts.
      normalized <- expr * (scale_factor / 1e6)
      normalized@x <- log1p(normalized@x)
      normalized
    },
    counts = {
      normalized <- expr %*% Matrix::Diagonal(x = scale_factor / Matrix::colSums(expr))
      normalized@x <- log1p(normalized@x)
      normalized
    },
    counts_linear = expr %*% Matrix::Diagonal(x = scale_factor / Matrix::colSums(expr)),
    counts_sqrt = {
      normalized <- expr %*% Matrix::Diagonal(x = scale_factor / Matrix::colSums(expr))
      normalized@x <- sqrt(normalized@x)
      normalized
    },
    counts_positive_reference = {
      normalized <- expr %*% Matrix::Diagonal(x = scale_factor / Matrix::colSums(expr))
      positive_mean <- Matrix::rowSums(normalized)/pmax(1,Matrix::rowSums(expr>0))
      reference_factor <- log1p(positive_mean)/pmax(positive_mean,.Machine$double.xmin)
      normalized <- Matrix::Diagonal(x=reference_factor)%*%normalized
      dimnames(normalized) <- dimnames(expr)
      normalized
    },
    provided_expression = expr,
    stop("Unknown LR input scale: ", input_scale)
  )
}

BuildLRDB <- function(expr, species = "human") {
  if (species == "celegans") {
    native <- data.table::fread(CONFIG$lr_resources$celegans)
    native[, `:=`(ligand_row = match(Ligand_WB, rownames(expr)),
      receptor_row = match(Receptor_WB, rownames(expr)))]
    native <- native[!is.na(ligand_row) & !is.na(receptor_row)]
    pair_map <- native[, .(pair_idx = seq_len(.N),
      ligand_genes = Ligand_WB, receptor_genes = Receptor_WB,
      ligand_id = Ligand_symbol, receptor_id = Receptor_symbol,
      interaction_name = paste(Ligand_symbol, Receptor_symbol, sep = "_"),
      pathway_name = `LR Function`,
      interaction_aliases = paste(Ligand_symbol, Receptor_symbol, sep = "_"),
      source_doi, evidence_level = Pair_confidence)]
    return(list(ligand_rows = as.list(native$ligand_row),
      receptor_rows = as.list(native$receptor_row), pair_gene_map = pair_map))
  }
  if (!species %in% c("human", "mouse", "zebrafish")) {
    return(list(ligand_rows = list(), receptor_rows = list(), pair_gene_map = NULL))
  }
  cc_db <- switch(species, mouse = CellChat::CellChatDB.mouse,
    human = CellChat::CellChatDB.human, zebrafish = CellChat::CellChatDB.zebrafish)
  interactions <- data.table::as.data.table(cc_db$interaction)
  complex_table <- data.table::as.data.table(cc_db$complex, keep.rownames = "partner")
  subunits <- data.table::melt(complex_table, id.vars = "partner",
    variable.name = "subunit_slot", value.name = "gene", variable.factor = FALSE)
  subunits <- subunits[!is.na(gene) & nzchar(gene), .(partner, gene)]
  partners <- unique(c(interactions$ligand, interactions$receptor))
  monomers <- setdiff(partners, complex_table$partner)
  subunits <- data.table::rbindlist(list(subunits,
    data.table::data.table(partner = monomers, gene = monomers)))
  input_keys <- rownames(expr)
  switch(as.character(species == "human" &&
    all(grepl("^ENSG[0-9]+$", input_keys))),
    "TRUE" = {
      symbol_map <- data.table::as.data.table(AnnotationDbi::select(
        org.Hs.eg.db::org.Hs.eg.db,
        keys = input_keys, columns = "SYMBOL", keytype = "ENSEMBL"
      ))
      symbol_map <- unique(symbol_map[
        !is.na(SYMBOL), .(ENSEMBL, SYMBOL)
      ])
      symbol_map[, n_symbols := data.table::uniqueN(SYMBOL), by = ENSEMBL]
      symbol_map[, n_ids := data.table::uniqueN(ENSEMBL), by = SYMBOL]
      symbol_map <- symbol_map[n_symbols == 1L & n_ids == 1L]
      matched_symbols <- symbol_map$SYMBOL[
        match(input_keys, symbol_map$ENSEMBL)
      ]
      input_keys[!is.na(matched_symbols)] <-
        matched_symbols[!is.na(matched_symbols)]
    }
  )
  subunits[, gene_row := match(gene, input_keys)]
  if (species == "zebrafish" ||
      (species == "mouse" && all(input_keys == toupper(input_keys)))) {
    # Upper-case deposited symbols match only unique database gene symbols.
    input_keys <- toupper(rownames(expr))
    unique_rows <- which(!duplicated(input_keys) & !duplicated(input_keys, fromLast = TRUE))
    subunits[, gene_row := unique_rows[match(toupper(gene), input_keys[unique_rows])]]
  }
  subunits[, input_gene := rownames(expr)[gene_row]]
  resolved <- subunits[, .(rows = list(gene_row),
    genes = paste(input_gene, collapse = "_"), complete = all(!is.na(gene_row))), by = partner]
  resolved <- resolved[complete == TRUE]
  interactions[, `:=`(ligand_index = match(ligand, resolved$partner),
    receptor_index = match(receptor, resolved$partner))]
  interactions <- interactions[!is.na(ligand_index) & !is.na(receptor_index)]
  # A database synonym (for example CXCL4/PF4) is one molecular pair.
  interactions[, interaction_aliases := paste(unique(interaction_name), collapse = ";"),
    by = .(ligand_index, receptor_index)]
  interactions <- unique(interactions, by = c("ligand_index", "receptor_index"))
  pair_map <- interactions[, .(pair_idx = seq_len(.N),
    ligand_genes = resolved$genes[ligand_index],
    receptor_genes = resolved$genes[receptor_index],
    ligand_id = ligand, receptor_id = receptor, interaction_name, pathway_name,
    interaction_aliases)]
  list(ligand_rows = resolved$rows[interactions$ligand_index],
    receptor_rows = resolved$rows[interactions$receptor_index], pair_gene_map = pair_map)
}
