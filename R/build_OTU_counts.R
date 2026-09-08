## =============================================================================
## R/build_OTU_counts.R
## =============================================================================

#' Build and filter a phyloseq object from microbiome data
#'
#' Two entry points in v2:
#' \enumerate{
#'   \item Pass an already-built phyloseq object via `physeq =` and use this
#'     function purely for the filtering cascade (abundance / prevalence /
#'     rarity / variance / agglomeration).
#'   \item Pass file paths via `biom = ` and `sample_table = ` (plus optional
#'     `tax_tables = `) to import from disk and then filter.
#' }
#' Exactly one of `physeq` or (`biom` + `sample_table`) must be supplied.
#'
#' @param physeq A prebuilt phyloseq object. If supplied, `biom` /
#'   `sample_table` / `tax_tables` are ignored and the object goes straight
#'   into the filtering cascade. Default NULL.
#' @param biom Character path to a BIOM file (JSON or HDF5). Used only when
#'   `physeq` is NULL. Default NULL.
#' @param sample_table Character path to a QIIME-format sample metadata file.
#'   Used only when `physeq` is NULL. Default NULL.
#' @param tax_tables Optional taxonomy table (path or object accepted by
#'   phyloseq::tax_table). Used only when `physeq` is NULL. Default NULL.
#' @param taxa_level Character taxonomic rank to agglomerate to (e.g. "Genus"),
#'   or NULL for no agglomeration.
#' @param include_taxonomy Logical, whether to attach taxonomy when importing
#'   from files. Default TRUE.
#' @param abundance_threshold Numeric, minimum total abundance to retain an
#'   OTU, or NULL.
#' @param prevalence_threshold Numeric in (0,1), minimum fraction of samples an
#'   OTU must appear in, or NULL.
#' @param rarity_threshold Numeric, remove OTUs with total abundance below this,
#'   or NULL.
#' @param variance_threshold Numeric, minimum across-sample variance to retain,
#'   or NULL.
#' @param force_build Logical, currently unused. Default FALSE.
#' @param verbose Logical, print filtering progress. Default FALSE.
#' @return A filtered phyloseq object.
#' @importFrom phyloseq import_biom import_qiime_sample_data merge_phyloseq
#' @importFrom phyloseq tax_table tax_glom prune_taxa taxa_sums
#' @importFrom phyloseq taxa_names ntaxa otu_table sample_data
#' @export
build_OTU_counts <- function(physeq               = NULL,
                             biom                 = NULL,
                             sample_table         = NULL,
                             tax_tables           = NULL,
                             taxa_level           = NULL,
                             include_taxonomy     = TRUE,
                             abundance_threshold  = NULL,
                             prevalence_threshold = NULL,
                             rarity_threshold     = NULL,
                             variance_threshold   = NULL,
                             force_build          = FALSE,
                             verbose              = FALSE) {
  
  ## ---- resolve the input: prebuilt object OR import from files -------------
  if (!is.null(physeq)) {
    if (!inherits(physeq, "phyloseq"))
      stop("`physeq` must be a phyloseq object; got class ", class(physeq)[1])
    if (!is.null(biom) || !is.null(sample_table))
      message("`physeq` supplied - ignoring `biom` / `sample_table` / `tax_tables`.")
    phylo <- physeq
    if (verbose) message("Using supplied phyloseq object: ",
                         phyloseq::ntaxa(phylo), " taxa x ",
                         phyloseq::nsamples(phylo), " samples")
    
  } else {
    if (is.null(biom) || is.null(sample_table))
      stop("Provide either `physeq`, or both `biom` and `sample_table`.")
    
    biom_obj <- phyloseq::import_biom(biom)
    samples  <- phyloseq::import_qiime_sample_data(sample_table)
    
    ## tax_tables may be: a file path to a tab-delimited tax file (OTU IDs in
    ## column 1), an already-built matrix, or a phyloseq tax_table object.
    ## phyloseq::tax_table() does NOT read a path, so detect and read it here.
    tt2_tax_test <- if (is.null(tax_tables)) {
      NULL
    } else if (is.character(tax_tables) && length(tax_tables) == 1 && file.exists(tax_tables)) {
      tax_mat <- as.matrix(utils::read.delim(tax_tables, header = TRUE, row.names = 1,
                                             sep = "\t", stringsAsFactors = FALSE,
                                             check.names = FALSE))
      phyloseq::tax_table(tax_mat)
    } else {
      phyloseq::tax_table(tax_tables)
    }
    
    if (include_taxonomy && !is.null(tt2_tax_test)) {
      phylo <- phyloseq::merge_phyloseq(biom_obj, samples, tt2_tax_test)
      tax_col <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
      ## only rename if the tax_table has the expected 7 columns
      if (ncol(phyloseq::tax_table(phylo)) == length(tax_col))
        colnames(phyloseq::tax_table(phylo)) <- tax_col
    } else {
      phylo <- phyloseq::merge_phyloseq(biom_obj, samples)
    }
    if (verbose) message("Imported from files: ",
                         phyloseq::ntaxa(phylo), " taxa x ",
                         phyloseq::nsamples(phylo), " samples")
  }
  
  ## ---- filtering cascade (identical for both entry points) ----------------
  ## Guard: the prevalence/variance filters below use apply(otu_table, 1, ...),
  ## which is per-taxon ONLY when taxa are rows. A samples-as-rows object would
  ## silently filter on the wrong dimension (wrong-length logical -> prune_taxa
  ## error). Force taxa-as-rows so every filter is orientation-safe.
  if (!phyloseq::taxa_are_rows(phylo)) {
    if (verbose) message("  taxa_are_rows is FALSE - transposing for filtering")
    phylo <- phyloseq::t(phylo)
  }
  
  if (!is.null(taxa_level) &&
      !is.null(phyloseq::tax_table(phylo, errorIfNULL = FALSE))) {
    phylo <- phyloseq::tax_glom(phylo, taxa_level)
    if (verbose) message("  after tax_glom(", taxa_level, "): ", phyloseq::ntaxa(phylo), " taxa")
  }
  
  if (!is.null(abundance_threshold)) {
    phylo <- phyloseq::prune_taxa(phyloseq::taxa_sums(phylo) > abundance_threshold, phylo)
    if (verbose) message("  after abundance > ", abundance_threshold, ": ", phyloseq::ntaxa(phylo), " taxa")
  }
  
  if (!is.null(prevalence_threshold)) {
    ## taxa are rows here (guaranteed by the guard above); one value per taxon
    otu <- as(phyloseq::otu_table(phylo), "matrix")
    prevalence <- rowSums(otu > 0) / ncol(otu)
    phylo <- phyloseq::prune_taxa(prevalence > prevalence_threshold, phylo)
    if (verbose) message("  after prevalence > ", prevalence_threshold, ": ", phyloseq::ntaxa(phylo), " taxa")
  }
  
  if (!is.null(rarity_threshold)) {
    total_abundance <- phyloseq::taxa_sums(phylo)
    rare_taxa <- names(total_abundance[total_abundance < rarity_threshold])
    phylo <- phyloseq::prune_taxa(!phyloseq::taxa_names(phylo) %in% rare_taxa, phylo)
    if (verbose) message("  after rarity >= ", rarity_threshold, ": ", phyloseq::ntaxa(phylo), " taxa")
  }
  
  if (!is.null(variance_threshold)) {
    otu <- as(phyloseq::otu_table(phylo), "matrix")
    var_filter <- apply(otu, 1, stats::var)     # taxa are rows (guard above)
    phylo <- phyloseq::prune_taxa(var_filter > variance_threshold, phylo)
    if (verbose) message("  after variance > ", variance_threshold, ": ", phyloseq::ntaxa(phylo), " taxa")
  }
  
  if (verbose) message("Final OTU count after filtering: ", phyloseq::ntaxa(phylo))
  phylo
}
