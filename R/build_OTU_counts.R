
#' Batch - build_OTU_counts analysis of many comparisons
#
#' @description The build_OTU_counts() function combines biom, taxa and sample tables
#' to create different microbiome plots and differential abundance analyses.
#
#' @param biom A "phyloseq" object with included groups
#'  to be analysed. For format specifications see ?phyloseq E.g.
#'  accessible as "build_OTU_counts$otu". Groups are used to automate colouring of
#'  samples. Default = NULL
#'
#' @param sample_table A sample table to map groups to corresponding samples.  Default = NULL
#'
#' @param taxa taxanomic rank to match ranks. Default = NULL
#'
#' @param verbose Verbosity ON/OFF. Default=FALSE
#'
#' @return A list of all OTU comparisons conducted.
#'
#' @export build_OTU_counts
#'
#' @importFrom phyloseq import_biom import_qiime_sample_data merge_phyloseq tax_table tax_table<-
# @import biom file, sample table and tax 


build_OTU_counts <- function(biom = NULL,
                              sample_table = NULL,
                              taxa = NULL,
                              abundance_threshold = NULL,
                              prevalence_threshold = NULL,
                              rarity_threshold = NULL,
                              variance_threshold = NULL,
                              taxa_level = NULL,
                              include_taxonomy = TRUE, # Option to include/exclude taxonomy
                              force_build = FALSE,
                              verbose = FALSE){
  
  ####///---- Check Inputs ----\\###
  if(is.null(biom) | is.null(sample_table)) {
    stop("EITHER a biom or sample_table is not provided. Please provide filenames with full path and rerun.")
  }
  
  # Import biom and sample table data
  biom <- import_biom(biom)
  samples <- import_qiime_sample_data(sample_table)
  
  # Merge into a phyloseq object
  phylo <- merge_phyloseq(biom, samples)
  
  # Include taxonomy if provided
  if (include_taxonomy && !is.null(taxa)) {
    tax <- import_biom(taxa)
    phylo <- merge_phyloseq(phylo, tax)
  }
  
  # Filter by taxa level (e.g., Genus, Family)
  if (!is.null(taxa_level) && !is.null(tax_table(phylo, errorIfNULL = FALSE))) {
    phylo <- tax_glom(phylo, taxa_level)
  }
  
  # Filter by abundance
  if (!is.null(abundance_threshold)) {
    phylo <- prune_taxa(taxa_sums(phylo) > abundance_threshold, phylo)
  }
  
  # Filter by prevalence (percentage of samples an OTU is present in)
  if (!is.null(prevalence_threshold)) {
    prevalence <- apply(otu_table(phylo), 1, function(x) sum(x > 0) / length(x))
    phylo <- prune_taxa(prevalence > prevalence_threshold, phylo)
  }
  
  # Filter by rarity (low-abundance taxa)
  if (!is.null(rarity_threshold)) {
    total_abundance <- taxa_sums(phylo)
    rare_taxa <- names(total_abundance[total_abundance < rarity_threshold])
    phylo <- prune_taxa(!taxa_names(phylo) %in% rare_taxa, phylo)
  }
  
  # Filter by variance (removes low variance features)
  if (!is.null(variance_threshold)) {
    var_filter <- apply(otu_table(phylo), 1, var)
    phylo <- prune_taxa(var_filter > variance_threshold, phylo)
  }
  
  if (verbose) {
    cat("Final OTU count after filtering:", ntaxa(phylo), "\n")
  }
  
  return(phylo)
}


