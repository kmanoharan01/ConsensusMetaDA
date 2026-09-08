# Load the packages
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(ALDEx2)
library(edgeR)
library(metagenomeSeq)
library(ADAPT)
library(dplyr)
library(tibble)
library(vegan)
library(ggplot2)
library(ConsensusMetaDA)

build_OTU_counts_NoTaxa <- function(
    biom = NULL,
    sample_table = NULL,
    tax_tables = NULL,
    taxa_level = NULL,
    include_taxonomy = TRUE,
    abundance_threshold = NULL,
    prevalence_threshold = NULL,
    rarity_threshold = NULL,
    variance_threshold = NULL,
    force_build = FALSE,
    verbose = FALSE
) {
  
  #### Check Inputs ####
  if (is.null(biom) || is.null(sample_table)) {
    stop("Both `biom` and `sample_table` must be provided.")
  }
  
  # Import biom and sample table
  biom <- import_biom(biom)
  samples <- import_qiime_sample_data(sample_table)
  
  # Taxonomy (optional)
  if (!is.null(tax_tables) && include_taxonomy) {
    tt2_tax_test <- tax_table(tax_tables)
    phylo <- merge_phyloseq(biom, samples, tt2_tax_test)
    
    # assign taxonomy names (optional but safe)
    tax_col <- c("Kingdom","Phylum","Class","Order",
                 "Family","Genus","Species")
    if (ncol(tax_table(phylo)) == length(tax_col)) {
      colnames(tax_table(phylo)) <- tax_col
    }
    
  } else {
    phylo <- merge_phyloseq(biom, samples)
  }
  
  
  #### Taxonomic Glom by level ####
  if (!is.null(taxa_level) && !is.null(tax_table(phylo, errorIfNULL = FALSE))) {
    phylo <- tax_glom(phylo, taxa_level)
  }
  
  
  #### Abundance filter ####
  if (!is.null(abundance_threshold)) {
    phylo <- prune_taxa(taxa_sums(phylo) > abundance_threshold, phylo)
  }
  
  
  #### Prevalence filter ####
  if (!is.null(prevalence_threshold)) {
    prevalence <- apply(otu_table(phylo), 1, function(x) sum(x > 0) / length(x))
    phylo <- prune_taxa(prevalence > prevalence_threshold, phylo)
  }
  
  
  #### Rarity filter ####
  if (!is.null(rarity_threshold)) {
    total_abundance <- taxa_sums(phylo)
    rare_taxa <- names(total_abundance[total_abundance < rarity_threshold])
    phylo <- prune_taxa(!taxa_names(phylo) %in% rare_taxa, phylo)
  }
  
  
  #### Variance filter ####
  if (!is.null(variance_threshold)) {
    var_filter <- apply(otu_table(phylo), 1, var)
    phylo <- prune_taxa(var_filter > variance_threshold, phylo)
  }
  
  
  #### Verbose output ####
  if (verbose) {
    cat("Final OTU count after filtering:", ntaxa(phylo), "\n")
    cat("Samples:", nsamples(phylo), "\n")
    if (!is.null(tax_table(phylo, errorIfNULL = FALSE))) {
      cat("Taxonomy included\n")
    } else {
      cat("No taxonomy\n")
    }
  }
  
  return(phylo)
}
################### Figure 3 #####################
wget ftp.microbio.me/emp/release1/otu_tables/deblur/emp_deblur_150bp.subset_2k.rare_5000.biom or 

Download via https://ftp.microbio.me/emp/release1/otu_tables/deblur/ 

file name: emp_deblur_150bp.subset_2k.rare_5000.biom

or
unzip ../Data/emp_deblur_150bp.subset_2k.rare_5000.biom

## emp
## Load data
biome_file <- "./Data/emp_deblur_150bp.subset_2k.rare_5000.biom"

sample_table_file <-  "./Data/samples_table_emp.txt"


## Load data to create phyloseq object
emp <- build_OTU_counts_NoTaxa(biom = biome_file, sample_table = sample_table_file)

## Plot using phyloseq object
emp <- OTU_plots(emp)


# An example vignette is provided (https://github.com/kmanoharan01/ConsensusMetaDA/blob/main/inst/doc/ConsensusMetaDA_Manual.html)
