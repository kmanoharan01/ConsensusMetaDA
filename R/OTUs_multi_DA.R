#' Multi-Method Differential Abundance Analysis for OTU Data
#'
#' @title OTUs_multi_DA: Comprehensive Differential Abundance Analysis Suite
#'
#' @description
#' This function performs differential abundance analysis using five different
#' statistical methods: EdgeR, DESeq2, ALDEx2, ADAPT, and metagenomeSeq. It
#' automatically compares all pairs of groups within the Age_Group variable
#' and generates comprehensive results tables for each comparison. The function
#' is designed for robust microbiome differential abundance analysis with
#' multiple validation approaches.
#'
#' @param build_OTU_counts_output (Required).
#'
#'
#' @return A list containing results from all differential abundance methods
#'   and comparisons.
#'
#' @importFrom ggplot2 ggplot aes geom_bar
#' @importFrom phyloseq get_variable nsamples otu_table tax_table psmelt
#' @importFrom DESeq2 DESeqDataSet DESeq results
#' @importFrom edgeR DGEList exactTest estimateCommonDisp estimateTagwiseDisp topTags estimateDisp estimateGLMCommonDisp estimateGLMTagwiseDisp calcNormFactors glmFit
#' @importFrom ALDEx2 aldex
#' @importFrom ADAPT adapt
#' @importFrom metagenomeSeq newMRexperiment cumNorm fitFeatureModel MRfulltable
#' @importFrom metagenomeSeq newMRexperiment cumNorm cumNormStat fitFeatureModel MRfulltable
#' @importFrom Biobase pData
#' @importFrom phyloseq otu_table sample_data tax_table taxa_are_rows nsamples
#' @importFrom edgeR DGEList calcNormFactors estimateDisp exactTest topTags
#' @importFrom ALDEx2 aldex
#' @importFrom ADAPT adapt summary
#' @importFrom Biobase pData
#' @importFrom stats model.matrix
#'
#'
#'
#' @param physeq (Required).  A \code{\link{phyloseq-class}} or
#'  an \code{\link{otu_table-class}} object.
#'  The latter is only appropriate if \code{group} argument is also a
#'  vector or factor with length equal to \code{nsamples(physeq)}.
#'
#' @param group (Required). A character vector or factor giving the experimental
#'  group/condition for each sample/library. Alternatively, you may provide
#'  the name of a sample variable. This name should be among the output of
#'  \code{sample_variables(physeq)}, in which case
#'  \code{get_variable(physeq, group)} would return either a character vector or factor.
#'  This is passed on to \code{\link[edgeR]{DGEList}},
#'  and you may find further details or examples in its documentation.
#'
#' @param method (Optional). The label of the edgeR-implemented normalization to use.
#'  See \code{\link[edgeR]{calcNormFactors}} for supported options and details.
#'  The default option is \code{"RLE"}, which is a scaling factor method
#'  proposed by Anders and Huber (2010).
#'  At time of writing, the \link[edgeR]{edgeR} package supported
#'  the following options to the \code{method} argument:
#'
#'  \code{c("TMM", "RLE", "upperquartile", "none")}.
#'
#' @param ... Additional arguments passed on to \code{\link[edgeR]{DGEList}}
#'

OTUs_multi_DA <- function(build_OTU_counts_output,
                          force_build = FALSE,
                          verbose = FALSE){


################# R Version 4.4 #############
phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  # Check required packages
  if(!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package 'edgeR' is required but not installed")
  }
  if(!requireNamespace("phyloseq", quietly = TRUE)) {
    stop("Package 'phyloseq' is required but not installed")
  }
  
  # Enforce orientation
  if(!phyloseq::taxa_are_rows(physeq)) {
    physeq <- phyloseq::t(physeq) 
  }
  
  # Convert OTU table to matrix
  x = as(phyloseq::otu_table(physeq), "matrix")
  
  # Add one to protect against overflow, log(0) issues
  x = x + 1
  
  # Check `group` argument
  if(length(group) == 1 && phyloseq::nsamples(physeq) > 1) {
    # Assume that group was a sample variable name (must be categorical)
    group = phyloseq::sample_data(physeq)[[group]]
    if(is.null(group)) {
      stop(paste("Sample variable", group, "not found in sample_data"))
    }
  }
  
  # Define gene annotations (`genes`) as tax_table
  taxonomy = phyloseq::tax_table(physeq, errorIfNULL=FALSE)
  if(!is.null(taxonomy)) {
    taxonomy = data.frame(as(taxonomy, "matrix"))
    rownames(taxonomy) = rownames(x)
  } 
  
  # Now turn into a DGEList
  y = edgeR::DGEList(
    counts = x, 
    group = group, 
    genes = taxonomy, 
    remove.zeros = TRUE, 
    ...
  )
  
  # Calculate the normalization factors
  z = edgeR::calcNormFactors(y, method = method)
  
  # Check for division by zero inside `calcNormFactors`
  if(!all(is.finite(z$samples$norm.factors))) {
    stop("Non-finite normalization factors detected. Consider changing the 'method' argument")
  }
  
  # Estimate dispersions - updated to use the recommended workflow in edgeR 4.4
  z = edgeR::estimateDisp(z)
  
  return(z)
}

##########################################################

   ####///---- check inputs ----\\\###
  if(is.null(build_OTU_counts_output)) {
    stop("A build_OTU_counts_output object is not provided. Please provide filenames with full path and rerun.")
  }  
  
  
  treat_list <- unique(sample_data(build_OTU_counts_output)$Age_Group)
  
  # treat_list <- unique(sample_data(build_OTU_counts_output)$Age_Group)
  # Initialize a list to store the results of all comparisons
  all_comparisons_results <- list()
  
  # Prepare contrast list
  contrast_list <- list()
  DESeq2_OTU_DE_results <- list()
  for (treatment in treat_list) {
    for (other_treatment in treat_list) {
      if (treatment != other_treatment) {
        
        contrast_list <- c(contrast_list, list(c("Age_Group", treatment, other_treatment)))
        
        contrast_pair <- paste0(treatment, "_vs_", other_treatment)
        cat("Running comparison:", contrast_pair, "\n")
        
        # Initialize a list to store the results of each method for this comparison
        comparison_results <- list()
        
        ######## EdgeR ###########
        
        print(treatment)
        print(other_treatment)
        
        
        # Subset the samples based on the pair of treatments
        filtered_sample_data <- filtered_sample_data <- sample_data(build_OTU_counts_output)[
          sample_data(build_OTU_counts_output)$Age_Group %in% c(treatment, other_treatment), ]
        
        # Filter OTU table and taxonomy based on filtered sample data
        filtered_OTU_table <- otu_table(build_OTU_counts_output)[, rownames(filtered_sample_data)]
        #  filtered_tax_table <- tax_table(build_OTU_counts_output)
        
        # Create a new phyloseq object with the filtered data
        # build_OTU_subset <- phyloseq(filtered_OTU_table, filtered_tax_table, filtered_sample_data)
        
        build_OTU_subset <- phyloseq(filtered_OTU_table,  filtered_sample_data)
        
        # Convert the subsetted phyloseq object to an edgeR-compatible object
        test_phylo_reads_edgeR <- phyloseq_to_edgeR(build_OTU_subset, group = "Age_Group")
        
        # Perform differential abundance analysis
        et <- edgeR::exactTest(test_phylo_reads_edgeR) 
        
        # Extract top differentially abundant features
        tt <- edgeR::topTags(et, n = nrow(test_phylo_reads_edgeR$table), adjust.method = "BH", sort.by = "PValue")
        
        #comparison
        print(tt[1,1])
        print("EdgeR: ")
        print( tt$comparison)
        
        # Extract results
        EdgeR_OTU_DE_results <- tt$table
        
        comparison_results$EdgeR <- EdgeR_OTU_DE_results         
        
        write.table(EdgeR_OTU_DE_results, file=paste0(treatment,"_vs_",other_treatment,"_edgeR_DE_results.txt"), quote=F, sep="\t", col.names = TRUE)
        
        
        ######## DESeq2 ########
        print(treatment)
        print(other_treatment)
        
        
        phylo_reads_collapsed_deseq <- phyloseq_to_deseq2(build_OTU_counts_output, ~Age_Group)
        

        #         phylo_reads_collapsed_deseq <- DESeq(phylo_reads_collapsed_deseq, test="Wald",fitType="parametric")
        phylo_reads_collapsed_deseq <- DESeq(phylo_reads_collapsed_deseq, sfType = "poscounts")
        
        DESeq2_OTU_DE_results = results(phylo_reads_collapsed_deseq , contrast=c("Age_Group", treatment, other_treatment), tidy=T, format="DataFrame")
        
        comparison_results$DESeq2 <- DESeq2_OTU_DE_results
        
        write.table(DESeq2_OTU_DE_results, file=paste0(treatment,"_vs_",other_treatment,"_DESeq_DE_results.txt"), quote=F, sep="\t", col.names = TRUE)
        
        ############## ALDEx2 #####
        
        print(treatment)
        print(other_treatment)
        
        otu_table <- as.matrix(otu_table(build_OTU_counts_output))
        
        group_info <- sample_data(build_OTU_counts_output)$Age_Group
        
        # Filter OTU table and group info for the specified conditions
        condition_indices <- which(group_info %in% c(treatment, other_treatment))
        
        count_treatment <- sum(sample_data(build_OTU_counts_output)$Age_Group == treatment)
        
        count_other_treatment <- sum(sample_data(build_OTU_counts_output)$Age_Group == other_treatment)
        
        filtered_otu <- otu_table[, condition_indices]
        
        #print(filtered_otu)
        
        filtered_groups <- group_info[condition_indices]
        
        #     conds <-  c(rep(treatment, count_treatment), rep(other_treatment, count_other_treatment))
        
        print(group_info)
        
        print(filtered_groups) 
        
        ALDEx2_OTU_DE_results <- aldex(reads = filtered_otu,  conditions = filtered_groups, mc.samples = 128, denom = "all", verbose = TRUE, useMC = TRUE, cores = 28, test="t", effect=TRUE, paired.test=FALSE)
        
        comparison_results$ALDEx2 <- ALDEx2_OTU_DE_results
        
        #print(ALDEx2_OTU_DE_results$wi.eBH < 0.05)
        
        write.table(ALDEx2_OTU_DE_results, file=paste0(treatment,"_vs_",other_treatment,"_aldex_DE_results.txt"), quote=FALSE, sep='\t', col.names = NA)
        
        ######## ADAPT ########
        print("Running ADAPT Analysis...")
        
        # Convert phyloseq object to a suitable format for ADAPT
        #   adapt_input <- phyloseq_to_edgeR(build_OTU_subset, group = "Age_Group")
        
        # Run ADAPT (modify parameters as needed)
        adapt_results <- ADAPT::adapt(build_OTU_subset, cond.var = "Age_Group")
        
        DAtaxa_result <- ADAPT::summary(adapt_results, select="all")
        
        # Store results
        comparison_results$ADAPT <- DAtaxa_result
        
        # Save results to file
        write.table(DAtaxa_result, file=paste0(treatment,"_vs_",other_treatment,"_ADAPT_DE_results.txt"), 
                    quote=FALSE, sep='\t', col.names=NA)
        
        ######## metagenomeSeq ########
        print("Running metagenomeSeq Analysis...")
        
        # Convert phyloseq object to metagenomeSeq MRexperiment format
        otu_mat <- as(otu_table(build_OTU_counts_output), "matrix")
        sample_data_df <- as(sample_data(build_OTU_counts_output), "data.frame")
        
        # Create a new MRexperiment object
        pheno_data <- AnnotatedDataFrame(sample_data_df)
        feature_data <- AnnotatedDataFrame(data.frame(OTU_ID = rownames(otu_mat), row.names = rownames(otu_mat)))
        meta_obj <- newMRexperiment(otu_mat, phenoData = pheno_data, featureData = feature_data)
        
        # Normalize data using cumulative sum scaling (CSS)
        meta_obj <- cumNorm(meta_obj, p = cumNormStat(meta_obj))
        
        # Define model matrix for differential analysis
        mod <- model.matrix(~ Age_Group, data = pData(pheno_data))
        
        # Fit metagenomeSeq model
        fit_meta <- fitFeatureModel(meta_obj, mod)
        
        # Extract results
        #metagenomeSeq_OTU_DE_results <- MRcoefs(fit_meta)
        
        metagenomeSeq_OTU_DE_results <- MRfulltable(fit_meta, number = nrow(otu_mat))
        
        # Store results in the comparison list
        comparison_results$metagenomeSeq <- metagenomeSeq_OTU_DE_results
        
        # Save results to file
        write.table(metagenomeSeq_OTU_DE_results, file=paste0(treatment,"_vs_",other_treatment,"_metagenomeSeq_DE_results.txt"), 
                    quote=FALSE, sep='\t', col.names=NA)
        
        
        
        # Fit metagenomeSeq model
       # fit_meta <- fitFeatureModel(meta_obj, mod)
        
        # Extract results
      #  metagenomeSeq_OTU_DE_results <- MRcoefs(fit_meta)
        
      #  metagenomeSeq_OTU_DE_results <- MRfulltable(fit_meta)
        
        
        #############################
        
        all_comparisons_results[[contrast_pair]] <- comparison_results
        
      }
    }
  }
  
  return(all_comparisons_results)

}

# Assumes `OTUs_multi_DA` function and all methods are loaded
# This script adds the logic for extracting significant IDs and plotting

run_and_plot_DE_analysis <- function(build_OTU_counts_output, pval_threshold = 0.05) {
  de_results <- OTUs_multi_DE(build_OTU_counts_output)
  
  for (contrast_pair in names(de_results)) {
    res <- de_results[[contrast_pair]]
    
    # Extract significant IDs (adjusted p-value < threshold)
    significant_results_edgeR <- subset(res$EdgeR, FDR < pval_threshold)
    significant_results_DESeq2 <- subset(res$DESeq2, padj < pval_threshold)
    significant_results_ALDEx2 <- subset(res$ALDEx2, wi.eBH < pval_threshold)
    significant_results_ADAPT <- subset(res$ADAPT, padj < pval_threshold)
    significant_results_metagenomeSeq <- subset(res$metagenomeSeq, adjPvalues < pval_threshold)
    
    print(length(significant_results_edgeR$ID))
    print(length(significant_results_DESeq2$ID))
    print(length(significant_results_ALDEx2$ID)) 
    print(length(significant_results_ADAPT$ID)) 
    print(length(significant_results_metagenomeSeq$ID)) 
    
    edgeR_ids <- significant_results_edgeR$ID
    DESeq2_ids <- significant_results_DESeq2$ID
    ALDEx2_ids <- significant_results_ALDEx2$ID
    ADAPT_ids <- significant_results_ADAPT$ID
    metagenomeSeq_ids <- significant_results_metagenomeSeq$ID
    
    list_input <- list(edgeR = edgeR_ids, DESeq2 = DESeq2_ids, ALDEx2 = ALDEx2_ids, ADAPT = ADAPT_ids, metagenomeSeq = metagenomeSeq_ids)
    
    venn1 <- euler::euler(list_input)
    treatments <- strsplit(contrast_pair, "_vs_")[[1]]
    treatment <- treatments[1]
    other_treatment <- treatments[2]
    
    png(paste0(treatment, "_vs_", other_treatment, "_ProportionalVennDiagram.png"))
    plot(venn1, fills = c("orange", "green", "purple", "pink", "lightblue"), col = "transparent", labels = TRUE, quantities = TRUE)
    dev.off()
    
    venn1_plot <- plot(venn1, fills = c("orange", "green", "purple", "pink", "lightblue"), col = "transparent", labels = TRUE, quantities = TRUE)
    ggsave(filename = paste0(treatment, "_vs_", other_treatment, "_ProportionalVennDiagram.pdf"), plot = venn1_plot, device = "pdf")
    
    upset_plot <- UpSetR::upset(UpSetR::fromList(list_input),  mainbar.y.label = "Adj. Pvalue DE Genes", main.bar.color = "brown", sets.x.label = "DE Gene Counts", sets.bar.color = "red")
    pdf(file=paste0(treatment,"_vs_",other_treatment,"_UpSet_plot.pdf"))
    print(upset_plot)
    dev.off()
  }
}

