
#' Consolidate Differential Abundance Results
#' @param edgeR_res EdgeR results
#' @param DESeq2_res DESeq2 results
#' @param ALDEx2_res ALDEx2 results
#' @param metagenomeSeq_res metagenomeSeq results
#' @param ADAPT_res ADAPT results
#' @export

consolidate_DA_results <- function(edgeR_res = NULL,
                                   DESeq2_res = NULL,
                                   ALDEx2_res = NULL,
                                   metagenomeSeq_res = NULL,
                                   ADAPT_res = NULL) {

  res_list <- list()

  # edgeR
  if (!is.null(edgeR_res)) {
    edgeR_df <- data.frame(
      Taxa = rownames(edgeR_res),
      logFC = edgeR_res$logFC,
      pval = edgeR_res$PValue,
      padj = edgeR_res$FDR
    ) %>% rename_with(~ paste0(., "_edgeR"), -Taxa)
    res_list <- append(res_list, list(edgeR_df))
  }

  # DESeq2
  if (!is.null(DESeq2_res)) {

    DESeq2_df <- data.frame(
      Taxa = DESeq2_res$Taxa,
      logFC = DESeq2_res$log2FoldChange,
      pval = DESeq2_res$pvalue,
      padj = DESeq2_res$padj
    ) %>% rename_with(~ paste0(., "_DESeq2"), -Taxa)
    res_list <- append(res_list, list(DESeq2_df))
  }

  # ALDEx2
  if (!is.null(ALDEx2_res)) {
    ALDEx2_df <- data.frame(
      Taxa = rownames(ALDEx2_res),
      logFC = ALDEx2_res$effect,
      pval = ALDEx2_res$we.ep,
      padj = ALDEx2_res$we.eBH
    ) %>% rename_with(~ paste0(., "_ALDEx2"), -Taxa)
    res_list <- append(res_list, list(ALDEx2_df))
  }

  # metagenomeSeq
  if (!is.null(metagenomeSeq_res)) {
    meta_df <- data.frame(
      Taxa = rownames(metagenomeSeq_res),
      logFC = metagenomeSeq_res$logFC,
      pval = metagenomeSeq_res$pvalues,
      padj = metagenomeSeq_res$adjPvalues
    ) %>% rename_with(~ paste0(., "_metaSeq"), -Taxa)
    res_list <- append(res_list, list(meta_df))
  }

  # ADAPT
  if (!is.null(ADAPT_res)) {
    adapt_df <- data.frame(
      Taxa = ADAPT_res$Taxa,
      logFC = ADAPT_res$log10foldchange,
      pval = ADAPT_res$pval,
      padj = ADAPT_res$adjusted_pval
    ) %>% rename_with(~ paste0(., "_ADAPT"), -Taxa)
    res_list <- append(res_list, list(adapt_df))
  }

  # Merge all
 # merged_res <- reduce(res_list, full_join, by = "Taxa")
  # # Extract each data frame from the list
  # edgeR_df2 <- res_list[[1]]
  # DESeq2_df2 <- res_list[[2]]
  # ALDEx2_df2 <- res_list[[3]]
  # metaSeq_df2 <- res_list[[4]]
  # ADAPT_df <- res_list[[5]]

 # print(ADAPT_df)

  # Now properly merge them by Taxa
    merged_res <- edgeR_df %>%
    full_join(DESeq2_df, by = "Taxa") %>%
    full_join(ALDEx2_df, by = "Taxa") %>%
    full_join(meta_df, by = "Taxa") %>%
    full_join(adapt_df, by = "Taxa")


  return(merged_res)

}


################# R Version 4.4 #############
phyloseq_to_edgeR <-  function(physeq, group, method="RLE", ...){
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


###################################################################################
###############################  OTUs_multi_DA ####################################
###################################################################################

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
#' @import phyloseq
#' @import ggplot2
#' @import dplyr
#' @importFrom VennDiagram venn.diagram
#' @importFrom grid grid.draw
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom vegan specnumber rarefy
#' @importFrom utils write.table
#' @importFrom grDevices dev.new dev.off pdf
#' @importFrom graphics abline
#' @importFrom stats reorder sd setNames var
#'
#'
#' @param build_OTU_counts_output (Required).  A \code{\link[phyloseq]{phyloseq-class}} or
#'  an \code{\link[phyloseq]{otu_table-class}} object.
#'  The latter is only appropriate if \code{group} argument and taxa table.
#'@param force_build Logical. Whether to force rebuilding even if processed
#'   data exists. Currently not implemented. Default: FALSE
#' @param verbose Verbosity ON/OFF. Default=FALSE
#'
#' @export

OTUs_multi_DA <- function(build_OTU_counts_output,
                          force_build = FALSE,
                          verbose = FALSE){

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

        #print(treatment)
        #print(other_treatment)
        print("Running EdgeR Analysis...")


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
        #print(treatment)
        #print(other_treatment)
        print("Running DESeq2 Analysis...")


        phylo_reads_collapsed_deseq <- phyloseq_to_deseq2(build_OTU_counts_output, ~Age_Group)


        #         phylo_reads_collapsed_deseq <- DESeq(phylo_reads_collapsed_deseq, test="Wald",fitType="parametric")
        phylo_reads_collapsed_deseq <- DESeq(phylo_reads_collapsed_deseq, sfType = "poscounts")

        DESeq2_OTU_DE_results = results(phylo_reads_collapsed_deseq , contrast=c("Age_Group", treatment, other_treatment), tidy=T, format="DataFrame")

        colnames(DESeq2_OTU_DE_results)[colnames(DESeq2_OTU_DE_results) == "row"] <- "Taxa"

        comparison_results$DESeq2 <- DESeq2_OTU_DE_results

        write.table(DESeq2_OTU_DE_results, file=paste0(treatment,"_vs_",other_treatment,"_DESeq_DE_results.txt"), quote=F, sep="\t", col.names = TRUE)

        ############## ALDEx2 #####

        #print(treatment)
        #print(other_treatment)
        print("Running ALDEx2 Analysis...")

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

        ##################### DA consolidation ###############


        OTUs_multi_DA_consolidated <-  consolidate_DA_results( edgeR_res = EdgeR_OTU_DE_results,
                                                               DESeq2_res = DESeq2_OTU_DE_results,
                                                               ALDEx2_res = ALDEx2_OTU_DE_results,
                                                               metagenomeSeq_res = metagenomeSeq_OTU_DE_results,
                                                               ADAPT_res = DAtaxa_result
        )

        print(OTUs_multi_DA_consolidated)

        #OTUs_multi_DA_consolidated
        write.table(OTUs_multi_DA_consolidated, file=paste0(treatment,"_vs_",other_treatment,"_ALLDE_results.txt"), quote=F, sep="\t", row.names = FALSE, col.names = TRUE)

        ############################ Venn Diagram Plots ##################

        pval_threshold = 0.05

        significant_results_edgeR <- subset(OTUs_multi_DA_consolidated, padj_edgeR < pval_threshold)
        significant_results_DESeq2 <- subset(OTUs_multi_DA_consolidated, padj_DESeq2 < pval_threshold)
        significant_results_ALDEx2 <- subset(OTUs_multi_DA_consolidated, padj_ALDEx2 < pval_threshold)
        significant_results_ADAPT <- subset(OTUs_multi_DA_consolidated, padj_metaSeq < pval_threshold)
        significant_results_metagenomeSeq <- subset(OTUs_multi_DA_consolidated, padj_ADAPT < pval_threshold)

        print(length(significant_results_edgeR$Taxa))
        print(length(significant_results_DESeq2$Taxa))
        print(length(significant_results_ALDEx2$Taxa))
        print(length(significant_results_ADAPT$Taxa))
        print(length(significant_results_metagenomeSeq$Taxa))

        edgeR_ids <- significant_results_edgeR$Taxa
        DESeq2_ids <- significant_results_DESeq2$Taxa
        ALDEx2_ids <- significant_results_ALDEx2$Taxa
        ADAPT_ids <- significant_results_ADAPT$Taxa
        metagenomeSeq_ids <- significant_results_metagenomeSeq$Taxa

        list_input <- list(edgeR = edgeR_ids, DESeq2 = DESeq2_ids, ALDEx2 = ALDEx2_ids, ADAPT = ADAPT_ids, metagenomeSeq = metagenomeSeq_ids)

        print(list_input)
        #treatments <- strsplit(contrast_pair, "_vs_")[[1]]
        #treatment <- treatments[1]
        #other_treatment <- treatments[2]


        venn1_plot <- venn.diagram(
          x = list_input,
          category.names = names(list_input),
          filename = NULL,
          output = TRUE,
          fill = c("red", "blue", "green", "yellow", "purple"),
          alpha = 0.1,
          cex = 1.5,
          fontfamily = "serif",
          fontface = "bold",
          cat.cex = 1.5,
          cat.fontface = "bold",
          margin = 0.2
        )
        dev.new()  # or plot.new()

        comparison_results$VennDiagram <- venn1_plot

        grid.draw(venn1_plot)

        ggsave(filename = paste0(treatment, "_vs_", other_treatment, "_ProportionalVennDiagram.pdf"), plot = venn1_plot, device = "pdf", height = 10, width = 8)

        ############################ UpSetR Plots ##################

        upset_plot <- UpSetR::upset(UpSetR::fromList(list_input),  mainbar.y.label = "Adj. Pvalue DE Genes", main.bar.color = "brown", sets.x.label = "DE Gene Counts", sets.bar.color = "red")
        print(upset_plot)

        comparison_results$upset_plot <- upset_plot

        pdf(file=paste0(treatment,"_vs_",other_treatment,"_UpSet_plot.pdf"))
        #print(upset_plot)
        dev.off()


        #############################

        all_comparisons_results[[contrast_pair]] <- comparison_results

      }
    }

  }
  return(all_comparisons_results)
}



