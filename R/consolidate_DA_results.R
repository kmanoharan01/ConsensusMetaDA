## =============================================================================
## R/consolidate_DA_results.R
## =============================================================================

#' Consolidate per-tool differential abundance results into one table
#'
#' The single, correct definition. Merges each tool's Taxa/logFC/pval/padj
#' into one wide table (suffixed `_<tool>`), flags significance per tool at
#' `fdr_cut`, and adds `n_sig` (how many tools called each taxon significant)
#' and `n_tools` (how many tools ran) — the two columns every consensus
#' function downstream depends on.
#'
#' New in v2: when `add_combined = TRUE` (the default), the three p-value
#' combination methods (Fisher, Stouffer, dependence-robust Cauchy) are
#' appended directly, so the consolidated table — and therefore the
#' `*_ALLDE_results.txt` file written by [OTUs_multi_DA()] — carries
#' `p_<method>`, `padj_<method>`, and `sig_<method>` for each method without a
#' separate call. Set `add_combined = FALSE` to reproduce the v1 output.
#'
#' @param edgeR_res,DESeq2_res,ALDEx2_res,metagenomeSeq_res,ADAPT_res,
#'   ANCOMBC2_res,MaAsLin3_res Per-tool raw result data.frames, or NULL to
#'   omit a tool. Pass whatever OTUs_multi_DA() produced per tool.
#' @param fdr_cut Numeric, adjusted p-value threshold defining significance.
#' @param add_combined Logical. If TRUE (default) append Fisher/Stouffer/Cauchy
#'   combined p-values, their BH-adjusted values, and significance flags.
#' @param combine_methods Character vector, subset of
#'   c("fisher","stouffer","cauchy"). Which combination methods to append when
#'   `add_combined = TRUE`.
#' @return data.frame, one row per taxon, or NULL if no tool results supplied.
#' @export
consolidate_DA_results <- function(edgeR_res         = NULL,
                                   DESeq2_res        = NULL,
                                   ALDEx2_res        = NULL,
                                   metagenomeSeq_res = NULL,
                                   ADAPT_res         = NULL,
                                   ANCOMBC2_res      = NULL,
                                   MaAsLin3_res      = NULL,
                                   fdr_cut           = 0.05,
                                   add_combined      = TRUE,
                                   combine_methods   = c("fisher", "stouffer", "cauchy")) {

  res_list <- list()
  add <- function(df, tool) {
    df <- dplyr::mutate(df, sig = !is.na(padj) & padj < fdr_cut)
    df <- dplyr::rename_with(df, ~ paste0(., "_", tool), -Taxa)
    res_list[[length(res_list) + 1]] <<- df
  }

  if (!is.null(edgeR_res))
    add(data.frame(Taxa = rownames(edgeR_res), logFC = edgeR_res$logFC,
                   pval = edgeR_res$PValue, padj = edgeR_res$FDR), "edgeR")

  if (!is.null(DESeq2_res))
    add(data.frame(Taxa = DESeq2_res$Taxa, logFC = DESeq2_res$log2FoldChange,
                   pval = DESeq2_res$pvalue, padj = DESeq2_res$padj), "DESeq2")

  ## ALDEx2's `effect` is a standardized effect size (median group difference
  ## / max within-group dispersion), NOT log2(fold-change) — kept here for
  ## reference/transparency but never folded into a composite fold-change.
  if (!is.null(ALDEx2_res))
    add(data.frame(Taxa = rownames(ALDEx2_res), logFC = ALDEx2_res$effect,
                   pval = ALDEx2_res$we.ep, padj = ALDEx2_res$we.eBH), "ALDEx2")

  if (!is.null(metagenomeSeq_res))
    add(data.frame(Taxa = rownames(metagenomeSeq_res), logFC = metagenomeSeq_res$logFC,
                   pval = metagenomeSeq_res$pvalues,
                   padj = metagenomeSeq_res$adjPvalues), "metaSeq")

  ## FIX: ADAPT reports log10 fold-change; rescale to log2 (* log2(10) ==
  ## * 3.32193) so it is on the same footing as edgeR/DESeq2.
  if (!is.null(ADAPT_res))
    add(data.frame(Taxa = ADAPT_res$Taxa,
                   logFC = ADAPT_res$log10foldchange * log2(10),
                   pval = ADAPT_res$pval, padj = ADAPT_res$adjusted_pval), "ADAPT")

  if (!is.null(ANCOMBC2_res))
    add(data.frame(Taxa = ANCOMBC2_res$taxon, logFC = ANCOMBC2_res$lfc,
                   pval = ANCOMBC2_res$p_val, padj = ANCOMBC2_res$q_val), "ANCOMBC2")

  if (!is.null(MaAsLin3_res))
    add(data.frame(Taxa = MaAsLin3_res$feature, logFC = MaAsLin3_res$coef,
                   pval = MaAsLin3_res$pval, padj = MaAsLin3_res$qval), "MaAsLin3")

  if (!length(res_list)) return(NULL)

  merged <- Reduce(function(a, b) dplyr::full_join(a, b, by = "Taxa"), res_list)
  sig_cols <- grep("^sig_", names(merged), value = TRUE)
  merged$n_sig   <- rowSums(merged[, sig_cols, drop = FALSE], na.rm = TRUE)
  merged$n_tools <- length(sig_cols)

  ## v2: append combined p-values so the consolidated table (and the ALLDE
  ## file) carries them. Needs >= 2 pval_ columns; silently skipped otherwise.
  if (isTRUE(add_combined)) {
    pcols <- grep("^pval_", names(merged), value = TRUE)
    if (length(pcols) >= 2) {
      merged <- add_combined_pvalues(merged, methods = combine_methods,
                                     fdr_cut = fdr_cut, verbose = FALSE)
    }
  }

  merged
}
