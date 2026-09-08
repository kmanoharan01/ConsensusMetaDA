## =============================================================================
## R/build_summary_table.R
## =============================================================================

#' Compact consensus summary table
#'
#' Per-tool adjusted p-value, per-tool rank, rank_sum, and the order-statistic
#' family p_union / p_atleast2 / p_intersect_all — no per-tool sig_/vote_
#' boolean columns, no logFC composite.
#'
#' @section Why there is no LogFC/LogFC_sd column:
#' p-values are comparable across tools by construction (all on [0,1]); logFC
#' is NOT (ALDEx2's `effect` is a standardized effect size, and
#' ANCOMBC2/MaAsLin3/metaSeq each carry their own normalization with no
#' confirmed constant rescaling). Read `logFC_<tool>` directly from
#' `consolidated` if needed.
#'
#' @param consolidated Output of consolidate_DA_results(fdr_cut = ...) — needs
#'   `pval_<tool>`, `padj_<tool>` per tool.
#' @param rank_on "padj" (default) or "pval".
#' @return data.frame, one row per taxon, sorted by rank_sum ascending, with a
#'   hidden `padj_matrix` attribute used by p_atleast_k().
#' @export
build_summary_table <- function(consolidated, rank_on = c("padj", "pval")) {
  rank_on <- match.arg(rank_on)
  ## per-TOOL columns only: exclude combined/ranksum padj columns
  all_pval <- grep("^pval_", names(consolidated), value = TRUE)
  all_padj <- grep("^padj_", names(consolidated), value = TRUE)
  combined <- c("padj_fisher","padj_stouffer","padj_cauchy","padj_combined","padj_ranksum",
                "pval_fisher","pval_stouffer","pval_cauchy","pval_combined")
  pval_cols <- setdiff(all_pval, combined)
  padj_cols <- setdiff(all_padj, combined)
  tools <- sub("^padj_", "", padj_cols)
  if (!length(tools)) stop("No per-tool padj_<tool> columns found in `consolidated`")

  out <- data.frame(Taxa = consolidated$Taxa, stringsAsFactors = FALSE)

  for (tl in tools) {
    pc <- paste0("padj_", tl)
    out[[paste0(tl, "_adj_p")]] <- if (pc %in% names(consolidated)) consolidated[[pc]] else NA_real_
  }

  rank_src <- if (rank_on == "padj") padj_cols else pval_cols
  R <- sapply(rank_src, function(cl) {
    x <- consolidated[[cl]]
    r <- rep(NA_real_, length(x)); ok <- !is.na(x)
    r[ok] <- rank(x[ok], ties.method = "average")
    r[!ok] <- sum(ok) + 1
    r
  })
  for (i in seq_along(tools)) out[[paste0(tools[i], "_rank")]] <- R[, i]
  out$rank_sum <- rowSums(R, na.rm = TRUE)

  P <- as.matrix(consolidated[, padj_cols, drop = FALSE])
  kth_smallest <- function(x, k) { x <- sort(x[!is.na(x)]); if (length(x) < k) NA_real_ else x[k] }
  out$p_union         <- apply(P, 1, kth_smallest, k = 1)
  out$p_atleast2      <- apply(P, 1, kth_smallest, k = 2)
  out$p_intersect_all <- apply(P, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))

  ord <- order(out$rank_sum)
  out <- out[ord, ]
  P   <- P[ord, , drop = FALSE]
  attr(out, "padj_matrix") <- P
  out
}

#' p-value for "at least k tools agree", for any k
#'
#' @param summary_tab Output of build_summary_table().
#' @param k Integer.
#' @return Numeric vector, same length/order as summary_tab's rows.
#' @export
p_atleast_k <- function(summary_tab, k) {
  P <- attr(summary_tab, "padj_matrix")
  if (is.null(P)) {
    adj_cols <- grep("_adj_p$", names(summary_tab), value = TRUE)
    if (!length(adj_cols))
      stop("No padj_matrix attribute and no <tool>_adj_p columns found - ",
           "rebuild with build_summary_table() on the original consolidated table")
    message("padj_matrix attribute not present - reconstructing from ", length(adj_cols),
            " <tool>_adj_p column(s)")
    P <- as.matrix(summary_tab[, adj_cols, drop = FALSE])
  }
  apply(P, 1, function(x) { x <- sort(x[!is.na(x)]); if (length(x) < k) NA_real_ else x[k] })
}

#' Column dictionary for build_summary_table() output
#' @export
describe_summary_table <- function() {
  data.frame(
    column = c("<tool>_adj_p", "<tool>_rank", "rank_sum", "p_union", "p_atleast2", "p_intersect_all"),
    description = c("Per-tool p-value", "Rank of the p-value reported by that tool", "Sum of ranks",
                    "Smallest adjusted p-value observed", "2nd-smallest adjusted p-value observed",
                    "Largest adjusted p-value observed"),
    detail = c(
      "BH-adjusted p-value from that tool",
      "Smallest rank = most significant reported by that tool",
      "Sum of every <tool>_rank column. Table is sorted by this",
      "Threshold for the UNION: p_union < cutoff means >=1 tool called it",
      "Threshold for AT LEAST 2 tools. Use p_atleast_k(summary_tab, k) for other k",
      "Threshold for the INTERSECT: every tool called it"),
    stringsAsFactors = FALSE)
}
