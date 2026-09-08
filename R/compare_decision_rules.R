## =============================================================================
## R/compare_decision_rules.R
## =============================================================================

#' Compare all decision rules side by side
#'
#' Union, intersection, every min-support threshold, Fisher, Stouffer,
#' Cauchy, and rank-sum aggregation on the same consolidated table.
#'
#' @param consolidated Output of consolidate_DA_results(fdr_cut = ...).
#' @param truth Character vector of true positive taxa (simulation only).
#' @param fdr_cut Numeric.
#' @return data.frame, one row per decision rule.
#' @export
compare_decision_rules <- function(consolidated, truth = NULL, fdr_cut = 0.05) {
  vr <- consensus_vote_rules(consolidated, truth = truth)
  out <- vr$summary
  names(out)[names(out) == "n_taxa"] <- "n_called"
  out$method <- "vote"

  add_row <- function(tab, name, called) {
    row <- data.frame(rule = name, votes = NA, n_called = sum(called, na.rm = TRUE),
                      method = "combination", stringsAsFactors = FALSE)
    if (!is.null(truth)) {
      is_true <- consolidated$Taxa %in% truth
      row$precision <- if (any(called, na.rm = TRUE)) mean(is_true[which(called)]) else NA_real_
      row$recall    <- sum(is_true & called, na.rm = TRUE) / max(1, sum(is_true))
    }
    rbind(tab, row)
  }

  ## use pre-computed sig_ columns if present (v2 consolidated), else compute
  get_sig <- function(cons, method) {
    col <- paste0("sig_", method)
    if (col %in% names(cons)) cons[[col]]
    else combine_pvalues(cons, method = method, fdr_cut = fdr_cut, verbose = FALSE)$sig_combined
  }

  out <- add_row(out, "Fisher combined",                    get_sig(consolidated, "fisher"))
  out <- add_row(out, "Stouffer combined",                  get_sig(consolidated, "stouffer"))
  out <- add_row(out, "Cauchy combined (dependence-robust)", get_sig(consolidated, "cauchy"))

  d_r <- rank_sum_aggregate(consolidated, fdr_cut = fdr_cut)
  out <- add_row(out, "Rank-sum aggregation", d_r$sig_ranksum)

  out[, c("rule", "method", "votes", "n_called", intersect(c("precision", "recall"), names(out)))]
}
