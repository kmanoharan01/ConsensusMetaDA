## =============================================================================
## R/print_consensus_summary.R
## =============================================================================

#' Print a plain-text consensus summary for one contrast
#'
#' v2: if the consolidated table carries combined-p columns
#' (`sig_fisher` / `sig_stouffer` / `sig_cauchy`), their significant-taxa
#' counts are printed too, so the console summary reflects the same columns
#' now written to the ALLDE file.
#'
#' @param consolidated Output of consolidate_DA_results(fdr_cut = ...).
#' @param contrast_pair Character, label for the printed header.
#' @param fdr_cut Numeric, threshold used (for the header only — significance
#'   was already applied in consolidate_DA_results()).
#' @return Invisibly, a list of the summary counts.
#' @export
print_consensus_summary <- function(consolidated, contrast_pair, fdr_cut = 0.05) {
  if (is.null(consolidated) || nrow(consolidated) == 0) {
    cat(sprintf("\n  No consolidated results for %s\n", contrast_pair))
    return(invisible(NULL))
  }

  n_tools <- if ("n_tools" %in% names(consolidated)) consolidated$n_tools[1] else
    length(grep("^sig_[A-Za-z]", names(consolidated)))

  n_sig <- consolidated$n_sig
  n_sig[is.na(n_sig)] <- 0
  total_taxa <- nrow(consolidated)

  cat(sprintf("\n%s\n", strrep("-", 50)))
  cat(sprintf("  Consensus summary - %s  (FDR < %g)\n", contrast_pair, fdr_cut))
  cat(sprintf("%s\n", strrep("-", 50)))
  cat(sprintf("  Tools that ran      : %d\n", n_tools))
  cat(sprintf("  Total taxa tested   : %d\n", total_taxa))

  cat("\n  Taxa significant by N tools:\n")
  for (k in seq_len(n_tools))
    cat(sprintf("    %d/%d tools : %d\n", k, n_tools, sum(n_sig == k)))

  n_any       <- sum(n_sig >= 1)
  n_majority  <- sum(n_sig >= ceiling(n_tools / 2))
  n_strict    <- sum(n_sig >= max(1, n_tools - 1))
  n_unanimous <- sum(n_sig == n_tools)

  cat("\n  Consensus tiers (vote rules):\n")
  cat(sprintf("    Any tool (>=1)          : %d\n", n_any))
  cat(sprintf("    Majority (>=%d)          : %d\n", ceiling(n_tools / 2), n_majority))
  cat(sprintf("    Strict (>=%d, all-but-1) : %d\n", max(1, n_tools - 1), n_strict))
  cat(sprintf("    Unanimous (%d/%d)        : %d\n", n_tools, n_tools, n_unanimous))

  ## v2: combined-p methods, if present
  combo_map <- c(fisher = "sig_fisher", stouffer = "sig_stouffer", cauchy = "sig_cauchy")
  present <- combo_map[combo_map %in% names(consolidated)]
  combo_counts <- NULL
  if (length(present)) {
    cat("\n  Combined-p methods:\n")
    combo_counts <- list()
    for (nm in names(present)) {
      cnt <- sum(consolidated[[present[[nm]]]], na.rm = TRUE)
      combo_counts[[nm]] <- cnt
      lbl <- if (nm == "cauchy") "Cauchy (dependence-robust)" else
        paste0(toupper(substring(nm, 1, 1)), substring(nm, 2))
      cat(sprintf("    %-26s: %d\n", lbl, cnt))
    }
  }
  cat(sprintf("%s\n\n", strrep("-", 50)))

  invisible(list(contrast = contrast_pair, n_tools = n_tools, n_any = n_any,
                 n_majority = n_majority, n_strict = n_strict, n_unanimous = n_unanimous,
                 combined = combo_counts))
}
