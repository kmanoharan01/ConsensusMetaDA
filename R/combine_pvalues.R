## =============================================================================
## R/combine_pvalues.R
## =============================================================================

#' Combine per-tool p-values into one significance call per taxon
#'
#' One pass over `pval_<tool>` columns — no rerunning tools, unlike
#' estimate_fdr_by_permutation(). Fisher/Stouffer require the tools'
#' p-values to be independent under the null; this dataset's tools are
#' documented to be correlated (e.g. co-firing false positives), which makes
#' Fisher/Stouffer anti-conservative. Prefer `method = "cauchy"`, which is
#' valid under arbitrary dependence (Liu & Xie 2020) and is the only method
#' here with an actual theoretical FDR guarantee under dependence.
#'
#' @param consolidated Output of consolidate_DA_results() — needs `pval_<tool>`.
#' @param method One of "fisher", "stouffer", "cauchy".
#' @param weights Named numeric vector of per-tool weights, Stouffer only.
#' @param fdr_cut Numeric, BH threshold on the combined p-value.
#' @param verbose Logical, print the significance count message.
#' @return `consolidated` + `p_combined`, `padj_combined`, `sig_combined`.
#' @export
combine_pvalues <- function(consolidated, method = c("fisher", "stouffer", "cauchy"),
                            weights = NULL, fdr_cut = 0.05, verbose = TRUE) {
  method <- match.arg(method)
  res <- .combine_one(consolidated, method = method, weights = weights, fdr_cut = fdr_cut)

  consolidated$p_combined    <- res$p
  consolidated$padj_combined <- res$padj
  consolidated$sig_combined  <- res$sig

  if (isTRUE(verbose)) {
    message(sprintf("combine_pvalues(%s): %d / %d taxa significant at FDR < %.2f",
                    method, sum(consolidated$sig_combined, na.rm = TRUE),
                    sum(!is.na(res$p)), fdr_cut))
    if (method %in% c("fisher", "stouffer"))
      message("  NOTE: assumes independent p-values across tools. Validate with ",
              "estimate_fdr_by_permutation() or compare against method = 'cauchy' first.")
  }
  consolidated
}

#' Append several p-value combination methods at once
#'
#' Like [combine_pvalues()] but writes method-suffixed columns
#' (`p_<method>`, `padj_<method>`, `sig_<method>`) for every method in
#' `methods`, so a single consolidated table can carry Fisher, Stouffer, and
#' Cauchy side by side. This is what [consolidate_DA_results()] calls when
#' `add_combined = TRUE`, and what puts the combined columns into the
#' `*_ALLDE_results.txt` file.
#'
#' @param consolidated Output of consolidate_DA_results() — needs `pval_<tool>`.
#' @param methods Character vector, subset of c("fisher","stouffer","cauchy").
#' @param weights Named numeric vector of per-tool weights (Stouffer only).
#' @param fdr_cut Numeric, BH threshold on each combined p-value.
#' @param verbose Logical.
#' @return `consolidated` with `p_<method>`, `padj_<method>`, `sig_<method>`
#'   appended for each requested method.
#' @export
add_combined_pvalues <- function(consolidated,
                                 methods = c("fisher", "stouffer", "cauchy"),
                                 weights = NULL, fdr_cut = 0.05, verbose = FALSE) {
  methods <- match.arg(methods, several.ok = TRUE)
  pcols <- grep("^pval_", names(consolidated), value = TRUE)
  if (length(pcols) < 2) {
    warning("add_combined_pvalues: need >= 2 pval_<tool> columns; found ",
            length(pcols), " - returning table unchanged.")
    return(consolidated)
  }

  for (m in methods) {
    res <- .combine_one(consolidated, method = m, weights = weights, fdr_cut = fdr_cut)
    consolidated[[paste0("p_",    m)]] <- res$p
    consolidated[[paste0("padj_", m)]] <- res$padj
    consolidated[[paste0("sig_",  m)]] <- res$sig
    if (isTRUE(verbose))
      message(sprintf("  %-8s: %d significant at FDR < %.2f",
                      m, sum(res$sig, na.rm = TRUE), fdr_cut))
  }
  consolidated
}

#' @keywords internal
#' @noRd
## Core combiner: returns list(p, padj, sig) for ONE method. Shared by
## combine_pvalues() and add_combined_pvalues() so the maths lives in one place.
.combine_one <- function(consolidated, method, weights = NULL, fdr_cut = 0.05) {
  pcols <- grep("^pval_", names(consolidated), value = TRUE)
  if (length(pcols) < 2)
    stop("Need >= 2 pval_<tool> columns; found: ", paste(pcols, collapse = ", "))

  P <- as.matrix(consolidated[, pcols, drop = FALSE])
  ## clamp away from 0 and 1: Fisher's log() and Cauchy's tan() both blow up at
  ## the boundaries (edgeR routinely reports p ~ 1e-60, and exact zeros/ones
  ## occur). 1e-300 / 1-1e-16 keeps everything finite without materially
  ## changing any p-value that matters.
  P <- pmin(pmax(P, 1e-300), 1 - 1e-16)
  k_i <- rowSums(!is.na(P))

  if (method == "fisher") {
    X <- -2 * rowSums(log(P), na.rm = TRUE)
    p_combined <- ifelse(k_i >= 2, stats::pchisq(X, df = 2 * k_i, lower.tail = FALSE), NA_real_)

  } else if (method == "stouffer") {
    tool_names <- sub("^pval_", "", pcols)
    w <- if (is.null(weights)) rep(1, length(pcols)) else {
      ww <- weights[tool_names]
      if (anyNA(ww)) stop("weights missing for: ", paste(tool_names[is.na(ww)], collapse = ", "))
      ww
    }
    Z  <- t(t(stats::qnorm(1 - P)) * w)
    Zc <- rowSums(Z, na.rm = TRUE) /
      sqrt(rowSums(matrix(w, nrow(P), length(w), byrow = TRUE)^2 * !is.na(P)))
    p_combined <- ifelse(k_i >= 2, stats::pnorm(Zc, lower.tail = FALSE), NA_real_)

  } else { # cauchy / ACAT
    Tc <- rowMeans(tan((0.5 - P) * pi), na.rm = TRUE)
    p_combined <- stats::pcauchy(Tc, lower.tail = FALSE)
  }

  padj <- stats::p.adjust(p_combined, method = "BH")
  list(p = p_combined, padj = padj, sig = !is.na(padj) & padj < fdr_cut)
}
