## =============================================================================
## R/rank_aggregate.R
## =============================================================================

#' Rank-sum aggregation across tools
#'
#' Sums each taxon's within-tool p-value rank across the K tools (rank 1 =
#' most significant). Analytic null assumes independent ranks; same
#' correlated-tools caveat as Fisher/Stouffer, checkable with `permute = TRUE`.
#'
#' @param consolidated Output of consolidate_DA_results() — needs `pval_<tool>`.
#' @param permute Logical, also compute a permutation-calibrated p-value.
#' @param n_perm Integer, permutations if `permute = TRUE`.
#' @param fdr_cut Numeric.
#' @return `consolidated` + `rank_sum`, `p_ranksum`, `padj_ranksum`,
#'   `sig_ranksum`, and `p_ranksum_perm` if `permute`.
#' @export
rank_sum_aggregate <- function(consolidated, permute = FALSE, n_perm = 1000, fdr_cut = 0.05) {
  pcols <- grep("^pval_", names(consolidated), value = TRUE)
  if (length(pcols) < 2) stop("Need >= 2 pval_<tool> columns")

  P <- as.matrix(consolidated[, pcols, drop = FALSE])
  N <- nrow(P)

  R <- apply(P, 2, function(x) {
    r <- rep(NA_real_, length(x)); ok <- !is.na(x)
    r[ok] <- rank(x[ok], ties.method = "average")
    r[!ok] <- sum(ok) + 1
    r
  })

  k_i      <- rowSums(!is.na(P))
  rank_sum <- rowSums(R, na.rm = TRUE)
  mu       <- k_i * (N + 1) / 2
  sigma2   <- k_i * (N^2 - 1) / 12
  z        <- (mu - rank_sum) / sqrt(sigma2)
  p_analytic <- stats::pnorm(z, lower.tail = FALSE)

  consolidated$rank_sum     <- rank_sum
  consolidated$p_ranksum    <- p_analytic
  consolidated$padj_ranksum <- stats::p.adjust(p_analytic, method = "BH")
  consolidated$sig_ranksum  <- !is.na(consolidated$padj_ranksum) & consolidated$padj_ranksum < fdr_cut

  if (permute) {
    message("Permuting rank-sum null (", n_perm, " reps)...")
    perm_null <- matrix(NA_real_, N, n_perm)
    for (b in seq_len(n_perm)) {
      Rb <- apply(R, 2, sample)
      perm_null[, b] <- rowSums(Rb, na.rm = TRUE)
    }
    p_perm <- vapply(seq_len(N), function(i)
      (sum(perm_null[i, ] <= rank_sum[i]) + 1) / (n_perm + 1), numeric(1))
    consolidated$p_ranksum_perm <- p_perm
    infl <- median(p_perm / pmax(p_analytic, 1e-12), na.rm = TRUE)
    message(sprintf("  median(p_permuted / p_analytic) = %.2f - ", infl),
            if (infl > 1.5) "substantial inflation: use p_ranksum_perm, not the analytic version"
            else "modest inflation: analytic version reasonably trustworthy here")
  }
  consolidated
}
