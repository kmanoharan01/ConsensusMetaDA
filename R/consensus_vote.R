## =============================================================================
## R/consensus_vote.R
## =============================================================================

#' Build vote-rule consensus calls from a consolidated DA table
#'
#' No closed-form FDR/power guarantee exists for union/intersect/min-support
#' rules — see estimate_fdr_by_permutation() for per-dataset calibration.
#'
#' @param consolidated Output of consolidate_DA_results() — must carry
#'   `sig_<tool>` columns and `n_sig`/`n_tools`.
#' @param truth Character vector of true-positive taxa (simulation only).
#' @return list(table = consolidated + vote_1..vote_K logical columns,
#'   summary = data.frame of taxa called per rule, with precision/recall
#'   if `truth` supplied).
#' @export
consensus_vote_rules <- function(consolidated, truth = NULL) {
  stopifnot(is.data.frame(consolidated), "n_sig" %in% names(consolidated))
  k <- consolidated$n_tools[1]
  if (is.na(k) || k < 1) stop("n_tools missing or zero - was fdr_cut passed to consolidate_DA_results()?")

  n_sig <- consolidated$n_sig
  n_sig[is.na(n_sig)] <- 0

  for (v in seq_len(k)) consolidated[[paste0("vote_", v)]] <- n_sig >= v

  rule_name <- function(v)
    if (v == 1) sprintf("Union (1 of %d)", k) else
      if (v == k) sprintf("Intersect (%d of %d)", k, k) else sprintf("%d of %d", v, k)

  summ <- data.frame(rule = vapply(seq_len(k), rule_name, character(1)), votes = seq_len(k),
                     n_taxa = vapply(seq_len(k), function(v) sum(n_sig >= v), integer(1)),
                     stringsAsFactors = FALSE)

  if (!is.null(truth)) {
    is_true <- consolidated$Taxa %in% truth
    summ$precision <- vapply(seq_len(k), function(v) {
      called <- n_sig >= v
      if (!any(called)) return(NA_real_)
      mean(is_true[called])
    }, numeric(1))
    summ$recall <- vapply(seq_len(k), function(v)
      sum(is_true & n_sig >= v) / max(1, sum(is_true)), numeric(1))
  }

  list(table = consolidated, summary = summ)
}

#' Empirical FDR of vote rules via permutation
#'
#' Estimates each vote rule's false discovery proportion for THIS dataset by
#' permuting group labels (null true by construction), rerunning the full
#' pipeline, and comparing null call counts to the observed. Expensive:
#' reruns all `tools` x `n_perm` times.
#'
#' @param physeq,group_var,covariate,group_order,fdr_cut,tools,cores As in
#'   OTUs_multi_DA().
#' @param n_perm Integer, number of permutations.
#' @param seed Integer.
#' @param outdir As in OTUs_multi_DA(); permutation runs use subdirectories
#'   and skip plots.
#' @return data.frame, one row per vote rule per contrast.
#' @export
estimate_fdr_by_permutation <- function(physeq, group_var = "Group", covariate = NULL,
                                        group_order = NULL, fdr_cut = 0.05,
                                        tools = c("edgeR","DESeq2","ALDEx2","metaSeq",
                                                  "ADAPT","ANCOMBC2","MaAsLin3"),
                                        cores = 1L, n_perm = 10L, seed = 1L,
                                        outdir = tempdir()) {

  message("Observed run...")
  obs <- OTUs_multi_DA(physeq, group_var = group_var, covariate = covariate,
                       group_order = group_order, fdr_cut = fdr_cut, tools = tools,
                       cores = cores, make_plots = FALSE, add_combined = FALSE,
                       outdir = file.path(outdir, "observed"))

  meta <- methods::as(phyloseq::sample_data(physeq), "data.frame")
  set.seed(seed)
  perm_rows <- list()

  for (contrast in names(obs)) {
    obs_cons  <- obs[[contrast]]$consolidated
    obs_votes <- consensus_vote_rules(obs_cons)$summary
    k <- obs_cons$n_tools[1]
    null_counts <- matrix(NA_integer_, nrow = n_perm, ncol = k)

    for (p in seq_len(n_perm)) {
      message(sprintf("  [%s] permutation %d/%d", contrast, p, n_perm))
      perm_meta <- meta
      perm_meta[[group_var]] <- sample(perm_meta[[group_var]])
      perm_physeq <- physeq
      phyloseq::sample_data(perm_physeq) <-
        phyloseq::sample_data(perm_meta[rownames(meta), , drop = FALSE])

      perm_out <- tryCatch(
        OTUs_multi_DA(perm_physeq, group_var = group_var, covariate = covariate,
                      group_order = group_order, fdr_cut = fdr_cut, tools = tools,
                      cores = cores, make_plots = FALSE, add_combined = FALSE,
                      outdir = file.path(outdir, sprintf("perm_%d", p))),
        error = function(e) { warning("permutation ", p, " failed: ", conditionMessage(e)); NULL })

      if (!is.null(perm_out) && contrast %in% names(perm_out)) {
        pv <- consensus_vote_rules(perm_out[[contrast]]$consolidated)$summary
        null_counts[p, ] <- pv$n_taxa[match(seq_len(k), pv$votes)]
      }
    }

    for (v in seq_len(k)) {
      obs_n  <- obs_votes$n_taxa[obs_votes$votes == v]
      null_v <- null_counts[, v]; null_v <- null_v[!is.na(null_v)]
      perm_rows[[length(perm_rows) + 1]] <- data.frame(
        contrast = contrast, rule = obs_votes$rule[obs_votes$votes == v], votes = v,
        n_perm_valid = length(null_v), observed_n = obs_n,
        null_mean = if (length(null_v)) mean(null_v) else NA_real_,
        null_sd   = if (length(null_v)) stats::sd(null_v) else NA_real_,
        est_FDP   = if (obs_n > 0 && length(null_v)) min(1, mean(null_v) / obs_n) else NA_real_,
        p_empirical = if (length(null_v)) (sum(null_v >= obs_n) + 1) / (length(null_v) + 1) else NA_real_,
        stringsAsFactors = FALSE)
    }
  }

  res <- do.call(rbind, perm_rows)
  res <- res[order(res$contrast, res$votes), ]
  utils::write.table(res, file.path(outdir, "permutation_FDR_estimate.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE)
  res
}

#' Recommend a vote threshold from a permutation FDR table
#'
#' @param perm_result Output of estimate_fdr_by_permutation().
#' @param target_fdr Numeric, default 0.05.
#' @return data.frame, one row per contrast: recommended vote count and its
#'   estimated FDP. If none meets target, recommends Intersect and flags it.
#' @export
recommend_vote_threshold <- function(perm_result, target_fdr = 0.05) {
  do.call(rbind, lapply(split(perm_result, perm_result$contrast), function(d) {
    ok <- d[!is.na(d$est_FDP) & d$est_FDP <= target_fdr, ]
    if (nrow(ok)) {
      best <- ok[which.min(ok$votes), ]
      data.frame(contrast = d$contrast[1], recommended_votes = best$votes, rule = best$rule,
                 est_FDP = best$est_FDP, note = "meets target", stringsAsFactors = FALSE)
    } else {
      worst <- d[which.max(d$votes), ]
      data.frame(contrast = d$contrast[1], recommended_votes = worst$votes, rule = worst$rule,
                 est_FDP = worst$est_FDP,
                 note = sprintf("NO rule met target_fdr=%.2f; even Intersect estimated at %.3f - increase n_perm or treat all calls with caution",
                                target_fdr, worst$est_FDP),
                 stringsAsFactors = FALSE)
    }
  }))
}
