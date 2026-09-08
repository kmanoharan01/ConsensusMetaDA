## =============================================================================
## R/reproducibility.R
## =============================================================================

#' Reproducibility of a decision rule across independent cohorts
#'
#' Splits `physeq` by `cohort_var`, runs the full pipeline + every decision
#' rule independently per cohort, and reports pairwise Jaccard overlap of
#' the significant-taxa sets between cohorts, per rule.
#'
#' @param physeq A phyloseq object spanning >= 2 cohorts.
#' @param cohort_var Metadata column identifying cohort membership.
#' @param group_var,covariate,fdr_cut,tools As in OTUs_multi_DA().
#' @return data.frame: rule x cohort-pair Jaccard overlap.
#' @export
reproducibility_across_cohorts <- function(physeq, cohort_var, group_var = "Group",
                                           covariate = NULL, fdr_cut = 0.05,
                                           tools = c("edgeR","DESeq2","ALDEx2","metaSeq",
                                                     "ADAPT","ANCOMBC2","MaAsLin3")) {
  meta <- methods::as(phyloseq::sample_data(physeq), "data.frame")
  cohorts <- unique(as.character(meta[[cohort_var]]))
  if (length(cohorts) < 2) stop("Need >= 2 cohorts in ", cohort_var)

  per_cohort <- list()
  for (co in cohorts) {
    message("Cohort: ", co)
    sub <- phyloseq::prune_samples(meta[[cohort_var]] == co, physeq)
    res <- OTUs_multi_DA(sub, group_var = group_var, covariate = covariate,
                         fdr_cut = fdr_cut, tools = tools, make_plots = FALSE)
    cons <- res[[1]]$consolidated
    cmp  <- compare_decision_rules(cons, fdr_cut = fdr_cut)

    padj_cols <- grep("^padj_", names(cons), value = TRUE)
    padj_cols <- setdiff(padj_cols, c("padj_combined", "padj_ranksum",
                                      "padj_fisher", "padj_stouffer", "padj_cauchy"))
    for (pc in padj_cols) {
      tool <- sub("^padj_", "", pc)
      called <- !is.na(cons[[pc]]) & cons[[pc]] < fdr_cut
      cmp <- rbind(cmp, data.frame(rule = tool, method = "single_tool", votes = NA,
                                   n_called = sum(called), precision = NA, recall = NA)[names(cmp)])
    }
    per_cohort[[co]] <- list(consolidated = cons, comparison = cmp)
  }

  jaccard <- function(a, b) if (length(union(a, b)) == 0) NA_real_ else
    length(intersect(a, b)) / length(union(a, b))

  pairs <- utils::combn(cohorts, 2, simplify = FALSE)
  rows  <- list()
  for (pr in pairs) {
    c1 <- per_cohort[[pr[1]]]$consolidated
    c2 <- per_cohort[[pr[2]]]$consolidated

    get_set <- function(cons, rule_name) {
      if (rule_name == "Union (1 of 7)") cons$Taxa[cons$n_sig >= 1]
      else if (rule_name == "Fisher combined") cons$Taxa[!is.na(cons$sig_fisher) & cons$sig_fisher]
      else if (rule_name == "Cauchy combined (dependence-robust)") cons$Taxa[!is.na(cons$sig_cauchy) & cons$sig_cauchy]
      else if (rule_name == "Rank-sum aggregation") { d <- rank_sum_aggregate(cons, fdr_cut = fdr_cut); d$Taxa[d$sig_ranksum] }
      else if (grepl("^padj_", rule_name)) cons$Taxa[!is.na(cons[[rule_name]]) & cons[[rule_name]] < fdr_cut]
      else character(0)
    }

    for (rule_name in unique(per_cohort[[1]]$comparison$rule)) {
      s1 <- get_set(c1, rule_name); s2 <- get_set(c2, rule_name)
      rows[[length(rows) + 1]] <- data.frame(
        cohort_1 = pr[1], cohort_2 = pr[2], rule = rule_name,
        n1 = length(s1), n2 = length(s2), overlap = length(intersect(s1, s2)),
        jaccard = jaccard(s1, s2), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}
