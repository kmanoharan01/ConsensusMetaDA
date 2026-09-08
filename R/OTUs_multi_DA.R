## =============================================================================
## R/OTUs_multi_DA.R
## =============================================================================

utils::globalVariables(c("metadata", "qval_individual", "pval_individual",
                         "feature", "coef", "padj"))

#' @keywords internal
#' @noRd
.safe <- function(label, expr) {
  message("  running ", label, " ...")
  tryCatch(force(expr), error = function(e) {
    warning(sprintf("%s failed: %s", label, conditionMessage(e)), call. = FALSE)
    NULL
  })
}

#' @keywords internal
#' @noRd
.da_plots <- function(cons, contrast_pair, fdr_cut, outdir) {
  cols <- c(edgeR = "#E41A1C", DESeq2 = "#377EB8", ALDEx2 = "#4DAF4A",
            metaSeq = "#FF7F00", ADAPT = "#984EA3", ANCOMBC2 = "#A65628",
            MaAsLin3 = "#F781BF")
  ## only per-tool padj columns - exclude combined ones so the UpSet stays tools-only
  padj_cols <- grep("^padj_", names(cons), value = TRUE)
  padj_cols <- setdiff(padj_cols, c("padj_fisher", "padj_stouffer", "padj_cauchy",
                                    "padj_combined", "padj_ranksum"))
  if (!length(padj_cols)) return(invisible(NULL))
  sets <- lapply(padj_cols, function(p) cons$Taxa[!is.na(cons[[p]]) & cons[[p]] < fdr_cut])
  names(sets) <- sub("^padj_", "", padj_cols)
  sets <- Filter(function(x) length(x) > 0, sets)
  if (length(sets) < 2) { message("  <2 tools with DA taxa - no UpSet plot"); return(invisible(NULL)) }

  if (requireNamespace("UpSetR", quietly = TRUE)) {
    grDevices::pdf(file.path(outdir, paste0(contrast_pair, "_UpSet.pdf")), width = 12, height = 8)
    print(UpSetR::upset(UpSetR::fromList(sets), nsets = length(sets),
                        mainbar.y.label = "Shared DA taxa", sets.x.label = "DA taxa per method",
                        sets.bar.color = unname(cols[names(sets)]), order.by = "freq"))
    grDevices::dev.off()
  } else message("  UpSetR not installed - skipping UpSet plot")
  invisible(NULL)
}

#' Multi-method differential abundance with consensus
#'
#' Runs up to seven DA tools on each pairwise contrast of a grouping variable
#' and consolidates the results into one table carrying per-taxon agreement
#' counts (`n_sig`, `n_tools`) and, when `add_combined = TRUE`, the
#' Fisher/Stouffer/Cauchy combined p-values.
#'
#' @param physeq A phyloseq object.
#' @param group_var Character. Metadata column holding the grouping variable.
#' @param covariate Character or NULL. Optional second metadata column to
#'   adjust for; model becomes `~ group_var + covariate` and edgeR switches
#'   from exactTest to a quasi-likelihood GLM.
#' @param group_order Character vector giving level order; first element of
#'   each pair is the reference.
#' @param fdr_cut Numeric. Adjusted p-value threshold defining significance.
#' @param outdir Directory for result files. Defaults to tempdir().
#' @param tools Character vector of tools to run.
#' @param cores Integer, cores for ALDEx2's Monte Carlo sampling.
#' @param make_plots Logical, write an UpSet diagram per contrast.
#' @param add_combined Logical, append Fisher/Stouffer/Cauchy combined p-values
#'   to the consolidated table (and the ALLDE file). Default TRUE.
#' @param verbose Logical.
#' @return Invisibly, a named list with one element per contrast, each
#'   holding the per-tool raw results and `$consolidated`.
#' @importFrom rlang .data
#' @export
OTUs_multi_DA <- function(physeq,
                          group_var    = "Group",
                          covariate    = NULL,
                          group_order  = NULL,
                          fdr_cut      = 0.05,
                          outdir       = tempdir(),
                          tools        = c("edgeR","DESeq2","ALDEx2","metaSeq",
                                           "ADAPT","ANCOMBC2","MaAsLin3"),
                          cores        = 1L,
                          make_plots   = TRUE,
                          add_combined = TRUE,
                          verbose      = FALSE) {

  stopifnot(inherits(physeq, "phyloseq"))
  tools <- match.arg(tools, several.ok = TRUE)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  if (!phyloseq::taxa_are_rows(physeq)) {
    message("taxa_are_rows(physeq) is FALSE - transposing")
    physeq <- phyloseq::t(physeq)
  }

  meta <- methods::as(phyloseq::sample_data(physeq), "data.frame")

  if (!group_var %in% colnames(meta))
    stop(sprintf("group_var '%s' not in sample metadata. Available: %s",
                 group_var, paste(colnames(meta), collapse = ", ")))
  if (!is.null(covariate)) {
    if (!covariate %in% colnames(meta))
      stop(sprintf("covariate '%s' not in sample metadata. Available: %s",
                   covariate, paste(colnames(meta), collapse = ", ")))
    if (identical(covariate, group_var)) stop("covariate and group_var are the same column")
  }

  rhs  <- if (is.null(covariate)) group_var else paste(group_var, "+", covariate)
  form <- stats::as.formula(paste("~", rhs))

  observed <- unique(as.character(meta[[group_var]]))
  if (!is.null(group_order)) {
    miss <- setdiff(group_order, observed)
    if (length(miss))
      stop(sprintf("group_order has levels absent from %s: %s",
                   group_var, paste(miss, collapse = ", ")))
    treat_list <- group_order
  } else treat_list <- observed
  if (length(treat_list) < 2)
    stop(sprintf("Need >= 2 levels in %s; found %d.", group_var, length(treat_list)))

  pair_mat <- utils::combn(treat_list, 2)   # each unordered pair ONCE

  message(sprintf("ConsensusMetaDA | model: %s | %d contrast(s) | tools: %s",
                  paste(deparse(form), collapse = ""), ncol(pair_mat),
                  paste(tools, collapse = ", ")))
  if (!is.null(covariate)) {
    message("Samples per group x covariate:")
    print(table(meta[[group_var]], meta[[covariate]]))
  }

  out <- list()

  for (k in seq_len(ncol(pair_mat))) {
    ref  <- pair_mat[1, k]
    test <- pair_mat[2, k]
    contrast_pair <- paste0(test, "_vs_", ref)
    message("\n== ", contrast_pair, " ==")

    keep     <- meta[[group_var]] %in% c(test, ref)
    sd_sub   <- phyloseq::sample_data(physeq)[keep, ]
    otu_sub  <- phyloseq::otu_table(physeq)[, rownames(sd_sub)]
    sub      <- phyloseq::phyloseq(otu_sub, sd_sub)
    sub_meta <- methods::as(phyloseq::sample_data(sub), "data.frame")

    r <- list()

    if ("edgeR" %in% tools) r$edgeR <- .safe("edgeR", {
      dge <- phyloseq_to_edgeR(sub, group = group_var)
      if (is.null(covariate)) {
        tt <- edgeR::topTags(edgeR::exactTest(dge), n = Inf, adjust.method = "BH")
      } else {
        design <- stats::model.matrix(form, data = sub_meta)
        dge <- edgeR::estimateDisp(dge, design)
        fit <- edgeR::glmQLFit(dge, design)
        cf <- grep(paste0("^", group_var), colnames(design), value = TRUE)   # anchored
        if (!length(cf)) stop("no design column starts with '", group_var, "'")
        tt <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = cf), n = Inf, adjust.method = "BH")
      }
      tt$table
    })

    if ("DESeq2" %in% tools) r$DESeq2 <- .safe("DESeq2", {
      dds <- phyloseq::phyloseq_to_deseq2(sub, form)
      dds <- DESeq2::DESeq(dds, sfType = "poscounts", quiet = TRUE)
      res <- DESeq2::results(dds, contrast = c(group_var, test, ref), tidy = TRUE, format = "DataFrame")
      colnames(res)[colnames(res) == "row"] <- "Taxa"
      res
    })

    if ("ALDEx2" %in% tools) r$ALDEx2 <- .safe("ALDEx2", {
      ALDEx2::aldex(reads = as.matrix(phyloseq::otu_table(sub)),          # pair subset, not full object
                    conditions = as.character(sub_meta[[group_var]]),
                    mc.samples = 128, denom = "all", verbose = verbose,
                    useMC = cores > 1L, cores = cores,
                    test = "t", effect = TRUE, paired.test = FALSE)
    })

    if ("ADAPT" %in% tools) r$ADAPT <- .safe("ADAPT", {
      ADAPT::summary(ADAPT::adapt(sub, cond.var = group_var), select = "all")
    })

    if ("metaSeq" %in% tools) r$metaSeq <- .safe("metagenomeSeq", {
      m  <- methods::as(phyloseq::otu_table(sub), "matrix")
      pd <- Biobase::AnnotatedDataFrame(sub_meta)
      fd <- Biobase::AnnotatedDataFrame(data.frame(OTU_ID = rownames(m), row.names = rownames(m)))
      mg <- metagenomeSeq::newMRexperiment(m, phenoData = pd, featureData = fd)
      mg <- metagenomeSeq::cumNorm(mg, p = metagenomeSeq::cumNormStat(mg))
      mod <- stats::model.matrix(form, data = Biobase::pData(pd))
      metagenomeSeq::MRfulltable(metagenomeSeq::fitFeatureModel(mg, mod), number = nrow(m))
    })

    if ("ANCOMBC2" %in% tools) r$ANCOMBC2 <- .safe("ANCOMBC2", {
      a <- ANCOMBC::ancombc2(data = sub, fix_formula = rhs, p_adj_method = "BH",
                             prv_cut = 0.10, lib_cut = 0, group = group_var,
                             struc_zero = TRUE, neg_lb = TRUE, verbose = FALSE)$res
      lfc <- grep(paste0("^lfc_", group_var), names(a), value = TRUE)[1]
      pv  <- grep(paste0("^p_",   group_var), names(a), value = TRUE)[1]
      qv  <- grep(paste0("^q_",   group_var), names(a), value = TRUE)[1]
      if (any(is.na(c(lfc, pv, qv))))
        stop("could not locate ", group_var, " columns in ancombc2 output; saw: ",
             paste(head(names(a), 12), collapse = ", "))
      data.frame(taxon = a$taxon, lfc = a[[lfc]], p_val = a[[pv]], q_val = a[[qv]])
    })

    if ("MaAsLin3" %in% tools) r$MaAsLin3 <- .safe("MaAsLin3", {
      inp <- as.data.frame(t(as.matrix(phyloseq::otu_table(sub))))
      md  <- data.frame(row.names = phyloseq::sample_names(sub), stringsAsFactors = FALSE)
      md[[group_var]] <- as.character(sub_meta[[group_var]])
      if (!is.null(covariate)) md[[covariate]] <- as.character(sub_meta[[covariate]])
      stopifnot(identical(rownames(inp), rownames(md)))
      fx <- if (is.null(covariate)) group_var else c(group_var, covariate)
      fit <- maaslin3::maaslin3(
        input_data = inp, input_metadata = md,
        output = file.path(outdir, paste0("maaslin3_", contrast_pair)),
        fixed_effects = fx, reference = paste0(group_var, ",", ref),
        normalization = "TSS", transform = "LOG", min_prevalence = 0.1,
        min_abundance = 0, max_significance = 0.1, correction = "BH", cores = 1,
        plot_summary_plot = FALSE, plot_associations = FALSE, max_pngs = 0, verbosity = "ERROR")
      res <- dplyr::filter(fit$fit_data_abundance$results,
                           .data$metadata == group_var, !is.na(.data$qval_individual))
      dplyr::select(res, "feature", "coef", pval = "pval_individual", qval = "qval_individual")
    })

    for (nm in names(r)) if (!is.null(r[[nm]]))
      utils::write.table(r[[nm]],
                         file.path(outdir, sprintf("%s_%s_DE_results.txt", contrast_pair, nm)),
                         quote = FALSE, sep = "\t", col.names = NA)

    cons <- consolidate_DA_results(
      edgeR_res = r$edgeR, DESeq2_res = r$DESeq2, ALDEx2_res = r$ALDEx2,
      metagenomeSeq_res = r$metaSeq, ADAPT_res = r$ADAPT,
      ANCOMBC2_res = r$ANCOMBC2, MaAsLin3_res = r$MaAsLin3,
      fdr_cut = fdr_cut, add_combined = add_combined)

    utils::write.table(cons, file.path(outdir, paste0(contrast_pair, "_ALLDE_results.txt")),
                       quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

    print_consensus_summary(cons, contrast_pair, fdr_cut = fdr_cut)
    if (make_plots) .da_plots(cons, contrast_pair, fdr_cut, outdir)

    r$consolidated <- cons
    out[[contrast_pair]] <- r
  }

  invisible(out)
}
