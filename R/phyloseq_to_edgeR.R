## =============================================================================
## R/phyloseq_to_edgeR.R
## =============================================================================

#' Convert a phyloseq object to an edgeR DGEList
#'
#' @param physeq A phyloseq object.
#' @param group Character (a sample_data column name) or a vector of group
#'   labels, one per sample.
#' @param method Normalization method passed to edgeR::calcNormFactors().
#' @param ... Passed to edgeR::DGEList().
#' @return An edgeR DGEList with dispersions estimated.
#' @export
phyloseq_to_edgeR <- function(physeq, group, method = "RLE", ...) {
  if (!requireNamespace("edgeR", quietly = TRUE)) stop("Package 'edgeR' is required")
  if (!requireNamespace("phyloseq", quietly = TRUE)) stop("Package 'phyloseq' is required")

  if (!phyloseq::taxa_are_rows(physeq)) physeq <- phyloseq::t(physeq)

  x <- as(phyloseq::otu_table(physeq), "matrix")
  x <- x + 1   # protect against log(0) / overflow downstream

  if (length(group) == 1 && phyloseq::nsamples(physeq) > 1) {
    group <- phyloseq::sample_data(physeq)[[group]]
    if (is.null(group)) stop("Sample variable '", group, "' not found in sample_data")
  }

  taxonomy <- phyloseq::tax_table(physeq, errorIfNULL = FALSE)
  if (!is.null(taxonomy)) {
    taxonomy <- data.frame(as(taxonomy, "matrix"))
    rownames(taxonomy) <- rownames(x)
  }

  y <- edgeR::DGEList(counts = x, group = group, genes = taxonomy,
                      remove.zeros = TRUE, ...)
  z <- edgeR::calcNormFactors(y, method = method)
  if (!all(is.finite(z$samples$norm.factors)))
    stop("Non-finite normalization factors detected. Consider changing 'method'.")
  edgeR::estimateDisp(z)
}
