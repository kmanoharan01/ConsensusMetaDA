
#' Package imports and global variables
#'
#' @name ConsensusMetaDA-imports
#' @import ggplot2
#' @import dplyr
#' @import phyloseq
#' @importFrom Biobase AnnotatedDataFrame
#' @importFrom VennDiagram venn.diagram
#' @importFrom grid grid.draw
#' @importFrom vegan specnumber rarefy
#' @importFrom magrittr %>%
NULL

utils::globalVariables(c(
  "Abundance", "Age_Group", "Class", "Family", "Genus", "Kingdom",
  "Mean_Abundance", "Order", "Phylum", "Sample", "Taxa",
  "padj_ADAPT", "padj_ALDEx2", "padj_DESeq2", "padj_edgeR", "padj_metaSeq"
))
