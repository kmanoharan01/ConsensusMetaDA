## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
# install supporting packages

# ADAPT R package requires R version 4.4 hence need all R updated to R v4.4.2

# One line per package
# if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# if (!require("phyloseq", quietly = TRUE)) BiocManager::install("phyloseq")
# if (!require("DESeq2", quietly = TRUE)) BiocManager::install("DESeq2")
# if (!require("ALDEx2", quietly = TRUE)) BiocManager::install("ALDEx2")
# if (!require("edgeR", quietly = TRUE)) BiocManager::install("edgeR")
# if (!require("ADAPT", quietly = TRUE)) BiocManager::install("ADAPT")
# if (!require("metagenomeSeq", quietly = TRUE)) BiocManager::install("metagenomeSeq")
# 
# if (!require("ConsensusMetaDA", quietly = TRUE)) {
#   if (!require("remotes", quietly = TRUE)) install.packages("remotes")
#   remotes::install_github("kmanoharan01/ConsensusMetaDA")
# }


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------

# Load the packages
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(ALDEx2)
library(edgeR)
library(metagenomeSeq)
library(ADAPT)


# load ConsensusMetaDA

library(ConsensusMetaDA)


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
library(ADAPT)

data("ecc_saliva")

otu_mat <- as(otu_table(ecc_saliva), "matrix")
sample_df <- as(sample_data(ecc_saliva), "data.frame")
tax_mat <- as(tax_table(ecc_saliva), "matrix")

#otu_final <- data.frame(sampleID = rownames(otu_mat), otu_mat)

#write.table(otu_mat, file = "ecc_saliva_Otu.txt",  col.names = TRUE, sep = "\t", quote = FALSE)

#write.table(sample_df, file = "ecc_saliva_sample_table.txt", col.names = TRUE, sep = "\t", quote = FALSE)

#write.table(tax_mat, file = "ecc_saliva_tax_table.txt", col.names = TRUE, sep = "\t", quote = FALSE)

################

#system("module load python3/3.10.4")

#system("~/.local/bin/biom convert -i ecc_saliva_Otu.txt   -o ecc_saliva_Otu.biom --table-type='OTU table' --to-json")

ecc_saliva_biom <- "../inst//extdata/ecc_saliva_Otu.biom"

ecc_saliva_sample_table <- "../inst//extdata/ecc_saliva_sample_table.txt"

ecc_saliva_tax <- tax_mat

taxa_names(tax_mat)
#tax_mat_test <- tax_table(tax_mat)


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
## build
ecc_saliva_build_OTU_counts <- build_OTU_counts(ecc_saliva_biom, 
                                                 ecc_saliva_sample_table, 
                                                 ecc_saliva_tax, 
                                                 abundance_threshold = 10, 
                                                 prevalence_threshold = 0.1, 
                                                 rarity_threshold = 0, 
                                                 variance_threshold = 0 
)

ecc_saliva_build_OTU_counts



## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------

ecc_saliva_build_OTU_counts <- build_OTU_counts(ecc_saliva_biom, 
                                                 ecc_saliva_sample_table, 
                                                 ecc_saliva_tax, 
                                                 abundance_threshold = 0, 
                                                 prevalence_threshold = 0, 
                                                 rarity_threshold = 0.1, 
                                                 variance_threshold = 0 
)

ecc_saliva_build_OTU_counts


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------

ecc_saliva_build_OTU_counts <- build_OTU_counts(ecc_saliva_biom, 
                                                 ecc_saliva_sample_table, 
                                                 ecc_saliva_tax, 
                                                 abundance_threshold = 10, 
                                                 prevalence_threshold = 0, 
                                                 rarity_threshold = 0, 
                                                 variance_threshold = 0 
)

ecc_saliva_build_OTU_counts


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------

ecc_saliva_build_OTU_counts <- build_OTU_counts(ecc_saliva_biom, 
                                                 ecc_saliva_sample_table, 
                                                 ecc_saliva_tax, 
                                                 abundance_threshold = 10, 
                                                 prevalence_threshold = 0.1, 
                                                 rarity_threshold = 0, 
                                                 variance_threshold = 0 
)

ecc_saliva_build_OTU_counts


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------

ecc_saliva_build_OTU_counts <- build_OTU_counts(ecc_saliva_biom, 
                                                 ecc_saliva_sample_table, 
                                                 ecc_saliva_tax, 
                                                 abundance_threshold = 0, 
                                                 prevalence_threshold = 0, 
                                                 rarity_threshold = 0, 
                                                 variance_threshold = 0.2 
)

ecc_saliva_build_OTU_counts


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
## build
ecc_saliva_build_OTU_counts <- build_OTU_counts(ecc_saliva_biom, 
                                                 ecc_saliva_sample_table, 
                                                 ecc_saliva_tax, 
                                                 abundance_threshold = 10, 
                                                 prevalence_threshold = 0.1, 
                                                 rarity_threshold = 0, 
                                                 variance_threshold = 0 
)

ecc_saliva_build_OTU_counts

## Perform DA

ecc_saliva_DA <- OTUs_multi_DA(ecc_saliva_build_OTU_counts)



## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------

print(ecc_saliva_DA$Control_vs_Case$upset_plot)



## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
library(grid)
grid.newpage()
grid.draw(ecc_saliva_DA$Control_vs_Case$VennDiagram)


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------

ecc_saliva_plots <- OTUs_plots(ecc_saliva_build_OTU_counts)


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------

ecc_saliva_plots <- OTUs_plots(ecc_saliva_build_OTU_counts)


## ----echo=TRUE, message=FALSE-------------------------------------------------

print(ecc_saliva_plots$AlphaDiv$tax_level)


## ----echo=TRUE, message=FALSE-------------------------------------------------

print(ecc_saliva_plots$BetaDiv)


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------

print(ecc_saliva_plots$all_comparisons_results$Control_vs_Case$Genus_plot)


## ----echo=TRUE, message=FALSE, warning=FALSE----------------------------------
#citation("ConsesnsusMetaDA")

