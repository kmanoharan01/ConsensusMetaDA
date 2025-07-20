# ConsensusMetaDA

ConsensusMetaDA is an R package for microbiome analysis using multiple algorithms - reaching consensus. No one tool perform best DA anlayes on metaganome data. Hence it needs a consensus approach to obtain more robust differential abundant microbiome. ConsensusMetaDA helps achieve this using five popularly used tools. ConsensusMetaDA uses popular biom format and samples table as input to perform DA along with it generates standard microbiome visualization plots such as Alpha Diversity, Beta Diversity, Rarefaction curve, Bidirectional plot and Scale plot at various taxa levels.

## Installing ConsensusMetaDA

To obtain the original version from github, install devtools in R and use the following:

```R

required packages

# Install ggplot2 
install.packages("ggplot2")

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("phyloseq")
BiocManager::install("DESeq2")
BiocManager::install("ALDEx2")
BiocManager::install("edgeR")
BiocManager::install("ADAPT")
BiocManager::install("edgeR")
BiocManager::install("metagenomeSeq")


# Load the packages
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(ALDEx2)
library(edgeR)
library(ADAPT)
library(metagenomeSeq)


# Install ConsensusMetaDA
devtools::install_github("kmanoharan01/ConsensusMetaDA")
# or
remotes::install_github("kmanoharan01/ConsensusMetaDA")

```

## Examples

To run ConsensusMetaDA, load the library and follow the examples in the vignette.

```R
library(ConsensusMetaDA)

GWMC_HOT_COLD <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)

DEs <- OTUs_multi_DE(GWMC_HOT_COLD)

OTU_plots(GWMC_HOT_COLD)

```

## Contact

For more details, contact Manoharan Kumar:
manoharan.kumar@jcu.edu.au
# ConsensusMetaDA
