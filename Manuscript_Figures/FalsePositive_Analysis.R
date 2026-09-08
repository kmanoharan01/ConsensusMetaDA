# Load the packages
library(ggplot2)
library(phyloseq)
library(DESeq2)
library(ALDEx2)
library(edgeR)
library(metagenomeSeq)
library(ADAPT)
library(dplyr)
library(tibble)
library(vegan)
library(ggplot2)
library(ConsensusMetaDA)

setwd("./Manuscript_Figures/")


biome_file <- "./Data/GWMC_HOT_COLD_genus_table.biom"

sample_table_file <- "./Data/GWMC_HOT_COLD_metadata.csv"

GWMC_HOT_COLD <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)

GWMC_HOT_COLD_DA <- OTUs_multi_DA(GWMC_HOT_COLD)


#### Sample 1: GWMC_HOT_COLD ########

metadata <- read.table("./Data/GWMC_HOT_COLD_metadata.csv", sep = "\t", header = TRUE)

head(metadata)

# Step 1: Identify the most frequent group
most_frequent_group <- metadata %>%
  dplyr::count(Group) %>%             # Replace 'Group' with the column containing "Lake" or "Watershed"
  dplyr::arrange(desc(n)) %>%
  dplyr::slice(1) %>%
  dplyr::pull(Group)

# Step 2: Filter for the most frequent group
filtered_metadata <- metadata %>%
  dplyr::filter(Group == most_frequent_group)


# Loop to generate 100 replicates
for (j in 1:10) {
  setwd("./Data/")
  
  # Create a new directory for each replicate
  dir_name <- paste0("Replicate_", j)
  dir.create(dir_name, showWarnings = FALSE)  # Create the directory, suppress warning if it exists
  setwd(paste0("./Data/", dir_name))
  
  # Randomly assign "Case" or "Control" to Group
  replicate_data <- filtered_metadata %>% dplyr::mutate(Group = sample(c("Case", "Control"), n(), replace = TRUE))
  
  # Print the replicate number and data
  print(paste("Replicate", j))
  print(replicate_data)
  
  # Save the replicate data into the corresponding directory
  write.table(replicate_data, file = paste0( "filtered_labeled_metadata_replicate_", j, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  #dir_name2 = paste0(dir_name)
  
  
  # biome_file <- "./Data/GWMC_HOT_COLD_genus_table.biom"
  
  #  sample_table_file <-  paste0("filtered_labeled_metadata_replicate_", j, ".txt")
  
  # GWMC_HOT_COLD <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)
  
  #  DEs <- OTUs_multi_DA(GWMC_HOT_COLD)
  
}

for (j in 1:10) {
  dir_name <- paste0("Replicate_", j)
  
  setwd(paste0("./Data/", dir_name))
  
  
  biome_file <- "./Data/GWMC_HOT_COLD_genus_table.biom"
  
  sample_table_file <-  paste0("filtered_labeled_metadata_replicate_", j, ".txt")
  
  GWMC_HOT_COLD <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)
  
  DAs <- OTUs_multi_DA(GWMC_HOT_COLD)
  
}



#### Sample 2: Office ########

setwd("./Office/")


biome_file <- "./Data/Office_genus_table.biom"

sample_table_file <- "./Data/Office_metadata_sample_table.txt"

Office <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)

Office_DAs <- OTUs_multi_DA(Office)


metadata <- read.table("Office_metadata_sample_table.txt", sep = "\t", header = TRUE)

head(metadata)

# Step 1: Identify the most frequent group
most_frequent_group <- metadata %>%
  dplyr::count(Group) %>%             # Replace 'Group' with the column containing "Lake" or "Watershed"
  dplyr::arrange(desc(n)) %>%
  dplyr::slice(1) %>%
  dplyr::pull(Group)

# Step 2: Filter for the most frequent group
filtered_metadata <- metadata %>%
  dplyr::filter(Group == most_frequent_group)


# Loop to generate 100 replicates
for (j in 1:10) {
  setwd("./Data/")
  
  # Create a new directory for each replicate
  dir_name <- paste0("Replicate_", j)
  dir.create(dir_name, showWarnings = FALSE)  # Create the directory, suppress warning if it exists
  setwd(paste0("./Data/", dir_name))
  
  # Randomly assign "Case" or "Control" to Group
  replicate_data <- filtered_metadata %>% dplyr::mutate(Group = sample(c("Case", "Control"), n(), replace = TRUE))
  
  # Print the replicate number and data
  print(paste("Replicate", j))
  print(replicate_data)
  
  # Save the replicate data into the corresponding directory
  write.table(replicate_data, file = paste0( "filtered_labeled_metadata_replicate_", j, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  #dir_name2 = paste0(dir_name)
  
  
  
  
}

for (j in 1:10) {
  dir_name <- paste0("Replicate_", j)
  
  setwd(paste0("./Data/", dir_name))
  
  biome_file <- "./Data/Office_genus_table.biom"
  
  sample_table_file <-  paste0("filtered_labeled_metadata_replicate_", j, ".txt")
  
  Office <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)
  
  Office_DAs <- OTUs_multi_DA(Office)
  
}
