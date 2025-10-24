setwd("./GWMC_HOT_COLD/")

print(i)
print(filtered_metadata)

biome_file <- "./GWMC_HOT_COLD/GWMC_HOT_COLD_genus_table.biom"

sample_table_file <- "./GWMC_HOT_COLD_metadata.csv"

GWMC_HOT_COLD <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)

GWMC_HOT_COLD_DA <- OTUs_multi_DA(GWMC_HOT_COLD)


#### Sample 1: GWMC_HOT_COLD ########

metadata <- read.table("./GWMC_HOT_COLD_metadata.csv", sep = "\t", header = TRUE)

head(metadata)

# Step 1: Identify the most frequent group
most_frequent_group <- metadata %>%
  dplyr::count(Age_Group) %>%             # Replace 'Group' with the column containing "Lake" or "Watershed"
  dplyr::arrange(desc(n)) %>%
  dplyr::slice(1) %>%
  dplyr::pull(Age_Group)

# Step 2: Filter for the most frequent group
filtered_metadata <- metadata %>%
  dplyr::filter(Age_Group == most_frequent_group)


# Loop to generate 100 replicates
for (j in 1:10) {
  setwd("./GWMC_HOT_COLD/")
  
  # Create a new directory for each replicate
  dir_name <- paste0("Replicate_", j)
  dir.create(dir_name, showWarnings = FALSE)  # Create the directory, suppress warning if it exists
  setwd(paste0("./GWMC_HOT_COLD/", dir_name))
  
  # Randomly assign "Case" or "Control" to Age_Group
  replicate_data <- filtered_metadata %>% dplyr::mutate(Age_Group = sample(c("Case", "Control"), n(), replace = TRUE))
  
  # Print the replicate number and data
  print(paste("Replicate", j))
  print(replicate_data)
  
  # Save the replicate data into the corresponding directory
  write.table(replicate_data, file = paste0( "filtered_labeled_metadata_replicate_", j, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  #dir_name2 = paste0(dir_name)
  
  
  # biome_file <- "./GWMC_HOT_COLD/GWMC_HOT_COLD_genus_table.biom"
  
  #  sample_table_file <-  paste0("filtered_labeled_metadata_replicate_", j, ".txt")
  
  # GWMC_HOT_COLD <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)
  
  #  DEs <- OTUs_multi_DA(GWMC_HOT_COLD)
  
}

for (j in 1:10) {
  dir_name <- paste0("Replicate_", j)
  
  setwd(paste0("./GWMC_HOT_COLD/", dir_name))
  
  
  biome_file <- "./GWMC_HOT_COLD/GWMC_HOT_COLD_genus_table.biom"
  
  sample_table_file <-  paste0("filtered_labeled_metadata_replicate_", j, ".txt")
  
  GWMC_HOT_COLD <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)
  
  DAs <- OTUs_multi_DA(GWMC_HOT_COLD)
  
}



#### Sample 2: Office ########

setwd("./Office/")


biome_file <- "Office_genus_table.biom"

sample_table_file <- "Office_metadata_sample_table.txt"

Office <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)

Office_DAs <- OTUs_multi_DA(Office)


metadata <- read.table("Office_metadata_sample_table.txt", sep = "\t", header = TRUE)

head(metadata)

# Step 1: Identify the most frequent group
most_frequent_group <- metadata %>%
  dplyr::count(Age_Group) %>%             # Replace 'Group' with the column containing "Lake" or "Watershed"
  dplyr::arrange(desc(n)) %>%
  dplyr::slice(1) %>%
  dplyr::pull(Age_Group)

# Step 2: Filter for the most frequent group
filtered_metadata <- metadata %>%
  dplyr::filter(Age_Group == most_frequent_group)


# Loop to generate 100 replicates
for (j in 1:10) {
  setwd("./Office/")
  
  # Create a new directory for each replicate
  dir_name <- paste0("Replicate_", j)
  dir.create(dir_name, showWarnings = FALSE)  # Create the directory, suppress warning if it exists
  setwd(paste0("./Office/", dir_name))
  
  # Randomly assign "Case" or "Control" to Age_Group
  replicate_data <- filtered_metadata %>% dplyr::mutate(Age_Group = sample(c("Case", "Control"), n(), replace = TRUE))
  
  # Print the replicate number and data
  print(paste("Replicate", j))
  print(replicate_data)
  
  # Save the replicate data into the corresponding directory
  write.table(replicate_data, file = paste0( "filtered_labeled_metadata_replicate_", j, ".txt"), sep = "\t", row.names = FALSE, quote = FALSE)
  
  #dir_name2 = paste0(dir_name)
  
  
  
  
}

for (j in 1:10) {
  dir_name <- paste0("Replicate_", j)
  
  setwd(paste0("./Office/", dir_name))
  
  biome_file <- "./Office/Office_genus_table.biom"
  
  sample_table_file <-  paste0("filtered_labeled_metadata_replicate_", j, ".txt")
  
  Office <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)
  
  Office_DAs <- OTUs_multi_DA(Office)
  
}
