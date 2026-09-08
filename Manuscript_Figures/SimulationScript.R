##Load library

library(sparseDOSSA)

setwd("./Manuscript_Figures/Data/")
# Define the number of repetitions
num_reps <- 100

#set.seed(123)  # Ensures different runs get different seeds but are still reproducible
#seeds <- sample(1:10000, num_reps, replace = FALSE) 


################### Looping #########

# Define parameter ranges
spike <- c(10)
dpercent <- c("0.05", "0.10", "0.15","0.25","0.50","1")

# Loop through each spike strength
for (s in spike) {
  # Loop through each percent spiked
  for (p in dpercent) {
    # Loop through each repetition
    for (i in 1:num_reps) {
      # Define folder name with spike and percent parameters
      main_folder <- paste0("./Manuscript_Figures/spike", s, "_percent", p)
      rep_folder <- paste0(main_folder, "/reps", i)
      
      # Create directories if they do not exist
      if (!dir.exists(main_folder)) {
        dir.create(main_folder, recursive = TRUE)
      }
      if (!dir.exists(rep_folder)) {
        dir.create(rep_folder, recursive = TRUE)
      }
      
      setwd(rep_folder)
      
      #### Simulation code ####
      #?sparseDOSSA()
      rep2 <- sparseDOSSA(
        strNormalizedFileName = "rep2-Normalized.pcl",      # Output: Normalized abundance data
        strCountFileName = "rep2-Counts.pcl",   # Output: Raw counts data
        parameter_filename = "rep2-SyntheticMicrobiomeParameterFile.txt",  # Output: Simulation parameters
        # bugs_to_spike = 0,     # Introduce 10 spiked-in differentially abundant features
        datasetCount = 1,       # Generate one dataset
        read_depth = 20000,     # Read depth per sample
        number_features = 500,  # Number of microbial species (OTUs)
        number_samples = 100,    # Number of samples to simulate
        percent_spiked = as.numeric(p),  # Use percent from loop
        minLevelPercent = 0.5,
        spikeStrength = as.character(s),  # Use spike strength from loop
        seed = i,             # Set a seed for reproducibility
        number_metadata = 1,
        verbose = TRUE,          # Enable logging
        noZeroInflate =  TRUE
      )
      
      #?sparseDOSSA()
      
      #SampleID
      
      rep2$OTU_count[[1]][1] <- "SampleID"
      
      rep2$OTU_count[[1]][c(1,4),]
      
      rep2$OTU_count[[1]][c(1,6:1505),]
      
      write.table(t(rep2$OTU_count[[1]][c(1,4),]), "sample_table.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      #metadata <- read.table("samples_data.txt", sep = "\t", header = TRUE) # Simulated groups
      #rownames(metadata) <- metadata$SampleID
      
      write.table(rep2$OTU_count[[1]][c(1,6:1505),], "simulation_counts_data.txt", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
      
      
      
      # Read the normalized abundance data
      synthetic_data <- read.table("simulation_counts_data.txt", header = TRUE, sep = "\t")
      
      row.names(synthetic_data) <- synthetic_data$SampleID
      
      synthetic_data <- synthetic_data[, -1]
      
      # View the first few rows
      head(synthetic_data)
      
      #rep1$OTU_count
      metadata <- read.table("temp_sample_table.txt", sep = "\t", header = TRUE) # Simulated groups
      
      row.names(metadata) <- metadata$SampleID
      metadata$Group <- metadata$Metadata3
      
      
      write.table(metadata, "sample_table.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      
      identical(colnames(synthetic_data), rownames(metadata))
      
      # Print progress
      cat("Completed: spike =", s, ", percent =", p, ", rep =", i, "\n")
    }
  }
}

######################## Simulation data to biom and DA analyses ##############

##Load libraries

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

setwd(paste0("./Manuscript_Figures/Data/spike", s, "_percent", p, "/reps", i, "/"))

num_reps <- 100

# Define parameter ranges
#spike <- c(1,2,3,4,5,10)

spike <- c(10)

#Decimal Percentage

dpercent <- c("0.05", "0.10", "0.15","0.25","0.50","1")

# Loop through each spike strength
for (s in spike) {
  # Loop through each percent spiked
  for (p in dpercent) {
    # Loop through each repetition
    for (i in 1:num_reps) {
      
      setwd(paste0("./Manuscript_Figures/Data/spike", s, "_percent", p, "/reps", i, "/"))
      
      system("module load python3/3.10.4")
      
      system("~/.local/bin/biom convert -i simulation_counts_data.txt   -o simulation_counts_data.biom --table-type='OTU table' --to-json")
      
      sample_t1 <- "sample_table2.txt"
      
      test_biom <- "./simulation_counts_data.biom"
      
      sim_Rep_build <- build_OTU_counts(test_biom, sample_t1)
      
      OTUs_multi_DA(sim_Rep_build)
    }
  }
}

