

library(phyloseq)
library(ggplot2)
library(DESeq2)
library(ALDEx2)
library(dplyr)
library(tibble)
library(vegan)




OTU_plots <- function( build_OTU_counts_output = NULL,
                       force_build = FALSE,
                       verbose = FALSE){
  
  ####///---- check inputs ----\\\###
  if(is.null(build_OTU_counts_output)) {
    stop("A build_OTU_counts_output object is not provided. Please provide filenames with full path and rerun.")
  }  
  
  
  ###### color function ##########
  
  # Original set of named colors
  color_palette <- c(
    "royalblue4", "deepskyblue", "blue", "cyan2", "darkorchid", "gold1", "forestgreen", "firebrick", "mediumspringgreen",
    "darkorange1", "saddlebrown", "deeppink", "slategray2", "seagreen", "bisque"
  )
  
  # Function to generate random hexadecimal color codes
  generate_random_colors <- function(n) {
    replicate(n, paste0("#", paste(sample(c(0:9, letters[1:6]), 6, replace = TRUE), collapse = "")))
  }
  
  # Generate additional colors to make up a total of 200 unique colors
  additional_colors <- generate_random_colors(200 - length(color_palette))
  
  # Combine both sets
  all_colors <- c(color_palette, additional_colors)
  
  # Ensure we have exactly 200 colors
  if (length(all_colors) < 200) {
    stop("Not enough colors defined.")
  }
  
  # Use the first 200 colors for your plot
  colors_to_use <- all_colors[1:200]
  
  ##################### bidirectional plot ########
  
  bidirectional_plot <- function(physeq2, tax_level) {
    
    cat("Processing:", tax_level, "\n")
    
    # Define taxonomic levels
    taxonomic_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
    
    # Validate taxonomic level
    if (!(tax_level %in% taxonomic_levels)) {
      stop("Invalid taxonomic level! Choose from: ", paste(taxonomic_levels, collapse = ", "))
    }
    
    # Agglomerate taxa at specified level
    physeq_taxa <- tax_glom(physeq2, taxrank = tax_level)
    
    # Normalize abundance (relative abundance)
    physeq_taxa_norm <- transform_sample_counts(physeq_taxa, function(x) x / sum(x))
    
    # Convert to data frame
    tax_data <- psmelt(physeq_taxa_norm)
    tax_data[, tax_level] <- as.character(tax_data[, tax_level])
    
    # Replace low-abundance taxa
    tax_data[tax_data$Abundance < 0.01, tax_level] <- "<1%_Abundance"
    
    # Write filtered data table
    write.table(tax_data, file = paste0("filtered_physeq_", tax_level, "_data.txt"), 
                col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    # Summarize abundance & prevalence
    tax_summary <- tax_data %>%
      group_by(Age_Group, !!sym(tax_level)) %>%
      summarize(
        Mean_Abundance = mean(Abundance * 100, na.rm = TRUE),
        SE_Abundance = sd(Abundance * 100, na.rm = TRUE) / sqrt(n()),
        Prevalence = sum(Abundance > 0)
      ) %>%
      ungroup()
    
    write.table(tax_summary, file = paste0(tax_level, "_summary_table.txt"), 
                col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    
    # Generate bidirectional plot
    adults_negative_summary <- tax_summary %>% filter(Age_Group == "SKM")
    adults_positive_summary <- tax_summary %>% filter(Age_Group == "SKB")
    
    plot <- ggplot() +
      geom_bar(data = adults_negative_summary, aes(x = reorder(!!sym(tax_level), -Mean_Abundance), 
                                                   y = -Mean_Abundance, fill = "SKM"), 
               stat = "identity", alpha = 0.7, width = 0.5) +
      geom_bar(data = adults_positive_summary, aes(x = reorder(!!sym(tax_level), -Mean_Abundance), 
                                                   y = Mean_Abundance, fill = "SKB"), 
               stat = "identity", alpha = 0.7, width = 0.5) +
      coord_flip() +
      labs(title = paste("Bidirectional Plot of", tax_level, "Prevalence and Abundance"),
           x = tax_level,
           y = "Mean Abundance / Prevalence") +
      scale_y_continuous(labels = abs, breaks = scales::pretty_breaks(n = 10)) +
      theme_minimal()
    
    # Save plot
    ggsave(paste0(tax_level, "_bidirectional_plot.pdf"), plot = plot, width = 15, height = 10)
    
    cat("Completed:", tax_level, "\n")
    
  }
  
  
  ####### Rarefaction Cure #######
  
  
  # Build phyloseq object from BIOM file
  build_OTU_counts_output <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)
  
  # Extract OTU table from phyloseq object
  otu_mat <- as(otu_table(build_OTU_counts_output), "matrix")
  
  # Transpose OTU table for rarefaction
  Data_t <- t(otu_mat)
  
  # Compute species count per sample
  S <- specnumber(Data_t)
  
  # Determine minimum read count per sample (for rarefaction)
  raremax <- min(rowSums(Data_t))
  
  # Perform rarefaction
  Srare <- rarefy(Data_t, raremax)
  
  # Plot observed vs. rarefied species
  pdf("Rarefaction_curve.pdf")
  plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
  abline(0, 1)
  
  # Generate rarefaction curves
  rarecurve(Data_t, step = 1, sample = raremax, col = "blue", cex = 0.4)
  dev.off()
  
  
  ####### PCoA / MDS #########
  
  
  pcoa <- ordinate(build_OTU_counts_output,"PCoA")
  
  pcoa_plot <- plot_ordination(build_OTU_counts_output,ordinate(build_OTU_counts_output, "PCoA"),color="Age_Group", label = "SampleID") + theme_bw() 
  
  pcoa_plot2 <- pcoa_plot + theme(panel.grid = element_blank(), panel.border = element_blank(),
                                  axis.line = element_line(color = "black")) + 
    guides(colour=guide_legend(override.aes=list(size=5)))  + 
    scale_color_manual(values = c("royalblue4", "deepskyblue", "blue","cyan2", "darkorchid","gold1", "forestgreen", "firebrick", "mediumspringgreen", "darkorange1", "deeppink","grey","slategray2", "dodgerblue", "orange2", "maroon","navy"))
  
  ggsave("Both_groups_mds_plot.pdf", plot = pcoa_plot2, width = 8, height = 6)
  
  
   ###### Alpha Diversity : taxa levels ##########
  
  
  # Loop through each taxonomic level
  for (tax_level in taxonomic_levels) {
    # Glom the phyloseq object at the current taxonomic level
    tax_glom_physeq <- tax_glom(build_OTU_counts_output, taxrank = tax_level)
    
    # Estimate Shannon diversity at the current taxonomic level
    shannon_diversity <- estimate_richness(tax_glom_physeq, measures = "Shannon")
    
    write.table(shannon_diversity, file = paste0("Age_Group_alpha_Diversity_",tax_level,".txt"), col.names = TRUE,  row.names = T, sep = "\t", quote = FALSE)
    
  }
  
  
  # Loop through each taxonomic level and plot richness
  for (tax_level in taxonomic_levels) {
    # Glom the phyloseq object at the current taxonomic level
    tax_glom_physeq <- tax_glom(build_OTU_counts_output, taxrank = tax_level)
    
    # Create the plot for Shannon and Fisher diversity indices
    richness_plot <- plot_richness(tax_glom_physeq, measures = c("Shannon"), x = "Age_Group", color = "Age_Group") +
      scale_color_manual(values = c("dodgerblue", "orange2", "maroon","navy")) + 
      theme_minimal() + 
      theme(panel.grid = element_blank(), panel.border = element_blank(), axis.line = element_line(color = "black")) +
      ggtitle(paste("Richness at", tax_level, "level")) +
      theme(  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    
    # Print the plot
    print(richness_plot)
    
    ggsave(paste0("Age_Group_alpha_Diversity_plot_shannon_",tax_level,".pdf"), plot = richness_plot, width = 8, height = 6)
    
  }
  

  
  ####### scale plots ########

  tax_glom_physeq_kingdom <- tax_glom(build_OTU_counts_output, taxrank = "Kingdom")
  tax_glom_physeq_phylum <- tax_glom(build_OTU_counts_output, taxrank = "Phylum")
  tax_glom_physeq_class <- tax_glom(build_OTU_counts_output, taxrank = "Class")
  tax_glom_physeq_order <- tax_glom(build_OTU_counts_output, taxrank = "Order")
  tax_glom_physeq_family <- tax_glom(build_OTU_counts_output, taxrank = "Family")
  tax_glom_physeq_genus <- tax_glom(build_OTU_counts_output, taxrank = "Genus")
  tax_glom_physeq_species <- tax_glom(build_OTU_counts_output, taxrank = "Species")
  
  
  build_OTU_counts_output_kingdom_tmp <- transform_sample_counts(tax_glom_physeq_kingdom,function(x) x / sum(x))
  build_OTU_counts_output_phylum_tmp <- transform_sample_counts(tax_glom_physeq_phylum,function(x) x / sum(x))
  build_OTU_counts_output_class_tmp <- transform_sample_counts(tax_glom_physeq_class,function(x) x / sum(x))
  build_OTU_counts_output_order_tmp <- transform_sample_counts(tax_glom_physeq_order,function(x) x / sum(x))
  build_OTU_counts_output_family_tmp <- transform_sample_counts(tax_glom_physeq_family,function(x) x / sum(x))
  build_OTU_counts_output_genus_tmp <- transform_sample_counts(tax_glom_physeq_genus,function(x) x / sum(x))
  build_OTU_counts_output_species_tmp <- transform_sample_counts(tax_glom_physeq_species,function(x) x / sum(x))
  
  build_OTU_counts_output_kingdom <- tax_glom(build_OTU_counts_output_kingdom_tmp,taxrank='Kingdom')
  build_OTU_counts_output_phylum <- tax_glom(build_OTU_counts_output_phylum_tmp,taxrank='Phylum')
  build_OTU_counts_output_class <- tax_glom(build_OTU_counts_output_class_tmp,taxrank='Class')
  build_OTU_counts_output_order <- tax_glom(build_OTU_counts_output_order_tmp,taxrank='Order')
  build_OTU_counts_output_family <- tax_glom(build_OTU_counts_output_family_tmp,taxrank='Family')
  build_OTU_counts_output_genus <- tax_glom(build_OTU_counts_output_genus_tmp,taxrank='Genus')
  build_OTU_counts_output_species <- tax_glom(build_OTU_counts_output_species_tmp,taxrank='Species')
  
  build_OTU_counts_output_kingdom_data <- psmelt(build_OTU_counts_output_kingdom)
  build_OTU_counts_output_phylum_data <- psmelt(build_OTU_counts_output_phylum)
  build_OTU_counts_output_class_data <- psmelt(build_OTU_counts_output_class)
  build_OTU_counts_output_order_data <- psmelt(build_OTU_counts_output_order)
  build_OTU_counts_output_family_data <- psmelt(build_OTU_counts_output_family)
  build_OTU_counts_output_genus_data <- psmelt(build_OTU_counts_output_genus)
  build_OTU_counts_output_species_data <- psmelt(build_OTU_counts_output_species)
  
  build_OTU_counts_output_kingdom_data$Kingdom <-  as.character(build_OTU_counts_output_kingdom_data$Kingdom)
  build_OTU_counts_output_phylum_data$Phylum <-  as.character(build_OTU_counts_output_phylum_data$Phylum)
  build_OTU_counts_output_class_data$Class <-  as.character(build_OTU_counts_output_class_data$Class)
  build_OTU_counts_output_order_data$Order <-  as.character(build_OTU_counts_output_order_data$Order)
  build_OTU_counts_output_family_data$Family <-  as.character(build_OTU_counts_output_family_data$Family)
  build_OTU_counts_output_genus_data$Genus <-  as.character(build_OTU_counts_output_genus_data$Genus)
  build_OTU_counts_output_species_data$Species <-  as.character(build_OTU_counts_output_species_data$Species)
  
  build_OTU_counts_output_kingdom_data$Kingdom[build_OTU_counts_output_kingdom_data$Abundance < 0.01] <- "<1%_Abundance"
  build_OTU_counts_output_phylum_data$Phylum[build_OTU_counts_output_phylum_data$Abundance < 0.01] <- "<1%_Abundance"
  build_OTU_counts_output_class_data$Class[build_OTU_counts_output_class_data$Abundance < 0.01] <- "<1%_Abundance"
  build_OTU_counts_output_order_data$Order[build_OTU_counts_output_order_data$Abundance < 0.01] <- "<1%_Abundance"
  build_OTU_counts_output_family_data$Family[build_OTU_counts_output_family_data$Abundance < 0.01] <- "<1%_Abundance"
  build_OTU_counts_output_genus_data$Genus[build_OTU_counts_output_genus_data$Abundance < 0.01] <- "<1%_Abundance"
  build_OTU_counts_output_species_data$Species[build_OTU_counts_output_species_data$Abundance < 0.01] <- "<1%_Abundance"
  
  build_OTU_counts_output_kingdom_data[build_OTU_counts_output_kingdom_data == 0] <- NA
  build_OTU_counts_output_phylum_data[build_OTU_counts_output_phylum_data == 0] <- NA
  build_OTU_counts_output_class_data[build_OTU_counts_output_class_data == 0] <- NA
  build_OTU_counts_output_order_data[build_OTU_counts_output_order_data == 0] <- NA
  build_OTU_counts_output_family_data[build_OTU_counts_output_family_data == 0] <- NA
  build_OTU_counts_output_genus_data[build_OTU_counts_output_genus_data == 0] <- NA
  build_OTU_counts_output_species_data[build_OTU_counts_output_species_data == 0] <- NA
  
  
  unique(build_OTU_counts_output_kingdom_data$Kingdom)
  unique(build_OTU_counts_output_phylum_data$Phylum)
  unique(build_OTU_counts_output_class_data$Class)
  unique(build_OTU_counts_output_order_data$Order)
  unique(build_OTU_counts_output_family_data$Family)
  unique(build_OTU_counts_output_genus_data$Genus)
  unique(build_OTU_counts_output_species_data$Species)
  
  length(unique(build_OTU_counts_output_kingdom_data$Kingdom))
  length(unique(build_OTU_counts_output_phylum_data$Phylum))
  length(unique(build_OTU_counts_output_class_data$Class))
  length(unique(build_OTU_counts_output_order_data$Order))
  length(unique(build_OTU_counts_output_family_data$Family))
  length(unique(build_OTU_counts_output_genus_data$Genus))
  length(unique(build_OTU_counts_output_species_data$Species))
  
  
  
  level_kingdom <- unique(build_OTU_counts_output_kingdom_data$Kingdom)
  level_phylum <- unique(build_OTU_counts_output_phylum_data$Phylum)
  level_class <- unique(build_OTU_counts_output_class_data$Class)
  level_order <- unique(build_OTU_counts_output_order_data$Order)
  level_family <- unique(build_OTU_counts_output_family_data$Family)
  level_genus <- unique(build_OTU_counts_output_genus_data$Genus)
  level_species <- unique(build_OTU_counts_output_species_data$Species)
  
  #tmp as before
  build_OTU_counts_output_kingdom_data$Kingdom <- factor(build_OTU_counts_output_kingdom_data$Kingdom, levels = level_kingdom)
  build_OTU_counts_output_phylum_data$Phylum <- factor(build_OTU_counts_output_phylum_data$Phylum, levels = level_phylum)
  build_OTU_counts_output_class_data$Class <- factor(build_OTU_counts_output_class_data$Class, levels = level_class)
  build_OTU_counts_output_order_data$Order <- factor(build_OTU_counts_output_order_data$Order, levels = level_order)
  build_OTU_counts_output_family_data$Family <- factor(build_OTU_counts_output_family_data$Family, levels = level_family)
  build_OTU_counts_output_genus_data$Genus <- factor(build_OTU_counts_output_genus_data$Genus, levels = level_genus)
  build_OTU_counts_output_species_data$Species <- factor(build_OTU_counts_output_species_data$Species, levels = level_species)
  
  
  scale_plot_kingdom <- ggplot(data=build_OTU_counts_output_kingdom_data,aes(x=Sample,y=Abundance,fill=Kingdom)) + facet_grid(~Age_Group,scales="free") + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values = colors_to_use) + theme(  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  scale_plot_phylum <- ggplot(data=build_OTU_counts_output_phylum_data,aes(x=Sample,y=Abundance,fill=Phylum)) + facet_grid(~Age_Group,scales="free") + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values = colors_to_use) + theme(  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  scale_plot_class <- ggplot(data=build_OTU_counts_output_class_data,aes(x=Sample,y=Abundance,fill=Class)) + facet_grid(~Age_Group,scales="free") + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values = colors_to_use) + theme(  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  scale_plot_order <- ggplot(data=build_OTU_counts_output_order_data,aes(x=Sample,y=Abundance,fill=Order)) + facet_grid(~Age_Group,scales="free") + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values = colors_to_use) + theme(  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  scale_plot_family <- ggplot(data=build_OTU_counts_output_family_data,aes(x=Sample,y=Abundance,fill=Family)) + facet_grid(~Age_Group,scales="free") + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values = colors_to_use) + theme(  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  scale_plot_genus <- ggplot(data=build_OTU_counts_output_genus_data,aes(x=Sample,y=Abundance,fill=Genus)) + facet_grid(~Age_Group,scales="free") + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values = colors_to_use) + theme(  axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #scale_plot_species <- ggplot(data=build_OTU_counts_output_species_data,aes(x=Sample,y=Abundance,fill=Species)) + facet_grid(~Age_Group,scales="free") + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values = colors_to_use) + theme( legend.position = "", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  #scale_plot_species2 <- ggplot(data=build_OTU_counts_output_species_data,aes(x=Sample,y=Abundance,fill=Species)) + facet_grid(~Age_Group,scales="free") + geom_bar(aes(), stat="identity", position="stack") + scale_fill_manual(values = colors_to_use) + theme( legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
  ggsave("scale_plot_kingdom.pdf", plot = scale_plot_kingdom, width = 15, height = 10)
  ggsave("sscale_plot_phylum.pdf", plot = scale_plot_phylum, width = 15, height = 10)
  ggsave("scale_plot_class.pdf", plot = scale_plot_class, width = 15, height = 10)
  ggsave("scale_plot_order.pdf", plot = scale_plot_order, width = 15, height = 10)
  ggsave("scale_plot_family.pdf", plot = scale_plot_family, width = 15, height = 10)
  ggsave("scale_plot_genus.pdf", plot = scale_plot_genus, width = 15, height = 10)
  #ggsave("scale_plot_species.pdf", plot = scale_plot_species, width = 15, height = 10)
  #ggsave("scale_plot_species2.pdf", plot = scale_plot_species2, width = 15, height = 10)
  
  
  
  ########## Bidirectional Plot #######
  
  bidirectional_plot(build_OTU_counts_output, "Kingdom")
  bidirectional_plot(build_OTU_counts_output, "Phylum")
  bidirectional_plot(build_OTU_counts_output, "Class")
  bidirectional_plot(build_OTU_counts_output, "Order")
  bidirectional_plot(build_OTU_counts_output, "Family")
  bidirectional_plot(build_OTU_counts_output, "Genus")
  bidirectional_plot(build_OTU_counts_output, "Species")
  
  
  ####################
  
  OTU_plot_tables <- list(
    
  #  scale_plot_table <- data_phylo2_class,
  #  AlphaDiv_plot_table <- species_shannon_diversity,
   # PCoA_table <- pcoa$vectors
  )
  
  
  return(OTU_plot_tables)
}






