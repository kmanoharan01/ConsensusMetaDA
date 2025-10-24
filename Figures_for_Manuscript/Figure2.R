
############### Figure 2: Differential abundance UpSet and Venn Diagram ##################


#1 Marine plastics
biome_file <- "/Users/jd116080/myprojects/Matt/MicrobiomeConsensusDE/MicrobiomeConsensusDE/TestData/Hackathon/Studies/sed_plastic_rosato/sed_plastic_rosato_ASVs_table.biom"

sample_table_file <-  "/Users/jd116080/myprojects/Matt/MicrobiomeConsensusDE/MicrobiomeConsensusDE/TestData/Hackathon/Studies/sed_plastic_rosato/sed_plastic_rosato_metadata.tsv"

## Marine plastics
sediment <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)

## Marine plastics Differential abundance
sediment_DA <- OTUs_multi_DA(sediment)

# An example vignette is provided (https://github.com/kmanoharan01/ConsensusMetaDA/blob/main/inst/doc/ConsensusMetaDA_Manual.html)
