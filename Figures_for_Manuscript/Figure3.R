################### Figure 3 #####################


## emp
biome_file <- "./emp_deblur_150bp.subset_2k.rare_5000.biom"

sample_table_file <-  "./samples_table_emp.txt"


## emp
emp <- build_OTU_counts(biom = biome_file, sample_table = sample_table_file)

emp <- OTU_plots(emp)


# An example vignette is provided (https://github.com/kmanoharan01/ConsensusMetaDA/blob/main/inst/doc/ConsensusMetaDA_Manual.html)
