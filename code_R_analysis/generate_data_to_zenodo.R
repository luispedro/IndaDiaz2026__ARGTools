library(dplyr)
library(tidyverse)


#setwd("~/Documents/GitHub/arg_compare/code_R_analysis") 
options(dplyr.summarise.inform = FALSE)
# source("helper.R")


# results per individual tool 
lst <- readRDS("code_R_analysis/output_abundance_diversity_resistome/results_tools_all_GMGC.rds")

# metagenomes' metadata
metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")
metadata <- metadata %>% select(sample_id, habitat, insertsHQ)

# abundance <- readRDS("output_abundance_diversity_resistome/abundance_diversity.rds")
abundance <- readRDS("code_R_analysis/output_abundance_diversity_resistome/abundance_diversity.rds")
abundance <- abundance %>% rename(geneclass = gene, pipeline = tool)

conversion_file <- readRDS("code_R_analysis/output_abundance_diversity_resistome/conversion_ARO_parent_new_level.rds")
conversion_file <- conversion_file %>% rename(geneclass = new_level)


lst2 <- lapply(lst, function(x) x <- x %>% rename(geneclass_argcompare = new_level,
                                                  pipeline = tool))
lst2 <- lst2[1:15]

lapply(names(lst2), function(nm) {
  write.csv(lst2[[nm]], paste0("code_R_analysis/data_to_Zenodo/pipelines_output/",nm,".csv"), row.names = FALSE)
})


write.csv(abundance, file = gzfile("code_R_analysis/data_to_Zenodo/abundance_richness.csv.gz"), row.names = F)
write.csv(metadata, file = "code_R_analysis/data_to_Zenodo/metagenomes_metadata.csv", row.names = F)
write.csv(conversion_file, file = "code_R_analysis/data_to_Zenodo/conversion_aro_geneclass.csv", row.names = F)


reported_unigenes_as_ARG_per_habitat <- read.delim("code_R_analysis/output_abundance_diversity_resistome/reported_unigenes_as_ARG_per_habitat.csv", sep = ",")
reported_unigenes_as_ARG_per_habitat <- reported_unigenes_as_ARG_per_habitat %>% rename(unigene = X)

write.csv(reported_unigenes_as_ARG_per_habitat, file = "code_R_analysis/data_to_Zenodo/reported_unigenes_as_ARG_per_habitat.csv", row.names = F)



