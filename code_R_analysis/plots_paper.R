library(dplyr)
library(ggplot2)
#library(plotly)
library(ggVennDiagram)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
#library(ggpattern)
library(grid)
library(eulerr)
library(reactable)
library(Cairo)


setwd("~/Documents/GitHub/ARG_tools_shiny/app") 
options(dplyr.summarise.inform = FALSE)

df2 <- readRDS(file = "data/conversion_ARO_parent_new_level.rds")
lst <- readRDS("data/results_tools.rds")
metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")
abundance <- readRDS("data/abundance_diversity.rds")
abundance <- abundance %>% mutate(unigenes = unigenes_rarefied)
core <- readRDS(file = "~/Documents/GitHub/arg_compare/code_R_analysis/output_abundance_diversity_resistome/core_resistome.rds")
core <- core %>% rename(new_level = new_level_centroid, X = centroid)
pan <- readRDS(file = "data/pan_resistome.rds")

sum(table(abundance$habitat[!duplicated(abundance$sample)]))
table(abundance$habitat[!duplicated(abundance$sample)])

write.csv(abundance %>% ungroup() %>% select(sample, habitat) %>% filter(!habitat %in% c("amplicon","built-environment","isolate")) %>%
  rename(accession = sample) %>% group_by(accession) %>% slice_head(n = 1), "~/Documents/plots_project3/accession_metagenomes.csv", row.names = FALSE)

write.csv(abundance %>% ungroup() %>% select(sample, habitat) %>% filter(!habitat %in% c("amplicon","built-environment","isolate")) %>%
            rename(accession = sample) %>% group_by(accession) %>% slice_head(n = 1) %>% 
            ungroup() %>% group_by(habitat) %>% summarise(samples = n()), "~/Documents/plots_project3/summary_metagenomes.csv", row.names = FALSE)


# name of gene classes to show
gene_classes = data.frame(rbind(
  c("glycopeptide resistance (van)" , "van"), 
  c("protein(s) and two-component regulatory system modulating antibiotic efflux" , "efflux pump"),
  c("gene altering cell wall charge" , "cell wall charge"), 
  c("rifamycin-resistant beta-subunit of RNA polymerase (rpoB)" , "rpoB"),
  c("class A beta-lactamase" , "class A beta-lactamase"), 
  c("class B beta-lactamase" , "class B beta-lactamase"), 
  c("class C beta-lactamase" , "class C beta-lactamase"), 
  c("class D beta-lactamase" , "class D beta-lactamase"),
  c("gene modulating beta-lactam resistance" , "beta-lactam modulation resistance"),
  c("aminoglycoside acetyltransferase (aac)" , "aac"), 
  c("aminoglycoside phosphotransferase (aph)" , "aph"), 
  c("aminoglycoside nucleotidyltransferase (ant)" , "ant"), 
  c("aminoglycoside bifunctional resistance protein" , "bifunctional aminoglycoside"),
  c("tetracycline-resistant ribosomal protection protein (tet RPG)" , "tet RPG"), 
  c("tetracycline inactivation enzyme (tet enzyme)" , "tet enzyme"),
  c("major facilitator superfamily antibiotic efflux pump (MFS efflux pump)" , "MFS efflux pump"), 
  c("erm 23S ribosomal RNA methyltransferase (erm)" , "erm"), 
  c("macrolide phosphotransferase (mph)"  , "mph"), 
  c("quinolone resistance protein (qnr)" , "qnr"),
  c("beta-lactam resistant penicillin-binding proteins (PBP)" , "PBP"), 
  c("antibiotic target modifying enzyme" , "target-modifying enzyme"),
  c("ciprofloxacin phosphotransferase (crpP)" , "crpP"), 
  c("rifampin-resistant RNA polymerase-binding protein (rifampin Rbp)" , "rifampin Rbp"),
  c("chloramphenicol phosphotransferase (cpt)" , "cpt"), 
  c("chloramphenicol acetyltransferase (cat)" , "cat"),
  c("sulfonamide resistant (sul)" , "sul"), 
  c("viomycin phosphotransferase (vph)" , "vph"),
  c("antibiotic resistant gene variant or mutant (variant or mutant)" , "variant or mutant"), 
  c("fosfomycin inactivation enzyme (fos)" , "fos"),
  c("antibiotic inactivation enzyme" , "antibiotic inactivation enzyme"), 
  c("antibiotic resistant dihydrofolate reductase (dfr)" , "dfr"),
  c("streptothricin acetyltransferase (sat)" , "sat"), 
  c("nitroimidazole reductase (nim)" , "nim"), 
  c("ABC-F ATP-binding cassette ribosomal protection protein (abcF)" , "abcF"), 
  c("protein modulating permeability to antibiotic" , "permeability modulation"),
  c("protein(s) conferring resistance via host-dependent nutrient acquisition" , "host-dependent nutrient acquisition"), 
  c("gene involved in antibiotic sequestration" , "antibiotic sequestration"),
  c("lincosamide nucleotidyltransferase (lnu)" , "lnu"), 
  c("macrolide glycosyltransferase (mgt)" , "mgt"),
  c("fusidic acid inactivation enzyme (fai)" , "fai"), 
  c("target protecting FusB-type protein conferring resistance to Fusidic acid (fusB-type)" , "fusB-type"), 
  c("macrolide esterase (mel)" , "mel"), 
  c("bah amidohydrolase (bah)" , "bah"),
  c("cpa acetyltransferase (cpa)" , "cpa"), 
  c("gene involved in self-resistance to antibiotic" , "self-resistance"),
  c("streptogramin inactivation enzym (vat)" , "vat"), 
  c("gene conferring resistance via absence" , "resistance by absence"),
  c("capreomycin phosphotransferase (cph)" , "cph"), 
  c("edeine acetyltransferase (edeQ)" , "edeQ"), 
  c("protein(s) conferring antibiotic resistance via molecular bypass" , "molecular bypass"),
  c("rifampin inactivation enzyme", "rifampin inactivation enzyme")))
colnames(gene_classes) <- c("old", "new")


gene_classes <- setNames(as.list(gene_classes$new), gene_classes$old)

# FORMAT PLOTS

pal_6 <- brewer.pal(6, "Dark2")
pal_10_q <- brewer.pal(10, "Paired")
pal_10_q <- pal_10_q[c(6,2,8,7,4,3,10,9,1,5)]
pal_6 <- pal_10_q[c(1, 3, 5, 9, 10, 8)]
general_size <- 5

pal8 <- brewer.pal(8, "Dark2")
pal8 <- pal8[c(1:4,6:8)]
pal8 <- rev(pal8)
pal8 <- c(pal8[1:5],"#7570B399", pal8[6], "#D95F0299", pal8[7], "#1B9E7799")

pal_10_q<- pal8

# HABITATS
EN <- c("human gut", "human oral",  "human skin", "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine", 
        "freshwater", "soil" , "amplicon", "isolate",  "built-environment" )

# SOURCE FOR EACH HABITAT
SO <- c(rep("humans", 5), rep("mammals", 4),  
        "wastewater", "marine", "freshwater", 
        "soil", rep("other", 2), "built-environment")

names(SO) <- EN


tools_levels_2 <- c("DeepARG", "fARGene",
                    "ABRicate-ARGANNOT", "ABRicate-MEGARes",
                    "RGI-DIAMOND", "ABRicate-CARD",
                    "AMRFinderPlus", "ABRicate-NCBI",
                    "ResFinder", "ABRicate-ResFinder")

# environments that we are not interested in
not_env <- c("amplicon", "isolate", "built-environment" )

EN2  <- EN[!EN %in% not_env]
h2 <- c("humans","mammals","wastewater", "freshwater", "soil", "marine")

source("helper.R")


abundance <- fix_abundance(abundance, tools_levels_2, not_env, EN, SO, h2)
abundance_tool_sample <- complete_per_sample_abundance_diversity(abundance, tools_levels_2, EN, general_size, pal_10_q, tool_label, tools_levels_2)
abundance_class <- fix_abundance_class(abundance, not_env, tools_levels_2, metadata) 
unigenes <- fix_unigenes(lst)

unigenes <- as_tibble(unigenes) %>% 
  filter(tool %in% tools_levels_2) %>%
  mutate(tool = factor(tool, levels = tools_levels_2))


core <- fix_core(core, not_env, tools_levels_2, EN2)

core <- core %>% 
  mutate(tool = factor(tool, levels = tools_levels_2))

pan <- fix_pan(pan, not_env, tools_levels_2, EN2)

pan <- pan %>% 
  mutate(tool = factor(tool, levels = tools_levels_2))


sumpan2 <- summarize_pan(pan)
recall_fnr <- create_recall_fnr(unigenes)

JI_all <- return_overlap_tools(unigenes)

top20 <- c("van", "efflux pump", "cell wall charge", "rpoB", "tet RPG", "class A beta-lactamase", "class B beta-lactamase", "class C beta-lactamase", "class D beta-lactamase", 
           "aac", "aph", "MFS efflux pump", "tet enzyme", "erm", "mph", "target-modifying enzyme", "qnr")



sumpan2 %>% ungroup() %>% group_by(habitat) %>% mutate(deep = mn[tool == "DeepARG"][1], rgi = mn[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(!tool %in% c("DeepARG", "RGI-DIAMOND")) %>%
  mutate(over_deepARG = deep / mn, over_rgi = rgi / mn) %>%
  group_by(habitat) %>%
  summarise(min_deep = min(over_deepARG), max_deep = max(over_deepARG), 
            min_rgi = min(over_rgi), max_rgi = max(over_rgi)) %>% 
  ungroup() %>%
  summarise(min_deep = min(min_deep), max_deep = max(max_deep), 
            min_rgi = min(min_rgi), max_rgi = max(max_rgi))


habitat %in% c("wastewater", "marine", "freshwater", "soil")

sumpan2 %>% ungroup() %>% group_by(habitat) %>% mutate(deep = mn[tool == "AMRFinderPlus"][1]) %>% 
  ungroup() %>%
  filter(tool %in% c("ABRicate-CARD")) %>%
  mutate(over_deepARG = deep / mn) %>%
  summarise(mn = mean(over_deepARG))

sumpan2 %>% ungroup() %>% group_by(habitat) %>% mutate(deep = mn[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(tool %in% c("ABRicate-CARD")) %>%
  mutate(over_deepARG = deep / mn) %>%
  summarise(mn = mean(over_deepARG))


sumpan2 %>% ungroup() %>% group_by(habitat) %>% mutate(deep = mn[tool == "ResFinder"][1]) %>% 
  ungroup() %>%
  filter(tool %in% c("ABRicate-ResFinder")) %>%
  mutate(over_deepARG = deep / mn) %>%
  summarise(mn = mean(over_deepARG))



sumcore <- sum_core_adjust(core, 250, 0.2)
sumcore %>% 
  ungroup() %>%
  group_by(tool, habitat) %>% 
  summarise(unigenes = sum(unigenes)) %>%
  group_by(habitat) %>% mutate(deep = unigenes[tool == "DeepARG"][1], rgi = unigenes[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(!tool %in% c("DeepARG", "RGI-DIAMOND")) %>%
  filter(!habitat %in% c("wastewater", "marine", "freshwater", "soil", "human oral", "human skin", "human vagina", "human nose", "human oral")) %>% 
  mutate(over_deepARG = deep / unigenes, over_rgi = rgi / unigenes) %>%
  group_by(habitat) %>%
  summarise(min_deep = min(over_deepARG, na.rm = T), max_deep = max(over_deepARG, na.rm = T), 
            min_rgi = min(over_rgi, na.rm = T), max_rgi = max(over_rgi, na.rm = T)) %>% 
  ungroup() %>%
  summarise(min_deep = min(min_deep, na.rm = T), max_deep = max(max_deep, na.rm = T), 
            min_rgi = min(min_rgi, na.rm = T), max_rgi = max(max_rgi, na.rm = T))


sumcore %>% 
  ungroup() %>%
  group_by(tool, habitat) %>% 
  summarise(unigenes = sum(unigenes)) %>%
  group_by(habitat) %>% mutate(deep = unigenes[tool == "DeepARG"][1], rgi = unigenes[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(!tool %in% c("DeepARG", "RGI-DIAMOND")) %>%
  filter(habitat %in% c("human gut")) %>% 
  mutate(over_deepARG = deep / unigenes, over_rgi = rgi / unigenes) 


sumcore %>% 
  ungroup() %>%
  group_by(habitat) %>% mutate(deep = unigenes[tool == "DeepARG"][1], rgi = unigenes[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(habitat %in% c("pig gut"), new_level %in% "van")

sumcore %>% 
  ungroup() %>%
  group_by(tool, habitat) %>% 
  summarise(unigenes = sum(unigenes)) %>%
  group_by(habitat) %>% mutate(deep = unigenes[tool == "DeepARG"][1], rgi = unigenes[tool == "RGI-DIAMOND"][1]) %>% 
  ungroup() %>%
  filter(habitat %in% c("pig gut")) %>% 
  filter(!tool %in% c("DeepARG", "RGI-DIAMOND")) %>%
  mutate(over_deepARG = deep / unigenes, over_rgi = rgi / unigenes) 

sumcore %>% 
  ungroup() %>%
  group_by(tool, habitat) %>% 
  summarise(unigenes = sum(unigenes)) %>% 
  filter(habitat %in% "marine")

################################################################################################
# Fig 1

p1a <- plot_count_genes_tool(unigenes , tools_levels_2, general_size, pal_10_q, tools_levels_2, tools_levels_2) + 
  ggtitle("A") + theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold"))

p1b <- return_heatmap_overalp(JI_all, as.character(unique(JI_all$tool_ref)), general_size, tools_levels_2, unique(as.character(JI_all$tool_ref))) + ggtitle("B") + theme(title = element_text(size = general_size + 2, face = "bold"))

p2a <- plot_total_abundance_diversity_new_version(abundance_tool_sample, tools_levels_2, h2, general_size, pal_10_q , "abundance", sd = 2025, obs = 200)
p2b <- plot_total_abundance_diversity_new_version(abundance_tool_sample, tools_levels_2, h2, general_size, pal_10_q , "diversity", sd = 2025, obs = 200)


ggsave("~/Documents/plots_project5/fig1c.svg", p1a, width = 60, height = 75, unit = "mm")
ggsave("~/Documents/plots_project5/figS19.svg", p1b + ggtitle(""), width = 90, height = 75, unit = "mm")


ggsave("~/Documents/plots_project5/fig2.svg", p2a, width = 180, height = 120, unit = "mm")
ggsave("~/Documents/plots_project5/figS1.svg", p2b, width = 180, height = 120, unit = "mm")

  ################################################################################################
  # Fig 2
  
pan_core_env <- c("human gut", "pig gut",  "marine", "soil")
sumcore <- sum_core_adjust(core, 250, 0.2)

p2 <- plot_fig2(pan_resistome_plot(sumpan2, tools_levels_2, pan_core_env , general_size, pal_10_q, tools_levels_2, h2, tools_levels_2) + ggtitle("A"),
          core_resistome_plot(sumcore, tools_levels_2, pan_core_env, general_size, pal_10_q, tools_levels_2, h2, tools_levels_2) + ggtitle("B"))

ggsave("~/Documents/plots_project5/fig3.svg", p2, width = 180, height = 150, unit = "mm")

  ################################################################################################   
  # Fig 3

abundance_env_class <- abundance_medians_env_class(abundance_class, "human gut", tools_levels_2, top20)

unigenes_class <- get_unigenes_class(unigenes, tools_levels_2,  top20)

p3a1 <- plot_unigenes_env_tool_class(unigenes_class,  general_size, pal_10_q, tools_levels_2, tools_levels_2, tools_levels_2, 0, 15000, 5, T) + 
  ggtitle("") + theme(title = element_text(size = general_size + 2, face = "bold"))

p3a2 <- plot_unigenes_env_tool_class(unigenes_class,  general_size, pal_10_q, tools_levels_2, tools_levels_2, tools_levels_2, 15000, max(unigenes_class$n), 2, T) + 
  ggtitle("A") + ylab("") + theme(title = element_text(size = general_size + 2, face = "bold"), axis.text.x  = element_blank())

p3a <- grid.arrange(p3a2 + theme(plot.margin = margin(0, 0, -10, 0)), 
                    p3a1 + theme(plot.margin = margin(-10, 0, 0, 0)), 
                    layout_matrix = rbind(c(1, 1), c(2, 2), c(2, 2), c(2, 2)))
  
p3b0 <- plot_abundance_env_tool_class(abundance_env_class, tools_levels_2, "human gut", top20, general_size, pal_10_q, tools_levels_2, tools_levels_2) +
  ggtitle("B") + theme(title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank())


p3b <- p3b0 + theme(legend.position = "none")
p3c <- plot_core_env_tool_class(sumcore, tools_levels_2, "human gut", top20, general_size, pal_10_q, tools_levels_2, tools_levels_2) +
  ggtitle("C") + theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank()) +
  ylab("Number of ARGs")

p3 <- grid.arrange(p3a, p3b, p3c, g_legend(p3b0), 
                   layout_matrix = rbind(c(1, 1), c(1, 1), c(1, 1), c(1, 1), c(1, 1), 
                                         c(2, 2), c(2, 2), c(2, 2), c(2, 2), c(2, 2), 
                                         c(3, 3), c(3, 3), c(3, 3), c(3, 3), c(3, 3), c(4,4)))

ggsave("~/Documents/plots_project5/fig4.svg", p3, width = 180, height = 220, unit = "mm")


microbiomes <- c("human oral",  "human skin", "human nose", "human vagina", 
                 "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine", 
                 "freshwater", "soil")

for(microbiome in 1:length(microbiomes)){
  abundance_env_class <- abundance_medians_env_class(abundance_class, microbiomes[microbiome], tools_levels_2, top20)
  
  p3b0 <- plot_abundance_env_tool_class(abundance_env_class, tools_levels_2, microbiomes[microbiome], top20, general_size, pal_10_q, tools_levels_2, tools_levels_2) +
    ggtitle(paste0(microbiomes[microbiome], "\nA")) + theme(title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank())
  
  p3b <- p3b0 + theme(legend.position = "none")
  p3c <- plot_core_env_tool_class(sumcore, tools_levels_2, microbiomes[microbiome], top20, general_size, pal_10_q, tools_levels_2, tools_levels_2) +
    ggtitle("B") + theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank()) +
    ylab("Number of ARGs")
  
  p3 <- grid.arrange(p3b, p3c, g_legend(p3b0), layout_matrix = rbind(c(1, 1), c(1, 1), c(1, 1), c(1, 1), c(1, 1), 
                                                                     c(2, 2), c(2, 2), c(2, 2), c(2, 2), c(2, 2),
                                                                     c(3,3)))
  ggsave(paste0(paste0("~/Documents/plots_project5/figS", 6 + microbiome),".svg"), p3, width = 180, height = 150, unit = "mm")
}




  ################################################################################################
  # Fig 4
  
p4 <- plot_recall_fnr(recall_fnr, tools_levels_2, unique(unigenes$new_level), tools_levels_2, pal_10_q, general_size, tools_levels_2)

ggsave("~/Documents/plots_project5/fig5.svg", p4, width = 180, height = 80, unit = "mm")



################################################################################################
# Fig 5

#library(ggstep)

p5a_1 <- plot_id_levels(unigenes %>% filter(tool %in% c("fARGene", "AMRFinderPlus", "RGI-DIAMOND", "DeepARG", "ResFinder")), tools_levels_2, pal_10_q, tools_levels_2, general_size) + 
  ggtitle("A") + theme(legend.position = "bottom", title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank())
p5a_2 <- plot_id_levels(unigenes %>% filter(!tool %in% c("fARGene", "AMRFinderPlus", "RGI-DIAMOND", "DeepARG", "ResFinder")), tools_levels_2, pal_10_q, tools_levels_2, general_size) + 
  ggtitle("") + theme(legend.position = "bottom", title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank())

p5a <- grid.arrange(p5a_1, p5a_2 + labs(y = NULL), nrow = 1, padding = unit(0, "npc"))

p5b10 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "DeepARG", pal_10_q[3:5], general_size, tool_label, tools_levels_2) + ggtitle("B")
p5b1 <- p5b10 + theme(legend.position = "none")
p5b2 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "RGI-DIAMOND", pal_10_q[3:5], general_size, tool_label, tools_levels_2) + 
  ggtitle("C") + theme(legend.position = "none")

p5 <- grid.arrange(p5a, p5b1, p5b2, g_legend(p5b10), layout_matrix = rbind(c(1, 1), c(1, 1), c(1, 1),
                                                                        c(2, 3), c(2, 3), c(2, 3), c(2, 3), c(2, 3), c(2, 3), c(2, 3), c(2, 3), 
                                                                        c(4,4)))
ggsave("~/Documents/plots_project5/fig6.svg", p5, width = 180, height = 220, unit = "mm")
  



################################################################################################
# Fig S1

rgi.ven = ggVennDiagram(list("RGI-DIAMOND" = lst$rgi.diamond$query,
                             "RGI-BLAST" = lst$rgi.blast$query,
                             "RGI-DIAMOND (aa)" = lst$rgi.diamond.prot$query),
                        color = 1, lwd = 0.7, label_size = 2, set_size = 2) + 
  scale_fill_gradient(low = "#7570B366", high = pal_10_q[5]) +
  theme(legend.position = "none") +
  ggtitle("A") +
  theme(
    legend.text = element_text(size = general_size),
    legend.title = element_text(size = general_size + 1),
    text = element_text(size = general_size),
    axis.text  = element_blank(),
    title = element_text(size = general_size + 2, , face = "bold"))

rgi.ven


deeparg.ven = ggVennDiagram(list("DeepARG" = lst$deeparg.norm$query,
                                 "DeepARG (aa)" = lst$deeparg.norm.prot$query),
                            label_percent_digit = 1,
                            color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
  scale_fill_gradient(low = "#66666666", high = pal_10_q[1]) +
  theme(legend.position = "none") +
  scale_x_continuous(labels = ) +
  ggtitle("B") +
  theme(
    legend.text = element_text(size = general_size),
    legend.title = element_text(size = general_size + 1),
    text = element_text(size = general_size),
    axis.text  = element_blank(),
    title = element_text(size = general_size + 2, , face = "bold"))

deeparg.ven


fargene.ven = ggVennDiagram(list("fARGene" = lst$fargene$query,
                                 "fARGene (aa)" = lst$fargene.prot$query),
                            label_percent_digit = 2,
                            color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
  scale_fill_gradient(low = "#A6761D66", high = pal_10_q[2]) +  theme(legend.position = "none") +
  ggtitle("C") +
  theme(
    legend.text = element_text(size = general_size),
    legend.title = element_text(size = general_size + 1),
    text = element_text(size = general_size),
    axis.text  = element_blank(),
    title = element_text(size = general_size + 2, , face = "bold"))

fargene.ven

amrfinder.ven = ggVennDiagram(list("AMRFinderPlus" = lst$amrfinder.norm$query,
                                   "AMRFinderPlus (aa)" = lst$amrfinder.norm.prot$query),
                              label_percent_digit = 1,
                              color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
  scale_fill_gradient(low = "#D95F0266", high = pal_10_q[7]) +  theme(legend.position = "none") +
  ggtitle("D") +
  theme(
    legend.text = element_text(size = general_size),
    legend.title = element_text(size = general_size + 1),
    text = element_text(size = general_size),
    axis.text  = element_blank(),
    title = element_text(size = general_size + 2, , face = "bold"))

amrfinder.ven


ggsave("~/Documents/plots_project5/figS2.svg", grid.arrange(rgi.ven, deeparg.ven, fargene.ven, amrfinder.ven, nrow = 2), width = 140, height = 140, unit = "mm")



################################################################################################
# Fig S2



pan_core_env <- c("human skin", "human oral", "human nose", "human vagina")
sumcore <- sum_core_adjust(core, 250, 0.2)
ps2 <- plot_fig2(pan_resistome_plot(sumpan2, tools_levels_2, pan_core_env, general_size, pal_10_q, tools_levels_2, h2, tools_levels_2) + ggtitle("A"),
                 core_resistome_plot(sumcore, tools_levels_2, pan_core_env, general_size, pal_10_q, tools_levels_2, h2, tools_levels_2) + ggtitle("B"))


ggsave("~/Documents/plots_project5/figS3.svg", ps2, width = 180, height = 150, unit = "mm")


pan_core_env <- c("dog gut", "cat gut", "mouse gut")
sumcore <- sum_core_adjust(core, 250, 0.2)
ps2 <- plot_fig2(pan_resistome_plot(sumpan2, tools_levels_2, pan_core_env, general_size, pal_10_q, tools_levels_2, h2, tools_levels_2) + ggtitle("A"),
                 core_resistome_plot(sumcore, tools_levels_2, pan_core_env, general_size, pal_10_q, tools_levels_2, h2, tools_levels_2) + ggtitle("B"))

ggsave("~/Documents/plots_project5/figS4.svg", ps2, width = 180, height = 150, unit = "mm")


pan_core_env <- c("wastewater", "freshwater" )
sumcore <- sum_core_adjust(core, 250, 0.2)
ps2 <- plot_fig2(pan_resistome_plot(sumpan2, tools_levels_2, pan_core_env, general_size, pal_10_q, tools_levels_2, h2, tools_levels_2) + ggtitle("A"),
                 core_resistome_plot(sumcore, tools_levels_2, pan_core_env, general_size, pal_10_q, tools_levels_2, h2, tools_levels_2) + ggtitle("B"))


ggsave("~/Documents/plots_project5/figS5.svg", ps2, width = 180, height = 150, unit = "mm")


################################################################################################
# Fig S3


sets0 <- core %>% filter(cut %in% 0.2 & cnt > 250, habitat %in% c( "human gut" )) %>%
  group_by(tool) %>%
  summarise(query = list(X), .groups = "drop")

sets1 <- core %>% filter(cut %in% 0.2 & cnt > 250, habitat %in% c( "human gut" )) %>%
  group_by(new_level, tool) %>%
  summarise(query = list(X), .groups = "drop")

# Pairwise combinations 
pairwise <- sets1 %>%
  group_by(new_level) %>%
  summarise(pairs = list(expand_grid(tool_ref = tool, tool_comp = tool)), .groups = "drop") %>%
  unnest(pairs)

JI_core_class <- pairwise %>%
  left_join(sets1, by = c("new_level", "tool_ref" = "tool")) %>%
  rename(qc_ref = query) %>%
  left_join(sets1, by = c("new_level", "tool_comp" = "tool")) %>%
  rename(qc_comp = query) %>%
  left_join(sets0, by = c( "tool_ref" = "tool")) %>%
  rename(q_ref = query) %>% 
  left_join(sets0, by = c( "tool_comp" = "tool")) %>%
  rename(q_comp = query)


new_intersect <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <-  unlist(qc_ref)
  B <- unlist(q_ref)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  x <- intersect(A, setdiff(D, C))
  y <- setdiff(C, intersect(C, setdiff(B, A)))
  z <- union(x, y)
  d <- intersect(A, z)
  r <- ifelse(length(z) == 0, NA, length(d) / length(z))
  return(r)
}

new_difference <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <-  unlist(qc_ref)
  B <- unlist(q_ref)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  x <- intersect(A, setdiff(D, C))
  y <- setdiff(C, intersect(C, setdiff(B, A)))
  z <- union(x, y)
  d <- setdiff(A, z)
  r <- ifelse(length(A) == 0, NA, length(d) / length(A))
  return(r)
}

new_union <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <-  unlist(qc_ref)
  B <- unlist(q_ref)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  x <- intersect(A, setdiff(D, C))
  y <- setdiff(C, intersect(C, setdiff(B, A)))
  z <- union(x, y)
  d <- union(A, z)
  r <- ifelse(length(d) == 0, 0, length(d))
  return(r)
}

new_intersect2 <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <-  unlist(qc_ref)
  B <- unlist(q_ref)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  x <- intersect(A, setdiff(D, C))
  y <- setdiff(C, intersect(C, setdiff(B, A)))
  z <- union(x, y)
  d <- intersect(A, z)
  r <- ifelse(length(z) == 0, NA, length(d))
  return(r)
}


JI_core_class <- pairwise %>%
  left_join(sets1, by = c("new_level", "tool_ref" = "tool")) %>%
  rename(qc_ref = query) %>%
  left_join(sets1, by = c("new_level", "tool_comp" = "tool")) %>%
  rename(qc_comp = query) %>%
  left_join(sets0, by = c( "tool_ref" = "tool")) %>%
  rename(q_ref = query) %>% 
  left_join(sets0, by = c( "tool_comp" = "tool")) %>%
  rename(q_comp = query)

JI_core_class <- JI_core_class %>% rowwise() %>% mutate(recall = new_intersect(qc_ref, q_ref, qc_comp, q_comp))
JI_core_class <- JI_core_class %>% rowwise() %>% mutate(fnr = new_difference(qc_ref, q_ref, qc_comp, q_comp))
JI_core_class <- JI_core_class %>% rowwise() %>% mutate(union = new_union(qc_ref, q_ref, qc_comp, q_comp))
JI_core_class <- JI_core_class %>% rowwise() %>% mutate(intersect = new_intersect2(qc_ref, q_ref, qc_comp, q_comp))
JI_core_class <- JI_core_class %>% 
  mutate(ref_n_class = length(qc_ref), comp_n_class = length(qc_comp), ref_n_all = length(q_ref), comp_n_all = length(q_comp))
JI_core_class <- JI_core_class %>% 
  mutate(jaccard = intersect / union)
JI_core_class <- JI_core_class %>% 
  ungroup() %>% filter(tool_ref != tool_comp) %>% 
  select(-c(qc_ref, qc_comp, q_ref, q_comp))


data.frame(JI_core_class %>% ungroup() %>% group_by(new_level) %>% summarise(m = sum(jaccard), n = n()) %>% mutate(m = ifelse(is.na(m), 0, m / 90)) %>% arrange(desc(m)))

length(unique(paste(JI_core_class$tool_ref, JI_core_class$tool_comp)))

core_jaccard_human <- JI_core_class %>% ungroup() %>% #complete(tool, new_level, fill = list(unigenes = 0, total = 0)) %>%
  ungroup() %>% 
  filter(as.numeric(tool_ref) != as.numeric(tool_comp)) %>% 
  ggplot(aes(x = tool_ref, y = tool_comp, fill = jaccard)) +
  geom_tile() +
  scale_x_discrete(labels = tools_levels_2) +
  scale_y_discrete(labels = tools_levels_2) +
  facet_wrap(new_level ~ .) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(9, "YlOrBr"))) +
  theme_minimal() +
  labs(fill = "") +
  xlab("") +
  ylab("") +
  ggtitle("") + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = general_size ),
    panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
    #panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    strip.text = element_text(size = general_size + 1, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size))

core_jaccard_human

ggsave("~/Documents/plots_project5/figS6.svg", core_jaccard_human, width = 300, height = 300, unit = "mm") 


################################################################################################
# Fig S4




