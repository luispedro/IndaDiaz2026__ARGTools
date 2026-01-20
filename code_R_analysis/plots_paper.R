library(dplyr)
library(ggplot2)
#library(plotly)
library(ggVennDiagram)
library(gridExtra)
library(tidyverse)
library(RColorBrewer)
library(ggpattern)
library(grid)
library(eulerr)
library(reactable)
library(Cairo)
library(ggalluvial)

setwd("~/Documents/GitHub/arg_compare/code_R_analysis") 
options(dplyr.summarise.inform = FALSE)

source("helper.R")

# Sourced gene classes 
# gene_classes_list is a vector for labeling the plots 
# gene_classes is the list with the names from the prepared datasets and used for shiny


# FORMAT PLOTS OLD COLORS
pal_7 <- brewer.pal(8, "BrBG")
pal_7 <- pal_7[c(4,5,3,6,2,7,1)]
pal_10_q <- pal_7[c(1,2,3,4,5,5,6,6,7,7)]

pal_10_complete <- brewer.pal(10, "BrBG")
pal_10_complete <- pal_10_complete[c(-1,-10)]

# all HABITATS
EN <- c("human gut", "human oral",  "human skin", "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine", 
        "freshwater", "soil" , "amplicon", "isolate",  "built-environment" )

# SOURCE FOR EACH HABITAT
SO <- c(rep("humans", 5), rep("mammals", 4),  
        "wastewater", "marine", "freshwater", 
        "soil", rep("other", 2), "built-environment")

names(SO) <- EN

tools_levels <- c("DeepARG", "fARGene",
                    "ABRicate-ARGANNOT", "ABRicate-MEGARes",
                    "RGI-DIAMOND", "ABRicate-CARD",
                    "AMRFinderPlus", "ABRicate-NCBI",
                    "ResFinder", "ABRicate-ResFinder")

tools_labels <- c("DeepARG", "fARGene",
                    "ABRicate-ARGANNOT", "ABRicate-MEGARes",
                    "RGI-CARD", "ABRicate-CARD",
                    "AMRFinderPlus-NCBI", "ABRicate-NCBI",
                    "ResFinder", "ABRicate-ResFinder")

# we repeat the color but we add texture to the plots
tools_texture <- c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD", "ABRicate-NCBI", "ABRicate-ResFinder")

# environments that we are not interested in
not_env <- c("amplicon", "isolate", "built-environment" )

EN2  <- EN[!EN %in% not_env]
h2 <- c("humans","mammals","wastewater", "freshwater", "soil", "marine")


################################################################################

## DATASETS 

# df2 <- readRDS(file = "data/conversion_ARO_parent_new_level.rds")
# results per individual tool 
lst <- readRDS("output_abundance_diversity_resistome/results_tools.rds")

# metagenomes' metadata
metadata <- read.delim("../data/metadata_GMGC10.sample.meta.tsv")

# abundance aggreatated per ARO and per gene class considering all the detected genes per tool
# abundance <- readRDS("output_abundance_diversity_resistome/abundance_diversity.rds")
abundance <- bind_rows(readRDS("output_abundance_diversity_resistome/abundance_diversity_part1.rds"),
                       readRDS("output_abundance_diversity_resistome/abundance_diversity_part2.rds"),
                       readRDS("output_abundance_diversity_resistome/abundance_diversity_part3.rds"))

# needs to be loaded and cocatenated 3 times if the file is split
# abundance1 <- abundance[1:5327637,]
# abundance2 <- abundance[(5327637+1):10655274,]
# abundance3 <- abundance[(10655274+1):15982911,]

# saveRDS(abundance1, "output_abundance_diversity_resistome/abundance_diversity_part1.rds", compress = T)
# saveRDS(abundance2, "output_abundance_diversity_resistome/abundance_diversity_part2.rds", compress = T)
# saveRDS(abundance3, "output_abundance_diversity_resistome/abundance_diversity_part3.rds", compress = T)

# the column distinct_unigenes_rarefied <- alpha diversity (number of different genes after rarefaction) 
# no need to complete information here
abundance <- abundance %>% 
  mutate(unigenes = distinct_unigenes_rarefied) %>%  # alpha diversity (number of different genes after rarefaction) 
  mutate(habitat = factor(habitat, levels = EN), # convert to factors for ordering the plots
         habitat2 = factor(SO[habitat], levels = h2), # convert to factors for ordering the plots
         tool = factor(tool, levels = tools_levels)) %>%  # convert to factors for ordering the plots
  mutate(location = ifelse(habitat2 %in% c("humans","mammals","wastewater","built-environment"), "human-related","external")) %>% 
  mutate(location = factor(location, levels = c("human-related","external"))) %>% 
  filter(tool %in% tools_levels & !habitat %in% not_env) %>%
  filter(aggregation %in% "new_level") %>%  # take only gene class aggregation
  mutate(habitat2 = factor(as.character(habitat2), h2))

# normed10m is abundance, unigenes is the diversity (number of different genes )
# we need to complete information, we need to make sure all samples appear in all tools 

abundance_tool_sample <- abundance %>%
  group_by(tool, sample, habitat, habitat2) %>%  
  summarise(normed10m = sum(normed10m), unigenes = sum(unigenes)) %>% # sum the abundance and diversity
  ungroup() %>% 
  complete(sample, tool) %>% # complete with NAs
  left_join(abundance %>% select(sample, habitat, habitat2) %>% 
              distinct(), by = "sample") %>% # get habitat and habitat2 
  mutate(habitat  = coalesce(habitat.x, habitat.y), 
         habitat2 = coalesce(habitat2.x, habitat2.y)) %>%
  select(-habitat.x, -habitat.y, -habitat2.x, -habitat2.y) %>% 
  mutate(normed10m = replace_na(normed10m, 0)) %>%  # change NAs to 0
  mutate(unigenes = replace_na(unigenes, 0)) %>% # change NAs to 0
  arrange(tool, sample)
  
# normed10m is abundance, unigenes is the diversity (number of different genes )
# we need to complete information, we need to make sure all samples have all gene classes and appear in all tools 

abundance_class <- abundance %>% 
  ungroup() %>% 
  filter(tool %in% tools_levels, 
         !habitat %in% not_env, aggregation %in% "new_level") %>% 
  mutate(gene = factor(gene), sample = factor(sample), 
         tool = factor(tool, levels = tools_levels)) %>%
  select(!c(raw, raw_unique, scaled, distinct_unigenes_rarefied, 
            distinct_unigenes_raw, distinct_unigenes_raw_unique)) %>% 
  complete( sample, gene, tool, 
            fill = list(normed10m = 0,  unigenes = 0)) %>% ## complete for all samples, tools, and gene classes
  mutate(aggregation = "new_level", 
         gene = as.character(gene),
         habitat = metadata$habitat[match(sample, metadata$sample_id)],
         habitat2 = metadata$habitat2[match(sample, metadata$sample_id)],
         location = metadata$location[match(sample, metadata$sample_id)])  # not so dplyr way to fetch habitat and environment


## CORE RESISTOME 

core <- readRDS(file = "output_abundance_diversity_resistome/core_resistome.rds")
core <- core %>% 
  rename(new_level = new_level_centroid, 
         X = centroid) %>% 
  filter(tool %in% tools_levels, 
         !habitat %in% not_env) %>% 
  mutate(habitat = factor(habitat, levels = EN2), 
         tool = factor(tool, levels =  tools_levels)) %>% 
  mutate(tool = factor(tool, levels = tools_levels))

sumcore <- sum_core_adjust(core, 450, 0.5)

## PAN RESISTOME

pan <- readRDS(file = "output_abundance_diversity_resistome/pan_resistome.rds")
pan <- pan %>% 
  filter(tool %in% tools_levels, 
         !habitat %in% not_env, 
         aggregation %in% "new_level_centroid") %>% 
  mutate(habitat = factor(habitat, levels = EN2), 
         tool = factor(tool, levels =  tools_levels)) %>% 
  mutate(tool = factor(tool, levels = tools_levels))

# summarise pan-resistomes by tool, gene class, 
# take the median, mean, and standard deviation across subsamples

sumpan2 <- pan %>% ungroup() %>% 
  group_by(tool, habitat, aggregation, epoch) %>% 
  summarise(s = sum(unigenes)) %>%
  ungroup() %>% 
  group_by(tool, habitat, aggregation) %>% 
  summarise(md = median(s), mn = mean(s), sd = sd(s))

pan_core <- sumpan2 %>% 
  left_join((sumcore %>% 
               ungroup() %>% 
               group_by(tool, habitat) %>% 
               summarise(core = sum(unigenes))), by = c("tool", "habitat")) %>%
  mutate(core = ifelse(is.na(core), 0, core)) %>% 
  mutate(prop = core / md) %>%
  ungroup() %>% group_by(tool, habitat) %>% 
  mutate(texture = ifelse(tool %in% tools_texture, "yes", "no"))


## unigenes identified by tool

unigenes <- as_tibble(
  do.call(rbind, 
  lapply(lst, function(x) 
    x[,c("query","tool", "ARO", "parent", 
         "parent_description", "new_level", "id")]))) %>%
  filter(tool %in% tools_levels)  %>%
  mutate(tool = factor(tool, levels = tools_levels))


# Jaccard index, recall/class specific concordance, fnr/class specific non-overlap
# per class and tool
recall_fnr <- create_class_overlaps(unigenes)


# per tool
JI_all <- return_overlap_tools(unigenes)

# Jaccard index, recall/class specific concordance, fnr/class specific non-overlap
# JI_core_class <- create_class_overlaps_core(core, 0.5, 450, "human gut")

JI_core_class <- return_overlap_tools(core %>% 
  filter(cut %in% 0.5 & cnt > 450, habitat %in% "human gut") %>% 
  rename(query = X))


top20 <- c("van", "efflux pump", "cell wall charge", "rpoB", "tet RPG", "class A beta-lactamase", "class B beta-lactamase", "class C beta-lactamase", "class D beta-lactamase", 
           "aac", "aph", "MFS efflux pump", "erm", "target-modifying enzyme", "tet enzyme", "mph", "fos", "cat", "abcF")


pan_core_env <- c("human gut", "human skin", "pig gut", "dog gut", "wastewater",  "marine", "freshwater", "soil")

list_plot_core_class <- list()
list_plot_pan_class <- list()

levels_unigenes <- unigenes %>% 
  ungroup() %>% 
  group_by(query) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  group_by(new_level) %>% 
  summarise( n = n()) %>% 
  arrange(desc(n)) %>%
  ungroup() %>% 
  pull(new_level)


top20 <- c("van", "efflux pump", "cell wall charge", "rpoB", "tet RPG", "class A beta-lactamase", 
           "aph", "MFS efflux pump", "erm", "target-modifying enzyme", "tet enzyme")



################################################################################################
# PLOTS


# Fig 1

square_legends <- g_legend(
    plot_count_genes_tool_legends(
      unigenes = unigenes , tools_for_figure = tools_levels, 
      general_size = general_size, # font size
      pal_10_q = pal_10_q, # pallet
      tool_label = tools_labels,  # labels for the tools
      tools_levels = tools_levels) +  # factor names of the tools
    ggtitle("A") + theme( title = element_text(size = general_size + 2, face = "bold")))

p_counts <- plot_count_genes_tool(
  unigenes = unigenes, 
  tools_for_figure = tools_levels, 
  general_size = general_size, 
  pal_10_q = pal_10_q, 
  tool_label = tools_labels , 
  tools_levels = tools_levels, 
  texture = tools_texture) # which tools have stripes

square_texture_legends <- g_legend( p_counts + guides(pattern = "none"))

p1a <- p_counts + 
  ggtitle("A") + theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold")) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) +
  theme(panel.background = element_rect(colour = "black", fill = NA))

p1b <- return_heatmap_overalp(
  JI_all = JI_all, # overlap between tools
  tools_selected = as.character(unique(JI_all$tool_ref)), # tools selected to be plotted 
  general_size = general_size, # font size
  tool_label = tools_labels, # tool names in the plots
  tools_levels = unique(as.character(JI_all$tool_ref))) + # factor names 
  ggtitle("B") + theme(title = element_text(size = general_size + 2, face = "bold"))


p2a <- plot_total_abundance_diversity_new_version(
  dataset = abundance_tool_sample, # 
  tools_labels = tools_labels,  #
  tools_to_plot = tools_levels,  #
  environments_plot = c("human gut", "human skin", "dog gut", "pig gut", "wastewater", "marine", "freshwater", "soil"), # habitats to plot (aggregated humans and mammals)
  general_size = general_size, # font size
  pal_10_q = pal_10_q , # pallet
  metric = "abundance", # metric (abundance or diversity)
  sd = 2025, # seed to plot random samples in the distribution 
  obs = 200,  # number of samples to plot as dots per environment
  texture = tools_texture) # texture for repeated color 

p2a2 <- plot_total_abundance_diversity_new_version(
  dataset = abundance_tool_sample, # 
  tools_labels = tools_labels,  #
  tools_to_plot = tools_levels,  #
  environments_plot = EN[1:13], # habitats to plot (aggregated humans and mammals)
  general_size = general_size, # font size
  pal_10_q = pal_10_q , # pallet
  metric = "abundance", # metric (abundance or diversity)
  sd = 2025, # seed to plot random samples in the distribution 
  obs = 200,  # number of samples to plot as dots per environment
  texture = tools_texture) # texture for repeated color 

p2diversity <- plot_total_abundance_diversity_new_version(
  dataset = abundance_tool_sample, # 
  tools_labels = tools_labels,  #
  tools_to_plot = tools_levels,  #
  environments_plot = EN[1:13], # habitats to plot (aggregated humans and mammals)
  general_size = general_size, # font size
  pal_10_q = pal_10_q , # pallet
  metric = "diversity", # metric (abundance or diversity)
  sd = 2025, # seed to plot random samples in the distribution 
  obs = 200,  # number of samples to plot as dots per environment
  texture = tools_texture) # texture for repeated color 

p2b <- plot_total_abundance_diversity_new_version(
  abundance_tool_sample, 
  tools_labels, 
  tools_levels, 
  c("human gut", "human skin", "dog gut", "pig gut", "wastewater", "marine", "freshwater", "soil"), 
  general_size, 
  pal_10_q , 
  "diversity", 
  sd = 2025, 
  obs = 200, 
  tools_texture)

################################################################################################
# Fig 2
  
shape_tools <- rep(21, length(tools_labels))
shape_tools[tools_labels %in% tools_texture] <- 24

circle_triangle_legends <- g_legend(pan_core %>% 
                             ggplot(aes(x = md, y = core +1, fill = tool,  shape = texture)) +
                             geom_point(size = 3, color = "black") +
                             scale_fill_manual(values = pal_10_q, labels = lab_fn(tools_labels), name = NULL) + 
                             scale_shape_manual(values = c(21, 24)) +
                             guides(
                                 fill = guide_legend(
                                   override.aes = list(
                                     shape = shape_tools,
                                     fill  = pal_10_q)), shape = "none") +
                             theme_minimal() +
                             theme(legend.position = "bottom",
                                   legend.text = element_text(size = general_size)))


circle_legends <- g_legend(pan_core %>% 
  ggplot(aes(x = md, y = core +1, fill = tool)) +
    geom_point(size = 3, color = "black", shape = 21) +
    scale_fill_manual(values = pal_10_q, labels = tools_labels) + 
    labs(fill = "") +
    theme(legend.position = "bottom",
      legend.text = element_text(size = general_size)))

triangle_legends <- g_legend(pan_core %>% 
  ggplot(aes(x = md, y = core +1, fill = tool)) +
    geom_point(size = 3, color = "black", shape = 24) +
    scale_fill_manual(values = pal_10_q, labels = tools_labels) + 
    labs(fill = "") +
    theme(legend.position = "bottom",
      legend.text = element_text(size = general_size)))


p_pan_core <- dot_plot_pan_core(
    pan_core, 
    pan_core_env, 
    tools_labels, 
    tools_levels, 
    tools_levels, 
    tools_texture, 
    general_size, 
    pal_10_q) + 
  ggtitle("A") + theme(legend.position = "none")


core_all <- pan_core %>% 
  filter(habitat %in% pan_core_env) %>% 
  ggplot(aes(x = tool, y =  core)) +
  geom_jitter(aes(color = tool, shape = texture), size = 2) + 
  facet_wrap(~ habitat, scales = "free_x", nrow = 1) +
  scale_color_manual(values = pal_10_q) +
  theme_minimal() +
  xlab("") +
  ylab("Core-resistome") + 
  theme(
    legend.position = "none",
    strip.text = element_text(size = general_size),
    legend.text = element_text(size = general_size ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = general_size),
    panel.border = element_blank(),   
    panel.background = element_rect(colour = "black", fill = NA))


pan_all <- pan_core %>% 
  filter(habitat %in% pan_core_env) %>% 
  ggplot(aes(x = tool, y =  mn)) +
  geom_jitter(aes(color = tool, shape = texture), size = 2) + 
  facet_wrap(~ habitat, scales = "free_x", nrow = 1) +
  scale_color_manual(values = pal_10_q) +
  theme_minimal() +
  xlab("") +
  ylab("Pan-resistome") + 
  theme(
    legend.position = "none",
    strip.text = element_text(size = general_size),
    legend.text = element_text(size = general_size ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = general_size),
    panel.border = element_blank(),   
    panel.background = element_rect(colour = "black", fill = NA))


pan_core_all <- pan_core %>% select(!c(md, sd)) %>% pivot_longer(cols = c(mn, core), names_to = "metric", values_to = "value") %>% 
  mutate(metric = ifelse(metric %in% "mn", "Pan-resistome", metric)) %>%
  mutate(metric = ifelse(metric %in% "core", "Core-resistome", metric)) %>%
  mutate(metric = factor(metric, levels = c("Pan-resistome", "Core-resistome"))) %>%
  filter(habitat %in% c("human gut", "pig gut", "wastewater", "marine", "freshwater", "soil")) %>% 
  ggplot(aes(x = habitat, y =  value)) +
  geom_jitter(aes(fill = tool, shape = texture),  color = "black", stroke = 0.3, size = 2.5, width = .5) + 
  facet_grid(metric ~ habitat, scales = "free") +
  scale_fill_manual(values = pal_10_q) +
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
  theme_minimal() +
  xlab("") +
  ylab("ARGs") + 
  theme(
    legend.position = "none",
    strip.text.x   = element_text(size = general_size),
    #strip.text.y   = element_blank(),
    legend.text = element_text(size = general_size),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = general_size),
    panel.border = element_blank(),   
    panel.background = element_rect(colour = "black", fill = NA))

pan_core_all

fig2_d <- plot_grid( pan_all, core_all, ncol = 1,
           align = "v",
           axis = "l")





p_pan_core_supp <- dot_plot_pan_core(
  pan_core, 
  EN2[!EN2 %in% pan_core_env], 
  tools_labels, 
  tools_levels, 
  tools_levels, 
  tools_texture, 
  general_size, 
  pal_10_q) + 
  ggtitle("A") + theme(legend.position = "none")

for( e in EN2){
  plot_list <- plot_pan_core_class(
    pan = pan, # pan resistome data
    sumcore = sumcore,  # core resistome data
    top20 = top20, # gene classes to plot, the rest go to Other
    tools_labels = tools_labels, #
    tools_to_plot = tools_levels, #
    tools_levels = tools_levels, #
    environments_plot = e,  #
    general_size = general_size, 
    pal_10_q = pal_10_q, 
    tools_texture = tools_texture)  
  
  list_plot_core_class[[e]] <- plot_list[[1]]
  list_plot_pan_class[[e]] <- plot_list[[2]]
}

p_alluvial  <- plot_alluvial_classes(unigenes, levels_unigenes, 0.95, 0.01, tools_levels, tools_labels, tools_levels, pal_10_q, general_size, gene_classes_list) 
p_alluvial  <- p_alluvial + theme(panel.background = element_rect(colour = "black", fill = NA)) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) 



##########
##########
##########

pan_core_human <- bind_rows(
  sumcore %>% filter(habitat %in% "human gut") %>% ungroup() %>% mutate(dt = "core"),
  pan %>% 
    filter(habitat %in% "human gut", aggregation %in% "new_level_centroid") %>% 
    ungroup() %>% 
    group_by(tool, habitat, epoch, gene_class) %>% 
    summarise(s = sum(unigenes)) %>%
    ungroup() %>% 
    group_by(tool, habitat, gene_class) %>% 
    summarise(unigenes = ceiling(mean(s))) %>% 
    rename(new_level  = gene_class) %>% mutate(dt = "pan")) %>%
  mutate(dt = factor(dt, levels = c("pan", "core")))

pan_core_human2 <- pan_core_human %>% ungroup() %>%  mutate(new_level = as.character(new_level)) %>% 
  #mutate(new_level = ifelse(new_level %in% top20, new_level, "Other")) %>%
  #mutate(new_level = factor(new_level, levels = c(top20, "Other"))) %>% 
  group_by(tool, new_level, dt) %>% summarise( n = sum(unigenes)) %>% 
  arrange(tool, desc(n))

pan_core_human3 <- pan_core_human2 %>% 
  ungroup() %>%
  group_by(tool, dt) %>% 
  mutate(proportion = n / sum(n)) %>% 
  arrange(tool, desc(proportion)) %>% 
  mutate(cum_p = cumsum(proportion)) %>% 
  ungroup() %>%
  mutate( new_level = factor(as.character(new_level), 
                             levels = c("Other", rev(levels_unigenes))))


pan_core_human4 <- pan_core_human3 %>%
  group_by(tool, dt) %>%
  arrange(proportion) %>% 
  { 
    if (any(.$cum_p == 0.95)) {
      filter(., cum_p <= 0.95)
    } else {
      lower_part <- filter(., cum_p < 0.99999)
      upper_part <- filter(., cum_p > 0.99999) %>% slice_tail(n = 1)
      bind_rows(lower_part, upper_part)
    }
  } %>% 
  ungroup() %>% 
  arrange(tool, cum_p)


pan_core_human5 <- pan_core_human4 %>% 
  filter(proportion > 0.01) %>% 
  ungroup() %>% 
  complete(new_level, tool, dt,
           fill = list(n = 0,  proportion = 0.000, cum_p = 0.000))  %>%  
  mutate(gene_name = gene_classes_list[as.character(new_level)]) %>%  
  mutate(gene_name = ifelse(new_level %in% "Other", "Other", gene_name)) %>% 
  mutate(gene_name = factor(gene_name, 
                            levels = c("Other",gene_classes_list[levels(pan_core_human4$new_level)])))

p_alluvial_pan_core <- pan_core_human5 %>% 
  filter(tool %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")) %>%
  mutate(tool = factor(as.character(tool), levels = c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder"))) %>%
  ggplot(
         aes(x = dt,
             stratum = gene_name,
             alluvium = gene_name,
             y = proportion,
             fill = gene_name,
             label = gene_name)) +
  geom_flow(alpha = 0.5) +
  xlab("") + 
  geom_stratum(color = "black") +
  geom_text(data = pan_core_human5[pan_core_human5$proportion > 0.01 & pan_core_human5$tool %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder"), ],  
            stat = "stratum",
            size = 2,color = "black",hjust = 0.5,
            position = position_jitter(width = 0, height = 0)) +
  scale_y_continuous(expand = c(.01, .01), 
                     name = "Proportion", 
                     limits = c(0, 1), 
                     breaks = seq(from = 0, to = 1, by = 0.1)) +
  facet_wrap(tool ~ ., nrow = 1) +
  scale_fill_manual(values = c(rep(pal_10_complete,10))) +
  xlab("Pan- vs core-resistome") + 
  #scale_x_discrete( labels =  labels_plot) +
  theme_minimal() +
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = general_size),
    strip.text = element_text(size = general_size, face = "bold"),
    panel.background = element_rect(colour = "black", fill = NA))


p_alluvial_pan_core

plot3_new <- plot_grid( p_alluvial + ggtitle("A") + scale_x_discrete(labels = lab_fn) ,
           p_alluvial_pan_core + ggtitle("B"),
           abu_class_human + ggtitle("C"),
           ncol = 1,
           align = "v",
           axis = "lr", 
           rel_heights = c(7,7,4))


ggsave("output_plots/fig3.svg", plot3_new, width = 180, height = 210, unit = "mm")




#########
#########
#########



p42 <- plot_recall_fnr_texture(recall_fnr, tools_levels, unique(unigenes$new_level), tools_labels, pal_10_q, general_size, tools_levels, tools_texture)


p5b1 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "DeepARG", pal_10_q[3:5], general_size, tool_label, tools_levels_2) + 
  ggtitle("DeepARG")

p5b10 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "DeepARG", pal_10_q[3:5], general_size, tool_label, tools_levels_2) + 
  ggtitle("DeepARG") + 
  theme(legend.position = "none")

p5b2 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "RGI-DIAMOND", pal_10_q[3:5], general_size, tool_label, tools_levels_2) +
  ggtitle("RGI-CARD") + 
  theme(legend.position = "none")

p5b3 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "fARGene", pal_10_q[3:5], general_size, tool_label, tools_levels_2) +
  ggtitle("fARGene") + 
  theme(legend.position = "none")

p5b4 <- one_tool_recall_fnr_id(unigenes, recall_fnr, "AMRFinderPlus", pal_10_q[3:5], general_size, tool_label, tools_levels_2) +
  ggtitle("AMRFinderPlus-NCBI") + 
  theme(legend.position = "none")

p5 <- grid.arrange(p5b10, p5b2, g_legend(p5b1), 
                   layout_matrix = rbind(c(1, 2), c(1, 2), c(1, 2), c(1, 2), 
                                         c(1, 2), c(1, 2), c(1, 2), c(1, 2),
                                         c(4,4)))


recall_fnr %>% filter(new_level %in% top20, !is.na(recall)) %>% 
  ggplot(aes(x = tool_ref, y = recall)) +
  geom_point() +
  facet_grid(.~ new_level)


recall_fnr %>% 
  filter(new_level %in% top20, !is.na(recall)) %>% 
  mutate(texture = ifelse(tool_ref %in% tools_texture, "yes", "no")) %>%
  ggplot(aes(x = tool_ref, y = recall, fill  = tool_ref, pattern = texture)) +
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot_pattern", 
               position = position_dodge2(preserve = "single"), 
               color = "black", linewidth = 0.2, 
               pattern_color = "white",
               pattern_density = 0.1, 
               pattern_spacing = 0.025, 
               pattern_fill = "white",
               pattern_key_scale_factor = 0.6, outlier.shape = NA, outlier.size = 0) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal_10_q, 
                    labels = tools_levels) +
  facet_grid(.~ new_level)



fargene_genes <- recall_fnr %>% filter(tool_ref %in% "fARGene") %>% ungroup() %>% select(new_level) %>% distinct() %>% pull()

pal2 <- pal_10_q[match(tool_second_order, tools_levels)]
tool_second_order <- c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder", 
                       "ABRicate-CARD", "ABRicate-ResFinder", "ABRicate-ARGANNOT", "ABRicate-NCBI")


p_class_recall_1 <- recall_fnr %>% 
  filter(new_level %in% top20) %>% 
  filter(tool_ref %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND")) %>%
  mutate(tool_ref = factor(as.character(tool_ref), levels = tool_second_order)) %>% 
  mutate(texture = ifelse(tool_comp %in% tools_texture, "yes", "no")) %>%
  mutate(new_level = factor(new_level, levels = top20)) %>%
  ggplot(aes(x = new_level, fill = tool_ref, y = recall)) +
  geom_violin() +
  #geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
  geom_jitter(aes(fill = tool_comp, shape = texture), color = "black", stroke = 0.3, size = 1.5, width = 0.4,  height = 0) + 
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal2, 
                    labels = tools_levels) +
  scale_color_manual(values = pal2, 
                     labels = tools_levels) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) + 
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
  facet_grid(tool_ref ~ new_level, scales = "free_x") +
  scale_y_continuous(limits = c(-0.2,1.2), 
                     breaks = seq(0, 1, length.out = 3),
                     labels = scales::label_number()) + 
  ylab("CSTC") + 
  xlab("") + 
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size),
    panel.background = element_rect(colour = "black", fill = NA),
    strip.text = element_blank(),
    strip.background = element_blank())

p_class_recall_1

p_class_recall_2 <- recall_fnr %>% 
  filter(!is.na(recall)) %>% 
  filter(tool_ref %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND")) %>%
  mutate(texture = ifelse(tool_ref %in% tools_texture, "yes", "no")) %>%
  ungroup() %>% 
  group_by(tool_ref, new_level, texture) %>% 
  summarise(recall = median(recall)) %>%
  ggplot(aes(x = "All classes medians", fill = tool_ref, y = recall)) +
  geom_violin() +
  #geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
  geom_jitter(aes(shape = texture), size = 1.5, width = 0.4, height = 0) + 
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal_10_q[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")], 
                    labels = tools_levels[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")]) +
  scale_color_manual(values = pal_10_q[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")], 
                     labels = tools_levels[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")]) +
  facet_grid(tool_ref ~ ., scales = "free_x") +
  scale_y_continuous(limits = c(-0.2,1.2), 
                     breaks = seq(0, 1, length.out = 3),
                     labels = scales::label_number()) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) + 
  ylab("") + 
  xlab("") + 
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA),
    strip.text = element_text(size = general_size, face = "bold"),
    strip.background = element_blank())


p_class_recall_2


p_class_fnr_1 <- recall_fnr %>% 
  filter(new_level %in% top20) %>% 
  filter(tool_ref %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND")) %>%
  mutate(tool_ref = factor(as.character(tool_ref), levels = tool_second_order)) %>% 
  mutate(texture = ifelse(tool_comp %in% tools_texture, "yes", "no")) %>%
  mutate(new_level = factor(new_level, levels = top20)) %>%
  ggplot(aes(x = new_level, fill = tool_ref, y = fnr)) +
  geom_violin() +
  #geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
  geom_jitter(aes(fill = tool_comp, shape = texture), color = "black", stroke = 0.3, size = 1.5, width = 0.4, height = 0) + 
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal2, 
                    labels = tools_levels) +
  scale_color_manual(values = pal2, 
                     labels = tools_levels) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) + 
  scale_shape_manual(values = c("no" = 21, "yes" = 24)) +
  facet_grid(tool_ref ~ new_level, scales = "free_x") +
  scale_y_continuous(limits = c(-0.2,1.2), 
                     breaks = seq(0, 1, length.out = 3),
                     labels = scales::label_number()) + 
  ylab("CSNO") + 
  xlab("") + 
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size),
    panel.background = element_rect(colour = "black", fill = NA),
    strip.text = element_blank(),
    strip.background = element_blank())

p_class_fnr_1

p_class_fnr_2 <- recall_fnr %>% 
  filter(!is.na(recall)) %>% 
  filter(tool_ref %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND")) %>%
  mutate(texture = ifelse(tool_ref %in% tools_texture, "yes", "no")) %>%
  ungroup() %>% 
  group_by(tool_ref, new_level, texture) %>% 
  summarise(fnr = median(fnr)) %>%
  ggplot(aes(x = "Class medians", fill = tool_ref, y = fnr)) +
  geom_violin() +
  #geom_boxplot(outlier.shape = NA, position = position_dodge2(preserve = "single")) + 
  geom_jitter(aes(shape = texture), size = 1.5, width = 0.4, height = 0) + 
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal_10_q[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")], 
                    labels = tools_levels[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")]) +
  scale_color_manual(values = pal_10_q[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")], 
                     labels = tools_levels[tools_levels %in% c("DeepARG", "fARGene", "ABRicate-MEGARes", "RGI-DIAMOND","AMRFinderPlus", "ResFinder")]) +
  facet_grid(tool_ref ~ ., scales = "free_x") +
  scale_y_continuous(limits = c(-0.2,1.2), 
                     breaks = seq(0, 1, length.out = 3),
                     labels = scales::label_number()) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) + 
  ylab("") + 
  xlab("") + 
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA),
    strip.text = element_text(size = general_size, face = "bold"),
    strip.background = element_blank())


p_class_fnr_2



p_class_id_2 <- unigenes %>% 
  filter(tool %in% c("DeepARG",  "RGI-DIAMOND"), new_level %in% top20) %>%
  mutate(new_level = factor(new_level, levels = top20)) %>% 
  ggplot(aes(x = new_level, fill = tool, y = id )) +
  geom_violin() +
  #geom_jitter(size = 0.01, width = 0.3, alpha = 0.02, color = "grey") + 
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal_10_q[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")], 
                    labels = tools_levels[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")]) +
  scale_color_manual(values = pal_10_q[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")], 
                     labels = tools_levels[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")]) +
  facet_grid(tool ~ new_level, scales = "free_x") +
  scale_y_continuous(limits = c(20, 101), 
                     breaks = seq(20, 100, length.out = 3),
                     labels = scales::label_number()) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) + 
  ylab("% Identity") + 
  xlab("") + 
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size ),
    panel.background = element_rect(colour = "black", fill = NA),
    strip.text = element_blank(),
    strip.background = element_blank())


p_class_id_2


p_class_id_3 <- unigenes %>% 
  filter(tool %in% c("DeepARG",  "RGI-DIAMOND")) %>%
  ungroup() %>% 
  group_by(tool, new_level) %>% 
  summarise(id = median(id)) %>% 
  ggplot(aes(x = "All class medians", fill = tool, y = id )) +
  geom_violin() +
  geom_jitter(size = 1.5, width = 0.4, height = 0) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal_10_q[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")], 
                    labels = tools_levels[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")]) +
  scale_color_manual(values = pal_10_q[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")], 
                     labels = tools_levels[tools_levels %in% c("DeepARG",  "RGI-DIAMOND")]) +
  facet_grid(tool ~ ., scales = "free_x") +
  scale_y_continuous(limits = c(20,101), 
                     breaks = seq(20, 100, length.out = 3),
                     labels = scales::label_number()) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) + 
  ylab("") + 
  xlab("") + 
  theme(
    legend.position = "none",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_blank(),
    panel.background = element_rect(colour = "black", fill = NA),
    strip.text = element_text(size = general_size, face = "bold"),
    strip.background = element_blank())

p_class_id_3

plot4a <- plot_grid( p_class_recall_1 + theme(axis.text.x = element_blank()), 
                     p_class_recall_2 + theme(axis.text.x = element_blank()), 
                     ncol = 2,
                     align = "h",
                     axis = "tb",
                     rel_widths = c(3, 1)
)

plot4b <- plot_grid( p_class_fnr_1 + theme(axis.text.x = element_blank()), 
                     p_class_fnr_2 + theme(axis.text.x = element_blank()), 
                     ncol = 2,
                     align = "h",
                     axis = "tb",
                     rel_widths = c(3, 1))



plot4c <- plot_grid( p_class_id_2, 
                     p_class_id_3, 
                     ncol = 2,
                     align = "h",
                     axis = "tb",
                     rel_widths = c(3, 1))

plot4_leg <- plot_grid( circle_triangle_legends, NULL,
                     ncol = 2,
                     align = "h",
                     axis = "tb",
                     rel_widths = c(6, 1))

plot4 <- plot_grid( plot4a, plot4b, plot4c, plot4_leg, ncol = 1,
                     align = "v",
                     axis = "lr",
                     rel_heights = c(14, 14, 9, 3))

ggsave("output_plots/fig4.svg", plot4, width = 180, height = 210, unit = "mm")
ggsave("output_plots/fig5_1.svg", plot4, width = 180, height = 210, unit = "mm")


#p5 <- grid.arrange(p5b10, p5b2, g_legend(p5b10), layout_matrix = rbind(c(1, 1), c(1, 1), c(1, 1),
#                                                                           c(2, 3), c(2, 3), c(2, 3), c(2, 3), c(2, 3), c(2, 3), c(2, 3), c(2, 3),
#                                                                           c(4,4)))




ggsave("output_plots/fig1c.svg", p1a, width = 60, height = 75, unit = "mm")
ggsave("output_plots/fig1c_legends_1.svg", p1a_legends_1, width = 120, height = 75, unit = "mm")
ggsave("output_plots/fig1c_legends.svg", p1a_legends, width = 120, height = 75, unit = "mm")
ggsave("output_plots/figS19.svg", p1b + ggtitle(""), width = 90, height = 75, unit = "mm")
#ggsave("output_plots/fig2.svg", p2a, width = 90, height = 140, unit = "mm")
ggsave("output_plots/figS1.svg", p2b, width = 90, height = 140, unit = "mm")
#ggsave("output_plots/fig3.svg", p2 + theme(legend.position = "none"), width = 180, height = 120, unit = "mm")
ggsave("output_plots/fig3_supp.svg", p2_other + theme(legend.position = "none"), width = 180, height = 120, unit = "mm")

for( e in EN2){
  ggsave(paste0("output_plots/pan_core_class_", e,".svg"), 
         grid.arrange(list_plot_core_class[[e]], list_plot_pan_class[[e]], nrow = 2) , 
         width = 180, height = 180, unit = "mm")
}

ggsave("output_plots/alluvial.svg", p_alluvial , width = 180, height = 180, unit = "mm")
ggsave("output_plots/fig5.svg", p42, width = 180, height = 80, unit = "mm")


ggsave("output_plots/fig6.svg", p5, width = 180, height = 220, unit = "mm")

################################################################################################   
  # Fig 3

# abundance_env_class <- abundance_medians_env_class(abundance_class, "human gut", tools_levels, unique(unigenes$new_level))


## 

# abundance_env_class <- abundance_medians_env_class(abundance_class, "human gut", tools_levels, top20)
# unigenes_class <- get_unigenes_class(unigenes, tools_levels,  top20) %>% 
#   ungroup() %>% 
#   group_by(tool) %>% 
#   mutate(proportion = n / sum(n)) %>% 
#   arrange(desc(proportion)) %>% 
#   mutate(cum_p = cumsum(proportion)) %>% 
#   ungroup() 
# 
# p3a1 <- plot_unigenes_env_tool_class(unigenes_class,  general_size, pal_10_q, tools_levels, tools_levels, tools_levels, 0, 15000, 5, T) + 
#   ggtitle("") + theme(title = element_text(size = general_size + 2, face = "bold"))
# 
# p3a2 <- plot_unigenes_env_tool_class(unigenes_class,  general_size, pal_10_q, tools_levels, tools_levels, tools_levels, 15000, max(unigenes_class$n), 2, T) + 
#   ggtitle("A") + ylab("") + theme(title = element_text(size = general_size + 2, face = "bold"), axis.text.x  = element_blank())
# 
# p3a <- grid.arrange(p3a2 + theme(plot.margin = margin(0, 0, -10, 0)), 
#                     p3a1 + theme(plot.margin = margin(-10, 0, 0, 0)), 
#                     layout_matrix = rbind(c(1, 1), c(2, 2), c(2, 2), c(2, 2)))
#   
# p3b0 <- plot_abundance_env_tool_class(abundance_env_class, tools_levels, "human gut", top20, general_size, pal_10_q, tools_levels, tools_levels) +
#   ggtitle("B") + theme(title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank())
# 
# 
# p3b <- p3b0 + theme(legend.position = "none")
# p3c <- plot_core_env_tool_class(sumcore, tools_levels, "human gut", top20, general_size, pal_10_q, tools_labels, tools_levels) +
#   ggtitle("C") + theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank()) +
#   ylab("Number of ARGs")
# 
# p3 <- grid.arrange(p3a, p3b, p3c, g_legend(p3b0), 
#                    layout_matrix = rbind(c(1, 1), c(1, 1), c(1, 1), c(1, 1), c(1, 1), 
#                                          c(2, 2), c(2, 2), c(2, 2), c(2, 2), c(2, 2), 
#                                          c(3, 3), c(3, 3), c(3, 3), c(3, 3), c(3, 3), c(4,4)))

# ggsave("~/Documents/plots_project6/fig4.svg", p3, width = 180, height = 220, unit = "mm")


# for(microbiome in 1:length(microbiomes)){
#   abundance_env_class <- abundance_medians_env_class(abundance_class, microbiomes[microbiome], tools_levels_2, top20)
#   
#   p3b0 <- plot_abundance_env_tool_class(abundance_env_class, tools_levels_2, microbiomes[microbiome], top20, general_size, pal_10_q, tools_levels_2, tools_levels_2) +
#     ggtitle(paste0(microbiomes[microbiome], "\nA")) + theme(title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank())
#   
#   p3b <- p3b0 + theme(legend.position = "none")
#   p3c <- plot_core_env_tool_class(sumcore, tools_levels_2, microbiomes[microbiome], top20, general_size, pal_10_q, tools_levels_2, tools_levels_2) +
#     ggtitle("B") + theme(legend.position = "none", title = element_text(size = general_size + 2, face = "bold"), strip.text = element_blank()) +
#     ylab("Number of ARGs")
#   
#   p3 <- grid.arrange(p3b, p3c, g_legend(p3b0), layout_matrix = rbind(c(1, 1), c(1, 1), c(1, 1), c(1, 1), c(1, 1), 
#                                                                      c(2, 2), c(2, 2), c(2, 2), c(2, 2), c(2, 2),
#                                                                      c(3,3)))
#   ggsave(paste0(paste0("~/Documents/plots_project6/figS", 6 + microbiome),".svg"), p3, width = 180, height = 150, unit = "mm")
# }




  ################################################################################################
  # Fig 4
general_size <- 10

d <- abundance_class %>% 
  filter(habitat %in% "human gut") %>% 
  mutate(total = normed10m) %>%
  mutate(gene = ifelse(gene %in% top20, gene, "Other")) %>% 
  mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
  mutate(gene = factor(gene, levels = c(top20, "Other"))) %>%
  ungroup() %>% 
  group_by(sample, tool, gene) %>%
  summarise(total = sum(total), texture = texture[1])

max_abun_class <- d %>% group_by(tool, gene) %>% summarise(q75 = quantile(total, .75)) %>% ungroup() %>% summarise(m = max(q75)) %>% pull()


abu_class_legend0 <- d %>%
  ggplot(aes(x = gene, y = total, fill = tool, pattern = texture)) + 
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot", 
               position = position_dodge2(preserve = "single"), 
               color = "black", linewidth = 0.2, 
               outlier.shape = NA, outlier.size = 0) +
  scale_fill_manual(values = pal_10_q, 
                    labels = lab_fn(tools_levels)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 3000, length.out = 4),
                     labels = scales::label_number()) + 
  theme_minimal() +
  facet_grid(. ~ gene, scales = "free_x", space = "free_x") + 
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x})  + 
  ylab("Relative abundance") +
  xlab("") +
  labs(fill = "") +
  ggtitle("B") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size),
    panel.background = element_rect(colour = "black", fill = NA),
    strip.text = element_blank())

pattern_density <- 0.001 
pattern_spacing <- 0.025
pattern_fill <- "white"
pattern_size <- 0.12

abu_class_human <- d %>%
  ggplot(aes(x = gene, y = total, fill = tool, pattern = texture)) + 
  stat_summary(fun.data = calc_boxplot_stat, geom="boxplot_pattern", 
               position = position_dodge2(preserve = "single"), 
               color = "black", linewidth = 0.2, 
               pattern_color = "white",
               pattern_density = pattern_density, 
               pattern_spacing = pattern_spacing, 
               pattern_fill = pattern_fill,
               pattern_size = pattern_size,
               pattern_key_scale_factor = 1.2, outlier.shape = NA, outlier.size = 0) +
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  scale_fill_manual(values = pal_10_q, 
                    labels = lab_fn(tools_levels)) +
  scale_y_continuous(expand = c(0, 0), breaks = seq(0, 3000, length.out = 4),
                     labels = scales::label_number()) + 
  theme_minimal() +
  facet_grid(. ~ gene, scales = "free_x", space = "free_x") + 
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x})  + 
  ylab("Relative abundance") +
  xlab("") +
  labs(fill = "") +
  ggtitle("B") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = general_size ),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
    legend.margin = margin(0, 0, 0, 0, unit = "pt"),
    panel.spacing = unit(0, "pt"),
    title = element_text(size = general_size + 2, face = "bold"),
    axis.title = element_text(size = general_size + 1, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
    axis.text.y = element_text(size = general_size),
    panel.background = element_rect(colour = "black", fill = NA),
    strip.text = element_blank())


ggsave("output_plots/fig_human_class_abundance.svg", abu_class_human, width = 180, height = 75, unit = "mm")

########################################################################################




df0 <- unigenes %>% ungroup() %>%
  mutate(tool_2 = ifelse(tool %in% c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD",
                                     "ABRicate-ResFinder", "ABRicate-NCBI", "ResFinder"), "ABRicate and ResFinder", as.character(tool))) %>% 
  filter(!((tool %in% "AMRFinderPlus") & (query %in% lst$amrfinder.norm.prot$query[lst$amrfinder.norm.prot$Method %in% "HMM"]))) %>% 
  filter(!tool %in% "fARGene") %>% 
  mutate(bin = cut_width(round(id), width = 1, boundary = min(round(id), na.rm = T))) %>%
  group_by(tool_2, bin) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  mutate(tool_2 = factor(tool_2, levels = c("DeepARG", "RGI-DIAMOND", "AMRFinderPlus", "ABRicate and ResFinder"))) %>%
  rename(tool = tool_2)


df_min <- unigenes %>% ungroup() %>%
  mutate(tool_2 = ifelse(tool %in% c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD",
                                     "ABRicate-ResFinder", "ABRicate-NCBI", "ResFinder"), "ABRicate and ResFinder", as.character(tool))) %>% 
  filter(!((tool %in% "AMRFinderPlus") & (query %in% lst$amrfinder.norm.prot$query[lst$amrfinder.norm.prot$Method %in% "HMM"]))) %>% 
  filter(!tool %in% "fARGene") %>% 
  select(-tool) %>% 
  rename(tool = tool_2) %>% 
  group_by(tool) %>% summarise(mi = min(id, na.rm = T), ma = max(id, na.rm = T))

df0 <- df0 %>% left_join(df_min) %>% 
  mutate(tool = factor(tool, levels = c("DeepARG", "RGI-DIAMOND", "AMRFinderPlus", "ABRicate and ResFinder"))) 


df0 <- df0 %>% mutate(texture = ifelse(tool %in% tools_texture, "yes", "no")) %>%
  ungroup() %>%
  mutate(id2 =  as.numeric(sub(".*,(.*)\\]", "\\1", bin))) %>%
  mutate(id1 =  as.numeric(sub("^[\\(\\[]([0-9]+),.*\\]$", "\\1", bin))) %>% 
  group_by(tool) %>%
  mutate(xmin = ifelse(id1 < mi, mi, id1)) %>%
  mutate(xmax = ifelse(id2 > ma, ma, id2))


lab_fn <- function(x) {
  x <- gsub("-", "-\n", x)
  x <- gsub(" ", "\n", x)
  x <- gsub("/", "/\n", x)
  x
}

id_plot <- unigenes %>% 
  ungroup() %>%
  mutate(tool_2 = ifelse(tool %in% c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD",
                                     "ABRicate-ResFinder", "ABRicate-NCBI", "ResFinder"), "ABRicate/ResFinder", as.character(tool))) %>% 
  filter(!((tool_2 %in% "AMRFinderPlus") & (query %in% lst$amrfinder.norm.prot$query[lst$amrfinder.norm.prot$Method %in% "HMM"]))) %>% 
  filter(!tool_2 %in% "fARGene") %>% 
  mutate(tool_2 = factor(tool_2, levels = c("DeepARG", "RGI-DIAMOND", "AMRFinderPlus", "ABRicate/ResFinder"))) %>% 
  mutate(data_type = ifelse(tool_2 %in% c("DeepARG", "RGI-DIAMOND", "AMRFinderPlus"), "amino acid", "nucleotide")) %>% 
  ggplot(aes(x = id , color = tool_2, fill = tool_2, linetype = data_type)) +
  stat_ecdf(geom = "step", linewidth = 1) + 
  scale_linetype_manual(values = c("solid","dotted")) +
  theme_minimal() + 
  #scale_fill_manual(values = pal_10_q[c(1,5,7,9)], labels = c(tools_labels[c(1,5,7)],"ABRicate/ResFinder"), guide = guide_legend(nrow = 2)) +
  scale_color_manual(values = pal_10_q[c(1,5,7,9)], 
                     labels = lab_fn(c(tools_labels[c(1,5,7)],"ABRicate/ResFinder")), 
                     guide = guide_legend(nrow = 2)) +
  ggtitle("B") + 
  ylab("ECDF") + 
  xlab("% Identity") + 
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) +
  labs(color = "", fill = "", linetype = "") +
  scale_y_continuous(limits = c(0,1), expand = c(0,0)) + 
  theme(legend.position = "bottom",
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text(size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        #panel.grid.major.y = element_blank(),
        #panel.grid.minor = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        panel.background = element_rect(colour = "black", fill = NA))

id_plot


id_plot2 <- df0 %>%
  ggplot(aes(x = id2 , color = tool, y = frac, fill = tool, pattern = texture)) +
  geom_rect_pattern(aes(xmin = xmin, xmax = xmax, 
                ymin = 0, ymax = frac), linewidth = 0,
                pattern_color = "white",
                pattern_density = 0.1, 
                pattern_spacing = 0.2, 
                pattern_fill = "white",
                pattern_key_scale_factor = 0.05) + 
  scale_pattern_manual(values = c('no' = 'none', 'yes' = 'stripe')) +
  theme_minimal() + 
  scale_fill_manual(values = pal_10_q[c(1,5,7,9)], labels = tools_labels[c(1,5,7,9)], guide = guide_legend(nrow = 2)) +
  scale_color_manual(values = pal_10_q[c(1,5,7,9)], labels = tools_labels[c(1,5,7,9)], guide = guide_legend(nrow = 2)) +
  #scale_x_continuous(limits = c(20, 100), breaks = seq(from = 0, to = 100, by = 20)) + 
  ggtitle("B") + 
  ylab("") + 
  xlab("") + 
  labs(color = "", fill = "") +
  scale_y_continuous(limits = c(0,max(df0$frac)+0.001), expand = c(0,0)) + 
  facet_grid(. ~ tool , scales = "free") +
  theme(legend.position = "none",
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_blank(),
        strip.text = element_text(size = general_size, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        panel.background = element_rect(colour = "black", fill = NA))

id_plot2



ggsave("output_plots/ids.svg", id_plot,  width = 90, height = 70, unit = "mm")



p1ab <- plot_grid( p1a + ggtitle("") + ylab(""), id_plot + ggtitle(""), ncol = 2,
           align = "h",
           axis = "tb",
           rel_widths = c(1, 1))

ggsave("output_plots/fig1ab.svg", p1ab, width = 180, height = 70, unit = "mm")

ggsave("output_plots/fig1d.svg", fig2_d, width = 90, height = 140, unit = "mm")
ggsave("output_plots/fig1c.svg", p2a + theme(legend.position = "none"), width = 90, height = 140, unit = "mm")


fig2_d <- plot_grid( pan_all, core_all, ncol = 1,
                     align = "v",
                     axis = "lr")


fig2_bottom <- plot_grid(
  p2a + theme(legend.position = " none") + ggtitle("C"), pan_core_all + ggtitle("D"),
  ncol = 2,
  align = "h",
  axis = "tb",
  rel_widths = c(5,3)
)

fig2_top <- plot_grid(
  p1a + ggtitle("A")  , id_plot,
  ncol = 2, align = "h", axis = "tb",
  rel_widths = c(7,4)
)

p2abcd <- plot_grid(fig2_top, fig2_bottom,
                   nrow = 2,
                   align = "v",
                   axis = "lr",
                   rel_heights = c(2, 5))


ggsave("output_plots/fig2.svg", p2abcd, width = 180, height = 210, unit = "mm")


p2_1 <- plot_grid(fig2_top,
                  plot_grid(p_alluvial + ggtitle("C") + scale_x_discrete(labels = lab_fn), ncol = 1, align = "h", axis = "tb"),
                  ncol = 1,
                  align = "hv",
                  axis = "tblr",
                  rel_heights = c(5, 9))


ggsave("output_plots/fig2_1.svg", p2_1, width = 180, height = 210, unit = "mm")


g_legend(abu_class_human)

fig_tot_abund_human_abund <- plot_grid(
  p2a2 + theme(legend.position = " none") + ggtitle("A"), 
  abu_class_human + theme(legend.position = " none") + ggtitle("B"),
  g_legend(abu_class_legend0),
  nrow = 3,
  align = "v",
  axis = "lr",
  rel_heights = c(14, 8, 2))


ggsave("output_plots/fig3_1.svg", fig_tot_abund_human_abund, width = 180, height = 210, unit = "mm")
ggsave("output_plots/fig3_1_alt_legend.svg", g_legend(abu_class_human), width = 180, height = 30, unit = "mm")

ggsave("output_plots/fig_diversity.svg", p2diversity, width = 180, height = 140, unit = "mm")

p4_1 <- plot_grid(pan_core_all + ggtitle("A"),
                  circle_triangle_legends,
                  p_alluvial_pan_core + ggtitle("B") + xlab("Pan- and core-resistomes"),
                  ncol = 1,
                  align = "hv",
                  axis = "tblr",
                  rel_heights = c(7,2,15))

ggsave("output_plots/fig4_1.svg", p4_1, width = 180, height = 210, unit = "mm")

################################################################################################
# Fig S1
# 
# rgi.ven = ggVennDiagram(list("RGI-DIAMOND" = lst$rgi.diamond$query,
#                              "RGI-BLAST" = lst$rgi.blast$query,
#                              "RGI-DIAMOND (aa)" = lst$rgi.diamond.prot$query),
#                         color = 1, lwd = 0.7, label_size = 2, set_size = 2) + 
#   scale_fill_gradient(low = "#F6E8C3", high = pal_10_q[5]) +
#   theme(legend.position = "none") +
#   ggtitle("A") +
#   theme(
#     legend.text = element_text(size = general_size),
#     legend.title = element_text(size = general_size + 1),
#     text = element_text(size = general_size),
#     axis.text  = element_blank(),
#     title = element_text(size = general_size + 2, , face = "bold"))
# 
# rgi.ven
# 
# 
# deeparg.ven = ggVennDiagram(list("DeepARG" = lst$deeparg.norm$query,
#                                  "DeepARG (aa)" = lst$deeparg.norm.prot$query),
#                             label_percent_digit = 1,
#                             color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
#   scale_fill_gradient(low = "#F6E8C3", high = pal_10_q[5]) +
#   theme(legend.position = "none") +
#   scale_x_continuous(labels = ) +
#   ggtitle("B") +
#   theme(
#     legend.text = element_text(size = general_size),
#     legend.title = element_text(size = general_size + 1),
#     text = element_text(size = general_size),
#     axis.text  = element_blank(),
#     title = element_text(size = general_size + 2, , face = "bold"))
# 
# deeparg.ven
# 
# 
# fargene.ven = ggVennDiagram(list("fARGene" = lst$fargene$query,
#                                  "fARGene (aa)" = lst$fargene.prot$query),
#                             label_percent_digit = 2,
#                             color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
#   scale_fill_gradient(low = "#F6E8C3", high = pal_10_q[5]) +  
#   theme(legend.position = "none") +
#   ggtitle("C") +
#   theme(
#     legend.text = element_text(size = general_size),
#     legend.title = element_text(size = general_size + 1),
#     text = element_text(size = general_size),
#     axis.text  = element_blank(),
#     title = element_text(size = general_size + 2, , face = "bold"))
# 
# fargene.ven
# 
# amrfinder.ven = ggVennDiagram(list("AMRFinderPlus" = lst$amrfinder.norm$query,
#                                    "AMRFinderPlus (aa)" = lst$amrfinder.norm.prot$query),
#                               label_percent_digit = 1,
#                               color = 1, lwd = 0.1, label_size = 2, set_size = 2) + 
#   scale_fill_gradient(low = "#F6E8C3", high = pal_10_q[5]) +  
#   theme(legend.position = "none") +
#   ggtitle("D") +
#   theme(
#     legend.text = element_text(size = general_size),
#     legend.title = element_text(size = general_size + 1),
#     text = element_text(size = general_size),
#     axis.text  = element_blank(),
#     title = element_text(size = general_size + 2, , face = "bold"))
# 
# amrfinder.ven

theme_overlap <- theme_minimal() + 
  theme(legend.position = "none",
        panel.border = element_blank(),
        axis.text.x = element_text(size = general_size),
        axis.text.y = element_text( size = general_size),
        strip.text = element_text(size = general_size, face = "bold"),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        panel.background = element_rect(colour = "black", fill = NA))


rgi_query <- data.frame(gene = unique(c(lst$rgi.blast$query, lst$rgi.diamond$query, lst$rgi.diamond.prot$query)))
rgi_query <- rgi_query %>% mutate(dataset = ifelse(gene %in% intersect(intersect(lst$rgi.blast$query, lst$rgi.diamond$query), lst$rgi.diamond.prot$query), "All",
                                            ifelse(gene %in% intersect(lst$rgi.blast$query, lst$rgi.diamond$query), "B-nt and D-nt",
                                            ifelse(gene %in% intersect(lst$rgi.diamond.prot$query, lst$rgi.diamond$query), "D-nt and D-aa",
                                            ifelse(gene %in% intersect(lst$rgi.diamond.prot$query, lst$rgi.blast$query), "B-nt and D-aa",
                                            ifelse(gene %in% lst$rgi.diamond.prot$query, "D-aa",
                                            ifelse(gene %in% lst$rgi.diamond$query, "D-nt",
                                            ifelse(gene %in% lst$rgi.blast$query, "B-nt", NA)))))))) %>%
  mutate(dataset = factor(dataset, levels = c("All", "B-nt and D-nt", "D-nt and D-aa", "B-nt and D-aa", "D-aa", "D-nt", "B-nt")))
                                            
rgi_over <- rgi_query %>% 
  group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
  ggplot(aes(x = dataset, y = n, fill = dataset)) + 
  geom_col() + 
  scale_fill_manual(values = pal_10_q[3:10])  +
  geom_text(aes(label = scales::percent(p, accuracy = 0.1)), vjust = -0.3, size = general_size / .pt) +
  ylab("ARGs") + 
  xlab("") +
  scale_x_discrete(labels = function(x) {
    x <- gsub(" ", "\n", x)
    x}) +
  theme_overlap


deeparg_query <- data.frame(gene = unique(c(lst$deeparg.norm$query, lst$deeparg.norm.prot$query)))
deeparg_query <- deeparg_query %>% mutate(dataset = ifelse(gene %in% intersect(lst$deeparg.norm$query, lst$deeparg.norm.prot$query), "Both DeepARG",
                                                   ifelse(gene %in% lst$deeparg.norm$query, "nt", "aa"))) %>% 
  mutate(dataset = factor(dataset, levels = c("Both DeepARG", "nt", "aa")))


deeparg_over <- deeparg_query %>% 
  group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
  ggplot(aes(x = dataset, y = n, fill = dataset)) + 
  geom_col() + 
  scale_fill_manual(values = pal_10_q[3:10])  +
  geom_text(aes(label = scales::percent(p, accuracy = 0.1)), vjust = -0.3, size = general_size / .pt) +
  ylab("ARGs") + 
  xlab("") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) +
  theme_overlap


fargene_query <- data.frame(gene = unique(c(lst$fargene$query, lst$fargene.prot$query)))
fargene_query <- fargene_query %>% mutate(dataset = ifelse(gene %in% intersect(lst$fargene$query, lst$fargene.prot$query), "Both fARGene",
                                                           ifelse(gene %in% lst$fargene$query, "nt", "aa"))) %>% 
  mutate(dataset = factor(dataset, levels = c("Both fARGene", "nt", "aa")))


fargene_over <- fargene_query %>% 
  group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
  ggplot(aes(x = dataset, y = n, fill = dataset)) + 
  geom_col() + 
  scale_fill_manual(values = pal_10_q[3:10])  +
  geom_text(aes(label = scales::percent(p, accuracy = 0.1)), vjust = -0.3, size = general_size / .pt) +
  ylab("ARGs") + 
  xlab("") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) +
  theme_overlap



amrfinder_query <- data.frame(gene = unique(c(lst$amrfinder.norm$query, lst$amrfinder.norm.prot$query)))
amrfinder_query <- amrfinder_query %>% mutate(dataset = ifelse(gene %in% intersect(lst$amrfinder.norm$query, lst$amrfinder.norm.prot$query), "Both AMRFinder-Plus",
                                                           ifelse(gene %in% lst$amrfinder.norm$query, "nt", "aa"))) %>% 
  mutate(dataset = factor(dataset, levels = c("Both AMRFinder-Plus", "nt", "aa")))


amrfinder_over <-amrfinder_query %>% 
  group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
  ggplot(aes(x = dataset, y = n, fill = dataset)) + 
  geom_col() + 
  geom_text(aes(label = scales::percent(p, accuracy = 0.1)), vjust = -0.3, size = general_size / .pt) +
  scale_fill_manual(values = pal_10_q[3:10])  +
  ylab("ARGs") + 
  xlab("") +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) +
  theme_overlap


plot_overlaps <- plot_grid(
  rgi_over + ggtitle("A"), 
  deeparg_over + ggtitle("B"), 
  fargene_over + ggtitle("C"), 
  amrfinder_over + ggtitle("D"),
  ncol = 2,
  align = "hv"
)




plot_overlaps <- function(x, y, nxy, nx, ny){
  db_query <- data.frame(gene = unique(c(x, y)))
  db_query <- db_query %>% mutate(dataset = ifelse(gene %in% intersect(x, y), nxy,
                                                          ifelse(gene %in% x, nx, ny))) %>% 
    mutate(dataset = factor(dataset, levels = c(nxy, nx, ny)))
  
  
  db_query_plot <- db_query %>% 
    group_by(dataset) %>% summarise(n = n()) %>% mutate(p = n/sum(n)) %>%
    ggplot(aes(x = dataset, y = n, fill = dataset)) + 
    geom_col() + 
    scale_fill_manual(values = pal_10_q[3:10])  +
    geom_text(aes(label = scales::percent(p, accuracy = 0.1)), vjust = -0.3, size = general_size / .pt) +
    ylab("ARGs") + 
    xlab("") +
    scale_x_discrete(labels = function(x) {
      x <- gsub("-", "-\n", x)
      x <- gsub(" ", "\n", x)
      x}) +
    theme_overlap
  
  return(db_query_plot)
}

card_over <- plot_overlaps(lst$rgi.diamond.prot$query, lst$abricate.card.norm$query, "Both", "RGI", "ABRicate")
ncbi_over <- plot_overlaps(lst$amrfinder.norm.prot$query, lst$abricate.ncbi.norm$query, "Both", "AMRFinderPlus", "ABRicate")
resfinder_over <- plot_overlaps(lst$resfinder.norm$query, lst$abricate.resfinder.norm$query, "Both", "ResFinder", "ABRicate")

plot_db <- plot_grid(
  card_over  + ggtitle("A"), 
  ncbi_over + ggtitle("B"), 
  resfinder_over + ggtitle("C"),
  ncol = 2,
  align = "hv"
)


abricate_meg1 <- plot_overlaps(lst$abricate.resfinder.norm$query, lst$abricate.megares.norm$query, "Both", "ResFinder", "MEGARes")
abricate_meg2 <- plot_overlaps(lst$abricate.argannot.norm$query, lst$abricate.megares.norm$query, "Both", "ARGANNOT", "MEGARes")
abricate_meg3 <- plot_overlaps(lst$abricate.card.norm$query, lst$abricate.megares.norm$query, "Both", "CARD", "MEGARes")
abricate_meg4 <- plot_overlaps(lst$abricate.ncbi.norm$query, lst$abricate.megares.norm$query, "Both", "NCBI", "MEGARes")

plot_megares <- plot_grid(
  abricate_meg2 + ggtitle("A"), 
  abricate_meg3 + ggtitle("B"), 
  abricate_meg4 + ggtitle("C"),
  abricate_meg1  + ggtitle("D"), 
  ncol = 2,
  align = "hv"
)


ggsave("output_plots/overlaps.svg", plot_overlaps, width = 180, height = 160, unit = "mm")
ggsave("output_plots/overlaps_db.svg", plot_db, width = 180, height = 160, unit = "mm")
ggsave("output_plots/overlaps_megares.svg", plot_megares, width = 160, height = 210, unit = "mm")

################################################################################################
# Fig S2




# ggsave("~/Documents/plots_project6/figS3.svg", ps2, width = 180, height = 150, unit = "mm")
# ggsave("~/Documents/plots_project6/figS4.svg", ps2, width = 180, height = 150, unit = "mm")
# ggsave("~/Documents/plots_project6/figS5.svg", ps2, width = 180, height = 150, unit = "mm")


################################################################################################
# Fig S3


# class overlap for core resistome 


# data.frame(JI_core_class %>% ungroup() %>% group_by(new_level) %>% summarise(m = sum(jaccard), n = n()) %>% mutate(m = ifelse(is.na(m), 0, m / 90)) %>% arrange(desc(m)))
# length(unique(paste(JI_core_class$tool_ref, JI_core_class$tool_comp)))
# 
# core_jaccard_human <- JI_core_class %>% ungroup() %>% 
#   filter(as.numeric(tool_comp) < as.numeric(tool_ref)) %>% #complete(tool, new_level, fill = list(unigenes = 0, total = 0)) %>%
#   ungroup() %>% 
#   filter(as.numeric(tool_ref) != as.numeric(tool_comp)) %>% 
#   ggplot(aes(x = tool_ref, y = tool_comp, fill = jaccard)) +
#   geom_tile(na.rm = T) +
#   scale_x_discrete(labels = tools_levels) +
#   scale_y_discrete(labels = tools_levels) +
#   facet_wrap(new_level ~ .) +
#   scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "YlOrBr"), na.value = "white") +
#   theme_minimal() +
#   labs(fill = "") +
#   xlab("") +
#   ylab("") +
#   ggtitle("") + 
#   theme(
#     legend.position = "bottom",
#     legend.text = element_text(size = general_size ),
#     panel.border = element_rect(fill = "transparent", color = "black", linewidth = 1),
#     #panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank(),
#     plot.margin = margin(0, 0, 0, 0, unit = "pt"),
#     legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
#     legend.margin = margin(0, 0, 0, 0, unit = "pt"),
#     panel.spacing = unit(0, "pt"),
#     strip.text = element_text(size = general_size + 1, face = "bold"),
#     axis.title = element_text(size = general_size + 1, face = "bold"),
#     axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
#     axis.text.y = element_text(size = general_size))
# 
# core_jaccard_human
# 
# ggsave("~/Documents/plots_project6/figS6.svg", core_jaccard_human, width = 300, height = 300, unit = "mm") 


################################################################################################
# Fig S4

##  Tables 

# sum(table(abundance$habitat[!duplicated(abundance$sample)]))
# table(abundance$habitat[!duplicated(abundance$sample)])

# write.csv(abundance %>% ungroup() %>% select(sample, habitat) %>% filter(!habitat %in% c("amplicon","built-environment","isolate")) %>%
#  rename(accession = sample) %>% group_by(accession) %>% slice_head(n = 1), "~/Documents/plots_project3/accession_metagenomes.csv", row.names = FALSE)

# write.csv(abundance %>% ungroup() %>% select(sample, habitat) %>% filter(!habitat %in% c("amplicon","built-environment","isolate")) %>%
#             rename(accession = sample) %>% group_by(accession) %>% slice_head(n = 1) %>% 
#            ungroup() %>% group_by(habitat) %>% summarise(samples = n()), "~/Documents/plots_project3/summary_metagenomes.csv", row.names = FALSE)




###






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



sumcore <- sum_core_adjust(core, 450, 0.5)
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








core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  group_by(tool, new_level) %>%
  summarise(n_tool = n_distinct(X)) %>%
  mutate(p = n_tool/sum(n_tool)) %>%
  arrange(desc(n_tool))%>%
  filter(new_level %in% "tet RPG")

28 / 12
core %>% filter(habitat %in% "pig gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  group_by(tool, new_level) %>%
  summarise(n_tool = n_distinct(X)) %>%
  mutate(p = n_tool/sum(n_tool)) %>%
  arrange(desc(n_tool))

core %>% filter(habitat %in% "pig gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  group_by(tool) %>%
  summarise(n_tool = n_distinct(X)) 
  arrange(desc(n_tool))

core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  group_by(X, new_level) %>%
  summarise(n_tool = n_distinct(tool)) %>% 
  mutate(n_tool50 = ifelse(n_tool>=5,"yes","no")) %>% 
  ungroup() %>% 
  group_by( new_level) %>% 
  mutate(N = n_distinct(X)) %>%
  ungroup() %>%
  group_by(new_level, n_tool) %>% 
  summarise(N = N[1], np = n()/N) %>%
  arrange(desc(N)) %>%
  print(n = 50) %>% filter(n_tool == 1)



core_1 <- core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  filter(!tool %in% c("DeepARG", "RGI-DIAMOND")) %>%
  group_by(X) %>%
  summarise(n_tool = n_distinct(tool)) 


#167 unique human gut genes 

core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  filter(tool %in% c("DeepARG", "RGI-DIAMOND"), !X %in% core_1$X) %>%
  ungroup() %>% 
  summarise(n_distinct(X))
  


only_rgi_deep_arg <- core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  filter(tool %in% c("DeepARG", "RGI-DIAMOND"), !X %in% core_1$X) %>%
  ungroup() 

core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  filter(tool %in% c("DeepARG", "RGI-DIAMOND"), !X %in% core_1$X) %>%
  ungroup() %>% 
  group_by(new_level) %>%
  summarise(n = n_distinct(X)) %>%
  arrange(desc(n))





fig_numbe_genes_core <- core %>% filter(habitat %in% "human gut", cut %in% 0.5, cnt > 450) %>% 
  ungroup() %>% 
  group_by(X, new_level) %>%
  summarise(n_tool = n_distinct(tool)) %>% 
  mutate(tools = ifelse(X %in% only_rgi_deep_arg$X, "Exlusive RGI/DeepARG", "Other")) %>%
  ggplot(aes (x = new_level, y = n_tool, color = tools)) +
  geom_jitter(height = 0.1) +
  ylab("Number of tools") +
  xlab("") + 
  labs(color = "") +
  scale_y_continuous(expand = c(.01, .01), 
                     limits = c(0, 10), 
                     breaks = seq(from = 0, to = 10, by = 2)) +
  scale_color_manual(values = pal_10_q[3:4]) +
  scale_x_discrete(labels = function(x) {
    x <- gsub("-", "-\n", x)
    x <- gsub(" ", "\n", x)
    x}) +
  facet_grid(. ~ new_level, scales = "free_x", space = "free_x") +
  theme_minimal() + 
  theme(legend.position = "bottom",
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = general_size),
        axis.text.y = element_text( size = general_size),
        strip.text = element_blank(),
        plot.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
        legend.margin = margin(0, 0, 0, 0, unit = "pt"),
        panel.spacing = unit(0, "pt"),
        legend.text = element_text(size = general_size),
        title = element_text(size = general_size + 2, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(size = general_size + 1, face = "bold"),
        panel.background = element_rect(colour = "black", fill = NA))




#######


alluvial_pan_core_env <- function(sumcore, pan, h){
  pan_core_env <- bind_rows(
    sumcore %>% filter(habitat %in% h) %>% ungroup() %>% mutate(dt = "core"),
    pan %>% 
      filter(habitat %in% h, aggregation %in% "new_level_centroid") %>% 
      ungroup() %>% 
      group_by(tool, habitat, epoch, gene_class) %>% 
      summarise(s = sum(unigenes)) %>%
      ungroup() %>% 
      group_by(tool, habitat, gene_class) %>% 
      summarise(unigenes = ceiling(mean(s))) %>% 
      rename(new_level  = gene_class) %>% mutate(dt = "pan")) %>%
    mutate(dt = factor(dt, levels = c("pan", "core")))
  
  pan_core_env <- pan_core_env %>% ungroup() %>%  mutate(new_level = as.character(new_level)) %>% 
    group_by(tool, new_level, dt) %>% summarise( n = sum(unigenes)) %>% 
    arrange(tool, desc(n))
  
  pan_core_env <- pan_core_env %>% 
    ungroup() %>%
    group_by(tool, dt) %>% 
    mutate(proportion = n / sum(n)) %>% 
    arrange(tool, desc(proportion)) %>% 
    mutate(cum_p = cumsum(proportion)) %>% 
    ungroup() %>%
    mutate( new_level = factor(as.character(new_level), 
                               levels = c("Other", rev(levels_unigenes))))
  
  
  pan_core_env <- pan_core_env %>%
    group_by(tool, dt) %>%
    arrange(proportion) %>% 
    { 
      if (any(.$cum_p == 0.95)) {
        filter(., cum_p <= 0.95)
      } else {
        lower_part <- filter(., cum_p < 0.99999)
        upper_part <- filter(., cum_p > 0.99999) %>% slice_tail(n = 1)
        bind_rows(lower_part, upper_part)
      }
    } %>% 
    ungroup() %>% 
    arrange(tool, cum_p)
  
  
  pan_core_env <- pan_core_env %>% 
    filter(proportion > 0.01) %>% 
    ungroup() %>% 
    complete(new_level, tool, dt,
             fill = list(n = 0,  proportion = 0.000, cum_p = 0.000))  %>%  
    mutate(gene_name = gene_classes_list[as.character(new_level)]) %>%  
    mutate(gene_name = ifelse(new_level %in% "Other", "Other", gene_name)) %>% 
    mutate(gene_name = factor(gene_name, 
                              levels = c("Other",gene_classes_list[levels(pan_core_env$new_level)])))
  
  
  
  p_alluvial_pan_core_env <- pan_core_env %>% 
    ggplot(
      aes(x = dt,
          stratum = gene_name,
          alluvium = gene_name,
          y = proportion,
          fill = gene_name,
          label = gene_name)) +
    geom_flow(alpha = 0.5) +
    xlab("") + 
    geom_stratum(color = "black") +
    geom_text(data = pan_core_env[pan_core_env$proportion > 0.01 , ],  
              stat = "stratum",
              size = 2,color = "black",hjust = 0.5,
              position = position_jitter(width = 0, height = 0)) +
    scale_y_continuous(expand = c(.01, .01), 
                       name = "Proportion", 
                       limits = c(0, 1), 
                       breaks = seq(from = 0, to = 1, by = 0.1)) +
    facet_wrap(tool ~ ., nrow = 2) +
    scale_fill_manual(values = c(rep(pal_10_complete,10))) +
    xlab("Pan- vs core-resistome") + 
    #scale_x_discrete( labels =  labels_plot) +
    theme_minimal() +
    theme(
      legend.position = "none",
      legend.text = element_text(size = general_size),
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      plot.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.box.margin = margin(0, 0, 0, 0, unit = "pt"),
      legend.margin = margin(0, 0, 0, 0, unit = "pt"),
      panel.spacing = unit(0, "pt"),
      title = element_text(size = general_size + 2, face = "bold"),
      axis.title = element_text(size = general_size + 1, face = "bold"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(size = general_size),
      strip.text = element_text(size = general_size, face = "bold"),
      panel.background = element_rect(colour = "black", fill = NA)) +
    ggtitle(h)
  
  return(p_alluvial_pan_core_env)
}

alluvial_all <- list()

for(e in EN[-c(14:16)]){
  alluvial_all[[e]] <- alluvial_pan_core_env(sumcore, pan, e)
}


for(e in EN[-c(14:16)]){
  ggsave(paste0(paste0("output_plots/s_",e),".svg"), alluvial_all[[e]], width = 180, height = 210, unit = "mm")
}


