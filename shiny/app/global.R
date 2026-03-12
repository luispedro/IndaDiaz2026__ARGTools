library(shiny)
library(bslib)
library(tidyverse)
library(shinyWidgets)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(RColorBrewer)
library(ggpattern)
library(grid)
library(reactable)
library(Cairo)
library(ggalluvial)
library(cowplot)
library(scales)
library(tidyr)
library(magrittr)
library(shinycssloaders)
library(ggrastr)
library(ragg)

options(dplyr.summarise.inform = FALSE)
options(shiny.useragg = TRUE)

# Defining the root and app directory 
APP_DIR  <- normalizePath(getwd(), mustWork = TRUE)                
ROOT_DIR <- normalizePath(file.path(APP_DIR, "..", ".."), mustWork = TRUE)

# Building paths from ROOT_DIR to ensure they are correct regardless of where the app is launched from
DATA_DIR      <- file.path(ROOT_DIR, "code_R_analysis", "output_abundance_diversity_resistome")
METADATA_PATH <- file.path(ROOT_DIR, "data", "metadata_GMGC10.sample.meta.tsv")

# Sourcing the scripts using the constructed paths
source(file.path(ROOT_DIR, "code_R_analysis", "helper.R"), local = TRUE)

if (!exists("plot_total_abundance_diversity_new_version_shiny", mode = "function") &&
    exists("plot_total_abundance_diversity_new_version", mode = "function")) {
  plot_total_abundance_diversity_new_version_shiny <- plot_total_abundance_diversity_new_version
}


# Constants & Plotting config
general_size <- 10
pal_7 <- brewer.pal(8, "BrBG")[c(4, 5, 3, 6, 2, 7, 1)]
pal_10_q <- pal_7[c(1, 2, 3, 4, 5, 5, 6, 6, 7, 7)]
pal_10_complete <- brewer.pal(10, "BrBG")
pal_10_complete <- pal_10_complete[c(-1,-10)]

# pattern for plots

pattern_density <- 0.001 
pattern_spacing <- 0.025
pattern_fill <- "white"
pattern_size <- 0.12

# All HABITATS
EN <- c("human gut", "human oral",  "human skin", "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine", 
        "freshwater", "soil" , "amplicon", "isolate",  "built-environment" )

# SOURCE FOR EACH HABITAT

SO <- c(rep("humans", 5), rep("mammals", 4),  
        "wastewater", "marine", "freshwater", 
        "soil", rep("other", 2), "built-environment")

names(SO) <- EN

not_env <- c("amplicon", "isolate", "built-environment" )
EN2  <- EN[!EN %in% not_env]
h2 <- c("humans","mammals","wastewater", "freshwater", "soil", "marine")

# Tool Definitions
tools_levels <- c("DeepARG", "fARGene", "ABRicate-ARGANNOT", "ABRicate-MEGARes", 
                  "RGI-DIAMOND", "ABRicate-CARD", "AMRFinderPlus", 
                  "ABRicate-NCBI", "ResFinder", "ABRicate-ResFinder")

tools_labels <- c("DeepARG", "fARGene", "ABRicate-ARGANNOT", "ABRicate-MEGARes", 
                  "RGI-CARD", "ABRicate-CARD", "AMRFinderPlus-NCBI", 
                  "ABRicate-NCBI", "ResFinder", "ABRicate-ResFinder")

tools_texture <- c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD", 
                   "ABRicate-NCBI", "ABRicate-ResFinder")


tool_choices <- c("DeepARG" = "DeepARG", "fARGene" = "fARGene", 
                  "ABRicate-ARGANNOT" = "ABRicate-ARGANNOT" , 
                  "ABRicate-MEGARes" = "ABRicate-MEGARes", 
                  "RGI-CARD" = "RGI-DIAMOND", "ABRicate-CARD" = "ABRicate-CARD", 
                  "AMRFinderPlus-NCBI" = "AMRFinderPlus", 
                  "ABRicate-NCBI" = "ABRicate-NCBI", "ResFinder" = "ResFinder", 
                  "ABRicate-ResFinder" = "ABRicate-ResFinder")

top20 <- c("van", "efflux pump", "cell wall charge", "rpoB", "tet RPG", "class A beta-lactamase", 
           "aph", "MFS efflux pump", "erm", "target-modifying enzyme", "tet enzyme")


# Load all pre-calculated data

data_args      <- qs::qread(file.path(ROOT_DIR, "data", "data_args.qs"))
data_abundance <- qs::qread(file.path(ROOT_DIR, "data", "data_abundance.qs"))
data_pan_core  <- qs::qread(file.path(ROOT_DIR, "data", "data_pan_core.qs"))
data_overlap   <- qs::qread(file.path(ROOT_DIR, "data", "data_overlap.qs"))


# Ensure same variable names as before where possible

# Shared variable
gene_classes    <- data_args$levels_unigenes

# ARGs tab — now pre-filtered per threshold
data_list <- list(
  unigenes_prepped = data_args$unigenes_prepped,   # named list: default/60/70/80
  levels_unigenes  = data_args$levels_unigenes
)

# Abundance tab
abundance_prepped       <- data_abundance$abundance_prepped
abundance_class_prepped <- data_abundance$abundance_class_prepped

# Pan & Core tab
pan_prepped      <- data_pan_core$pan_prepped
core_prepped     <- data_pan_core$core_prepped
sumpan2_prepped  <- data_pan_core$sumpan2_prepped
core_sum_prepped <- data_pan_core$core_sum_prepped  # all 196 combos
pan_core_joined_prepped <- data_pan_core$pan_core_joined_prepped #The precomputed join for sumpan2_prepped and core_sum_prepped

# Overlap tab
data_list$recall_fnr         <- data_overlap$recall_fnr
data_list$recall_fnr60       <- data_overlap$recall_fnr60
data_list$recall_fnr70       <- data_overlap$recall_fnr70
data_list$recall_fnr80       <- data_overlap$recall_fnr80
data_list$cstc_summary_prepped <- data_overlap$cstc_summary_prepped  
data_list$csno_summary_prepped <- data_overlap$csno_summary_prepped 

rm(data_args, data_abundance, data_pan_core, data_overlap)


# Helpers
pal_for_tools <- function(selected_tools, tools_levels, pal_10_q) {
  sel  <- intersect(tools_levels, selected_tools)
  cols <- pal_10_q[match(sel, tools_levels)]
  stats::setNames(cols, sel)
}

# Lookup key for core_sum_prepped: "<threshold>|<proportion>|<sample_thr>"
core_sum_key <- function(thr, prop, samp) {
  paste(thr, prop, samp, sep = "|")
}