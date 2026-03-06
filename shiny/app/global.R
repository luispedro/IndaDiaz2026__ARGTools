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
# options(shiny.usecairo = TRUE)
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


# Constants & pLotting config
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
# data_list <- readRDS(file.path(ROOT_DIR, "data", "precomputed_data.rds"))
data_list <- qs::qread(file.path(ROOT_DIR, "data", "precomputed_data.qs"))


# Map the precomputed levels to the variable the UI expects
gene_classes <- data_list$levels_unigenes

# Extract standard structures expected by server.R for convenience
abundance_prepped       <- data_list$abundance_prepped
abundance_class_prepped <- data_list$abundance_class_prepped


# If you included pan/core in data_prep.R:
pan_prepped             <- data_list$pan_prepped
core_prepped            <- data_list$core_prepped
sumpan2_prepped         <- data_list$sumpan2_prepped


#For median abundance and diversity plot
pal_for_tools <- function(selected_tools, tools_levels, pal_10_q) {
  sel <- intersect(tools_levels, selected_tools)
  cols <- pal_10_q[match(sel, tools_levels)]
  stats::setNames(cols, sel)   # named vector is key
}

