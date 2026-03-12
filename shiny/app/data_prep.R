library(tidyverse)
library(qs)

DATA_DIR <- "code_R_analysis/output_abundance_diversity_resistome"
PREP_DIR <- "data"

EN <- c("human gut", "human oral", "human skin", "human nose", "human vagina",
        "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine",
        "freshwater", "soil", "amplicon", "isolate", "built-environment")

SO <- c(rep("humans", 5), rep("mammals", 4), "wastewater", "marine", "freshwater",
        "soil", rep("other", 2), "built-environment")
names(SO) <- EN

not_env <- c("amplicon", "isolate", "built-environment")
EN2 <- EN[!EN %in% not_env]
h2  <- c("humans", "mammals", "wastewater", "freshwater", "soil", "marine")

tools_levels <- c("DeepARG", "fARGene", "ABRicate-ARGANNOT", "ABRicate-MEGARes",
                  "RGI-DIAMOND", "ABRicate-CARD", "AMRFinderPlus",
                  "ABRicate-NCBI", "ResFinder", "ABRicate-ResFinder")

# Added tools_texture for the pre-join mutation
tools_texture <- c("ABRicate-ARGANNOT", "ABRicate-MEGARes", "ABRicate-CARD", 
                   "ABRicate-NCBI", "ABRicate-ResFinder")

source("code_R_analysis/helper.R")

# All discrete UI parameter values for pan/core — used to precompute all combinations
THRESHOLDS   <- c("default", "60", "70", "80")
PROPORTIONS  <- c(0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)   # threshold_proportion choices
SAMPLE_THRS  <- c(200, 250, 300, 350, 400, 450, 500)      # threshold_samples choices

metadata <- read.delim("data/metadata_GMGC10.sample.meta.tsv")


# TAB 1 - ARGs
# Precompute all 4 threshold versions of unigenes upfront

message("Loading unigenes...")
lst      <- readRDS(file.path(DATA_DIR, "results_tools.rds"))
unigenes <- as_tibble(do.call(rbind, lapply(lst, function(x)
  x[, c("query","tool","ARO","parent","parent_description","new_level","id")]))) %>%
  filter(tool %in% tools_levels) %>%
  mutate(tool = factor(tool, levels = tools_levels))

levels_unigenes <- unigenes %>%
  group_by(query) %>% slice_head(n = 1) %>% ungroup() %>%
  group_by(new_level) %>% summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n)) %>% pull(new_level)

# Precompute all 4 threshold-filtered versions so server never filters at runtime
unigenes_prepped <- list(
  "default" = unigenes,
  "60"      = unigenes %>% filter(!(tool %in% c("DeepARG","RGI-DIAMOND") & id < 60)),
  "70"      = unigenes %>% filter(!(tool %in% c("DeepARG","RGI-DIAMOND") & id < 70)),
  "80"      = unigenes %>% filter(!(tool %in% c("DeepARG","RGI-DIAMOND") & id < 80))
)

message("Saving data_args.qs ...")
qs::qsave(
  list(
    unigenes_prepped = unigenes_prepped,
    levels_unigenes  = levels_unigenes
  ),
  file.path(PREP_DIR, "data_args.qs")
)


# TAB 4 - Overlap
# Precompute CSTC/CSNO summary medians per threshold (removes runtime summarise)

message("Computing overlap data...")

make_recall <- function(uni) create_class_overlaps(uni)

recall_list <- list(
  "default" = make_recall(unigenes_prepped[["default"]]),
  "60"      = make_recall(unigenes_prepped[["60"]]),
  "70"      = make_recall(unigenes_prepped[["70"]]),
  "80"      = make_recall(unigenes_prepped[["80"]])
)

# Precompute CSTC medians (recall) per threshold — used by plot_overlap_cstc_summary
cstc_summary_prepped <- lapply(recall_list, function(df) {
  df %>%
    filter(!is.na(recall)) %>%
    ungroup() %>%
    group_by(tool_ref, new_level) %>%
    summarise(recall = median(recall), .groups = "drop")
})

# Precompute CSNO medians (fnr) per threshold — used by plot_overlap_csno_summary
csno_summary_prepped <- lapply(recall_list, function(df) {
  df %>%
    filter(!is.na(recall)) %>%
    ungroup() %>%
    group_by(tool_ref, new_level) %>%
    summarise(fnr = median(fnr), .groups = "drop")
})

message("Saving data_overlap.qs ...")
qs::qsave(
  list(
    levels_unigenes      = levels_unigenes,
    recall_fnr           = recall_list[["default"]],
    recall_fnr60         = recall_list[["60"]],
    recall_fnr70         = recall_list[["70"]],
    recall_fnr80         = recall_list[["80"]],
    cstc_summary_prepped = cstc_summary_prepped,  # NEW: precomputed medians
    csno_summary_prepped = csno_summary_prepped   # NEW: precomputed medians
  ),
  file.path(PREP_DIR, "data_overlap.qs")
)

rm(lst, unigenes, unigenes_prepped, recall_list,
   cstc_summary_prepped, csno_summary_prepped)
gc()


# TAB 2 - Abundance & Diversity

message("Loading abundance data...")

process_abundance_file <- function(file_paths) {
  bind_rows(lapply(file_paths, function(x) readRDS(file.path(DATA_DIR, x)))) %>%
    mutate(
      unigenes = distinct_unigenes_rarefied,
      habitat  = factor(habitat, levels = EN),
      habitat2 = factor(SO[habitat], levels = h2),
      tool     = factor(tool, levels = tools_levels)
    ) %>%
    filter(tool %in% tools_levels & !habitat %in% not_env, aggregation == "new_level")
}

ab_base      <- process_abundance_file(c("abundance_diversity_part1.rds",
                                         "abundance_diversity_part2.rds",
                                         "abundance_diversity_part3.rds"))
ab_60        <- process_abundance_file("abundance_diversity_60.rds")
ab_70        <- process_abundance_file("abundance_diversity_70.rds")
ab_80        <- process_abundance_file("abundance_diversity_80.rds")
ab_base_excl <- ab_base %>% filter(!tool %in% c("DeepARG", "RGI-DIAMOND"))

prep_abundance_grid <- function(df) {
  df    <- df %>% mutate(sample = factor(sample), tool = factor(tool, levels = tools_levels))
  grid  <- tidyr::expand_grid(sample = levels(df$sample), tool = levels(df$tool))
  hmap  <- df %>% group_by(sample) %>% summarise(habitat = first(habitat[!is.na(habitat)]), .groups = "drop")
  summ  <- df %>% group_by(tool, sample) %>%
    summarise(normed10m = sum(normed10m, na.rm = TRUE),
              unigenes  = sum(unigenes,  na.rm = TRUE), .groups = "drop")
  grid %>%
    left_join(summ, by = c("sample","tool")) %>%
    left_join(hmap,  by = "sample") %>%
    mutate(normed10m = replace_na(normed10m, 0), unigenes = replace_na(unigenes, 0)) %>%
    arrange(tool, sample)
}

process_ab_class <- function(df) {
  df %>% ungroup() %>%
    mutate(gene = factor(gene), sample = factor(sample)) %>%
    select(sample, gene, tool, normed10m, unigenes) %>%
    complete(sample, gene, tool, fill = list(normed10m = 0, unigenes = 0)) %>%
    mutate(
      aggregation = "new_level",
      gene        = as.character(gene),
      habitat     = metadata$habitat[ match(sample, metadata$sample_id)],
      habitat2    = metadata$habitat2[match(sample, metadata$sample_id)],
      location    = metadata$location[ match(sample, metadata$sample_id)]
    )
}

class_base <- process_ab_class(ab_base)
class_excl <- class_base %>% filter(!tool %in% c("DeepARG", "RGI-DIAMOND"))

message("Saving data_abundance.qs ...")
qs::qsave(
  list(
    levels_unigenes         = levels_unigenes,
    abundance_prepped       = list(
      "default" = prep_abundance_grid(ab_base),
      "60"      = prep_abundance_grid(bind_rows(ab_base_excl, ab_60)),
      "70"      = prep_abundance_grid(bind_rows(ab_base_excl, ab_70)),
      "80"      = prep_abundance_grid(bind_rows(ab_base_excl, ab_80))
    ),
    abundance_class_prepped = list(
      "default" = class_base,
      "60"      = bind_rows(class_excl, process_ab_class(ab_60)),
      "70"      = bind_rows(class_excl, process_ab_class(ab_70)),
      "80"      = bind_rows(class_excl, process_ab_class(ab_80))
    )
  ),
  file.path(PREP_DIR, "data_abundance.qs")
)

rm(ab_base, ab_60, ab_70, ab_80, ab_base_excl, class_base, class_excl)
gc()


# TAB 3 - Pan & Core Resistome
# Precompute ALL 196 combinations of threshold × proportion × sample_threshold
# so sum_core_adjust() is never called at runtime

message("Loading pan/core data...")

process_core_file <- function(file_name) {
  readRDS(file.path(DATA_DIR, file_name)) %>%
    rename(new_level = new_level_centroid, X = centroid) %>%
    filter(tool %in% tools_levels, !habitat %in% not_env) %>%
    mutate(habitat = factor(habitat, levels = EN2),
           tool    = factor(tool,    levels = tools_levels))
}

process_pan_file <- function(file_name) {
  readRDS(file.path(DATA_DIR, file_name)) %>%
    filter(tool %in% tools_levels, !habitat %in% not_env,
           aggregation %in% "new_level_centroid") %>%
    mutate(habitat = factor(habitat, levels = EN2),
           tool    = factor(tool,    levels = tools_levels))
}

calculate_sumpan2 <- function(pan_df) {
  pan_df %>% ungroup() %>%
    group_by(tool, habitat, aggregation, epoch) %>%
    summarise(s = sum(unigenes), .groups = "drop") %>%
    group_by(tool, habitat, aggregation) %>%
    summarise(md = median(s), mn = mean(s), sd = sd(s), .groups = "drop")
}

core_base    <- process_core_file("core_resistome.rds")
pan_base     <- process_pan_file("pan_resistome.rds")
sumpan2_base <- calculate_sumpan2(pan_base)

core_60    <- process_core_file("core_resistome60.rds")
pan_60     <- process_pan_file("pan_resistome60.rds")
sumpan2_60 <- calculate_sumpan2(pan_60)

core_70    <- process_core_file("core_resistome70.rds")
pan_70     <- process_pan_file("pan_resistome70.rds")
sumpan2_70 <- calculate_sumpan2(pan_70)

core_80    <- process_core_file("core_resistome80.rds")
pan_80     <- process_pan_file("pan_resistome80.rds")
sumpan2_80 <- calculate_sumpan2(pan_80)

tools_excl        <- c("DeepARG", "RGI-DIAMOND")
core_base_excl    <- core_base    %>% filter(!tool %in% tools_excl)
pan_base_excl     <- pan_base     %>% filter(!tool %in% tools_excl)
sumpan2_base_excl <- sumpan2_base %>% filter(!tool %in% tools_excl)

core_prepped <- list(
  "default" = core_base,
  "60"      = bind_rows(core_base_excl, core_60),
  "70"      = bind_rows(core_base_excl, core_70),
  "80"      = bind_rows(core_base_excl, core_80)
)

pan_prepped <- list(
  "default" = pan_base,
  "60"      = bind_rows(pan_base_excl, pan_60),
  "70"      = bind_rows(pan_base_excl, pan_70),
  "80"      = bind_rows(pan_base_excl, pan_80)
)

sumpan2_prepped <- list(
  "default" = sumpan2_base,
  "60"      = bind_rows(sumpan2_base_excl, sumpan2_60),
  "70"      = bind_rows(sumpan2_base_excl, sumpan2_70),
  "80"      = bind_rows(sumpan2_base_excl, sumpan2_80)
)

# Precompute all 196 combinations of threshold × proportion × sample_threshold
# key format: "<threshold>|<proportion>|<sample_thr>"  e.g. "default|0.5|450"
message("Precomputing all 196 sum_core_adjust combinations (4 x 7 x 7)...")
core_sum_prepped <- list()

for (thr in THRESHOLDS) {
  core_df <- core_prepped[[thr]]
  for (prop in PROPORTIONS) {
    for (samp in SAMPLE_THRS) {
      key <- paste(thr, prop, samp, sep = "|")
      core_sum_prepped[[key]] <- sum_core_adjust(core_df, samp, prop) %>%
        dplyr::ungroup() %>%
        dplyr::group_by(tool, habitat) %>%
        dplyr::summarise(core = sum(unigenes, na.rm = TRUE), .groups = "drop")
    }
  }
  message("  Done threshold: ", thr)
}

message("Executing final join between sumpan2 and core_sum across all combinations...")
pan_core_joined_prepped <- list()

for (key in names(core_sum_prepped)) {
  # Extract the threshold (e.g., "default", "60") to grab the matching sumpan2 data
  thr_name <- strsplit(key, "\\|")[[1]][1]
  
  sumpan2_df <- sumpan2_prepped[[thr_name]]
  core_sum   <- core_sum_prepped[[key]]
  
  # Join and apply mutations previously done in server.R
  pan_core_joined_prepped[[key]] <- sumpan2_df %>%
    dplyr::left_join(core_sum, by = c("tool", "habitat")) %>%
    dplyr::mutate(
      core    = tidyr::replace_na(core, 0),
      prop    = core / md,
      texture = ifelse(tool %in% tools_texture, "yes", "no")
    )
}
message("Saving data_pan_core.qs ...")
qs::qsave(
  list(
    levels_unigenes  = levels_unigenes,
    core_prepped     = core_prepped,
    pan_prepped      = pan_prepped,
    sumpan2_prepped  = sumpan2_prepped,
    core_sum_prepped = core_sum_prepped  # NEW: all 196 combinations precomputed
  ),
  file.path(PREP_DIR, "data_pan_core.qs")
)

rm(core_base, core_60, core_70, core_80,
   pan_base,  pan_60,  pan_70,  pan_80,
   sumpan2_base, sumpan2_60, sumpan2_70, sumpan2_80,
   core_base_excl, pan_base_excl, sumpan2_base_excl,
   core_prepped, pan_prepped, sumpan2_prepped, core_sum_prepped)
gc()

message("\nDone! Files saved to '", PREP_DIR, "':")
message("  data_args.qs      -> ARGs tab        (unigenes prepped for all 4 thresholds)")
message("  data_abundance.qs -> Abundance tab")
message("  data_pan_core.qs  -> Pan/Core tab    (196 sum_core_adjust combinations precomputed)")
message("  data_overlap.qs   -> Overlap tab     (CSTC/CSNO medians precomputed per threshold)")