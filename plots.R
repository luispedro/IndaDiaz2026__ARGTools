#rgi.blast
#rgi.diamond
#deeparg.norm
#fargene
#amrfinder.norm.prot
#abricate.argannot.norm
#abricate.card.norm
#abricate.megares.norm
#abricate.ncbi.norm
#abricate.resfinder.norm
#resfinder.norm

library(dplyr)
library(ggplot2)
library(ggVennDiagram)
library(gridExtra)
library(ggradar)
library(tidyverse)

df2 <- readRDS(file = "~/df2.rds")
abundance_parent <- readRDS("GitHub/arg_compare/abundance_per_tool.rds")


#pal0 <- c("#8c510a","#bf812d","#dfc27d","#f6e8c3","#c7eae5","#80cdc1","#35978f","#01665e")
pal <- c("#543005", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#f5f5f5", "#c7eae5", "#80cdc1", "#35978f", "#01665e", "#003c30")

pal2a <- c("#00B5E2", "#013B67", "#CC0C01", "#985428", "#73006D", "#E08728", "#377E60", "#FD8CC0")
pal3a <- c("#00B5E2", "#377E60", "#CC0C01", "#985428", "#73006D", "#E08728",  "#FD8CC0")
pal4 <- c("#a6611a", "#dfc27d", "#80cdc1", "#018571")
pal5 <- c("#a6611a", "#dfc27d", "grey30", "#80cdc1", "#018571")

# HABITATS
EN <- c("human gut", "human oral",  "human skin", "human nose", "human vagina", 
        "dog gut", "cat gut", "mouse gut", "pig gut", "wastewater", "marine", "freshwater",  
        "soil" , "amplicon", "isolate",  "built-environment" )
# SOURCE FOR EACH HABITAT
SO <- c(rep("Humans", 5), rep("Mammals", 4),  "Wastewater", "Marine", "Freshwater", "Soil", rep("Other", 3))
names(SO) <- EN

# changing habitats and tools to factor
abundance_parent$habitat <- factor(abundance_parent$habitat, levels = EN)
abundance_parent$habitat2 <- factor(SO[abundance_parent$habitat], levels = unique(SO))
abundance_parent$tool <- factor(abundance_parent$tool, 
                                levels = c("rgi.diamond", "rgi.diamond.prot", "rgi.blast", "deeparg", "deeparg.prot", 
                                "fargene", "fargene.prot", "resfinder", "amrfinder", "amrfinder.prot", "abricate.argannot", 
                                "abricate.card", "abricate.megares", "abricate.ncbi", "abricate.resfinder"))

# add ontology description
abundance_parent <- abundance_parent %>% mutate(parent_label = df2$Parent_Label[match(parent,df2$Parent_ID)])
# add the condense ontology abbreviations
abundance_parent <- abundance_parent %>% mutate(new_level = new_level_df$new[match(parent_label, new_level_df$old)])

## factors for new_level, we take the highest abundance per ontology by tool
factor_new_level <- abundance_parent %>% ungroup() %>% 
  group_by(tool, new_level) %>% summarise(total = sum(normed10m)) %>%
  ungroup() %>% arrange(tool, desc(total)) %>% 
  group_by(tool, new_level) %>%  ungroup() %>% select(new_level) %>% distinct() %>% pull()

factor_new_level2 <- c(factor_new_level[seq(1, length(factor_new_level), by = 4)],
                       factor_new_level[seq(1, length(factor_new_level), by = 4) + 1],
                       factor_new_level[seq(1, length(factor_new_level), by = 4) + 2],
                       factor_new_level[seq(1, length(factor_new_level), by = 4) + 3])



abundance_parent$new_level <- factor(abundance_parent$new_level, levels = factor_new_level2)

# factor for habitat
abundance_parent$habitat <- factor(abundance_parent$habitat, levels = EN)
abundance_parent$habitat <- as.character(abundance_parent$habitat)

# environments that we are not interested in
not_env <- c("built-environment", "amplicon", "isolate")
total_abundance_sample <- abundance_parent %>% group_by(sample) %>% summarise(total = sum(normed10m)) %>%
  arrange(desc(total))

# sample outliers
extreme_samples <- total_abundance_sample %>% filter(total > quantile(total, probs = .9977)) %>%
  ungroup() %>% select(sample) %>% pull()

# tools selected to plot
tool_selected <- c("rgi.diamond", "rgi.blast", "deeparg", "fargene", "amrfinder.prot", 
                   "abricate.argannot", "abricate.card", "abricate.megares", "abricate.ncbi",
                   "abricate.resfinder", "resfinder")

# number of samples per habitat
abundance_parent <- abundance_parent %>% group_by(habitat) %>% mutate(n1 = n_distinct(sample))
abundance_parent <- abundance_parent %>% group_by(habitat2) %>% mutate(n2 = n_distinct(sample))

# number of genes per sample
abundance_parent %>% group_by(sample, habitat) %>% summarise(n = n_distinct(parent))

# Abundance boxplot  habitat 2

fig1 <- abundance_parent  %>% 
  group_by(habitat, habitat2, tool, sample) %>% summarise(total = sum(normed10m)+1) %>%
  filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ggplot(aes( x = habitat2)) +
  geom_boxplot(aes(y = total, fill = tool), outlier.shape = NA) +
  scale_fill_manual(values = pal) +
  theme_minimal() +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))  + 
  ylab("Normalized abundance") +
  xlab("Source") +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90))
  #scale_y_continuous(labels = scales::log10_trans())
fig1

selected_levels <- c("AAC", "ABC-F", "APH", "Class A", "Class B", "Efflux pump", "GPA", "MFS - efflux pump", "rpoB")





# box plots
abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ungroup() %>% group_by(tool, new_level, habitat2) %>% summarise(mn = sum(normed10m)/n2[1]+1) %>%  ungroup() %>%
  ggplot(aes( x = tool, y = mn)) +
  geom_boxplot(aes(fill = tool), position = position_dodge(preserve = "single"), na.rm = TRUE, outlier.shape = NA) +
  scale_fill_manual(values = pal) +
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))  + 
  facet_grid(. ~ new_level) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90)) 

abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ungroup() %>% group_by(tool, new_level) %>% summarise(mn = sum(normed10m)/n2[1]+1)  %>%
  mutate(p = log10(mn) / max(log10(mn))) %>%
  ungroup() %>% arrange(tool, desc(p)) %>% group_by(tool) %>% slice_head(n = 10) %>%
  ggplot(aes( x = new_level, y = p, group = tool)) +
  geom_point(aes(color = tool)) +
  scale_color_manual(values = pal) +
#  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
#                labels = scales::trans_format("log10", scales::math_format(10^.x)))  + 
  facet_grid(. ~ tool, sc) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90)) 

# heatmap

factor_new_level_heatmap <- abundance_parent %>% ungroup() %>% 
  group_by(tool, new_level) %>% summarise(total = sum(normed10m)) %>%
  ungroup() %>% arrange(tool, desc(total)) %>% 
  group_by(tool, new_level) %>%  ungroup() %>% select(new_level) %>% distinct() %>% pull()


hm1 <- abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ungroup() %>% mutate(N = n_distinct(sample)) %>% 
  ungroup() %>% group_by(tool, new_level) %>% summarise(mn = sum(normed10m)/N[1])  %>%
  #ungroup() %>% group_by(tool, new_level) %>% summarise(med = median(mn))  %>%
  ungroup() %>% group_by(tool) %>% mutate(p = log10(mn+1) / max(log10(mn+1))) %>%
  ungroup()  %>% mutate(tool = factor(as.character(tool), levels = tool_selected)) %>%
  ungroup()  %>% complete(tool, new_level, fill = list(p = 0)) %>%
  mutate(new_level = factor(as.character(new_level), levels = factor_new_level_heatmap)) %>% 
   ggplot(aes( x = tool, y = new_level, fill = p)) +
   geom_tile() +
   #scale_fill_gradient(low = "white", high = pal[length(pal)]) +
  scale_fill_viridis_c() + 
#   facet_grid(. ~ habitat) +
  xlab("Tool") + 
  ylab("Ontology") + 
   theme_minimal() +
  labs(fill = "") +
   theme(
     legend.position = "bottom",
     panel.grid.major = element_blank(), 
     panel.grid.minor = element_blank(),
     panel.background = element_blank(),
     axis.text.x = element_text(angle = 90)) 


hm2 <- abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ungroup() %>% mutate(N = n_distinct(sample)) %>% 
  ungroup() %>% group_by(tool, new_level) %>% summarise(mn = sum(normed10m)/N[1])  %>%
  #ungroup() %>% group_by(tool, new_level) %>% summarise(med = median(mn))  %>%
  ungroup() %>% group_by(tool) %>% mutate(p = log10(mn+1) / max(log10(mn+1))) %>%
  ungroup()  %>% mutate(tool = factor(as.character(tool), levels = tool_selected)) %>%
  ungroup()  %>% complete(tool, new_level, fill = list(p = 0)) %>%
  mutate(new_level = factor(as.character(new_level), levels = factor_new_level_heatmap)) %>% 
  ggplot(aes( x = new_level, y = tool, fill = p)) +
  geom_tile() +
  #scale_fill_gradient(low = "white", high = pal[length(pal)]) +
  scale_fill_viridis_c() + 
  #   facet_grid(. ~ habitat) +
  ylab("Tool") + 
  xlab("Ontology") + 
  theme_minimal() +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90)) 

hm3 <- abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  ungroup() %>% mutate(N = n_distinct(sample)) %>% 
  ungroup() %>% group_by(tool, new_level) %>% summarise(mn = sum(normed10m)/N[1])  %>%
  #ungroup() %>% group_by(tool, new_level) %>% summarise(med = median(mn))  %>%
  ungroup() %>% group_by(new_level) %>% mutate(p = log10(mn+1) / max(log10(mn+1))) %>%
  ungroup()  %>% mutate(tool = factor(as.character(tool), levels = tool_selected)) %>%
  ungroup()  %>% complete(tool, new_level, fill = list(p = 0)) %>%
  mutate(new_level = factor(as.character(new_level), levels = factor_new_level_heatmap)) %>% 
  ggplot(aes( x = new_level, y = tool, fill = p)) +
  geom_tile() +
  #scale_fill_gradient(low = "white", high = pal[length(pal)]) +
  scale_fill_viridis_c() + 
  #   facet_grid(. ~ habitat) +
  ylab("Tool") + 
  xlab("Ontology") + 
  theme_minimal() +
  labs(fill = "") +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    axis.text.x = element_text(angle = 90)) 


png("hm_by_tool.png",  width= 1600,  height = 800, res = 150)
hm2
dev.off()

png("hm_by_class.png" , width = 1600,  height = 800, res = 150)
hm3
dev.off()








##### Radial per tool

rad0 <- abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  group_by(new_level) %>% mutate(N = n_distinct(sample)) %>% 
  ungroup() %>% group_by(tool, new_level) %>% summarise(mn = log10(sum(normed10m)/N[1] +1))  %>%
  ungroup() %>% arrange(tool, desc(mn))


df_long <- rad0 %>%
  pivot_longer(-c(tool, new_level), names_to = "variable", values_to = "value") %>%
  ungroup() %>%
  complete(tool, new_level, variable, fill = list(value = 0)) %>%
  group_by(tool, variable) %>%
  mutate(angle = 2 * pi * (row_number() - 1) / n(),  # angle for each variable
         x = sin(angle) * value,
         y = cos(angle) * value)

# adding the first point to close the polygone
df_poly_mn <- df_long %>% filter(variable %in% "mn") %>%
  group_by(tool, variable) %>%
  do(rbind(., slice(., 1))) %>%
  mutate(group = rep(unique(df$group), each = nrow(df_long) + 1)) 

data.frame(df_poly_mn %>% filter(tool == "resfinder"))

make_circle <- function(radius, n = 500, center_x = 0, center_y = 0) {
  tibble(
    x = cos(seq(0, 2 * pi, length.out = n)) * radius + center_x,
    y = sin(seq(0, 2 * pi, length.out = n)) * radius + center_y,
    r = radius
  )
}

# Combine all circles into one dataframe
circles_df_7 <- bind_rows(lapply(c( 0.5, 1, 1.5, 2, 2.5, 3), make_circle))

df_radial_plot <- df_poly_mn %>% filter( tool %in% c("fargene","deeparg", "rgi.diamond", "resfinder"), value > 0) %>% 
  ungroup() %>% group_by(tool, new_level) %>% 
  mutate(label = if_else(value == max(value) & value > 1, new_level, "")) %>% 
  mutate(hjust_val = ifelse(x < 0, 1, 0)) %>%
  mutate(hjust_val = ifelse(angle %in% c(0, pi / 2, pi, pi*3/2), 0.5, hjust_val)) %>%
  mutate(angle_text = abs(atan2(x,y)* 180/pi) + 90) %>%
  mutate(angle_text = ifelse(x > 0 & y > 0, 180 - angle_text, 
                             ifelse(x > 0 & y < 0, 180 - angle_text, 
                                    ifelse( x < 0 & y > 0,  angle_text - 180, 
                                            ifelse(x < 0 & y < 0, angle_text - 180, angle_text )))))

df_radial_plot <- bind_rows(df_radial_plot, df_radial_plot %>% mutate(x = 0, y = 0, label = ""))



fig2_2 <- df_radial_plot %>%
  ggplot( aes(x = x, y = y)) +
  geom_polygon(aes(fill = tool), alpha = 0.4) +
  geom_point(aes(color = tool),  size = 2, alpha = 0.6) +
  geom_text(aes(label = label, hjust = hjust_val, angle = angle_text), vjust = -0.5, size = 2) +
  geom_path(data = circles_df_7, aes(x, y, group = r), color = "grey",  show.legend = F, alpha = 0.4) + 
  geom_line(aes(group = new_level), alpha = 0.4) +
  scale_fill_manual(values = pal4) +
  scale_color_manual(values = pal4) +
  coord_equal() +
  facet_grid(. ~ tool) +
  #geom_point(data = data.frame(x = 0, y = 0, tool="a"), aes(x = x, y = y), color = "black", size = 2, show.legend = F) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(fill = "", color="") +
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  labs(title = "")

fig2_2


metadata <- read.delim("metadata_GMGC10.sample.meta.tsv")

query_tool_new_level <- do.call(rbind, 
        lapply(seq_along(1:length(lst1)), 
               function(j) data.frame(tool = rep(names(lst1)[j], dim(lst[[j]])[1]), 
                                      query = lst1[[j]]$query, 
                                      ARO = lst1[[j]]$ARO,
                                      new_level = lst1[[j]]$new_level,
                                      sample = lst1[[j]]$query )))
tool_name = c("rgi.diamond", "rgi.diamond.prot", "rgi.blast", "deeparg", "deeparg.prot", 
  "fargene", "fargene.prot", "resfinder", "amrfinder", "amrfinder.prot", "abricate.argannot", 
  "abricate.card", "abricate.megares", "abricate.ncbi", "abricate.resfinder")

names(tool_name) <- c("rgi.diamond", "rgi.diamond.prot", "rgi.blast", "deeparg.norm", "deeparg.norm.prot",
                      "fargene", "fargene.prot", "resfinder.norm", "amrfinder.norm", "amrfinder.norm.prot", "abricate.argannot.norm", 
                      "abricate.card.norm", "abricate.megares.norm", "abricate.ncbi.norm", "abricate.resfinder.norm")

query_tool_new_level <- query_tool_new_level %>% mutate(tool = factor(tool_name[tool], levels = levels(df_radial_plot$tool)))

query_tool_new_level_count <- query_tool_new_level %>% group_by(tool, new_level) %>% 
  summarise(n = n_distinct(query))  %>% 
  arrange(tool, desc(n)) %>% left_join(df_radial_plot %>% select(tool, new_level, angle)) %>%
  mutate(x = sin(angle) * log10(n+1),
         y = cos(angle) * log10(n+1), 
         t = "diversity")


fig3_2 <- query_tool_new_level_count %>% filter( tool %in% c("fargene","deeparg", "rgi.diamond", "resfinder"), log10(n+1)>0) %>%
  ggplot( aes(x = x, y = y)) +
#  geom_polygon(aes(fill = tool), alpha = 0.4) +
  geom_point(aes(color = tool),  size = 2, alpha = 0.6) +
  #geom_text(aes(label = label, hjust = hjust_val, angle = angle_text), vjust = -0.5, size = 2) +
  geom_path(data = circles_df_7, aes(x, y, group = r), color = "grey",  show.legend = F, alpha = 0.4) + 
  geom_line(aes(group = new_level), alpha = 0.4) +
  scale_fill_manual(values = pal4) +
  scale_color_manual(values = pal4) +
  coord_equal() +
  facet_grid(. ~ tool) +
  #geom_point(data = data.frame(x = 0, y = 0, tool="a"), aes(x = x, y = y), color = "black", size = 2, show.legend = F) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(fill = "", color="") +
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    plot.title = element_text(face = "bold"),
    legend.position = "none",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  labs(title = "")
fig3_2



## Radial
rad0 <- abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  group_by(tool, new_level, habitat) %>% 
  summarise(l = sum(log(normed10m + 1)), median = median(log(normed10m + 1)), mn = sum(log(normed10m + 1)) / max(n1)) %>% 
  ungroup() %>% arrange(tool, desc(l))

# this is why we cannot use median or we need to complement it with zeros
# not all samples detect all genes, the median on the genes will be innacurate due to the incomplete data
abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  group_by(tool, new_level, habitat) %>% 
  summarise(n = n_distinct(sample)) %>% ungroup() %>% arrange(tool, habitat, new_level)

rad1 <- abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected, !sample %in% extreme_samples) %>% 
  group_by(tool, new_level, habitat) %>% summarise(l = sum(log(normed10m + 1)), median = median(log(normed10m + 1))) %>% 
  ungroup() %>% arrange(tool, desc(l)) %>% 
  group_by(tool, habitat) %>% mutate(m = max(l), p = l / max(l) * 100)


df_long <- rad0 %>%
  pivot_longer(-c(tool, new_level, habitat), names_to = "variable", values_to = "value") %>%
  ungroup() %>%
  complete(tool, new_level, habitat, variable, fill = list(value = 0)) %>%
  group_by(tool, variable, habitat) %>%
  mutate(angle = 2 * pi * (row_number() - 1) / n(),  # angle for each variable
         x = sin(angle) * value,
         y = cos(angle) * value)

df_poly <- df_long %>% filter(variable %in% "median") %>%
  group_by(tool, habitat, variable) %>%
  do(rbind(., slice(., 1))) %>%
  mutate(group = rep(unique(df$group), each = nrow(df_long) + 1)) 

df_poly_mn <- df_long %>% filter(variable %in% "mn") %>%
  group_by(tool, habitat, variable) %>%
  do(rbind(., slice(., 1))) %>%
  mutate(group = rep(unique(df$group), each = nrow(df_long) + 1)) 


make_circle <- function(radius, n = 500, center_x = 0, center_y = 0) {
  tibble(
    x = cos(seq(0, 2 * pi, length.out = n)) * radius + center_x,
    y = sin(seq(0, 2 * pi, length.out = n)) * radius + center_y,
    r = radius
  )
}



# Combine all circles into one dataframe
circles_df <- bind_rows(lapply(c( 1, 3, 5), make_circle))
circles_df$tool = as.character(circles_df$r)

circles_df_7 <- bind_rows(lapply(c( 1, 3, 5, 7, 9, 11), make_circle))
circles_df_7$tool = as.character(circles_df_7$r)


df_poly$habitat <- factor(df_poly$habitat, levels = EN)
df_poly_mn$habitat <- factor(df_poly_mn$habitat, levels = EN)



fig2 <- df_poly_mn %>% filter(habitat %in% c("human gut", "wastewater", "soil", "pig gut"), 
                           tool %in% c("fargene","deeparg", "rgi.diamond", "resfinder"), value > 0) %>% 
  ungroup() %>% group_by(new_level, habitat) %>% mutate(label = if_else(value == max(value) & value > 1.5, new_level, "")) %>% 
  mutate(hjust_val = ifelse(x < 0, 1, 0)) %>%
  mutate(hjust_val = ifelse(angle %in% c(0, pi / 2, pi, pi*3/2), 0.5, hjust_val)) %>%
  mutate(angle_text = abs(atan2(x,y)* 180/pi) + 90) %>%
  mutate(angle_text = ifelse(x > 0 & y > 0, 180 - angle_text, 
                             ifelse(x > 0 & y < 0, 180 - angle_text, 
                                    ifelse( x < 0 & y > 0,  angle_text - 180, 
                                            ifelse(x < 0 & y < 0, angle_text - 180, angle_text ))))) %>%
  ggplot( aes(x = x, y = y)) +
  geom_polygon(aes(fill = tool), alpha = 0.4) +
  geom_point(aes(color = tool),  size = 2, alpha = 0.6) +
  geom_text(aes(label = label, hjust = hjust_val, angle = angle_text), vjust = -0.5, size = 2) +
  geom_path(data = circles_df_7, aes(x, y), color = "grey",  show.legend = F, alpha = 0.4) + 
  geom_line(aes(group = new_level), alpha = 0.4) +
  scale_fill_manual(values = pal4) +
  scale_color_manual(values = pal4) +
  coord_equal() +
  facet_grid(. ~ habitat) +
  geom_point(data = data.frame(x = 0, y = 0, tool="a"), aes(x = x, y = y), color = "black", size = 2, show.legend = F) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(fill = "", color="") +
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  labs(title = "")


fig2

fig3 <- df_poly_mn %>% filter(habitat %in% c("dog gut","cat gut","mouse gut", "pig gut"), 
                           tool %in% c("fargene","deeparg", "rgi.diamond", "resfinder"), value > 0) %>% 
  ungroup() %>% group_by(new_level, habitat) %>% mutate(label = if_else(value == max(value) & value > 1.5, new_level, "")) %>% 
  mutate(hjust_val = ifelse(x < 0, 1, 0)) %>%
  mutate(hjust_val = ifelse(angle %in% c(0, pi / 2, pi, pi*3/4), 0.5, hjust_val)) %>%
  mutate(angle_text = abs(atan2(x,y)* 180/pi) + 90) %>%
  mutate(angle_text = ifelse(x > 0 & y > 0, 180 - angle_text, 
                             ifelse(x > 0 & y < 0, 180 - angle_text, 
                                    ifelse( x < 0 & y > 0,  angle_text - 180, 
                                            ifelse(x < 0 & y < 0, angle_text - 180, angle_text ))))) %>%
  ggplot( aes(x = x, y = y)) +
  geom_polygon(aes(fill = tool), alpha = 0.4) +
  geom_point(aes(color = tool),  size = 2, alpha = 0.6) +
  geom_text(aes(label = label, hjust = hjust_val, angle = angle_text), vjust = -0.5, size = 2) +
  geom_path(data = circles_df_7, aes(x, y), color = "grey",  show.legend = F, alpha = 0.4) + 
  geom_line(aes(group = new_level), alpha = 0.4) +
  scale_fill_manual(values = pal4) +
  scale_color_manual(values = pal4) +
  coord_equal() +
  facet_grid(. ~ habitat) +
  geom_point(data = data.frame(x = 0, y = 0, tool="a"), aes(x = x, y = y), color = "black", size = 2, show.legend = F) +
  theme_minimal() +
  xlab("") +
  ylab("") +
  labs(fill = "", color="") +
  theme(
    strip.text = element_text(face = "bold"),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()) +
  labs(title = "")
fig3

g_legend <- function(a.gplot){
  tmp <- ggplotGrob(a.gplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


fig2_3 <- grid.arrange(fig2 + theme(legend.position = "none"), 
             fig3 + theme(legend.position = "none"), 
             g_legend(fig2), layout_matrix = rbind(rep(1, 5), rep(1, 5), rep(1, 5), rep(1, 5), 
                                               rep(2, 5), rep(2, 5), rep(2, 5), rep(2, 5),
                                               c(NA, NA, 3, NA, NA)))

pdf("fig1.pdf", width = 14, height = 10)
fig1
dev.off()

png("fig1.png", , width = 1600, height = 800, res = 150)
fig1
dev.off()

png("fig2.png", , width = 1600, height = 800, res = 150)
fig2
dev.off()

png("fig3.png", , width = 1600, height = 800, res = 150)
fig3
dev.off()
# bar saggered plot 




########################################################################################
########################################################################################
########################################################################################
########################################################################################


lst <- readRDS("GitHub/arg_compare/results_tools.rds")

lst$deeparg.norm.prot

number_unigenes_per_tool <- sapply(lst, function(x) length(unique(x$query)))

names(lst)

unigenes_per_tool <- do.call(rbind, lapply(seq_along(1:length(lst)), function(j) data.frame(tool = rep(names(lst)[j], dim(lst[[j]])[1]), query = lst[[j]]$query, 
                                                                        ARO = lst[[j]]$ARO, parent = lst[[j]]$parent, 
                                                                        parent_description = lst[[j]]$parent_description)))

unigenes_per_tool <- unigenes_per_tool %>% filter(!is.na(parent)) %>% distinct()

several_parents <- unigenes_per_tool %>% group_by(query) %>% 
  summarise(n = n_distinct(parent)) %>% 
  ungroup() %>% arrange(desc(n)) %>% filter(n>1) %>% 
  select(query) %>% pull()


different_class <- unigenes_per_tool %>% filter(query %in% several_parents) %>% 
  arrange(query) %>% group_by(query)  %>%
  summarise(pairs = list(combn(parent_description, 2, simplify = FALSE))) %>%
  unnest(pairs) %>% group_by(query, pairs) %>% 
  mutate(
    item1 = map_chr(pairs, ~ .x[1]),
    item2 = map_chr(pairs, ~ .x[2])) %>% 
  filter(item1 != item2) %>% 
  mutate(i1 = min(item1, item2),
         i2 = max(item1, item2)) %>% ungroup() %>%
  select(query, i1, i2) %>%  mutate(pair = paste(i1, i2, sep = "@")) %>%
  group_by(query, pair) %>% distinct() %>% ungroup() %>% 
  group_by(pair) %>% summarise(n = n()) %>% 
  mutate(a1 = sapply(strsplit(pair, split = "@"), function(x) x[1]), 
         a2 = sapply(strsplit(pair, split = "@"), function(x) x[2])) %>% 
  ungroup() %>% arrange(desc(n)) %>% select(-pair)

data.frame(different_class)

tool_2 <- c("deeparg.norm", "rgi.diamond", "fargene", "amrfinder.norm.prot", "abricate.argannot.norm",
            "abricate.card.norm", "abricate.megares.norm", "abricate.ncbi.norm", "abricate.resfinder.norm",
            "resfinder.norm")

tools_per_unigene <- unigenes_per_tool %>% ungroup() %>% filter(tool %in% tool_2) %>% 
  arrange(query) %>% 
  group_by(query) %>% 
  mutate(n_tools = n()) %>% 
  mutate(single = (n_tools ==1)) 


tools_per_unigene$new_level <- new_level_df$new[match(tools_per_unigene$parent_description, new_level_df$old)]



plot_counts <- tools_per_unigene %>% group_by(tool, single) %>% summarise(n=n()) %>%
  ggplot(aes(x = tool, y=n)) +
  geom_col(aes (fill = tool, alpha = !single), color = "black") +
  ggtitle("") +
  theme_minimal() +
  ylab("Unigenes") +
  xlab("Tool") +
  scale_fill_manual(values = c(pal4,pal4,pal4)) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(labels = scales::comma)
plot_counts

dev.off()
png("counts.png", , width = 1600, height = 800, res = 150)
plot_counts
dev.off()




plot_counts_deeparg <- tools_per_unigene %>% filter(tool %in% "deeparg.norm") %>%
  group_by(tool, new_level, single) %>% summarise(n=n()) %>%
  ggplot(aes(x = new_level, y=n)) +
  geom_col(aes (fill = tool, alpha = !single), color = "black") +
  ggtitle("") +
  theme_minimal() +
  ylab("Unigenes") +
  xlab("Tool") +
  scale_fill_manual(values = c(pal4,pal4,pal4)) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(labels = scales::comma)
plot_counts_deeparg


plot_counts_rgi <- tools_per_unigene %>% filter(tool %in% "rgi.diamond") %>%
  group_by(tool, new_level, single) %>% summarise(n=n()) %>%
  ggplot(aes(x = new_level, y = n)) +
  geom_col(aes (fill = tool, alpha = !single), color = "black") +
  ggtitle("") +
  theme_minimal() +
  ylab("Unigenes") +
  xlab("Tool") +
  scale_fill_manual(values = c(pal4,pal4,pal4)) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(labels = scales::comma)
plot_counts_rgi


plot_counts_fargene <- tools_per_unigene %>% filter(tool %in% "fargene") %>%
  group_by(tool, new_level, single) %>% summarise(n=n()) %>%
  ggplot(aes(x = new_level, y = n)) +
  geom_col(aes (fill = tool, alpha = !single), color = "black") +
  ggtitle("") +
  theme_minimal() +
  ylab("Unigenes") +
  xlab("Tool") +
  scale_fill_manual(values = c(pal4,pal4,pal4)) +
  theme(legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 90)) +
  scale_y_continuous(labels = scales::comma)
plot_counts_fargene


#lst1 <- lapply(lst, function(x) cbind.data.frame(x, 
#                              data.frame(single = tools_per_unigene$single[match(x$query, tools_per_unigene$query)],
#                              new_level = tools_per_unigene$new_level[match(x$query, tools_per_unigene$query)],
#                              n_tools = tools_per_unigene$n_tools[match(x$query, tools_per_unigene$query)])))


lst1 <- lapply(lst, function(x) cbind.data.frame(x, 
                                                 data.frame(single = tools_per_unigene$single[match(x$query, tools_per_unigene$query)],
                                                            new_level = new_level_df$new[match(x$parent_description, new_level_df$old)],
                                                            n_tools = tools_per_unigene$n_tools[match(x$query, tools_per_unigene$query)])))

deeparg_highest <- lst1$deeparg.norm %>% filter(!is.na(single) & !is.na(new_level) & new_level %in% c("Efflux p.","MFS - efflux p.","rpoB", "GPA", "Cell wall charge")) %>% 
  ungroup() %>%
  mutate(X1 = cut(identity, breaks = 30), X2 = cut(alignment.bitscore, breaks = 30)) %>%
  mutate(X1 =  as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", X1))) %>% 
  mutate(X2 =  as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", X2))) %>% 
  mutate(X1 =  factor(X1, levels = sort(unique(X1)))) %>% 
  mutate(X2 =  factor(X2, levels = sort(unique(X2)))) %>% 
  mutate(new_level = new_level_df$new[match(parent_description, new_level_df$old)]) %>% 
  group_by(new_level, X1, X2, single) %>% summarise(n = n()) %>% 
  ungroup() %>% group_by(new_level) %>% mutate( p = n/max(n)) %>%
  ungroup()  %>% complete(new_level, single, X1, X2, fill = list(n = 0,  p= 0)) %>%
  mutate(single =ifelse(single, "unique", "shared")) %>%
  ggplot(aes(x = X1, y = X2, fill = p)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_viridis_c() + 
  #scale_fill_continuous(low = "#f6e8c3", high = "#543005", na.value = "white") +
  facet_grid(single ~ new_level) +
  xlab("Identity") +
  ylab("Bit score") +
  labs(fill = "") +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    plot.title = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 90, size = 5),
    axis.text.y  = element_text(size = 5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  labs(title = "")
deeparg_highest


rgidiamond_highest <- lst1$rgi.diamond %>% filter(!is.na(single) & !is.na(new_level) & new_level %in% c("Efflux p.","MFS - efflux p.","RIF", "GPA", "NTR", "ABC-F")) %>% 
  ungroup() %>%
  mutate(X1 = cut(Best_Identities, breaks = 30, ), X2 = cut(Best_Hit_Bitscore, breaks = 30)) %>%
  mutate(X1 =  as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", X1))) %>% 
  mutate(X2 =  as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", X2))) %>% 
  mutate(X1 =  factor(X1, levels = sort(unique(X1)))) %>% 
  mutate(X2 =  factor(X2, levels = sort(unique(X2)))) %>% 
  mutate(new_level = new_level_df$new[match(parent_description, new_level_df$old)]) %>% 
  group_by(new_level, X1, X2, single) %>% summarise(n = n()) %>% 
  ungroup() %>% group_by(new_level) %>% mutate( p = n/max(n)) %>%
  ungroup()  %>% complete(new_level, single, X1, X2, fill = list(n = 0,  p= 0)) %>%
  mutate(single =ifelse(single, "unique", "shared")) %>%
  ggplot(aes(x = X1, y = X2, fill = p)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_viridis_c() +
  #scale_fill_continuous(low = "white", high = "red", na.value = "white") +
  facet_grid(single ~ new_level) +
  xlab("Identity") +
  ylab("Bit score") +
  labs(fill = "") +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    plot.title = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 90, size = 5),
    axis.text.y  = element_text(size = 5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  labs(title = "")
rgidiamond_highest


png("deeparg.png", , width = 1600, height = 800, res = 150)
deeparg_highest
dev.off()

png("rgidiamond.png", , width = 1600, height = 800, res = 150)
rgidiamond_highest
dev.off()


rgi.fargene <- read.delim("GitHub/arg_compare/check_missing_annot/fargene_predicted_fna.txt")
rgi.fargene <- rgi.fargene[,-c(1,18,19,20)]
rgi.fargene$query <- gsub('tmp_', '', gsub('.{2}$', '', gsub(' ', '', rgi.fargene$Contig)))

fargene_highest <- lst1$fargene %>% filter(!is.na(single) & !is.na(new_level)) %>% # & new_level %in% c("Class A", "Class B", "AAC", "APH", "Class D")) %>% 
  ungroup() %>%
  mutate(Best_Identities = rgi.fargene$Best_Identities[match(query, rgi.fargene$query)]) %>%
  mutate(X1 = cut(Best_Identities, breaks = 30, ), X2 = cut(hmm, breaks = 30)) %>%
  mutate(X1 =  as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", X1))) %>% 
  mutate(X2 =  as.numeric(sub("[^,]*,([^]]*)\\]", "\\1", X2))) %>% 
  mutate(X1 =  factor(X1, levels = sort(unique(X1)))) %>% 
  mutate(X2 =  factor(X2, levels = sort(unique(X2)))) %>% 
  mutate(new_level = new_level_df$new[match(parent_description, new_level_df$old)]) %>% 
  group_by(new_level, X1, X2, single) %>% summarise(n = n()) %>% 
  ungroup() %>% group_by(new_level) %>% mutate( p = n/max(n)) %>%
  ungroup()  %>% complete(new_level, single, X1, X2, fill = list(n = 0,  p= 0)) %>%
  mutate(single =ifelse(single, "unique", "shared")) %>%
  ggplot(aes(x = X1, y = X2, fill = p)) +
  geom_tile() +
  theme_minimal() +
  scale_fill_viridis_c() + 
  #scale_fill_continuous(low = "#f6e8c3", high = "#543005", na.value = "white") +
  facet_grid(single ~ new_level) +
  xlab("Identity") +
  ylab("Bit score") +
  labs(fill = "") +
  theme(
    strip.text = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "lightgrey", color = NA),
    plot.title = element_text(face = "bold"),
    axis.text.x  = element_text(angle = 90, size = 5),
    axis.text.y  = element_text(size = 5),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1)) +
  labs(title = "")


fargene_highest
png("fargene.png", , width = 1600, height = 800, res = 150)
fargene_highest
dev.off()



##
##

new_level <- c("chloramphenicol phosphotransferase", "CHL ph.",
"ciprofloxacin phosphotransferase", "CIP ph.",
"viomycin phosphotransferase", "VIO ph. ",
"major facilitator superfamily (MFS) antibiotic efflux pump", "MFS - efflux p.",
"class C beta-lactamase", "Class C",
"aminoglycoside bifunctional resistance protein", "Aminoglycoside b.",
"streptothricin acetyltransferase (SAT)", "SAT",
"AAC(2')", "AAC",
"AAC(6')", "AAC",
"AAC(3)", "AAC",
"APH(3')", "APH",
"APH(6)" , "APH",
"class A beta-lactamase", "Class A",
"macrolide phosphotransferase (MPH)", "MPH",
"class B (metallo-) beta-lactamase", "Class B",
"class D beta-lactamase", "Class D",
"quinolone resistance protein (qnr)", "QNR",
"Erm 23S ribosomal RNA methyltransferase", "ERM",
"APH(2'')", "APH",
"tetracycline inactivation enzyme", "TET - enzyme",
"tetracycline-resistant ribosomal protection protein", "TET - RPG",
"gene(s) or protein(s) associated with a glycopeptide resistance cluster", "GPA",
"protein(s) and two-component regulatory system modulating antibiotic efflux", "Efflux p.",
"fosfomycin inactivation enzyme", "FO",
"antibiotic resistant dihydrofolate reductase", "DHFR",
"rifampin-resistant RNA polymerase-binding protein", "RIF",
"antibiotic resistant gene variant or mutant", "V/M",            
"gene involved in antibiotic sequestration", "Sequestration",
"gene modulating beta-lactam resistance", "Beta-lactam modulation",
"rifampin inactivation enzyme", "RIF",
"nitroimidazole reductase", "NTR",
"cpa acetyltransferase", "CPA ac.",
"protein(s) conferring resistance via host-dependent nutrient acquisition", "Host-dependent",
"protein modulating permeability to antibiotic", "Permeability modulation",
"chloramphenicol acetyltransferase (CAT)", "CAT",
"streptogramin inactivation enzyme", "Streptogramin enzyme",
"gene altering cell wall charge", "Cell wall charge",
"ANT(2'')", "ANT",
"ANT(9)", "ANT",
"ANT(4')", "ANT",
"APH(3'')", "APH",
"lincosamide nucleotidyltransferase (LNU)", "LNU",
"ANT(3'')", "ANT",
"gene involved in self-resistance to antibiotic", "Self-resistance",
"ANT(6)", "ANT",
"APH(9)", "APH",
"APH(7'')", "APH",
"sulfonamide resistant sul", "SUL",
"ABC-F ATP-binding cassette ribosomal protection protein", "ABC-F",
"macrolide glycosyltransferase", "MGT",
"macrolide esterase", "MEs",
"fusidic acid inactivation enzyme", "FUS enzyme",
"Bah amidohydrolase", "Bah amidohydrolase",
"Target protecting FusB-type protein conferring resistance to Fusidic acid", "FUS protection",
"APH(4)", "APH",
"gene conferring resistance via absence", "Absence",
"glycopeptide resistance gene cluster", "GPA",
"beta-lactam resistant penicillin-binding proteins", "PBP",
"Edeine acetyltransferase", "EdeQ",
"rifamycin-resistant beta-subunit of RNA polymerase (rpoB)", "rpoB",
"efflux pump complex or subunit conferring antibiotic resistance", "Efflux p.",
"protein(s) conferring antibiotic resistance via molecular bypass", "Molecular bypass",
"antibiotic target modifying enzyme", "Target-modifying enzyme",
"antibiotic inactivation enzyme", "Inactivation enzyme")

odd_vals  <- new_level[seq(1, length(new_level), by = 2)]
even_vals <- new_level[seq(2, length(new_level), by = 2)]
new_level_df <- df <- data.frame(old = odd_vals, new = even_vals)
rm(new_level, even_vals, odd_vals)


saveRDS(df2, file = "~/df2.rds")
unique(df2$Parent_Label)



#############################################
#############################################
# HEATMAP 

j = 4
tool_selected[j]

m <- abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected[j], !sample %in% extreme_samples) %>%
  select(sample, parent, habitat, tool, normed10m) %>% ungroup() %>% 
  group_by(sample, parent, tool, habitat) %>%
  mutate(normed10m = replace_na(normed10m, 0)) %>%
  pivot_wider(names_from = parent, values_from = normed10m, values_fill = 0) %>%
  ungroup() %>%  select(!c(sample, habitat, tool)) %>%
  select(where(~ is.numeric(.) && sum(.) > 0))



n <- abundance_parent %>% filter(!habitat %in% not_env, tool %in% tool_selected[j], !sample %in% extreme_samples) %>%
  select(sample, parent, habitat, habitat2, tool, normed10m) %>% ungroup() %>% 
  group_by(sample, parent, tool, habitat) %>%
  mutate(normed10m = replace_na(normed10m, 0)) %>%
  pivot_wider(names_from = parent, values_from = normed10m, values_fill = 0) %>%
  ungroup() %>%  select(c(sample, habitat, tool, habitat2))


mpca <- prcomp(log(m+1), center = TRUE, scale. = TRUE)
dfpca <- as.data.frame(mpca$x)
dfpca$habitat <- n$habitat
dfpca$habitat2 <- n$habitat2


set.seed(2025)

dfpca2 <- dfpca %>% group_by(habitat) %>%
  group_modify(~ slice_sample(.x, n = min(500, nrow(.x)), replace = F))

grid.arrange(
  ggplot(dfpca2, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = habitat)),
  ggplot(dfpca2, aes(x = PC1, y = PC3)) +
    geom_point(aes(color = habitat)) , 
  ggplot(dfpca2, aes(x = PC2, y = PC3)) +
    geom_point(aes(color = habitat)),
  layout_matrix = rbind(c(1, 2), c(NA, 3)))


