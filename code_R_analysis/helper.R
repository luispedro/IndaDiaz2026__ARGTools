
# function to extract the legend from ggplots
g_legend <- function(a.gplot){
  tmp <- ggplotGrob(a.gplot)
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# function to extract the core resistome based on:
# minimum proportion of metagenomic samples in a subsample needed to count an ARG as subsample core, and
# minimum number of subsample cores an ARG needs to be in to be. counted in the core-resistome

sum_core_adjust <- function(core, cnt_subset = 900, threshold_samples = 0.5){
  return(core %>% 
           filter(cut %in% threshold_samples & cnt > cnt_subset) %>% 
           ungroup() %>% 
           group_by(gene_class, tool, habitat) %>% 
           summarise(unigenes = n_distinct(X)))
}

# gene classes

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

#  gene_classes0 and gene_classes_list are the labels for the genes in the plots 
# conversion of gene classes
gene_classes0 <- gene_classes
gene_classes0 <- gene_classes0 %>% mutate(new = ifelse(new %in% "class A beta-lactamase", "A beta-lactamase", new)) %>% 
  mutate(new = ifelse(new %in% "class B beta-lactamase", "B beta-lactamase", new)) %>% 
  mutate(new = ifelse(new %in% "class C beta-lactamase", "C beta-lactamase", new)) %>%
  mutate(new = ifelse(new %in% "class D beta-lactamase", "D beta-lactamase", new)) %>%
  mutate(new = ifelse(new %in% "beta-lactam modulation resistance", "beta-lactam modulation", new))


gene_classes_list <- gene_classes0$new
rm(gene_classes0)
names(gene_classes_list) <- gene_classes$new

# gene class vector
gene_classes <- setNames(as.list(gene_classes$new), gene_classes$old)

# function to calculate the CSC 
# qc_ref genes reported in the gene class by the reference tool 
# q_ref all genes reported by the reference tool 
# qc_comp genes reported in the gene class by the tool we are comparing against
# q_comp all genes reported by the tool we are comparing against 
new_intersect_lists <- function(qc_ref, q_ref, qc_comp, q_comp){
  A <- unlist(qc_ref)
  B <- unlist(q_ref)
  A_complement <- setdiff(B, A)
  C <- unlist(qc_comp)
  D <- unlist(q_comp)
  C_complement <- setdiff(D, C)
  x1 <- intersect(C_complement, A) # genes in a different class in comparison group but right class in reference group
  x2 <- setdiff(C, A_complement) # genes in the right class in comparison group and not in any other class in reference group
  comp_class <- unique(union(x1, x2))
  r <- ifelse(length(A) == 0 & length(comp_class) == 0, NA, 
              ifelse(length(A) == 0 & length(comp_class) != 0, NA, 
              ifelse(length(A) != 0 & length(comp_class) == 0, NA, 
                     length(intersect(A, comp_class)) / length(comp_class))))
  return(r)
}

# Calculate the CSC  for all tools and gene classes in unigenes 
create_class_overlaps <- function(unigenes){
  tools_per_unigene <- unigenes %>% ungroup()  %>% 
    arrange(query) %>% 
    group_by(query) %>% 
    mutate(n_tools = n_distinct(tool)) %>% 
    mutate(single = (n_tools ==1)) 
  
  sets0 <- tools_per_unigene %>%
    group_by(tool) %>%
    summarise(query = list(query), .groups = "drop")
  
  sets1 <- tools_per_unigene %>%
    group_by(gene_class, tool) %>%
    summarise(query = list(query), .groups = "drop")  
  
  pairwise <- sets1 %>%
    group_by(gene_class) %>%
    summarise(pairs = list(expand_grid(tool_ref = tool, tool_comp = tool)), .groups = "drop") %>%
    unnest(pairs)
  
  JI_class_other <- pairwise %>%
    left_join(sets1, by = c("gene_class", "tool_ref" = "tool")) %>%
    rename(qc_ref = query) %>%
    left_join(sets1, by = c("gene_class", "tool_comp" = "tool")) %>%
    rename(qc_comp = query) %>%
    left_join(sets0, by = c( "tool_ref" = "tool")) %>%
    rename(q_ref = query) %>% 
    left_join(sets0, by = c( "tool_comp" = "tool")) %>%
    rename(q_comp = query)
  
  JI_class_other <- JI_class_other %>% 
    rowwise() %>%
    mutate(csc = new_intersect_lists(qc_ref, q_ref, qc_comp, q_comp)) %>% 
    rowwise() %>% 
    mutate(ref_n_class = length(qc_ref), comp_n_class = length(qc_comp), ref_n_all = length(q_ref), comp_n_all = length(q_comp)) %>% 
    ungroup() %>% 
    filter(tool_ref != tool_comp) %>% 
    select(-c(qc_ref, qc_comp, q_ref, q_comp))
  return(JI_class_other)
}

# calculates the Jaccard Index for all tools in unigenes
return_overlap_tools <- function(unigenes) { 
  sets <- unigenes %>% ungroup()  %>% 
    arrange(query) %>% 
    group_by(query) %>% 
    mutate(n_tools = n_distinct(tool)) %>% 
    mutate(single = (n_tools ==1))  %>%  
    group_by(tool) %>%
    summarise(query = list(query), .groups = "drop") # put every query in a list
  
  JI_all <- expand_grid(tool_ref = sets$tool, tool_comp = sets$tool)  %>%
    left_join(sets, by = c("tool_ref" = "tool")) %>%
    rename(values1 = query) %>%
    left_join(sets, by = c("tool_comp" = "tool")) %>%
    rename(values2 = query) %>%
    mutate(jaccard = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length(union(.x, .y)))) %>%
    mutate(csc = map2_dbl(values1, values2, ~ length(intersect(.x, .y)) / length( .y))) %>%
    select(tool_ref, tool_comp, jaccard, csc) 
  
  return(JI_all)
  
}


# calculates the CSC without assigning gene classes from reference tool to the tool we compare against 
create_class_overlaps_no_shuffling <- function(unigenes){
  tools_per_unigene <- unigenes %>% ungroup()  %>% 
    arrange(query) %>% 
    group_by(query) %>% 
    mutate(n_tools = n_distinct(tool)) %>% 
    mutate(single = (n_tools ==1)) 
  
  sets0 <- tools_per_unigene %>%
    group_by(tool) %>%
    summarise(query = list(query), .groups = "drop")
  
  sets1 <- tools_per_unigene %>%
    group_by(gene_class, tool) %>%
    summarise(query = list(query), .groups = "drop")  
  
  pairwise <- sets1 %>%
    group_by(gene_class) %>%
    summarise(pairs = list(expand_grid(tool_ref = tool, tool_comp = tool)), .groups = "drop") %>%
    unnest(pairs)
  
  JI_class_other <- pairwise %>%
    left_join(sets1, by = c("gene_class", "tool_ref" = "tool")) %>%
    rename(qc_ref = query) %>%
    left_join(sets1, by = c("gene_class", "tool_comp" = "tool")) %>%
    rename(qc_comp = query) %>%
    left_join(sets0, by = c( "tool_ref" = "tool")) %>%
    rename(q_ref = query) %>% 
    left_join(sets0, by = c( "tool_comp" = "tool")) %>%
    rename(q_comp = query)
  
  JI_class_other <- JI_class_other %>% 
    rowwise() %>%
    mutate(csc = length(intersect(qc_ref, qc_comp))/ length(qc_comp)) %>% 
    rowwise() %>% 
    mutate(ref_n_class = length(qc_ref), comp_n_class = length(qc_comp), ref_n_all = length(q_ref), comp_n_all = length(q_comp)) %>% 
    ungroup() %>% 
    filter(tool_ref != tool_comp) %>% 
    select(-c(qc_ref, qc_comp, q_ref, q_comp))
  return(JI_class_other)
}
