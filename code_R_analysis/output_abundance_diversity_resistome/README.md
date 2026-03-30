
# File descriptions

## gene class - conversion_ARO_parent_new_level.csv 
Alternatively, `conversion_ARO_parent_new_level.rds`.

- Term_ID: ARO assigned by RGI or argNorm  
- Term_Label: ARO's description 
- Parent_ID: Higher onthology
- Parent_Label: Description 
- new_label: Gene class assigned

## detected unigenes - unigenes_per_tool.csv.gz 

- query: unigene name 
- tool: tool + database
- new_level: gene class
- id: identity level of the hit


## Abundance and diversity - abundance_diversity.csv.gz

sample: metagenomic sample ID
gene: gene class (new_level) 
tool: tool + database 
abundance: normalized read by gene size scaled to 1 million reads per metagenome
richness: number of unique genes after rarefaction
richness_no_rarified: number of unique genes without rarefaction


## core_resistome.csv.gz and pan_resistome.csv.gz
Alternatively, `core_resistome.rds` and `pan_resistome.rds`.

The identified ARGs were clusted at 90% identity level and the centroids were assigned to each gene in the clusters. The pan- and core-resistomes were calculated for each environment by randomly taking 500 subsamples of size 100 metagenomic samples from each habitat. For each subsample, we stored (1) the unique number of centroids (unigene ID) present in any sample (alpha diversity) grouped by tool, habitat, and aggregation level (ARO, parent description, new_level) for the pan-resistome calculations, and (2) the centroids present in at least a proportion of the samples determined by the "cut" level (commonly encountered gene set) grouped by tool, habitat, and aggregation level (new_level). 

The size of the pan-resistome for each environment will be calculated as the average alpha diversity over the 500 subsamples, determined by "epoch", grouped by tool, habitat, and aggregation level, and the core-resistome as the number of centroids in at least an arbitraty number determined by "cnt", e.g. 450, of the commonly encountered gene sets grouped by tool, habitat, and aggregation level.

### Core-resistome

- centroid: in the context of core- and pan-resistome, the gene
- new level centroid: gene class of the centroid
- tool: tool + database
- habitat: type of habitat that the core-resistome is calculated for
- cut: the centroid was found in at least the proportion of metagenomes in the subsample indicated in the cut level
- cnt: number of subsamples where the centroid was found above the cut level



### Pan-resistome

- tool: tool + database
- habitat: type of habitat that the pan-resistome is calculated for
- gene class: gene class of the centroid
- unigenes: number of centroids in the pan-resistome at that epoch (subsample)
- aggregation: new_level centroid
- epoch: the susbample 

## results_tools.rds

Alternatively, all the processed outputs of each tool are present in the R list `results_tools.rds`. Each element of the list contains the output of each tool after argNorm, and manual curation of AROs, parent level and new level. 



