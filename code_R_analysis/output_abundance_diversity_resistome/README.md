
# File descriptions

## conversion_ARO_parent_new_level.csv 
Alternatively, `conversion_ARO_parent_new_level.rds`.

| Term_ID    | Term_Label | Parent_ID | Parent_Label | new_label  |
| -------- | ------- | -------- | ------- | -------- |
| ARO assigned by RGI or argNorm  | ARO's description    | Higher onthology | Description | Manual class assignment |


## abundance_diversity.csv 
Alternatively, `abundance_diversity.rds`.

The abundance scaled, raw, raw_unique, normed10m, are retrieved form GMGC.
The table reports the sums of the four abundance metrics and counts of the unique unigenes identifiers (alpha diversity) grouping the data by sample and "gene", at the aggregation level (ARO, parent description or new_level). 


| sample    | gene | aggregation | scaled | raw  | raw_unique  | normed10m  | unigenes  | habitat  | habitat2 
| -------- | ------- | -------- | :-------: | :--------: | :--------: | :--------: | -------- | -------- | -------- 
| metagenomic sample ID  | AR0, parent description, or new_label   | By ARO, or parent, or new_level | scaled  | raw | raw_unique | normed10m | aplpha diversity (number of unigenes) | 15 different habitats | source of habitats (human, mammals, etc.)


## core_resistome.csv and pan_resistome.csv
Alternatively, `core_resistome.rds` and `pan_resistome.rds`.

The pan- and core-resistomes were calculated for each environment by randomly taking 500 subsamples of size 100 metagenomic samples from each habitat. For each subsample, we stored (1) the unique number of ARGs (unigene ID) present in any sample (alpha diversity) grouped by tool, habitat, and aggregation level (ARO, parent description, new_level) for the pan-resistome calculations, and (2) the ARGs present in at least a proportion of the samples determined by the "cut" level (commonly encountered gene set) grouped by tool, habitat, and aggregation level (new_level). 

The size of the pan-resistome for each environment will be calculated as the average alpha diversity over the 500 subsamples, determined by "epoch", grouped by tool, habitat, and aggregation level, and the core-resistome as the number of ARGs present in at least an arbitraty number determined by "cnt", e.g. 450, of the commonly encountered gene sets grouped by tool, habitat, and aggregation level.


### Core-resistome

| X    | new_level | tool | habitat | cut  | cnt
| -------- | ------- | -------- | ------- | -------- | --------
| unigene ID  | assigned manual gene class  | tool used to detect ARGs | habitat | cut level for the core resistome | how many times the unigene cleared the cut level


### Pan-resistome

| tool | habitat | gene_class  | unigenes  | aggregation  | epoch 
| -------- | ------- | -------- | ------- | -------- | --------
| tool used to detect ARGs | habitat | AR0, parent description, or new_label | number of unique unigenes id | By ARO, or parent, or new_level | subsample number 

## processed_tool.X.csv
Alternatively, all the processed outputs of each tool are present in the R list `results_tools.rds`.

Each file contains the output of each tool after argNorm, and manual curation of AROs, parent level and new level. For fARGene outputs, even the identity levels and AROs after running the predicted genes with RGI/DIAMOND are included. 

Some key columns are: 
- query: unigene id
- id: identity level by each tool (for fARGene is the output of RGI on the predicted genes by fARGene)
- ARO: ARO id after argNorm or manual curation
- tool: DeepARG (nt), DeepARG (a.a.), RGI (BLAST  nt), RGI (DIAMOND  nt), RGI (DIAMOND  a.a.), fARGene (nt), fARGene (a.a.), AMRFinderPlus (nt), AMRFinderPlus (a.a.), ABRicate (ARG-ANNOT - nt), ABRicate (CARD - nt), ABRicate (MEGARES - nt), ABRicate (NCBI - nt), ABRicate (ResFinder - nt), ResFinder (nt)
- parent and parent_description: assigned higher ontology (ARO id and description)
- new_level: manual converstion to gene class

