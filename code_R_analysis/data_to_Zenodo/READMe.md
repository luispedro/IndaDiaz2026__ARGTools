
# File descriptions

## Software

### Tools and databases 

The tools and databases used to detect ARGs are listed below.

| **Tool** | **Availability** | 
| :---: | :---: | 
| fARGene (v0.1) | https://github.com/fannyhb/fargene | 
| DeepARG (v2) | https://github.com/gaarangoa/deeparg | 
| AMRFinderPlus (v4.0.15), database 2024-12-18.1 | https://github.com/ncbi/amr |
| RGI (v6.0.3), database CARD (v4.0.0) | https://github.com/arpcard/rgi | 
| ResFinder (v2.4.0) | https://github.com/cadms/resfinder | 
| ABRicate v1.0.1, databases: ARG-ANNOT, CARD, MEGARes v2.0, ResFinder, and NCBI (all updated 2025-01-14) | https://github.com/tseemann/ABRICATE |  

### Normalization

The outputs of DeepARG, AMRFinderPlus, ABRicate, and ResFinder were processed with [argNorm v1.0.0](https://github.com/BigDataBiology/argNorm).

## Datasets

### Output from ARG detection pipelines

Each file is located in the folder `pipelines_output`. The result of each individual pipeline, after processing them with argNorm (except for fARGene), and adding their manually curated class. We include all unigenes from the Global Microbial Gene Catalog v1.0 ([GMGC](https://gmgc.embl.de/)) reported as ARG by the tools. 

For each pipeline, we add the columns from argNorm (we direct the reader to argNorm):

- ARO
- ARO_name
- Cut_Off
- confers_resistance_to
- confers_resistance_to_names
- resistance_to_drug_classes
- resistance_to_drug_classes_names

And from the manual curation of gene classes:

- pipeline: internal pipeline name
- id: standardize name for the identity level column
- parent: first level of ARO aggregation
- parent_description: description of the first level of aggregation
- geneclass_argcompare: gene class - final level of aggregation, should perhaps rename it, in the abundance file it’s called gene class

Pipeline abbreviations:

- DeepARG: DeepARG tool in nucleotide format
- DeepARG-aa: DeepARG tool in amino acid format
- fARGene: fARGene tool in nucleotide format
- fARGene-aa: fARGene tool in amino acid format
- RGI-DIAMOND: RGI tool in nucleotide format using DIAMOND aligner
- RGI-DIAMOND-aa: RGI tool in amino acid format using DIAMOND aligner
- RGI-BLAST: RGI tool in nucleotide format using BLAST aligner
- AMRFinderPlus: AMRFinderPlus tool in amino acid format
- AMRFinderPlus-nt: AMRFinderPlus tool in nucleotide format
- ResFinder: ResFinder tool in nucleotide format
- ABRicate-ARGANNOT: ABRicate tool using ARGANNOT reference database
- ABRicate-CARD: ABRicate tool using CARD reference database
- ABRicate-MEGARes: ABRicate tool using MEGARes reference database
- ABRicate-NCBI: ABRicate tool using NCBI reference database
- ABRicate-ResFinder: ABRicate tool using ResFinder reference database

Note: the pipelines DeepARG70, DeepARG80, DeepARG90, RGI-DIAMOND70, RGI-DIAMOND80, and RGI-DIAMOND90 were added to the abundance file to calculate aggregated abundance and richness by gene class, first filtering genes at 70%, 80%, and 90% identity. 

The list of unigenes per habitat reported as ARG by any pipeline is found in the file `reported_unigenes_as_ARG_per_habitat.csv`. 

### Conversion and aggregation of AROs to gene classes

The mapping of the AROs after argNorm to their parent level and the gene classes used here is found in the file `conversion_aro_geneclass.csv`. AROs were manually added to fARGene.

The columns include: 

- Term_ID: ARO from argNorm
- Term_Label: description of the ARO term
- Parent_ID: RO of the first higher aggregation
- Parent_Label: description of the Parent_ID
- geneclass: manually curated gene class

### Abundance and richness

The abundance and richness of each gene (only from the habitats we are interested in) are found in the file `abundance_richness.csv.gz`. The columns include:

- sample: metagenomic sample ID from the Global Microbial Gene Catalog v1.0 ([GMGC](https://gmgc.embl.de/)).
- geneclass: manually curated gene class
- pipeline: tool and database identifier
- abundance: normalized read by gene size, scaled to 1 million reads per metagenome aggregated at the gene class level
- richness: number of unique genes detected in the metagenome (raw count >0) after rarefaction aggregated at the gene class level
- richness_no_rarified: number of unique genes detected in the metagenome (raw count >0) aggregated at the gene class level

### Abundance per gene 

Original information from ([GMGC](https://gmgc.embl.de/)) after filtering for those unigenes identified as ARGs found in file `args_abundances.tsv.gz`. The file has the same structure as GMGC. Note that the abundance here is scaled to 10 million per metagenome. 


## Metadata of the metagenomic samples

The metadata of the metagenomes is found in the file `metagenomes_metadata.csv`. The file contains the following information:
- sample_id: metagenome identifier
- habitat
- insertsHQ: high-quality reads



