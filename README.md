# Genomics and Epigenetic Guided Safe Harbor mapper (GEG-SH mapper)
This program is designed to identify candidate Genomic Safe Harbor (GSH) sites for safe genome editing. It integrates  polymorphic mobile element insertions (pMEIs) with epigenomic signatures and 3D chromatin organization information.

------
> **_NOTE:_** The program is based on hg19 version.
------

<br/>

# Contents
- [Prerequisites](#1-prerequisites)
- [Input file format](#2-input-file-format)
- [Running the program](#3-running-the-program)
- [Output file description](#4-output-file-description)
- [Data Description](#5-data-description)

<br/>

<br/>

## 1) Prerequisites
> **_NOTE:_** After installation of conda you have to add **bioconda channel** to your conda. You can do that by: `conda config --add channels bioconda`.

<br/>

All the required dependencies for this program can be loaded by the following steps

    `conda env create -f environment.yml`
    
After this you can activate the created environment using following command:

    `conda activate safe_harbor`

   If user does not have conda installed. You can download either of it from the link below:
   - [anaconda](https://www.anaconda.com/distribution/)
   - [miniconda](https://docs.conda.io/en/latest/miniconda.html)

<br/>

<br/>

## 2) Input file format

The input file of pMEIs should be in tab delimited format with the following column name: “**id, position, fdr, af**”. The minimum requirements for the input file is the **id** and **position**. The picture below shows example format of the data:

![input_data_example](https://user-images.githubusercontent.com/22225447/133011347-809044d0-dfa7-40a6-8b13-4ab7403741b3.png)

<br/>

> **_NOTE:_** The column name has to be same. The names can be upper or lower case.

<br/>

<br/>

## 3) Running the program
------
> **_NOTE:_** Before running the program, make sure all files are extracted. You can do that by running the script **extract.sh** by running following command in terminal in mac or linux system: `sh ./extract.sh` 
------

The simplest way to run the program is `python safe_harbor.py -i input_filename`, assuming the input file has at least `id and position` column. In this case, the program assumes that the input file has only `id and position` column, and even if there are other columns such as `FDR, AF and eQTL`, since the user did not pass the required parameters while running the program those columns will be ignored. But for other required files such as `TAD domain information, repressive region, chromatin interaction` default files will be used. To see all the parameters available you can run `python safe_harbor.py -h`, which will give following details:
```
usage: safe_harbor.py [-h] [-i INPUT] [-t THRESH] [-eqtl EQTL_GENES]
                      [-tad TAD_DOMAIN] [-rr REPRESSIVE_REGION]
                      [-ar ACTIVE_REGION] [-gd GENE_DENSITY] [-o OUTPUT]
                      [-br BLACKLIST_REGION] [-af ALLELE_FREQ]
                      [-hic HIC_INTERACTION] [-l NEARBY_CANCER_GENES]
                      [-fname FILE_NAME]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input file with mobile element id
  -t THRESH, --thresh THRESH
                        fdr threshold to filter out variants
  -eqtl EQTL_GENES, --eqtl_genes EQTL_GENES
                        eQTL genes
  -tad TAD_DOMAIN, --tad_domain TAD_DOMAIN
                        custom tad domain file, else will use default file
                        provided
  -rr REPRESSIVE_REGION, --repressive_region REPRESSIVE_REGION
                        bed file containing regions with repressive mark
  -ar ACTIVE_REGION, --active_region ACTIVE_REGION
                        bed file containing regions with active transcription
                        mark
  -gd GENE_DENSITY, --gene_density GENE_DENSITY
                        gene density in the tad domains, mean gene density
                        will be used as default
  -o OUTPUT, --output OUTPUT
                        ouput folder name
  -br BLACKLIST_REGION, --blacklist_region BLACKLIST_REGION
                        bed file with the coordinates that the user does not
                        want to include in output
  -af ALLELE_FREQ, --allele_freq ALLELE_FREQ
                        allele frequency threshold for the variant
  -hic HIC_INTERACTION, --hic_interaction HIC_INTERACTION
                        chromatin interaction bed file
  -l NEARBY_CANCER_GENES, --nearby_cancer_genes NEARBY_CANCER_GENES
                        default = 0, takes number as input representing the
                        distance user wants to check for oncogenes or tumor
                        repressor genes in upstream or downstream of the pMEI
  -fname FILE_NAME, --file_name FILE_NAME
                        output file name
```
Here is the description of the different parameters:
- **input** : this parameter takes the input file. This is the only required parameter without which program will not run. The input file has to be in tab delimited format with at least two columns `id and position`.
<br/>

- **thresh** : this parameter is used to filter the column FDR, keeping all those variants with FDR > threshold. If the input file doesn't have this column, no need to use this parameter, otherwise it will give an error. But if the input file has FDR, but this parameter is not passed, then FDR column will be ignored and would not be considered for filtration steps. If there are multiple FDR values, each value should be separated by `,`.
<br/>

- **eqtl_genes** : this parameter takes boolean value as `True or False`. If your input file has eQTL genes associated with MEIs, then user need to pass the parameter as True, by default it is set to ***False***. If there are multiple eQTL genes, each gene should be separated by `,`. The eQTL genes names are expected to be in **gene symbol format**.
<br/>

- **tad_domain** : this parameter takes the TAD domain information file in bed file format. The default value is set to `None`, so if user does not pass any file, it will look for the provided file in data folder (***merged_gm12878.bed***). This file is used for several purpose, to assign TAD domain information to the variants and genes, calculate gene density, and check the common TAD domain between variant and tumor repressor or oncogenes
<br/>

- **repressive_region** : this parameter takes repressive chromatin region information in bed file format. This file is used to remove the variants that overlap with any repressive chromatin region. If user does not use this parameter, default file (**blood_repressive_marks.bed**) will be used.
<br/>

- **active_region** : this parameter takes active chromatin region information in bed file format. This file is only used to tag the variants that overlap with any active chromatin region. In absence of this parameter the tagging step is ignored. If user does not use this parameter, default file (**blood_active_transcription_marks.bed**) will be used.
<br/>

- **gene_density** : this parameter takes the gene density value that the user wants to allow in a TAD domain. If the user does not pass any value mean gene density of all TAD domains are used. Gene density for each TAD domain is calculated by: `(total gene in particular TAD domain/ length of TAD domain) * 1000000` as `genes per million TAD`
<br/>

- **output** : this parameter takes output folder name. **results** is used as default folder name.
<br/>

- **blacklist_region**: this parameter takes the bed file with coordinates that the user want to exclude from the result.
<br/>

- **allele_freq** : this parameter takes the threshold value used for filtering based on allele frequency of certain variant. The value provided will be considered as lower limit threshold and upper limit threshold will be calculated based on given value as `1-af`, keeping all those variants with AF > threshold and AF < 1- threshold . If the input file doesn't have this column, no need to use this parameter, otherwise it will give an error. But if the input file has AF, but this parameter is not passed, then AF column will be ignored and would not be considered for filtration steps.
<br/>

- **hic_interaction** : this parameter takes the chromatin interaction file. The file should have 4 columns: `chr, start, stop and gene`. This file is used to check if any of the variant involved interaction involve tumor repressor or oncogenes and dosage sensitive genes.
<br/>

- **nearby_cancer_genes** : this parameter takes the value in base pair. The value represents how many base pairs upstream and downstream, the user wants to look for tumor repressor or oncogenes. The default value is `0 (0 kb)`.
<br/>

- **file_name** : this parameter takes output filename for result. Default file name is `result.csv`.Two main output files will be generated. One will contain all the data with raw information, another will contain just the filtered result with `id, position and UCSC link.`

> **_NOTE:_** all other intermediate files and log file is stored in folder `temp_files` inside the output folder.
<br/>

------
> **_NOTE:_** if the user wants to use default files provided, please make sure the `blood_repressive_active_region.tar.gz and blood_hic_interaction.bed.tar.gz` are unzipped and in the data folder.
------
<br/>

<br/>

## 4) Output file Description

|Column name|	Description|
|-----------|--------------|
|id|	unique variant id|
|position|	chromosome coordinates for the variant|
|fdr|	false discovery rate|
|af|	population allele frequency|
|chr|	chromosome number of variant|
|start|	chromosome start position|
|stop|	chromosome end position|
|extended_start|	(this will only be in output if linear distance is used for filtration) upstream extension of start, based on the extension parameter, by default it is extended by x kb|
|extended_stop|	(this will only be in output if linear distance is used for filtration) downstream extension of stop, based on the extension parameter, by default it is extended by x kb|
|tad_name|	tad domain associated with the variant|
|same_cancer_tad|	represents boolean value, True if variant and cancer genes or tumor repressor genes are in same TAD, else represented as False|
|gene_density (x, genes per million tad)|	number of genes per mb tad|
|hic_interacted_genes|	genes interacting with the variant associated locus|
|common_tad_count|	common TAD between the variant and the interacted genes|
|dosage_sensitive interactions|	number of interactions with dosage sentive genes|
|hic_interacted_genes (oncogenic or tumor repressor)|	boolean value representing if the interaction with genes comprises any oncogenic or tumor repressor genes|
|repressive_region|	boolean value representing if the variant overlaps with the repressive region or not|
|repressive_region_info| provides chromatin state information for repressive region (Het, ReprPC, ReprPCWk, Quies)
|nearby_cancer_genes (x kb)|	boolean value representing if there is present of oncogenic or tumor repressor genes based on the user provided distance, default value is 50 kb||
|nearby_cancer_genes_names|	(this will only be in output if linear distance is used for filtration) oncogenic or tumor repressor genes with in user provided distance|
|active region| boolean value representing if the variant overlaps with the active chromatin region or not|
|active region_info| provides chromatin state information for sctive chromatin region (TssA, TssAFlnk, Tx, TxWk)
|fdr_test| boolean value representing if FDR value is greater than provided threshold
|passed_all_filter|	represents boolean value, True if the variant passes all the filters used, else represented as False||
|ucsc_link|	UCSC browser link to the variant position|

<br/>

<br/>

## 5) Data Description
- ***merged_gm12878.bed :*** This represents the TAD domain file for blood. It is downloaded from [Rao et al, Cell, 2014](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525). **GSE63525_GM12878_primary+replicate_Arrowhead_domainlist.txt** is further processed to generate bed file with **chr, start, stop and length** of TAD domain. ***bedtools merge*** is then applied to the generated bed file. Bed file should be sorted before applying bedtools merge.

- **brain_hc_tad_100kb_domains.bed :** This represents the TAD domain file for brain. Raw Hi-C data for hippocampus (SRA: SRR4094699) was downloaded from [GSE86189](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86189). It was then processed locally using [Hic-pro](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0831-x) (version 2.11.1) and [HiCExplorer](https://hicexplorer.readthedocs.io/en/latest/content/list-of-tools.html) (version 3.6).

- ***oncogenes_and_tumor_suppressor_genes.bed :*** This is the list of [tumor suppressors genes](https://bioinfo.uth.edu/TSGene/?csrt=14252334967482590952) and [oncogenes](http://ongene.bioinfo-minzhao.org) downloaded from [Zhao et al, Nucleic Acids Research, 2015](https://academic.oup.com/nar/article/44/D1/D1023/2503080) and [Liu et al, Journal of Genetics and Genomics, 2017](https://www.sciencedirect.com/science/article/pii/S1673852716302053?via%3Dihub) respectively.

- ***blood_hic_interaction.bed :*** This file list chromatin interactions involving different blood cells. It is downloaded from [Javierre et al, Cell, 2016](https://pubmed.ncbi.nlm.nih.gov/27863249/)

- ***brain_hic_interaction.bed :*** This file list chromatin interactions involving different brain cells. It is downloaded from [Jung et al, Nature Genetics, 2019](https://www.nature.com/articles/s41588-019-0494-8)

- ***dosage_sensitive_genes.bed :*** This is the list of dosage sensitive genes. This is downloaded from (Exome Aggregation Consortium, Nature, 2016)[https://www.nature.com/articles/nature19057]. The genes were selected based on pLI (probability loss of function intolerant) value greater than 0.9.

- ***sorted_gene_annotation.bed:*** This is genes list with their start and stop position. It is downloaded using Biomart hg19 version.

- ***blood_repressive_marks_state.bed :*** This contains the chromatin regions that are assigned to be repressive/quiescent/heterochromatin in nature based on ChromHMM. The files associated with **blood** are filtered using column **ANATOMY** from [Roadmap, 2013](https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15), and then regions labelled as `Het, ReprPC, ReprPCWk, Quies` were extracted.

- ***blood_active_transcription_marks_state.bed :*** This contains the chromatin regions that are assigned to be active in nature based on ChromHMM. The files associated with **blood** are filtered using column **ANATOMY** from [Roadmap, 2013](https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15), and then regions labelled as `TssA, TssAFlnk, Tx, TxWk` were extracted.

- ***brain_repressive_marks_state.bed :*** This contains the chromatin regions that are assigned to be repressivequiescent/heterochromatin in nature based on ChromHMM. The files associated with **brain** are filtered using column **ANATOMY** from [Roadmap, 2013](https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15), and then regions labelled as `Het, ReprPC, ReprPCWk, Quies` were extracted.

- ***brain_active_transcription_marks_state.bed :*** This contains the chromatin regions that are assigned to be active in nature based on ChromHMM. The files associated with **brain** are filtered using column **ANATOMY** from [Roadmap, 2013](https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15), and then regions labelled as `TssA, TssAFlnk, Tx, TxWk` were extracted.

- ***genes_tad.bed :*** This file provides the TAD domain information for the genes in **sorted_gene_annotation.bed**. The file is generated using [bedtools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) on the TAD domain file downloaded from ***Rao et al, Cell, 2014*** and gene annotation file from ***Ensembl Biomart hg19***. If the user provides own tad domain file, new file will be generated, else this file will be used as default.

- ***gene_density_all_tad.csv :*** This is the pre-calculated gene density information for the TAD domains in **merged_gm12878.bed**. This file is generated using **genes_tad.bed**.
