# Safe Harbor
The script allows filtering of the given different variants on different parameters provided by the user such as allele frequency, false discovery rate (FDR), nearby tumor repressor genes or oncogenes etc. The main purpose is to narrow down the variant positions which can be then tested and used for gene therapy purpose

------
> **_NOTE:_** The program is based on hg19 version.
------

<br/>

# Contents
- [Prerequisites](#1-prerequisites)
- [Input file format](#2-input-file-format)
- [Running the program](#3-running-the-program)
- [Data Description](#4-data-description)

<br/>

<br/>

## 1) Prerequisites
> **_NOTE:_** If you decide to use conda you have to add **bioconda channel** to your conda. You can do that by: `conda config --add channels bioconda`.

<br/>

Before running the program user has to install all the required dependencies. User can follow either of following steps for that purpose:
- **Using Conda :** User can use ***environment.yml*** using conda to create a new environment with required dependensies. If you have conda installed, following code can be used to create separate environment with the dependencies installed:

    `conda env create -f environment.yml`

After this you can activate the created environment using following command:

    `conda activate safe_harbor`

  If user does not have conda installed. You can download either of it from the link below:
  - [anaconda](https://www.anaconda.com/distribution/)
  - [miniconda](https://docs.conda.io/en/latest/miniconda.html)

- **Using pip :** The user have to install all the dependensies one by one using ***pip***. Here are the commands for installation:
  - `pip install pandas`
  - `pip install loguru` 
  
  **bedtools** can not be installed using **pip**, you can follow the steps [here](https://bedtools.readthedocs.io/en/latest/content/installation.html)
<br/>

<br/>

## 2) Input file format

The input file should be in tab delimited format. The minimum requirement for the input file is the ***id*** and ***position***. The picture below shows example format of the data:

![input_data_example](https://user-images.githubusercontent.com/22225447/111481096-9f761d80-8700-11eb-9f3a-86a2d2bb1dea.png)

<br/>

> **_NOTE:_** The column name has to be same. The names can be upper or lower case.

<br/>

<br/>

## 3) Running the program
The simplest way to run the program is `python safe_harbor.py -i input_filename`, assuming the input file has at least `id and position` column. In this case, the program assumes that the input file has only `id and position` column, and even if there are other columns such as `FDR, AF and eQTL`, since the user did not pass the required parameters while running the program those columns will be ignored. But for other required files such as `TAD domain information, repressive region, chromatin interaction` default files will be used. To see all the parameters available you can run `python safe_harbor.py -h`, which will give following details:
```
usage: safe_harbor.py [-h] [-i INPUT] [-t THRESH] [-eqtl EQTL_GENES]
                      [-tad TAD_DOMAIN] [-rr REPRESSIVE_REGION]
                      [-ar ACTIVE_REGION] [-gd GENE_DENSITY] [-o OUTPUT]
                      [-af ALLELE_FREQ] [-hic HIC_INTERACTION]
                      [-l NEARBY_CANCER_GENES] [-fname FILE_NAME]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input file with mobile element id
  -t THRESH, --thresh THRESH
                        fdr threshold to filter out MEIs
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
  -af ALLELE_FREQ, --allele_freq ALLELE_FREQ
                        threshold allele frequency of the snp to filter
  -hic HIC_INTERACTION, --hic_interaction HIC_INTERACTION
                        hic-promoter interaction bed file
  -l NEARBY_CANCER_GENES, --nearby_cancer_genes NEARBY_CANCER_GENES
                        any MEI with oncogenes or tumor repressor genes 50kb
                        upstream or downstream will be removed
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

- **allele_frequency** : this parameter is used to filter the column AF, keeping all those variants with AF > threshold. If the input file doesn't have this column, no need to use this parameter, otherwise it will give an error. But if the input file has AF, but this parameter is not passed, then AF column will be ignored and would not be considered for filtration steps.
<br/>

- **hic_interaction** : this parameter takes the chromatin interaction file. The file should have 4 columns: `chr, start, stop and gene`. This file is used to check if any of the variant involved interaction involve tumor repressor or oncogenes and dosage sensitive genes.
<br/>

- **nearby_cancer_genes** : this parameter takes the value in base pair. The value represents how many base pairs upstream and downstream, the user wants to look for tumor repressor or oncogenes. The default value is `50000 (50 kb)`.
<br/>

- **file_name** : this parameter takes output filename for result. Default file name is `result.csv`.Two main output files will be generated. One will contain all the data with raw information, another will contain just the filtered result with `id, position and UCSC link.`

> **_NOTE:_** all other intermediate files and log file is stored in folder `temp_files` inside the output folder.
<br/>

------
> **_NOTE:_** if the user wants to use default files provided, please make sure the `blood_repressive_active_region.tar.gz and blood_hic_interaction.bed.tar.gz` are unzipped and in the data folder.
------
<br/>

<br/>

## 4) Data Description
- ***merged_gm12878.bed :*** This is the TAD domain file. It is downloaded from [Rao et al, Cell, 2014](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525)

- ***oncogenes_and_tumor_repressor_genes.bed :*** This is the list of tumor repressors and oncogenes downloaded from [tumor repressor link](https://bioinfo.uth.edu/TSGene/?csrt=14252334967482590952) and [oncogenes link](http://ongene.bioinfo-minzhao.org)

- ***blood_hic_interaction.bed :*** This file list chromatin interactions involving different blood cells. It is downloaded from (Javierre et al, Cell, 2016)[https://pubmed.ncbi.nlm.nih.gov/27863249/]

- ***dosage_sensitive_genes.txt :*** This is the list of dosage sensitive genes. This is downloaded from (Exome Aggregation Consortium, Nature, 2016)[https://www.nature.com/articles/nature19057]

- ***sorted_gene_annotation.bed:*** This is genes list with their start and stop position. It is downloaded using Biomart hg19 version.

- ***blood_repressive_marks :*** This contains the chromatin regions that are assigned to be repressive in nature based on ChromHMM. The files associated with **blood** are filtered using column **ANATOMY** from [Roadmap, 2013](https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15), and then regions labelled as `Het, ReprPC, ReprPCWk, Quies` were extracted.

- ***blood_active_transcription_marks :*** This contains the chromatin regions that are assigned to be active in nature based on ChromHMM. The files associated with **blood** are filtered using column **ANATOMY** from [Roadmap, 2013](https://docs.google.com/spreadsheets/d/1yikGx4MsO9Ei36b64yOy9Vb6oPC5IBGlFbYEt-N6gOM/edit#gid=15), and then regions labelled as `TssA, TssAFlnk, Tx, TxWk` were extracted.

- ***genes_tad.bed :*** This file provides the TAD domain information for the genes in **sorted_gene_annotation.bed**. The file is generated using [bedtools intersect](https://bedtools.readthedocs.io/en/latest/content/tools/intersect.html).

- ***gene_density_all_tad.csv :*** This is the pre-calculated gene density information for the TAD domains in **merged_gm12878.bed**. This file is generated using **genes_tad.bed**.
