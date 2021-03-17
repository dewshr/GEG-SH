# Safe Harbor
The script allows filtering of the given Mobile Element Insertion sites (MEIs) based on different parameters provided by the user such as allele frequency, false discovery rate (FDR), nearby tumor repressor genes or oncogenes etc.

<br/>

# Contents
- [Prerequisites](#1-prerequisites)
- [Input file format](#2-input-file-format)
- [Running the program](#3-running-the-program)

<br/>

<br/>

## 1) Prerequisites
The dependencies required for the program are listed in ***environment.yml*** file. The dependencies can be installed individually or if you have conda installed, following code can be used to create separate environment with the dependencies installed:
`conda env create -f environment.yml`

<br/>

<br/>

## 2) Input file format

The input file should be in tab delimited format. The minimum requirement for the input file is the ***id*** and ***position***. The picture below shows example format of the data:

![input_data_example](https://user-images.githubusercontent.com/22225447/111481096-9f761d80-8700-11eb-9f3a-86a2d2bb1dea.png)

<br/>

> **_NOTE_** The column name has to be in same exact format but no particular order.

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

- **thresh** : this parameter is used to filter the column FDR. If the input file doesn't have this column, no need to use this parameter, otherwise it will give an error. But if the input file has FDR, but this parameter is not passed, then FDR column will be ignored and would not be considered for filtration steps. If there are multiple FDR values, each value should be separated by `,`.
<br/>

- **eqtl_genes** : this parameter takes boolean value as `True or False`. If your input file has eQTL genes associated with MEIs, then user need to pass the parameter as True, by default it is set to ***False***. If there are multiple eQTL genes, each gene should be separated by `,`.
<br/>

- **tad_domain** : this parameter takes the TAD domain information file in bed file format. The default value is set to `None`, so if user does not pass any file, it will look for the provided file in data folder (***merged_gm12878.bed***)
<br/>

- **repressive_region** : this parameter takes repressive chromatin region information in bed file format. If user does not use this parameter, default file (**blood_repressive_marks.bed**) will be used.
<br/>

- **active_region** : this parameter takes active chromatin region information in bed file format. This file is only used to tag, if the MEIs overlap with any active chromatin region or not. In absence of this parameter the tagging step is ignored.
<br/>

- **gene_density** : this parameter takes the gene density value that the user wants to allow in a tad domain. If the user does not pass any value mean gene density of all TAD domains are used. Gene density for each TAD domain is calculated by: `(total gene in particular TAD domain/ length of TAD domain) * 1000000` as `genes per million TAD`
<br/>

- **output** : this parameter takes output folder name. **results** is used as default folder name.
<br/>

- **allele_frequency** : this parameter is used to filter the column AF. If the input file doesn't have this column, no need to use this parameter, otherwise it will give an error. But if the input file has AF, but this parameter is not passed, then AF column will be ignored and would not be considered for filtration steps.
<br/>

- **hic_interaction** : this parameter takes the chromatin interaction file. The file should have 4 columns: `chr, start, stop and gene`.
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
