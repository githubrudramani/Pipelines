# Pipelines

 ## A. Bulk-RNAseq 
This is bulk RNAseq data processing pipe lines.
#### For data preprocessing (raw fastq files and outputs count matrix):
I have desinged the workflow using simple bash script as well as using Snakemake.
- For bash script follow the folder: ***bash_script***
- For Snakemake follow the folder: ***snakemake***
#### For Differential Analysis:
The corresponding folder is ***DESeq2*** 
This folder implements the workflow for  differential gene expression analysis using DESeq2 bioconductor package. Follow the instructions in ***Bulk RNAseq Analysis.pdf*** file. It also contains few examples in file ***DESeq2Design.pdf*** to get the desired results form DESeq2 (dds) object.

#### PathwaysAnalysis
This folder contains a pdf file named ***Functional Pathways Analysis.pdf***, which includes a brief introducion to the package ***clusterProfiler*** to do GO and KEEG pathways analysis.
