# Pipelines

 ## A. Bulk-RNAseq 
This is bulk RNAseq data processing pipe lines.
#### For data preprocessing (raw fastq files and outputs count matrix):
I have desinged the workflow using simple bash script as well as using Snakemake.
- For bash script follow the folder: ***bash_script***
- For Snakemake follow the folder: ***snakemake***
#### For Differential Analysis:
The corresponding folder is ***DESeq2***.\
This folder implements the workflow for  differential gene expression analysis using DESeq2 bioconductor package. Follow the instructions in ***Bulk RNAseq Analysis.pdf*** file. The file ***DESeq2Design.pdf***  contains few examples to design and get the desired results from DESeq2 (dds) object.

#### PathwaysAnalysis
This folder contains a pdf file named ***Functional Pathways Analysis.pdf***, which includes a brief introducion to the package ***clusterProfiler*** to do GO and KEEG pathways analysis.

## B. scRNAseq
In this folder there are two files:
- **sc-RNAseq cellrnager.pdf** is **cellranger** workflow for demultiplexing and getting count matrix from scRNAseq BCL files.
- **scRNA_data_analysis_complete_workflow.pdf** is complete Seurat workflow for scRNAseq data analysis starting from cellranger output.

## C. Flow Cytometry Data
The pdf file **Flow Cytometry data analysis workflow** describes the complete workflow for Flow Cytometry data.


## D. Inferential statistics
In the stat folder there is jupyter notebook **stat.ipynp** which includes detailed explanation of z-test, t-test, ANOVA, correlation, regression and chi-square test with examples incorpotated from Udacity's free statistical inference course


## E. nsolver_mRNA
It is a Python Package  to QC, normalize and analyze nano string nsolver mRNA data.
Look at the jupyter notebook (main.ipynb) in this folder for more detailed information.
