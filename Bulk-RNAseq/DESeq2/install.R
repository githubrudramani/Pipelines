
r = getOption("repos")
r["CRAN"] = "http://cran.us.r-project.org"
options(repos = r)
# install Bioconductor

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(version = "3.14")
# install statistical packages
BiocManager::install(c("edgeR", "limma", "DESeq2"))
# install required tools
install.packages("tidyverse")
install.packages("pheatmap")
install.packages("ape")
# for metaboanalyst
BiocManager::install("mzR")
BiocManager::install("MSnbase")
install.packages("remotes")
remotes::install_github("xia-lab/MetaboAnalystR")


