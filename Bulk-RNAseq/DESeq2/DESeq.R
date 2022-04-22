# Importing Libraries ----
library(ggplot2)
library(DESeq2)
library(tidyverse)
library(ape)
library(pheatmap)
library(RColorBrewer)
library(MetaboAnalystR)
rm(list=ls())
# Setting directories and work and variables ----
setwd("/xdisk/rpokhrel/Test/R/bulkRNA/Deseq2")
count_dir <- "data/gene_count.csv"

# Remember it should be matched with columns of count data.
# If you have more complex study design you can input it here
sample_dir <- "data/sample.csv"
# design of you study
#examples:  design <- ~ 0 + variable1 + variable2 + variable3 + variable1:variable2+ variable1:variable3 + variable2:variable3
design <- ~0 + class
# In factor put the variables appeared in design
factor <- c("class")
# gene filtering criteria
minimum_count <- 5 # in a sample
at_least_in_samples <- 2  #  recommend put 1 less than number of replicates

# Sample depth check
threshold <- 5000000
# figure dimensions
width <- 5
height <- 5
# main variable to plots from sample
## import utility functions
source("utils.R")
# Starts DEseq2 ----
#dev.off()
call_DESeq2()





