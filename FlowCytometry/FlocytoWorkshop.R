setwd("/Users/rudramanipokhrel/work/Pipelines/Presentation")
# Install Packages
install.packages("BiocManager")

BiocManager::install("flowCore") # read and interpret .fcs files
BiocManager::install("flowViz") #basic visulisation
BiocManager::install("ggcyto") #advanced visulisation using the ggPlot nomenclature
BiocManager::install("openCyto") #Used to link various analysis methodologies
BiocManager::install("flowWorkspace") #used to build anaysis templates
BiocManager::install("CytoML") #imports FlowJo and DiVA workspaces


# Loading the library
library(flowCore)
library(flowAI)
library(ggcyto)
library(flowViz)
library(CytoTree)
library(openCyto)
# Read a single file
file <- "~/work/FlowCytometry/MJ191016 WNV 15 color/CTRL_Tube_001_026.fcs"
fs <- read.FCS(file)
fs

# check the sample
rownames(fs)

# Check the channels and markers
colnames(fs)
markernames(fs)


# Check the expression
head(exprs(fs))

# Visualize with ggcyto
autoplot(fs)
autoplot(fs, x="CD3")
autoplot(fs, x="CD3", y="Time", bins = 256)
autoplot(fs, x="CD3", y="CD4", bins = 256)
autoplot(fs, x="FSC-A", y="Time", bins=256)
autoplot(fs, x="FSC-A", y="SSC-A", bins=256)
# With flowViz
densityplot(~ ., fs)
densityplot(~ `CD3`, fs )
histogram(~ `CD3`, fs)



#Compensation
# Before proceeding with further analysis of the data 
# it is important to properly compensate the data for spectral overlap 
# between fluorescence channels.
comp <- spillover(fs)
head(comp$SPILL) #  splill over matrix
fs_comp <-compensate(fs, comp$SPILL)

# Automatic quality control of flow cytometry data by flowAI
# The flowAI package allows to perform quality control on flow cytometry data
# in order to warrant superior results for both manual and automated 
# downstream analysis.
fs_comp_clean <- flow_auto_qc(fs_comp)

# Transformation
trans.list <- estimateLogicle(fs_comp_clean, colnames(fs_comp_clean[,7:21]))
?estimateLogicle
fs_comp_clean_trans <- transform(fs_comp_clean, trans.list)
# lets visualize again
autoplot(fs, x="CD3")
autoplot(fs_comp_clean_trans, x="CD3")
autoplot(fs, x="CD3", y="CD4", bins = 256)
autoplot(fs_comp_clean_trans, x="CD3", y="CD4", bins = 256)

# Gating and filtering 
autoplot(fs_comp_clean_trans, x="FSC-A", y="SSC-A", bins=256)  
# Gating for lymphocytes
autoplot(fs_comp_clean_trans, x="FSC-A", y="SSC-A", bins=256) +
        geom_vline(xintercept = c(30000, 200000), linetype="dotted", 
                       color = "red", size=1.5) + 
        geom_hline(yintercept = c(0, 100000), linetype="dotted", 
             color = "red", size=1.5) 
                                                                            
  
rg<-rectangleGate("FSC-A"=c(30000, 200000), "SSC-A" = c(0, 100000),
                   filterId = "NoneDebris")
result <- filter(fs_comp_clean_trans, rg)
summary(result)
# Filtering out the debris
fs_nonDebris <- Subset(fs_comp_clean_trans, rg)

autoplot(fs_nonDebris, x="FSC-A", y="SSC-A", bins=256) +
  geom_vline(xintercept = c(30000, 200000), linetype="dotted", 
             color = "red", size=1.5) + 
  geom_hline(yintercept = c(0, 100000), linetype="dotted", 
             color = "red", size=1.5) 

# Filtering out the dublets using channels ("FSC-A", "FSC-H")\
autoplot(fs_nonDebris, x="FSC-A", y="FSC-H", bins=256)

# Imposible to do manually for large data sets
# Automatic gating 
# Load all files
myfiles <- list.files(path = "/Users/rudramanipokhrel/work/FlowCytometry/MJ191016 WNV 15 color/", pattern="*fcs",
                      full.names=TRUE)
fs <- read.flowSet(myfiles)
# sample names
pData(fs)
# number of cells
fsApply(fs, nrow)
# dimensions
dim(exprs(fs[[1]]))
# expression values
exprs(fs[[1]])[1:5, 1:22]
fs
fs[[1]]
names(fs[[1]])
colnames(fs)
colnames(fs[[1]])
head(exprs(fs[[1]]))

#Compensation for spillover 
comp <- fsApply(fs, function(x) spillover(x)[[1]], simplify=FALSE)
fs_comp <- compensate(fs, comp)
fs_comp_clean <- flow_auto_qc(fs_comp)

# Transformation
trans.list <- estimateLogicle(fs_comp_clean[[1]], 
                              channels = colnames(fs)[7:21])
fs_comp_clean_trans <- transform(fs_comp_clean, trans.list)

#Visualise file 
autoplot(fs[[1]]) #before transformation 
autoplot(fs_comp_clean_trans[[1]]) # after transformation.
densityplot(~ `CD3`, fs_comp_clean_trans)
#automatic gating----
# openCyto: https://bioconductor.org/packages/release/bioc/vignettes/openCyto/inst/doc/openCytoVignette.html
#create the empty gating set
auto_gs<-GatingSet(fs_comp_clean_trans)
#cell gate 
fs_data<- gs_pop_get_data(auto_gs)
#noneDebris_gate<- fsApply(fs_data, function(fr) gate_flowclust_2d(fr, xChannel = "FSC-A", yChannel = "SSC-A", K = 3))
noneDebris_gate<- fsApply(fs_data, function(fr) openCyto:::.flowClust.2d(fr, channels= c("FSC-A","SSC-A"), K = 3))
gs_pop_add(auto_gs, noneDebris_gate, parent = "root", name="noneDebris_gate")
recompute(auto_gs)

#autoplot(auto_gs, x="FSC-A", y="SSC-A", "noneDebris_gate", bins=256)
autoplot(auto_gs[[1]], x="FSC-A", y="SSC-A", "noneDebris_gate", bins=256)
# get non debris stat
gs_pop_get_stats(auto_gs, "noneDebris_gate", "percent")
stat <- gs_pop_get_stats(auto_gs, "noneDebris_gate", "percent")
stat <- data.frame(stat)
ggplot(stat, aes(x = sample,y =  percent)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Percentage of non-Debris cells")

#Singlet gate
fs_data <- gs_pop_get_data(auto_gs, "noneDebris_gate") #get parent data
singlet_gate <- fsApply(fs_data, function(fr) openCyto:::.singletGate(fr, channels =c("FSC-A", "FSC-H")))
gs_pop_add(auto_gs, singlet_gate, parent = "noneDebris_gate", name = "singlets")
recompute(auto_gs)
#autoplot(auto_gs, x = 'FSC-A', y = 'FSC-H', "singlets", bins = 256)
autoplot(auto_gs[[1]], x = 'FSC-A', y = 'FSC-H', "singlets", bins = 256)
gs_pop_get_stats(auto_gs, "singlets", "percent")

stat2 <- gs_pop_get_stats(auto_gs, "singlets", "percent")
stat2 <- data.frame(stat2)
ggplot(stat2, aes(x = sample,y =  percent)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Percentage of Singlets")
#Quad gate
markernames(auto_gs)
fs_data <- gs_pop_get_data(auto_gs, "singlets") #get parent data
chnl <- c('BV750-A', 'BV711-A') # CD3 vs CD4 channels
BGquad_gate <- fsApply(fs_data, 
                       function(fr) openCyto:::.quadGate.seq(fr,
                                                             chnl, channels = chnl, gFunc="mindensity"))

gs_pop_add(auto_gs, BGquad_gate, parent = "singlets")
recompute(auto_gs)
gs_get_pop_paths(auto_gs[[1]])
plot(auto_gs)
autoplot(auto_gs[1:4], x = 'BV750-A', y = 'BV711-A', gs_get_pop_paths(auto_gs)[4:7], bins = 256)

stat3 <- gs_pop_get_stats(auto_gs, "BV750-A+BV711-A-", "percent")
stat3 <- data.frame(stat3)
ggplot(stat3, aes(x = sample,y =  percent)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Percentage of CD3")

stat4 <- gs_pop_get_stats(auto_gs, "BV750-A-BV711-A+", "percent")
stat4 <- data.frame(stat4)
ggplot(stat4, aes(x = sample,y =  percent)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 90)) +
  ggtitle("Percentage of CD4")

# Save the filtered and normalized data.
dir.create("fs_filt")
fs_filt <- gs_pop_get_data(auto_gs, y = "singlets")
outDir <- file.path("fs_filt")
write.flowSet(fs_filt, outDir)

# Clustering and differential Analysis
library(diffcyt)
library(CATALYST)
library(SingleCellExperiment)
library(CytoTree)
library(ggplot2)

# read the data
files <- list.files(path = "fs_filt", pattern="*fcs",
                    full.names=TRUE)
fs <- read.flowSet(files, transformation = FALSE, 
                   truncate_max_range = FALSE)
# sample names
pData(fs)
# number of cells
fsApply(fs, nrow)
# dimensions
dim(exprs(fs[[1]]))
# expression values
exprs(fs[[1]])[1:5, 1:22]

#Setup metadata
file_name <- as.character(pData(fs)$name)
sample_id = gsub(".fcs$", "", file_name)
condition = gsub("_Tube.*$", "", file_name)
md <- data.frame( file_name , sample_id, condition,
                  stringsAsFactors = FALSE)
head(md)

# Meta-data: marker information
# column indices of all markers, lineage markers, and functional markers
fcs.data <- runExprsMerge(files, comp = FALSE, transformMethod = "none")
colnames(fcs.data)
col <- sub("<NA>", "", colnames(fcs.data))
col <- sub(".*<", "",col)
col <- sub(">", "",col)
col <- c(col, "Time")

# channel and marker names
channel_name <- colnames(fs)
marker_name <- col

# marker classes
# note: using lineage markers for 'cell type', and functional markers for  'cell state'
cols_lineage <- c(7:10, 12:21)
cols_func <- c(11)
marker_class <- rep("none", ncol(fs[[1]]))
marker_class[cols_lineage] <- "type"
marker_class[cols_func] <- "state"
marker_class <- factor(marker_class, 
                       levels = c("type", "state", "none"))
panel <- data.frame(
  channel_name, antigen = marker_name, marker_class, stringsAsFactors = FALSE
)
panel

# spot check that all panel columns are in the flowSet object
all(panel$channel_name %in% colnames(fs))
# specify levels for conditions & sample IDs to assure desired ordering
md$condition <- factor(md$condition, levels = c("CTRL", "WNV"))
md$sample_id <- factor(md$sample_id, 
                       levels = md$sample_id[order(md$condition)])

# construct SingleCellExperiment

panel_cols = list(channel = "channel_name", class =
                    "marker_class")
md_cols = list(file = "file_name", id = "sample_id",
               factors = c("condition"))
names(md)
sce <- prepData(fs, panel, md, features = panel$channel_name,
                panel_cols = panel_cols,
                md_cols = md_cols,
                cofactor = 150
                )
# cofactor 150 for flowcytometry, 5 for cytof

counts(sce)[1:5,1:5]
##### Diagnostic plots
p <- plotExprs(sce, color_by = "condition")
p$facet$params$ncol <- 3 # Number of columns in plot
p
#cell persamples
n_cells(sce) 
plotCounts(sce, group_by = "sample_id", color_by = "condition")

# Multi-dimensional scaling (MDS) plot
p <- pbMDS(sce, color_by = "condition", label_by = NULL,
           dims = c(1,2))
data <- p$data
ggplot(data, aes(x, y, color= condition, size = n_cells)) +
  geom_point() +
  #geom_text(label=rownames(data))
  labs(title="MDS plot",
       x ="Dim 1", y = "Dim 2")
# Heatmap plot
plotExprHeatmap(sce, scale = "last",
                hm_pal = rev(hcl.colors(100, "Reds")))

plotExprHeatmap(sce, scale = "last",
                row_clust = FALSE,
                hm_pal = rev(hcl.colors(100, "Reds")))

# Marker ranking based on the non-redundancy score
# Markers with higher score explain a larger portion of 
# variability present in a given sample.
plotNRS(sce, features = "type", color_by = "condition")

##### Cell population identification with FlowSOM and  ConsensusClusterPlus
sce <- cluster(sce, features = "type",
               xdim = 10, ydim = 10, maxK = 20, seed = 42)
# view all available clustering
names(cluster_codes(sce))
# access specific clustering resolution
table(cluster_ids(sce, "meta10"))
barplot(table(cluster_ids(sce, "meta10")))
barplot(table(cluster_ids(sce, "meta16")))
# view delta area plot
# The delta area represents the amount of extra cluster stability 
# gained when clustering into k groups as compared to k-1 groups. 
# It can be expected that high stability of clusters can be reached 
# when clustering into the number of groups that best fits the data. 
# The "natural" number of clusters present in the data should thus 
# corresponds to the value of k where there is no longer a 
# considerable increase in stability (pleateau onset).
delta_area(sce)

plotExprHeatmap(sce, features = "type", 
                by = "cluster_id", k = "meta10", 
                row_clust = FALSE,
                bars = TRUE, perc = TRUE)

plotClusterExprs(sce, k = "meta5", features = "type")

# run t-SNE/UMAP on at most 500/1000 cells per sample
set.seed(42)
sce <- runDR(sce, "TSNE", cells = 1000, features = "type")
sce <- runDR(sce, "UMAP", cells = 2500, features = "type")

plotDR(sce, "UMAP", color_by = "CD4")
plotDR(sce, "UMAP", color_by = "meta10")
plotDR(sce, "UMAP", color_by = "meta20")
plotDR(sce, "UMAP", color_by = "meta16")
plotDR(sce, "TSNE", color_by = "CD4")
plotDR(sce, "TSNE", color_by = "meta10")
plotDR(sce, "TSNE", color_by = "meta20")
plotDR(sce, "TSNE", color_by = "meta16")
plotDR(sce, "UMAP", color_by = "CD3")
plotDR(sce, "TSNE", color_by = "CD3")



# Differential analysis
# Differential analysis of cell population abundance compares the 
# proportions of cell types across experimental conditions and aims 
# to highlight populations that are present at different ratios.
plotAbundances(sce, k = "meta16", by = "sample_id")
plotAbundances(sce, k = "meta16", by = "cluster_id")
exprs(sce)
colData(sce)
rowData(sce)
counts(sce)[1:5,1:5]
md
p <- plotAbundances(sce, k = "meta16",
               by = "cluster_id") # shape_by = "sample_id")
#p$facet$params$ncol <- 2
p
# cluster 1
data <- p$data
head(data)
d1 <- data %>% subset(cluster_id==3)
head(d1)
ggplot

# T test 
library(ggpubr)
library(rstatix)
library(tidyverse)
stat.test <- data %>%
  group_by(cluster_id) %>%
  wilcox_test(Freq ~ condition) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance("p.adj")

boxPlot <- function(data, i ) {
  d1 <- data %>% subset(cluster_id==i)
  stat.test <- d1 %>%
    group_by(cluster_id) %>%
    wilcox_test(Freq ~ condition) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  # Create a box plot
  bxp <- ggboxplot(
    d1, x = "condition", y = "Freq", 
    color = "condition", palette = c("blue", "red"),
    add = "dotplot"
    )
    # Add p-values onto the box plots
  stat.test <- stat.test %>%
    add_xy_position(x = "condition", dodge = 0.01)
  # Add 10% spaces between the p-value labels and the plot border
  bxp + stat_pvalue_manual(
    stat.test,  label = "p", tip.length = 0.01
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.5)))
}
par(mfrow=c(2,2)) 
boxPlot(data,4)
boxPlot(data,1)
dev.off()
data


# Differnetial design matrix
design <- createDesignMatrix(
  md, cols_design = c("condition")
)
head(design)
contrast <- createContrast(c(1, 0)) # creating contrast only for second columns
# of design matris

da_res1 <- diffcyt(
  d_input = sce, 
  experiment_info = md, 
  marker_info = panel, 
  design = design, 
  contrast = contrast, 
  analysis_type = "DA", 
  clustering_to_use = "meta16"
)



names(da_res1)
FDR_cutoff = 0.1
table(rowData(da_res1$res)$p_adj < FDR_cutoff)
res_Da_all <- topTable(da_res1, all = TRUE)
# Heatmap plot with normalized frequency
plotDiffHeatmap(sce, rowData(da_res1$res), 
                all = TRUE, fdr = FDR_cutoff)

# Method in tutorial
ei <- metadata(sce)$experiment_info
(da_formula1 <- createFormula(ei, 
                              cols_fixed = "condition", 
                              cols_random = "sample_id"))
contrast <- createContrast(c(0, 1))
da_res1 <- diffcyt(sce, 
                   formula = da_formula1, contrast = contrast,
                   analysis_type = "DA", method_DA = "diffcyt-DA-GLMM",
                   clustering_to_use = "meta16", verbose = FALSE)
names(da_res1)
FDR_cutoff = 0.1
table(rowData(da_res1$res)$p_adj < FDR_cutoff)
res_Da_all <- topTable(da_res1, all = TRUE)
plotDiffHeatmap(sce, rowData(da_res1$res), 
                all = TRUE, fdr = FDR_cutoff)


