library(flowCore)
library(flowAI)
library(ggcyto)
library(openCyto)

#Load many fcs files into a flow set using flowCore ----
myfiles <- list.files(path = "/Users/rudramanipokhrel/work/FlowCytometry/MJ191016 WNV 15 color/", pattern="*fcs",
                      full.names=TRUE)
fs <- read.flowSet(myfiles)
exprs(fs[[1]])[1:5, 12:15]
fs1 <- read.flowSet(myfiles, transformation=FALSE)
exprs(fs1[[1]])[1:5, 12:15]

# getting names/columns of the file
names(fs[[1]])

# get gene marker name
markernames(fs)
colnames(fs)
# getting expression of the file
exprs(fs[[1]])

# getting median of each columns
each_col(fs[[1]], median)

# gettting file name
keyword(fs[[1]])$FILENAME

#Compensation for spillover ----
# Before proceeding with further analysis of the data in a flowFrame or flowSet,
# it is important to properly compensate the data for spectral overlap between fluorescence channels.
comp <- fsApply(fs, function(x) spillover(x)[[1]], simplify=FALSE)
fs_comp <- compensate(fs, comp)
fs_comp_clean <- flow_auto_qc(fs_comp)
colnames(fs_comp_clean)


trans <- estimateLogicle(fs_comp_clean[[1]],  channels = colnames(fs)[7:21])
fs_comp_clean_trans <- transform(fs_comp_clean, trans)
?flowCore::transform
??estimateLogicle

## Inverse 
invLogicle <- inverseLogicleTransform(trans)
before <- transform(fs_comp_clean_trans, invLogicle)
exprs(before[[1]])[1:10,7:12]
exprs(fs[[1]])[1:10,7:12]

#Visualise file 
autoplot(fs[[1]]) #before transformation 
autoplot(fs_comp_clean_trans[[1]]) # after transformation.

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

autoplot(auto_gs, x="FSC-A", y="SSC-A", "noneDebris_gate", bins=256)
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
autoplot(auto_gs, x = 'FSC-A', y = 'FSC-H', "singlets", bins = 256)
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

#statistics ----
gs_pop_get_stats(auto_gs)
gs_pop_get_stats(auto_gs, "noneDebris_gate", "percent")
gs_pop_get_stats(auto_gs, "noneDebris_gate", type = pop.MFI) #Mean Fluorescent Intensity

pop.quantiles <- function(fr){
  chnls <- colnames(fr)
  res <- matrixStats::colQuantiles(exprs(fr), probs = 0.75)
  names(res) <- chnls
  res
}
gs_pop_get_stats(auto_gs, gs_get_pop_paths(auto_gs), type = pop.quantiles)

pop.mean <- function(fr){
  chnls <- colnames(fr)
  res <- colMeans(exprs(fr))
  names(res) <- chnls
  res
}
gs_pop_get_stats(auto_gs, gs_get_pop_paths(auto_gs), type = pop.mean)



#Faceting ----
# Plot metadata featrues one vs another
pData(auto_gs)$name
grep("CTRL",pData(auto_gs)$name )
grep("WNV",pData(auto_gs)$name )
# Create Metadata
meta.data <- data.frame(name = pData(auto_gs)$name,
                        group = sapply(strsplit(pData(auto_gs)$name,"_"), `[`, 1))

#a <- auto_gs
pData(auto_gs)$group = meta.data$group
p2 <- ggcyto(auto_gs[1:4], aes(x = "CD3", y = "CD4"), subset="singlets") 
p2 <- p2 + geom_hex(bins=256)
myPars <- ggcyto_par_set(limits = list(y = c(0,5), x = c(0,5)))
p2 <- p2 + myPars
p2 + facet_grid(group ~ name)

p3 <- ggcyto(auto_gs[1:9], aes(x = "CD3"), subset="singlets") 
p3 <- p3 + geom_histogram(bins = 256) 
p3 + facet_grid(group~name)

p3 <- ggcyto(auto_gs[10:18], aes(x = "CD3"), subset="singlets") 
p3 <- p3 + geom_histogram(bins = 256) 
p3 + facet_grid(group~name)


#overlay plot

p3 <- ggcyto(auto_gs[1:3], aes(x = "CD3", y = "CD4"), subset = "singlets")
p3 <- p3 + geom_hex(bins=128)
myPars <- ggcyto_par_set(limits = list(y = c(0,4), x = c(0,4)))
p3 <- p3  + myPars
p3 + geom_overlay(data = "BV750-A+BV711-A-",size = 0.01, alpha = 0.05, color = "purple")

#Exporting image
ggsave("plot1.png", p3)


### additional plot using flowviz ----
library(flowViz)
densityplot(~ `CD3`, fs_comp_clean_trans )

# Clustering using CytoTree ----
#https://ytdai.github.io/CytoTree/ti.html
library(CytoTree)
fcs.files <- list.files(path = "/Users/rudramanipokhrel/work/FlowCytometry/MJ191016 WNV 15 color/", pattern="*fcs",
                        full.names=TRUE)
fcs.data <- runExprsMerge(fcs.files, comp = TRUE, transformMethod = "logicle")
?runExprsMerge
colnames(fcs.data)
fcs.data <- fcs.data[,7:21]
# Refining the columns:
col <- sub("<NA>", "", colnames(fcs.data), perl=T)
col <- sub(".*<", "",col)
col <- sub(">", "",col)
colnames(fcs.data) <- col


# Create Metadata
meta.data <- data.frame(cell = rownames(fcs.data),
                        stage = sapply(strsplit(rownames(fcs.data),"_"), `[`, 1))
cyt <- createCYT(raw.data = fcs.data,markers = col,
                 meta.data = meta.data,
                 normalization.method = "None",
                 verbose = TRUE)
# Cluster cells by SOM algorithm
# Set random seed to make results reproducible
set.seed(1)
cyt <- runCluster(cyt, cluster.method = "som")
# Do not perform downsampling, if cell number is too large downsampling is recommended
#cyt <- processingCluster(cyt, downsampling.size = 0.1
cyt <- processingCluster(cyt)
# Visualization for cluster
plotCluster(cyt, item.use = c("tSNE_1", "tSNE_2"), category = "categorical",
            size = 100, show.cluser.id = TRUE)

# You can set xdim and ydim to specify the number of clusters
# the cluster number is xdim * ydim
set.seed(1)
cyt <- runCluster(cyt, cluster.method = "som", xdim = 5, ydim = 5)
cyt <- processingCluster(cyt, downsampling.size = 0.1)
plotCluster(cyt, item.use = c("tSNE_1", "tSNE_2"), category = "categorical",
            size = 100, show.cluser.id = TRUE) 


## Dimensionality Reduction

set.seed(42)
cyt <- runCluster(cyt, cluster.method = "som", xdim = 10, ydim = 10)
cyt <- processingCluster(cyt, perplexity = 5, downsampling.size = 0.1, 
                         force.resample = TRUE)

fs_comp_clean_trans
auto_gs
# Four popular dimensionality reduction methods are integrated 
# in CytoTree, namely PCA, tSNE, diffusion maps and UMAP.
# These four steps are optional steps

# run Principal Component Analysis (PCA)
cyt <- runFastPCA(cyt)
# run t-Distributed Stochastic Neighbor Embedding (tSNE)
set.seed(1)
cyt <- runTSNE(cyt)
# run Diffusion map
cyt <- runDiffusionMap(cyt)
# run Uniform Manifold Approximation and Projection (UMAP)
cyt <- runUMAP(cyt)

# 2D plot
plot2D(cyt, item.use = c("PC_1", "PC_2"), color.by = "CD3", 
       alpha = 1, main = "PCA", category = "numeric") + 
  scale_colour_gradientn(colors = c("#00599F","#EEEEEE","#FF3222"))

plot2D(cyt, item.use = c("tSNE_1", "tSNE_2"), color.by = "CD3", 
       alpha = 1, main = "tSNE", category = "numeric") + 
  scale_colour_gradientn(colors = c("#00599F","#EEEEEE","#FF3222"))
# 2D plot
plot2D(cyt, item.use = c("UMAP_1", "UMAP_2"), color.by = "CD3", 
       alpha = 1, main = "UMAP", category = "numeric") + 
  scale_colour_gradientn(colors = c("#00599F","#EEEEEE","#FF3222"))


# Build Trajectory
# 1. Build tree using raw expression matrix
cyt <- buildTree(cyt, dim.type = "raw")
# Tree plot
plotTree(cyt, color.by = "CD3", show.node.name = F, cex.size = 1) + 
  scale_colour_gradientn(colors = c("#00599F", "#EEEEEE", "#FF3222"))

# 5. Build tree using UMAP
cyt <- buildTree(cyt, dim.type = "umap", dim.use = 1:2)
# Tree plot
plotTree(cyt, color.by = "CD3", show.node.name = FALSE, cex.size = 1) + 
  scale_colour_gradientn(colors = c("#00599F", "#EEEEEE", "#FF3222"))


# Run differential expressed markers of different branch
diff.info <- runDiff(cyt)
# plot heatmap of clusters and branches
plotClusterHeatmap(cyt)

plotBranchHeatmap(cyt, colorRampPalette(c("#00599F", "#FFFFFF", "#FF3222"))(100), 
                  clustering_method = "complete")
?runDiff
?plotBranchHeatmap

table(meta.data$stage)

## Saving the filtered Data
fs_filt <- gs_pop_get_data(auto_gs, y = "singlets", inverse.transform = TRUE)
outDir <- file.path("~/Downloads", "fs_filt")
write.flowSet(fs_filt, outDir)

