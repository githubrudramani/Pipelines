## function for DE analysis between two groups
de <- function(group1, group2, dds = dds){
  res1 <-results(dds, contrast = list(group1, group2) )
  res1 <- as.data.frame(res1)
  res1 <- res1[res1$pvalue <= 0.05,]
  res1 <- res1[order(-res1$log2FoldChange),]
  return(res1)
}

## function to plot box plot for gene expression
boxPlot <- function (data, x,y, color_list, fontsize = 12, xlabel, ylabel) {
  ggplot(data = data, aes(data[,x], data[,y])) +theme_classic() +
    geom_boxplot() +
    geom_boxplot(fill = color_list) +
    theme(aspect.ratio=1,
          axis.text.x = element_text(size = fontsize, face = "bold", angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = fontsize, face = "bold", angle = 0, vjust = 0, hjust=0)
          ,axis.title.x=element_text(size=fontsize,face="bold", vjust = 0.5 )
          ,axis.title.y=element_text(size=fontsize,face="bold", hjust = 0.5, vjust = 1.5 ),
          plot.title = element_text(size = fontsize, face = "bold"))+
    stat_summary(fun=mean, geom="point", shape=12, size=4) +
    xlab("Samples") +
    ylab("Normalized Expression") +
    ggtitle(paste("Expression of ", y))
}

PCAplot <- function(vsd=vsd, sample=sample, vars = vars){
  (data <- plotPCA(vsd, intgroup=colnames(sample), returnData=TRUE))
  (percentVar <- 100*round(attr(data, "percentVar"),2))
  groups <- sample[,vars[1]]
  shape <- sample[,vars[2]]
  ggplot(data, aes(PC1,PC2, col=groups, shape = shape)) + geom_point(size = 3) +
    ylab(paste0("PC2: ",percentVar[2], " % variance"))+
    xlab(paste0("PC1: ",percentVar[1], " % variance"))+
    theme(axis.text=element_text(size=12),
          axis.title=element_text(size=12,face="bold"),
          legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))+
    coord_fixed(ratio = 1)
}

plotHeatmap <- function(data ) {
  mSet<-InitDataObjects("conc", "stat", FALSE)
  mSet<-Read.TextData(mSet, data, "rowu", "disc");
  mSet<-SanityCheckData(mSet)
  mSet<-ReplaceMin(mSet);
  mSet<-PreparePrenormData(mSet)
  mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
  mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
  mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
  mSet<-PlotHeatMap(mSet, "Analysis/plots/heatmap_ all_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "overview",  T, T, NULL, T, F)
  mSet<-PlotHeatMap(mSet, "Analysis/plots/heatmap_avg_", "pdf", 72, width=NA, "norm", "row", "euclidean", "ward.D","bwm", "overview", T, T, NULL, T, T)
}

print(paste("Imported five functions: ", "de,", "boxpPlot,", "PCAplot,","plotHeatmap", " and call_DEseq2"))

call_DESeq2 <- function() {
    count <- read.csv(count_dir, row.names = 1)
    sample <-  read.csv(sample_dir , row.names = 1)
    anno <- sample
    # converting sample columns to factor
    #sapply(colnames(sample), FUN= function(x) sample[,x] = as.factor(sample[,x] ))
    for(x in colnames(sample)) {
      sample[,x] = as.factor(sample[,x])
      }
    
    if (mean(rownames(sample)!=colnames(count))){
      print("your sample order is not matching with columns of count")
    }else{
      
    
    
    low_depth_samples <- colnames(count[, colSums(count) <= threshold])
    if (length(low_depth_samples) >= 1) {
      print(paste0("Samples having less than ", threshold, " counts:" ))
      print(low_depth_samples)
    } else {
      print("All samples passed the threshold")
    }
    
    ### Filter low count genes
    print("filtering low count genes")
    print(paste("Total genes before filtering:", dim(count)[1]))
    keep <- rowSums(count>minimum_count)> at_least_in_samples
    f <- count[keep,]
    print(paste("Total genes after filtering:", dim(f)[1]))
    
    col <- colnames(sample)
    
    print("creating model matrix")
    ml <- model.matrix(design, sample)
    ml_df = as.data.frame(unname(ml)) # some of last combinatins may be zeros
    idx <- which(colSums(ml_df)!=0)
    ml <- ml[,idx]
    
    dds <- DESeqDataSetFromMatrix(countData =f ,
                                  colData = sample,
                                  design = ml)
    
                      
    ## Normalize the data
    print("Normalizing the data")
    vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
    #rld <- rlog(dds, blind=FALSE)
    dir.create("Analysis")
    dir.create("Analysis/data")
    write.csv(assay(vsd), "Analysis/data/vsd_normalized.csv")
    #write.csv(assay(rld), "data/rld_normalized.csv")
    
    ## Plot PCAs
    print("Plotting PCA and dendogram")
    vars <- colnames(sample)
    dir.create("Analysis/plots")

    PCAplot(vsd, sample, vars)
    ggsave(paste0("Analysis/plots/", "pca_vsd_with_two_variables",".pdf"), width = width, height = height)
    
    for (v in vars){
      plotPCA(vsd, intgroup=v)
      ggsave(paste0("Analysis/plots/", "pca_vsd_", v,".pdf"), width = width, height = height)
      
    }
    ## Plot dendogram
    ## Plot the coloring
    
    hc <- hclust(dist(t(assay(vsd))))
    pdf("Analysis/plots/dendogram_vsd.pdf")
    plot(as.phylo(hc), cex = 0.8,
         no.margin = TRUE)
    dev.off()
    print("Calling DEseq() function")
    ## DE analysis 
    dds <- DESeq(dds)
    norm <- data.frame(counts(dds, normalized = T))
    write.csv(norm, "Analysis/data/dd_normalized.csv")
    saveRDS(dds, "Analysis/dds.rds")
    dir.create("Analysis/plots/expression") 
    group <- resultsNames(dds)
    print("resultsNames in dds:")
    print(group)
    top_genes <- c()
    for (i in 1:(length(group)-1)){
      for (j  in 2:length(group) ) {
        if (i <j){
          name = paste0(group[i],"_vs_", group[j])
          compare =   de(group[i], group[j], dds = dds)
          compare <- drop_na(compare)
          write.csv(compare, paste0("Analysis/data/", name, ".csv"))
          print(paste("Plotting the expression of significant 10 genes in", name))
          top5 <- head(compare,5) %>% filter(log2FoldChange > 1)
          bottom5 <- tail(compare,5) %>% filter(log2FoldChange < -1)
          #top10 <- assay(vsd)[rownames(rbind(top5, bottom5)),]
          top10 <- norm[rownames(rbind(top5, bottom5)),]
          data <- cbind(sample, t(top10))

          genes <- row.names(top10)
          top_genes <- append(top_genes, genes)
          folder <- paste0("Analysis/plots/expression/", name)
          dir.create(folder)
          # plotting for all variables in sample
          for (f in factor) {
            n_colors <- length(levels(sample[,f]))
            palette <- rainbow(n_colors) 
            folder2 <- paste0(folder,"/", f)
            dir.create(folder2)
            for (k in genes){
              df <- data[,c(f,k)]
              boxPlot(df, x = f, y =  k, color_list = palette, 
                      xlabel = f, ylabel = "Normalized Expression")
              ggsave(paste0(folder2,"/", k , ".pdf"), width = 5, height = 5 )
              
            }}
        }}}
    
    top_genes <- unique(top_genes)
    vsd_data <- assay(vsd)
    vsd_data <- vsd_data[top_genes,]
    
    print("Plotting heatmap")
    rownames(anno) <- rownames(t(assay(vsd)))
    color <- colorRampPalette(c("darkgreen", "gray", "darkred"))(1000)
    pdf("Analysis/plots/heatmap.pdf")
    pheatmap(assay(vsd)[top_genes,], cluster_rows=TRUE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=anno, color = color)
    dev.off()
    
    
    metabo <- cbind(sample, t(assay(vsd)))
    metabo <- metabo[, c(vars[1], top_genes)]
    write.csv(metabo, "Analysis/data/tometabo.csv")
    plotHeatmap("Analysis/data/tometabo.csv")

    
    }
unlink("*png")
unlink("*qs")
print("Processed finished")
print("All the analysis are  in Analysis folder")
print("dds instance  is saved to dds.RDS, load it to do your desired analysis")
 }


