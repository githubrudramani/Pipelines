## Enrichment Analysis
library(tidyverse)
library(limma)
library(gplots) #for heatmaps
library(DT) #interactive and searchable tables of our GSEA results
library(GSEABase) #functions and methods for Gene Set Enrichment Analysis
library(Biobase) #base functions for bioconductor; required by GSEABase
library(GSVA) #Gene Set Variation Analysis, a non-parametric and unsupervised method for estimating variation of gene set enrichment across samples.
library(gprofiler2) #tools for accessing the GO enrichment results using g:Profiler web resources
library(clusterProfiler) # provides a suite of tools for functional enrichment analysis
library(msigdbr) # access to msigdb collections directly within R
library(enrichplot) # great for making the standard GSEA enrichment plots
library(AnnotationHub)
library(org.Hs.eg.db)
library(pathview)


## GO enrichment analyiss ####
# First plot overall pathways
# load the logfold change stat files
stat = read.csv("top.csv", row.names = 1)
head(stat)


gost.res <- gost(rownames(stat), organism = "hsapiens", correction_method = "fdr")

gostplot(gost.res, interactive = F, capped = F) #set interactive=FALSE to get plot for publications

ggsave(paste0("summary_pathways.pdf"))


# Extract data
#filter(gost.res_up$result, p_value < 0.05)

#Second plot individual GO terms and KEGG for up vs down genes

data <-  stat

# put geneID information
data$geneID = rownames(data)

#Use mapIds to convert the gene name to Entrez ID or ENSEMBL ID

data$entrez <- mapIds(x = org.Hs.eg.db,
                      keys = data$geneID,
                      column = "ENTREZID",
                      keytype = "SYMBOL",
                      multiVals = "first")

# Remove any genes that do not have any entrez identifiers
entrez <- subset(data, is.na(entrez) == FALSE)

# Create a matrix of gene log2 fold changes
gene_matrix <- entrez$logFC
names(gene_matrix) <- entrez$entrez


# Enrich genes using the KEGG database
kegg_enrich <- enrichKEGG(gene = names(gene_matrix),
                          organism = 'hsa',
                          pvalueCutoff = 0.05,
                          qvalueCutoff = 0.10)


dotplot(kegg_enrich, showCategory = 15, font.size = 11)
ggsave("kegg_enrich.pdf")



# Gene ontology plot
c("BP", "CC", "MF")

go = "MF" # you can change this for three of above
go_enrich <- enrichGO(gene = names(gene_matrix),
                      OrgDb = 'org.Hs.eg.db',
                      readable = T,
                      ont = go,
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.10)

dotplot(go_enrich, showCategory = 10, font.size = 11)
ggsave("go_MF_not_found.pdf")

## GSEA
# for gsea dont appy any filter to logfc, however you can select pval <0.05
hs_gsea_c2 <- msigdbr(species = "Homo sapiens", # change depending on species your data came from
                      category = "C2") %>% # choose your msigdb collection of interest
  dplyr::select(gs_name, gene_symbol) #just get the columns corresponding to signature name and gene symbols of genes in each signature 

data$geneID =  rownames(data)

data <- dplyr::select(data, geneID, logFC)
mydata.gsea <- data$logFC
names(mydata.gsea) <- as.character(data$geneID)
mydata.gsea <- sort(mydata.gsea, decreasing = TRUE)
head(mydata.gsea)

# run GSEA using the 'GSEA' function from clusterProfiler
myGSEA.res <- GSEA(mydata.gsea, TERM2GENE=hs_gsea_c2, verbose=FALSE)

myGSEA.df <- as_tibble(myGSEA.res@result)
myGSEA.df = filter(myGSEA.df, p.adjust < 0.05)
myGSEA.df = myGSEA.df[order(-myGSEA.df$NES),]


# add a variable to this result that matches enrichment direction with phenotype
myGSEA.df <- myGSEA.df %>%
  mutate(phenotype = case_when(
    NES > 0 ~ "Day4", # change the naming here (it can be treatment)
    NES < 0 ~ "Day1")) # it is control

dim = dim(myGSEA.df)

myGSEA.df$NES
f = c(dim[1]-20, dim[1])
f

colnames(myGSEA.df)
# create 'bubble plot' to summarize y signatures across x phenotypes
ggplot(myGSEA.df, aes(x=phenotype, y=ID)) + 
  geom_point(aes(size=setSize, color = NES, alpha=-log10(p.adjust))) +
  scale_color_gradient(low="blue", high="red") +
  theme(axis.text.x = element_text(size = 10, face = "bold", angle = 60, vjust = 0.5, hjust=0.5),
        axis.text.y = element_text(size = 10, face = "bold"))
ggsave("gsea.pdf")




# view results as an interactive table

myGSEA.df <- as_tibble(myGSEA.res@result)
datatable(myGSEA.df, 
          extensions = c('KeyTable', "FixedHeader"), 
          caption = 'Signatures enriched in Interaction Term',
          options = list(keys = TRUE, searchHighlight = TRUE, pageLength = 10, lengthMenu = c("10", "25", "50", "100"))) %>%
  formatRound(columns=c(2:10), digits=2)

write.csv(myGSEA.df, "GSEA_c2_.csv")
# We are mostly interested in performing the GSEA using the ranked gene list 
# from the interaction term, and to have this done separated by the three 
# different tissue types (epithelium, stroma, and adipose). 


# create enrichment plots using the enrichplot package
gseaplot2(myGSEA.res, 
          geneSetID = c(25,34), # it is your id of interest, you can search id from above data table
          pvalue_table = FALSE) #can set this to FALSE for a cleaner plot
#title = myGSEA.res$Description[47]) #can also turn off this title

ggsave("GSEA_rank_c2_Interaction_down.pdf")

gseaplot2(myGSEA.res, 
          geneSetID = c(1,2),
          #geneSetID = c(25,34),
          pvalue_table = FALSE) #can set this to FALSE for a cleaner plot
#title = myGSEA.res$Description[47]) #can also turn off this title

ggsave("GSEA_rank_c2_Interaction_up.pdf")
