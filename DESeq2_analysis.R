#Load the Libraries
library("DESeq2")
library('edgeR')
library("tidyverse")
library('variancePartition')
library('BiocParallel')
library("dplyr")
library("ggplot2")
library("pheatmap")
library("vsn")
library("RColorBrewer")
library("EnhancedVolcano")
library("clusterProfiler")
library("enrichplot")
library("writexl")
library("AnnotationHub")
require(DOSE)

#Set working directory
setwd(getwd())

## DATA LOADING AND PREPROCESSING

# Load the quantification dataset
diff_exp <- read.table(file = "data.csv", sep = '\t', header = TRUE, row.names = 1)
head(diff_exp, 2)
nrow(diff_exp)

#Remove all the genes with 0 count across all samples
diff_exp <- diff_exp %>%
  filter(rowSums(.) != 0)
nrow(diff_exp)

#Select the columns (samples) to be used for pairwise analysis
count_matrix <- diff_exp[c("control1","control2","control3","treated1","treated2","treated3")]
count_matrix <- as.matrix(count_matrix)
head(count_matrix,2)

#Prepare sample information dataframe
coldata <- data.frame(
  sample = c("control1","control2","control3","treated1","treated2","treated3"),
  condition = c("Control","Control","Control","treated", "treated", "treated"),
  row.names = "sample")
coldata$condition <- as.factor(coldata$condition)
coldata

#Check if row names from the sample info dataframe correspond to the column names in quantification dataset
all(rownames(coldata) == colnames(count_matrix))

#Convert data to integer
count_matrix <- round(count_matrix, digits = 0)

#Construct a DESeqDataSet (dds) using the count matrix and sample information 
dds <- DESeqDataSetFromMatrix(countData = count_matrix, colData = coldata, design = ~ condition)
dds

#Remove genes with low counts i.e count below 10
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]
dds

#Estimate size factor
dds <- estimateSizeFactors(dds)

##DATA VISUALIZATION

#Plot PCA of log2 transformed counts
selog <- SummarizedExperiment(log2(counts(dds, normalized=TRUE) + 1), colData=colData(dds))

#Use DESeqTransform() to trigger the plotPCA method.
#use the default i.e top 500 genes for PCA plot
plotPCA(DESeqTransform(selog)) 
#use all genes for PCA plot
plotPCA(DESeqTransform(selog),ntop = nrow(count_matrix))

#Use other transformation approaches for the PCA plot
#Use regularized  log transformation (rlog)
rld <- rlog(dds, blind=FALSE)
head(assay(rld), 2)
plotPCA(rlog(dds)) #default 
plotPCA(rlog(dds), intgroup=c("condition"), ntop = nrow(count_matrix)) #all genes

#Use variance stabilizing transformation (vst)
vsd <- vst(dds, blind=FALSE)
head(assay(vsd), 2)
plotPCA(vsd, intgroup=c("condition")) #default ntop=500
plotPCA(vsd, intgroup=c("condition"), ntop = nrow(count_matrix)) #all genes

#Use ggplot for vst transformed data PCA plot
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaData <- plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed()

#Visualize the transformed data using meanSdP to see the effect transformations
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

#Visualize the transformed data using PHeatmap
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:20]
df <- as.data.frame(colData(dds)[c("condition")])
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)

#Visualize sample-to-sample distances in Heatmap
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

#Visualize sample distribution in clusters using pheatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

#Visualize sample distribution using plotMDS
pch <- c(16,17)
pch <- pch[as.numeric(coldata$condition)]
colors <- rep(c("#9457EB", "#177245"),3)
colors <- colors[as.numeric(coldata$condition)]
plotMDS(dds, col=colors, pch=pch,cex=3.0)
legend("bottomright", legend=levels(coldata$condition), pch=c(16,17), col=rep(c("#9457EB", "#177245"),3),cex=1.5)

#Visualize data using plotMD
plotMD(dds, column=6)
abline(h=0, col="red", lty=2, lwd=2)

#Estimate size factor of the dds
sizeFactors(dds)
plot(sizeFactors(dds), colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))

#Estimate dispersion of the dds and plot
dds = estimateDispersions(dds)
plotDispEsts(dds)

##DIFFERENTIAL EXPRESSION ANALYSIS

#Conduct DESeq analysis
dds <- DESeq(dds)

#view result
res <- results(dds) #default
summary(res)
head(res)

#save all the result from the analysis
res_ordered <- res[order(res$padj),] #sorting data using the padj. value
write.csv(as.data.frame(res_ordered), file="DEG_treated_vs_control.csv")

#Filter and save only the significant (padj < 0.05) result
#log2FoldChange threshold can also be used in selecting significant DEG
resFilt <- res[which(res$padj < 0.05 & abs(res$log2FoldChange) >=0.0), ]
resOrdered <- resFilt[order(resFilt$pvalue),] #sorting data using the pvalue
write.csv(resOrdered, file="DEG_ctrl_vs_treated_0.05.csv") 

#Check the number of DEG using adjusted pvalues threshold 
sum(res$padj < 0.05, na.rm=T) 

#Check read count of a single gene between control and treated group
plotCounts(dds, gene=which.min(res$padj), intgroup="condition")

#Count the DEG
resSigind = res[ which(res$padj < 0.05 & res$log2FoldChange > 0), ] #upregulated
resSigrep = res[ which(res$padj < 0.05 & res$log2FoldChange < 0), ] #downregulated
resSig = rbind(resSigind, resSigrep)
dim(resSigind)
dim(resSigrep)
rownames(resSigind)
rownames(resSigrep)
head(resSigind)

#Plot the padj. value using histogram
hist(res$padj, 
     breaks = seq(0, 1, 0.05), 
     col = "purple1",
     main = "Histogram of Padj. from treated vs control results",
     xlab = "Padj-value")


#If necessary, map Ensembl to gene symbol and make dataframe with it
library(org.Dr.eg.db)
res.df <- as.data.frame(res)
res.df$gene_symbol <- mapIds(org.Dr.eg.db, key = row.names(res.df), keytype = "ENSEMBL", column = "SYMBOL")

#Make volcano plot of DEGs result
#Determine point colors based on significance and sign of the logFC and the Padj.
res.df <- res.df %>% 
  mutate(point_color = case_when(
    padj < 0.05 & log2FoldChange < 0 ~ "down", # significantly down
    padj < 0.05 & log2FoldChange > 0 ~ "up", # significantly up
    TRUE ~ "NS") # not significant
  )

#Create the Volcano plot
volcano_plot <- ggplot(data = res.df, aes(x = log2FoldChange, y = -log10(padj), col = point_color)) +
  geom_hline(yintercept = -log10(0.05), col = "green4", linetype = 'dashed')+
  geom_point(color = ifelse(res.df$log2FoldChange < -0.001 & res.df$padj<0.05, "blue", ifelse(res.df$log2FoldChange > 0.001 & res.df$padj<0.05, "red", "grey")), size = 1.5) +
  geom_text(aes(label = ifelse(res.df$log2FoldChange < -2 & -log10(padj) >20 | res.df$log2FoldChange > 2 & -log10(padj) >20, res.df$gene_symbol, ""), size = 3),
            vjust = -0.5, hjust = 0, size = 3)+
  scale_color_manual(values = c("blue", "grey", "red"),
                     labels = c("Down", "NS", "Up"))+
  coord_cartesian(ylim = c(0, 5), xlim = c(-5, 5))+
  labs(color = 'Annotation', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"Adjp-value"), size = 5)+
  ggtitle('Control vs Treated') +
  guides(shape = guide_legend(override.aes = list(color = NULL, size = 4)))+
  scale_shape_manual(values = c(16, 17, 18))+
  theme(axis.line = element_line(),panel.background = element_blank())
print(volcano_plot)
ggsave("contol_treated_volcano_plot.png", plot = volcano_plot, width = 8, height = 6, dpi = 300)


##CORRELATIONS BETWEEN SAMPLE PAIRS
treated_control <- diff_exp[c("control1","control2","control3","treated1","treated2","treated3")]

#Find the correlation between samples
correlation_matrix <- cor(treated_control)

# Create the correlation heatmap using pheatmap
htmap <- pheatmap(correlation_matrix, 
         color = colorRampPalette(c("blue", "white", "red"))(100), 
         main = "Treated vs Control Correlation Heatmap",
         fontsize = 12,
         cluster_rows = TRUE,
         cluster_cols = TRUE,)
print(htmap)
ggsave("treated_contol_heatmap.png", plot = htmap, width = 8, height = 8, dpi = 300)


## GSEA ENRICHMENT ANALYSIS

#Select the desired organism's annotation;
#Download the Zebrafish annotation file from http://bioconductor.org/packages/release/BiocViews.html#___OrgDb
BiocManager::install("org.Dr.eg.db", character.only = TRUE)
library("org.Dr.eg.db", character.only = TRUE)
keytypes(org.Dr.eg.db)

#Select the gene name and log2FoldChange columns from the DEG dataframe
new_df = subset(res.df, select = c(Gene_Name,log2FoldChange))
head(new_df)

new_df$log2FoldChange = as.numeric(as.character(new_df$log2FoldChange)) #convert LFC column to numeric
df <- new_df[new_df$log2FoldChange != 0, ] #keep rows with non zero value keep

#Select only the log2FoldChange 
original_gene_list <- df$log2FoldChange
#original_gene_list <- as.numeric(original_gene_list) #convert character to numeric, if needed

#Assign log2FoldChange to their respective genes
names(original_gene_list) <- df$Gene_Name 

#Drop NA values, if any 
gene_list<-na.omit(original_gene_list)

#Sort the list in decreasing order of log2FoldChange
gene_list = sort(gene_list, decreasing = TRUE)

#For keytype use GENENAME for gene symbol or ENSEMBL for ensembl
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL",  
             nPerm = 10000, 
             minGSSize = 2, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Dr.eg.db, 
             pAdjustMethod = "BH")

#convert enrichment object to dataframe
enrich_pathways<-as.data.frame(gse) 
write_xlsx(enrich_pathways, 'GSEA_treated_control.xlsx')

#Visualize enrichment result in dotplot
dotplot(gse, showCategory=5, split=".sign") + facet_grid(.~.sign)

#Visualize enrichment result in emapplot
x2 <- pairwise_termsim(gse)
emapplot(x2)

#Visualize enrichment result in cnetplot
cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 3)

#Visualize enrichment result in Ridge
ridgeplot(gse) + labs(x = "enrichment distribution")

#Visualize gene set in a given term
#Use the `Gene Set` param for the index in the title, and as the value for geneSetId.
#The first gene set is 1, second gene set is 2, etc.
gseaplot(gse, by = "all", title = gse$Description[15], geneSetID = 15) #both gse$Description and genesetID must be thesame