# Load packages
library(tidyverse)
library(openxlsx)
library(edgeR)

#Set working directory
setwd(getwd())

#1. Load gene expression data
counts <- read.table(file = "gene_expression_count.csv", sep = ',', header = TRUE, row.names = 1)
head(counts)

#select desired samples, if need be
counts <- counts[c("control1","control2","control3","PFBSHC1","PFBSHC2","PFBSHC3","PFBSLC1","PFBSLC2","PFBSLC3","PFOSHC1","PFOSHC2","PFOSHC3","PFOSLC1","PFOSLC2","PFOSLC3")]
head(counts,2)

#Make DGEList object
d0 <- DGEList(counts)

## PREPROCESS THE GENE EXPRESSION DATASET
#Calculate normalization factors
d0 <- calcNormFactors(d0)
d0

#Filter low-expressed genes
cutoff <- 1
drop <- which(apply(cpm(d0), 1, max) < cutoff)
d <- d0[-drop,] 
dim(d)  

#or continue without filtering
d <- d0
dim(d) 

#check the sample names
snames <- colnames(counts)
snames

##DATA VISUALIZATION

#Visualize the data using multidimensional scaling (MDS) plot
plotMDS(d, col = as.numeric(snames)) 

# Visualize the data using Mean difference (MD) plot
plotMD(d, column=8) #each column rep a sample
abline(h=0, col="red", lty=2, lwd=2)

# group samples by treatment, remove the number (1,2,3) from each sample name
group <- substr(snames, 1, nchar(snames) - 1) #remove the last letter
group

## VOOM TRANSFORMATION
#Specify the model to be fitted. 
mm <- model.matrix(~0 + group)
mm

#Conduct voom
y <- voom(d, mm, plot = T) #if the voom plot is not very good, then filter the features

# Fit the linear models of limma
fit <- lmFit(y, mm)
head(coef(fit))

#Get DEG against control
#Comparisons between groups (log fold-changes) are obtained as contrasts of these fitted linear models:
#Specify which groups to compare:e.g control and PFOSHC
contr <- makeContrasts(groupControl - groupPFOSLC, levels = colnames(coef(fit)))
contr

#Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

# Use the Empirical Bayes to smooth standard errors
tmp <- eBayes(tmp)

#Filter the genes that are most differentially expressed
top.table <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table, 5) 

#Check the number of DEGs at padj <0.05
length(which(top.table$adj.P.Val < 0.05))

#save the DEGs into file
write.table(top.table, file = "DEG_control_vs_PFOSLC.xlsx", row.names = T, sep = "\t")

## CORELATION ANALYSIS BETWEEN BEHAVIOURAL DATA AND GENE EXPRESSION COUNT
# Load the behavioural dataset
behav <- read.xlsx("behavioural_data.xlsx", rowNames = T)

#if required select the desired columns
behav <- behav[c("control1","control2","control3","PFBSHC1","PFBSHC2","PFBSHC3","PFBSLC1","PFBSLC2","PFBSLC3","PFOSHC1","PFOSHC2","PFOSHC3","PFOSLC1","PFOSLC2","PFOSLC3")]
head(behav, 8)

#Select the desired behavioural endpoint (LON1/LOFF/LON2) to the correlated with gene expression by their row indexes (e.g 1)
#Only the numeric values are required here as such index and column name must be removed
#Only one end point must be selected at once e.g [1,] refers to LON1 endpoints
beh_cor <-unlist(behav[1,], use.names = FALSE)
beh_cor

# Specify the model matrix
mm <- model.matrix(~beh_cor)
head(mm)

#Conduct voom
y <- voom(d, mm, plot = F)

# Fit the linear models of limma
fit <- lmFit(y, mm)
tmp <- contrasts.fit(fit, coef = 2) # test "behavioural" coefficient

# Use the Empirical Bayes to smooth standard errors
tmp <- eBayes(tmp)

#Filter the genes that are most differentially expressed and sort by pvalues
top.table1 <- topTable(tmp, sort.by = "P", n = Inf)
head(top.table1, 20)

#save all the behavioural correlation results
write.xlsx(top.table1, file = "LON1_cor_treated_expression_count.xlsx", rowNames = T)

#select the most correlated gene with behavioural endpoint (LON1) base on padj
length(which(top.table1$adj.P.Val < 0.05)) 
sig <- top.table1[which(top.table1$adj.P.Val < 0.05),]
write.xlsx(sig, file = "LON2_cor_adjpval_treated_expression_count.xlsx", rowNames = T)

#select the most correlated gene with  behavioural endpoint (LON1) base on pval
length(which(top.table1$P.Value < 0.05))
sig1 <- top.table1[which(top.table1$P.Value < 0.05),]
head(sig, 10)
write.xlsx(sig1, file = "LON2_cor_pval_treated_expression_count.xlsx", rowNames = T)

# Checking the correlation result with one gene
# Limma fit a linear regression model (i.e straight line) to each gene, with the slope and intercept defined by the model coefficients:
ENSDARG00000018542 <- y$E["ENSDARG00000018542",]
plot(ENSDARG00000018542 ~ beh_cor, ylim = c(0, 2.0))
intercept <- coef(fit)["ENSDARG00000018542", "(Intercept)"]
View(intercept)

slope <- coef(fit)["ENSDARG00000018542", "beh_cor"] 
#the log fold change logFC is the slope of the line. 
#Here, a logFC of -0.062 means a 0.062 log2 CPM decrease in gene expression for each unit increase in LON1 or  a 4% decrease on the CPM scale (2^0.062 = 1.04).
View(slope)
abline(a = intercept, b = slope)

#The log2 counts per million are more variable at lower expression levels. The variance weights calculated by Voom addresses this situation.