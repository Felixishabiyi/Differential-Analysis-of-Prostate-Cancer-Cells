#reading the count matrix file intp R studio
pros <- read.delim("C:/Users/Administrator/Documents/RNA-seq/Datasets/GSE220097_LNCaP_expression.tsv",
                   sep = "\t")
#checking the dimension
dim(pros)

#observing the columns names
colnames(pros)

#Cleaning the columns names
colnames(pros) <- gsub("_1nM_R1881_","",colnames(pros))
colnames(pros) <- gsub("_R1_001","",colnames(pros))

colnames(pros)
#loading the GEOquery library to extract the metadata
library(GEOquery)

#loading the metadata from the NCBI via GEOquery
prosgser <- getGEO("GSE220097",GSEMatrix = TRUE)
prosgser

#extracting the metadata file
meta <- pData(phenoData(prosgser[[1]]))

#loading the tidyverse library for data manipulation
library(tidyverse)

#cleaning the metadata to extract our columns for the metadata
modi <- meta %>%
  select(1,14,21) %>%
  rename(treatment = characteristics_ch1.4)

#head of the count data
head(pros)
#summary of the count data
summary(pros)
#structure of the count data
str(pros)
#View pros
view(pros)

#Data Manipulation
library(dplyr)

#Extracting genenames from count data
pros_names <- pros$Gene

#checking the gene names
pros_names

#Extracting the count from the count matrix
pros <- subset(pros,select = -Gene)

#assigning row names to the count matrix
rownames(pros) <- pros_names

#checking the count matrix
pros
#normalization - count per million
pros_cpm <- apply(pros,2,function(x) {x/sum(as.numeric(x)) * 10^6})

#the total of each columns
colSums(pros_cpm)

#calculating the variant genes
pros_v <- apply(pros_cpm, 1, var)

#extracting the top 100 variant genes
pros_v_names <- names(pros_v[order(pros_v,decreasing = T)][1:100])

#assigning rownames to the metadata 
rownames(modi) <- colnames(pros)

#extarting the treatment column from the metadata
modiiio <- modi %>%
  select(2)

#loading heatmap
library(pheatmap)

#visualizing heatmap
pheatmap(pros[pros_v_names,],scale = 'row', show_rownames = F,angle_col = 45,
         annotation_col = modiiio)

#Differential Expression Analysis

#Countdata
countData <- as.matrix(pros)

class(modi)
view(modi)

#Design formula for the differential expression
designFormula <- "~ treatment"

library(DESeq2)
library(stats)

#Creating a Deseq2 Dataset from the count matrix and the metaData
dds <- DESeqDataSetFromMatrix(countData = round(countData), 
                              colData = modi, 
                              design = as.formula(designFormula))

#Printing the dataset
print(dds)

row.names(dds)
colnames(dds)

#Removing genes that do not have least 1 count in the count matrix
dds <- dds[ rowSums(DESeq2::counts(dds)) > 1, ]

#Implementing the DESeq function on dds
dds <- DESeq(dds)

#Computig the contrast for the 'treatment'
DEresults = results(dds, contrast = c("treatment", 
                                      'treatment: R1881 1nM + Enzalutamide 2uM',
                                      'treatment: R1881 1nM'))

#sort by increasing Pvalue
DEresults <- DEresults[order(DEresults$pvalue),]
print(DEresults)

#Activating the enhanced volcano library
library(EnhancedVolcano)
#Plotting the Volcano
EnhancedVolcano(DEresults, x = 'log2FoldChange', y = 'pvalue', 
                lab = rownames(DEresults))

DESeq2::plotMA(object = dds, ylim = c(-5, 5))
