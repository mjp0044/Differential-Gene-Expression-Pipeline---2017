### You will need to set your working directory to the location you have your data.
#setwd("~/Documents/Functional Genomics/Counts")

#### Install the DESeq2 package if you have not already
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")

## Load the DESeq2 library 
library(DESeq2)

#####   1.3 Input data
##  Input the count data

countdata <-(read.csv("gene_count_matrix.csv", header=TRUE, row.names=1))
dim(countdata)
head(countdata)

##Import the metadata
colData <-(read.table("colData.txt", header=TRUE, row.names=1))
dim(colData)
head(colData)

## Creat the dataset and define model  (page 6 of the manual)
dds <- DESeqDataSetFromMatrix(countData = countdata, 
                              colData=colData,
                              design = ~treatment)
#look at it
dds

##### starting at  1.3.6 Prefiltering
# Here we perform a minimal pre-filtering to remove rows that have only 0 or 1 read.
dds <- dds[ rowSums(counts(dds)) > 1, ]

## set factors  ###### Note you need to change condition to treatment and levels to control and restricted
# example:
#dds$condition <- factor(dds$condition, levels=c("untreated","treatment"))
dds$condition <- factor(dds$treatment, levels=c("Control","Restricted"))


### 1.4 Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res


#We can order our results table by the smallest adjusted p value:
resOrdered <- res[order(res$padj),]
resOrdered
#We can summarize some basic tallies using the summary function.
summary(res)


#How many adjusted p-values were less than 0.05?
sum(res$padj < 0.05, na.rm=TRUE)

#If the adjusted p value cutoff will be a value other than 0:1, alpha should be set to that value:
res01 <- results(dds, alpha=0.05)
summary(res01)
sum(res01$padj < 0.05, na.rm=TRUE)


##1.5.1 MA-plot
##plotMA shows the log2 fold changes attributable to a given variable over the meanof normalized counts. Points will be colored red if the adjusted p value is less than 0.05. Points which fall out of the window are plotted as open triangles pointing either up or down
plotMA(res05, main="DESeq2", ylim=c(-2,2))

#After calling plotMA, one can use the function identify to interactively detect the row number of individualgenes by clicking on the plot. One can then recover the gene identiers by saving the resulting indices:

idx <- identify(res$baseMean, res$log2FoldChange)
rownames(res)[idx]

# A column lfcMLE with the unshrunken maximum likelihood estimate (MLE) for the log2 fold change will be added with an additional argument to results:
resMLE <- results(dds, addMLE=TRUE)
head(resMLE, 4)

#One can make an MA-plot of the unshrunken estimates
plotMA(resMLE, MLE=TRUE, main="unshrunken LFC", ylim=c(-2,2))


##  1.5.2 Plot counts
# You can select the gene to plot by rowname or by numeric index.
plotCounts(dds, gene=which.min(res$padj), intgroup="treatment")

##  
d <- plotCounts(dds, gene=which.min(res$padj), intgroup="treatment",returnData=TRUE)
library("ggplot2")

ggplot(d, aes(x=treatment, y=count)) +
  geom_point(position=position_jitter(w=0.1,h=0)) +
  scale_y_log10(breaks=c(25,100,400))
##
write.csv(as.data.frame(resOrdered),
          file="condition_treated_results.csv")  

## 2.1.2 Extracting transformed values
rld <- rlog(dds)
vsd <- varianceStabilizingTransformation(dds)
head(assay(rld), 3)

## 2.2.1 Heatmap of the count matrix
library("pheatmap")
select <- order(rowMeans(counts(dds,normalized=TRUE)),decreasing=TRUE)[1:50]
nt <- normTransform(dds) # defaults to log2(x+1)
log2.norm.counts <- assay(nt)[select,]
df <- as.data.frame(colData(dds)[,c("treatment","type")])
pheatmap(log2.norm.counts, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)
pheatmap(assay(rld)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)  
pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df)



#2.2.2 Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(rld)))
library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$treatment, rld$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)


# 2.2.3 Principal component plot of the samples

library(ggplot2)
Obj1 <- plotPCA(rld, intgroup=c("treatment", "type"), returnData=TRUE)
ggplot(data = Obj1, aes(x = PC1, y = PC2, col = treatment)) +
  geom_point(size = 2) +
  #geom_label()+
  geom_label(aes(label = Obj1$name), position = position_nudge(y = 5)) +
  theme_classic() +
  theme(legend.position = "none") +
  xlab ("PC1 (69%)") + ylab ("PC2 (12%)")
