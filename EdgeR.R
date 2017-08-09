#edgeR

library(limma)
library(edgeR)
library(RColorBrewer)
library(mixOmics)
library(HTSFilter)

directory <- "~/Desktop/"
dir(directory)

#read in data & pheno info
rawCountTable <- read.csv(paste0(directory,"gene_count_matrix.csv"),header=TRUE,row.names=1)
head(rawCountTable)
sampleInfo <- read.csv(paste0(directory,"ballgown/pheno_data.csv"), header=TRUE,
                         row.names=1)

#check number of reads
nrow(rawCountTable)

#replace with gene names like MAGIC
rawCountTable <- rawCountTable[,match(rownames(sampleInfo),
                                      colnames(rawCountTable))]
head(rawCountTable)

#create table with treatment & library sizes
dgeFull <- DGEList(rawCountTable, group=sampleInfo$treatment)
dgeFull

#add sample info into DGElist data (optional)
dgeFull$sampleInfo <- sampleInfo
dgeFull

#preparing data for analysis

#data exploration & quality assessment
pseudoCounts <- log2(dgeFull$counts+1)
head(pseudoCounts)
hist(pseudoCounts[,"daphE1"])
quartz()
boxplot(pseudoCounts, col=sampleInfo$treatment)

quartz()
par(mfrow=c(1,2))
## treated2 vs treated3
# A values
avalues <- (pseudoCounts[,1] + pseudoCounts[,2])/2
# M values
mvalues <- (pseudoCounts[,1] - pseudoCounts[,2])
plot(avalues, mvalues, xlab="A", ylab="M", pch=19, main="caloric restriction")
abline(h=0, col="red")
## untreated3 vs untreated4
# A values
avalues <- (pseudoCounts[,3] + pseudoCounts[,4])/2
# M values
mvalues <- (pseudoCounts[,3] - pseudoCounts[,4])
plot(avalues, mvalues, xlab="A", ylab="M", pch=19, main="ad lib")
abline(h=0, col="red")

#visualize treatments in 2D space
quartz()
plotMDS(pseudoCounts)

#HEATMAP FINALLY YAY on pseudoCounts
sampleDists <- as.matrix(dist(t(pseudoCounts)))
sampleDists
cimColor <- colorRampPalette(rev(brewer.pal(9, "Blues")))(16)
cimColor
quartz()
cim(sampleDists, color=cimColor, symkey = FALSE)

#remove genes with 0 counts for all samples
dge <- DGEList(dgeFull$counts[apply(dgeFull$counts, 1, sum) != 0, ],
               group = dgeFull$sampleInfo$treatment)
dge$sampleInfo <- dgeFull$sampleInfo
head(dge$counts)

#est nomalization factors
dge <- calcNormFactors(dge)
dge$samples

#est common & tagwise dispersion
dge <- estimateCommonDisp(dge)
dge <- estimateTagwiseDisp(dge)
dge

#perform exact test
dgeTest <- exactTest(dge)
dgeTest

#filtering
quartz()
filtData <- HTSFilter(dge)$filteredData

dgeTestFilt <- exactTest(filtData)
dgeTestFilt

#diagnostic plot
quartz()
hist(dgeTest$table[,"PValue"], breaks=50)
quartz()
hist(dgeTestFilt$table[,"PValue"], breaks=50)

#inspecting results
resNoFilt <- topTags(dgeTest, n=nrow(dgeTest$table))
head(resNoFilt)
resNoFilt
resFilt <- topTags(dgeTestFilt, n=nrow(dgeTest$table))
head(resFilt)

#number before filtering
sum(resNoFilt$table$FDR < 0.05)

#number after filtering
sum(resFilt$table$FDR < 0.05)

#extract & sort DE genes
sigDownReg <- resFilt$table[resFilt$table$FDR<0.05,]
sigDownReg <- sigDownReg[order(sigDownReg$logFC),]
head(sigDownReg)

sigUpReg <- sigDownReg[order(sigDownReg$logFC, decreasing=TRUE),]
head(sigUpReg)
nrow(sigUpReg)
write.csv(sigDownReg, file="sigDownReg.csv")
write.csv(sigUpReg, file="sigUpReg.csv")

quartz()
plotSmear(dgeTestFilt,
          de.tags = rownames(resFilt$table)[which(resFilt$table$FDR<0.01)])

#VOLCANOPLOTTT
volcanoData <- cbind(resFilt$table$logFC, -log10(resFilt$table$FDR))
colnames(volcanoData) <- c("logFC", "negLogPval")
head(volcanoData)
quartz()
plot(volcanoData, pch=19)




#transform norm counts in log counts per mill
y <- cpm(dge, log=TRUE, prior.count = 1)
head(y)
resFilt$table$FDR
#select 1% DE genes & make heatmappp
selY <- y[rownames(resFilt$table)[resFilt$table$FDR<0.0000000000001 & 
                                    abs(resFilt$table$logFC)>2.5],]
?cim
head(selY)
cimColor1 <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)[255:1]
cimColor1
quartz()
  cim(t(selY), color=cimColor1, symkey=FALSE, col.cex=.8, margins=c(11,5))
dev.off()

