setwd("~/Desktop")

####Ballgown Differential Expression Analysis Basline####
#Code compiled from Pertea et. al 2016 Transcript-level expression analysis of RNA-seq experiments with HISAT, StringTie, and Ballgown

#Load relevant R packages, downloaded previously
#Ballgown and genefilter DL from bioconductor, RSkittleBrewer from AlyssaFreeze github, and Devtools and dplyr installed using install.packages

library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

#Import the metadata summarizing the samples and treatments (must create this separately in a text editor;similar to colData file in garter snake practice)
pheno_data = read.csv("ballgown/pheno_data.csv")
pheno_data
#Import the countdata table made using the merge script.
# dataDir = the parent directory holding you StringTie output subdirectories/files
#samplePattern = pattern for ballgown to find in all your sample subdirectories
daphnia_countdata = ballgown(dataDir = "ballgown", samplePattern = "daph", pData =pheno_data)
daphnia_countdata
####Display data in tabular format####

#Filter low abundance genes out 
daphnia_countdata_filt = subset(daphnia_countdata, "rowVars(texpr(daphnia_countdata)) >1", genomesubset=TRUE)
daphnia_countdata_filt
#Identify transcripts that show sig difs between groups
results_transcripts = stattest(daphnia_countdata_filt, feature="transcript", covariate = "treatment", getFC = TRUE, meas="FPKM")
head(results_transcripts)
#Identify genes that show statistically sig diffs between groups
results_genes = stattest(daphnia_countdata_filt, feature="gene", covariate = "treatment", getFC = TRUE, meas="FPKM")
results_genes
#Add gene names and gene IDs to the results_trascripts data frame
results_transcripts = data.frame(geneNames=ballgown::geneNames(daphnia_countdata_filt), geneIDs=ballgown::geneIDs(daphnia_countdata_filt), results_transcripts)


#Sort the results from smallest P value to the largest
results_transcripts = arrange(results_transcripts,pval)
results_genes = arrange(results_genes,pval)

#Write results to a csv file that can be shared and distributed
write.csv(results_transcripts, "daphnia_transcript_results.csv", row.names = FALSE)
write.csv(results_genes, "daphnia_genes_results.csv", row.names = FALSE)

#Identify transcripts and genes with a q value <0.05
sig.transcripts <- subset(results_transcripts,results_transcripts$qval<0.05)
sig.genes <- subset(results_genes,results_genes$qval<0.05)
head(sig.transcripts)
count(sig.transcripts, header=TRUE)
count(sig.genes)

write.csv(sig.transcripts, "sigTranscriptsInitial.csv", row.names= FALSE)
write.csv(sig.genes, "sigGenesInitial.csv", row.names= FALSE)
nrow(sig.transcripts) -1
nrow(sig.genes) -1

number((unique(sig.transcripts$geneIDs))
####Display data in visual format####

#Make pretty plot colors
tropical = c('darkorange', 'dodgerblue', 'hotpink', 'limegreen', 'yellow')
palette(tropical)

#Show the distribution of gene abundances (measured as FPKM with ballgown) across samples, colored by treatment
fpkm = texpr(daphnia_countdata,meas="FPKM")
fpkm = log2(fpkm+1)
quartz()
boxplot(fpkm,col=as.numeric(pheno_data$treatment),las=2,ylab='log2(FPKM+1)')
?boxplot
#Plotting individual trascripts if you have one you want to target
ballgown::transcriptNames(daphnia_countdata)[58299
                                            ]
##  12
## "GTPBP6"
ballgown::geneIDs(daphnia_countdata)[hxAUG26us17g343t1]
##  12
## "GTPBP6

View(as.data.frame(daphnia_countdata))

quartz()
plot(fpkm[58648,] ~ pheno_data$treatment, border=c(1,2),main=paste(ballgown::geneNames(daphnia_countdata)[17777],' : ',ballgown::transcriptNames(daphnia_countdata)[17777]),pch=19,xlab="Treatments",ylab='log2(FPKM+1)')
points(fpkm[58648,] ~ jitter(as.numeric(pheno_data$treatment)),col=as.numeric(pheno_data$treatment))

#Plot structure and expression levels in a sample of all transcripts that share the same gene locus
plotTranscripts(ballgown::geneIDs(daphnia_countdata)[4], daphnia_countdata, main=c('Gene XIST in sample daphC3_CCGAAG'), sample=c('daphC3_CCGAAG'))

#Plot the average expression levels for all transcripts of a gene within different groups using the plotMeans function. Need to specify which gene to plot
quartz()
plotMeans('MSTRG.19923', daphnia_countdata_filt, groupvar = "treatment",legend=TRUE)

#MA-plots btwn control or treated
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
     
     
