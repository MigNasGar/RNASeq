##########################################
## Build count matrix with HTSeq counts ##
##########################################
#setwd("") to define the directory where the files are

library("scatterplot3d")
count.files<-c("T1.counts","T2.counts","T3.counts", "A1.counts", "A2.counts",'A3.counts')
#Define the files .counts to be used, for the control and affected factors

# Build the count matrix
count_matrix<-lapply(1:length(count.files), function(x){
  count.vec<-read.delim(file=paste(count.files[x]), header=FALSE, stringsAsFactor=FALSE)
  return(count.vec)
  }) 
#make a matrix by assigning count.files values to a function
count_matrix<-do.call(cbind, count_matrix) 
#merge the two columns of each file

##take a look on data
head(count_matrix)
tail(count_matrix)

# exclude the last 5 rows (they are aligment informaion)
count_matrix<-count_matrix[-c(1:4),] 
#exclue as colunas contendo N_unmapped, N_multimapping, N_noFeature, N_ambiguous uma vez que não são utilizadas na análise de RNAseq 
tail(count_matrix)
dim(count_matrix)

#the first column of each count.file has the gene name
#set gene names as row names
row.names(count_matrix)<-count_matrix[,1]
# keep counts only
count_matrix<-count_matrix[,c(2,4,6,8,10,12)]
#set names of count files as column names
colnames(count_matrix)<-c("T1","T2","T3", "A1", "A2", "A3")

#explore the count matrix
head(count_matrix)
##how meny trasncript do we have?
dim(count_matrix)

# build 2 data frames. The first one will contain sample information and the second gene annotation
sample_data<-data.frame(sample=names(count_matrix), condition=c(rep("T", times=3), rep("A", times=3)))


# We must ensure that all genes have 10 counts in at least one condition
summary(rowSums(count_matrix[,1:3]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
summary(rowSums(count_matrix[,4:6]))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 

# build a index of genes with good counts
#in this case good counts are defined as those satisfying count >= 10
index_good_counts<-sapply(1:nrow(count_matrix), function(x){
    idx<-all(count_matrix[x,1:3] > rep(9, 3)) | all(count_matrix[x,4:6] > rep(9,3)) 
    return(idx)
})
table(index_good_counts)
#index_good_counts
#FALSE  TRUE 

####################
## comment on this number
# keep those genes that have good counts
count_matrix<-count_matrix[index_good_counts,]
dim(count_matrix)
# [1] 8223    6

######################
## Data exploration ##
######################
# ATENTION!!
head(sample_data$condition)
#[1] "T"   "T"   "T"   "A" "A" "A"
#Levels: T A

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

# Build a DESeqDataset object
# source("http://bioconductor.org/biocLite.R")
# biocLite("DESeq2")
library(DESeq2)

# ?DESeqDataSetFromMatrix
# We use count_matrix as countData, sample_data as colData and we choose condition column of sample_data as principal factor to the Binomial Negative Model
y_DESeq<-DESeqDataSetFromMatrix(countData=count_matrix, colData=sample_data, design= ~condition)
#set size factor equal 1 (no library size effect) to avoind rlogTransformation does for you
sizeFactors(y_DESeq)<-rep(1,6)

# To the first exploration, we make a log2 conversion over count data
#this data does not have any factor correction
matrix_rlog_DESeq<-rlogTransformation(y_DESeq, blind=TRUE)

head(count_matrix_rlog_DESeq<-assays(matrix_rlog_DESeq)[[1]])

names(count_matrix_rlog_DESeq)<-names(count_matrix)

#takes a look on data distributions

###################################
##PCA over log2 transformed data ##
###################################
# The first exploration must be done over our count_matrix in order to determine if our samples are clustered togheter using count values
# Principal Component Analysis: explore samples separability  
pr_comp_y<-prcomp(t((count_matrix_rlog_DESeq)))
resumen <- summary(pr_comp_y)
labX <- signif((resumen$importance[2,1])*100, 3)
labY <- signif((resumen$importance[2,2])*100, 3)
pdf(file="PCA_DESeq_counts.pdf")
par(mfrow=c(1,2))
plot(x=pr_comp_y$x[,1], y=pr_comp_y$x[,2], col=c(rep("red",4), rep("blue",4)), xlab=paste("PC1 (",labX,"%)", sep=""), ylab=paste("PC2 (",labY,"%)", sep=""))
abline(h=0)
abline(v=0)
title(main="PCA DESeq counts")
biplot(pr_comp_y,cex=c(1,0.001),xlabs=c("T1","T2","T3","A1","A2", "A3"),ylabs=rep("",nrow(count_matrix_rlog_DESeq)),var.axes=FALSE)
title(main="Biplot DESeq counts")
dev.off()
# samples are separated according the condition

colors <- c("#999999", "#E69F00")
colors <- colors[as.numeric(iris$Species)]
scatterplot3d(x=pr_comp_y$x[,1],y=pr_comp_y$x[,2],z=pr_comp_y$x[,3], pch = 16, color = c(rep("red",3), rep("blue",3)), angle = 60)
legend(legend=c("T7","WT"),x=3, fill = c("red","blue"),bty="n", cex=0.8,x.intersp=0.5,y.intersp = 0.8)

#######################################
##Boxplot over log2 transformed data ##
#######################################
# The second exploration is over all counts for each sample in order to determine if samples are comparable
# Boxplot of the raw and normalized read counts
pdf(file="Comparisson_of_boxplot_distributions.pdf")
par(mfrow=c(1,2))
boxplot(count_matrix_rlog_DESeq, col=colors()[c(137, 134, 59,30,132,125)], main="DESeq", names=c("T1","T2","T3","A1","A2","A3"))
  FDP_DESeq<-sapply(1:ncol(count_matrix_rlog_DESeq), function(x){
  fdp<-density(count_matrix_rlog_DESeq[,x])
  return(fdp)})
FDP_DESeq[,1]$call<-"log2counts Densities"
plot(FDP_DESeq[,1], col=colors()[137], xlab="log2DESeqCounts", ylab="Density", xlim=c(0,16), ylim=c(0,0.500),type="l", lwd=4)
lines(FDP_DESeq[,2],col=colors()[134], lwd=4)
lines(FDP_DESeq[,3],col=colors()[59], lwd=4)
lines(FDP_DESeq[,4],col=colors()[30], lwd=4)
lines(FDP_DESeq[,5],col=colors()[132], lwd=4)
lines(FDP_DESeq[,6],col=colors()[125], lwd=4)
legend(legend=c("T1","T2","T3","A1","A2","A3"),x="topright", fill = colors()[c(137, 134, 59,30,132,125)])
dev.off()

################################
## Library size normalization ##
################################
# There are some diferences between density functions. We should remember that existing differences in library size
# In order to estimate those diferences, we use estimateSizeFactors function
y_DESeq<-estimateSizeFactors(y_DESeq)
sizeFactors(y_DESeq)

norm_rlog_DESeqData<-rlogTransformation(y_DESeq, blind=TRUE)
count_matrix_norm_rlog_DESeqData<-assays(norm_rlog_DESeqData)[[1]]
colnames(count_matrix_norm_rlog_DESeqData)<-sample_data$sample
head(count_matrix_norm_rlog_DESeqData)

pdf(file="boxplot_distributions_norm.pdf")
par(mfrow=c(1,2))

boxplot(count_matrix_norm_rlog_DESeqData, col=colors()[c(137, 134, 59,30,132,125)], main="Normalized Rlog Counts", names=c("T1","T2","T3","A1","A2","A3"))
FDP_DESeq<-sapply(1:ncol(count_matrix_norm_rlog_DESeqData), function(x){
  fdp<-density(count_matrix_norm_rlog_DESeqData[,x])
  return(fdp)})
  
FDP_DESeq[,1]$call<-"log2counts Densities"
plot(FDP_DESeq[,1], col=colors()[137], xlab="log2DESeqcounts", ylab="Density", xlim=c(0,16), ylim=c(0,0.50),type="l", lwd=4)
lines(FDP_DESeq[,2],col=colors()[134], lwd=4)
lines(FDP_DESeq[,3],col=colors()[59], lwd=4)
lines(FDP_DESeq[,4],col=colors()[30], lwd=4)
lines(FDP_DESeq[,5],col=colors()[132], lwd=4)
lines(FDP_DESeq[,6],col=colors()[124], lwd=4)
lines(FDP_DESeq[,7],col=colors()[53], lwd=4)
lines(FDP_DESeq[,8],col=colors()[65], lwd=4)

#legend(legend=c("c1","c2","c3","p1","p2","p3","b1","b2","b3"),x="topright", fill = colors()[c(137, 59,101,30,125,128,20,30,45)]) dev.off()
dev.off()
#########################################
## Binomial Negative Model and DE Test ##
#########################################
# Before DE testing is necessary estimate dispersion of the model
y_DESeq<-estimateDispersions(y_DESeq, fitType="local")
par(mfrow=c(1,1)) 
plotDispEsts(y_DESeq)
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates

#Test de DE
et_DESeq<-nbinomWaldTest(y_DESeq)

# Outlier's analysis. For outlier detection, we should determine which of all genes has its pvalue as NA
table(is.na(results(et_DESeq)$pvalue))

decide<-function(dds, p.adj=TRUE, value=0.05, FC_teste=1,cond1,cond2){
  if(class(dds) != "DESeqDataSet") stop("dds must be a DESeqDataSet")
  var<-"padj"
  if(!p.adj) var<-"pvalue"
  dds<-as.data.frame(results(dds,contrast = c("condition",cond1,cond2)))
  de_dds<-rep(0, nrow(dds))
  de_dds<-sapply(1:nrow(dds), function(x){
    
    if(!is.na(dds[x, var]) & dds[x, var] <= value & dds[x,"log2FoldChange"] >= FC_teste ) {
      aux<-1} else{
        if(!is.na(dds[x, var]) & dds[x, var] <= value & dds[x,"log2FoldChange"] <= -(FC_teste)) {
          aux<-(-1)
        }else {
          if(is.na(dds[x, var])) {
            aux<-NA}else aux<-0
        }}
    return(aux)
  })
  return(de_dds)
}

# We consider DE genes those that have a adjusted pavalue lower than 0.01
table(de_DESeq <-decide(et_DESeq, p.adj=TRUE, value=0.05,FC_teste=1,cond1 = "T", cond2="A" ))

head(results(et_DESeq)[order(results(et_DESeq,contrast = c("condition","T7","DAC1"))$padj),])


DEResults<-results(et_DESeq,contrast = c("condition","T","A"))
write.table(DEResults, file = "T-vs-A_expression_analysis.tab", sep= "\t")