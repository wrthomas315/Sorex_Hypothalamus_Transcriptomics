###Here is the R code for processing the temporal transcriptomics of Dehnel's phenomenon
#Reads were preprocessed on Noctillio server, where they were trimmed using fastp and aligned/quantified using Kalisto
#First I want to get the quantification files onto R and convert transcript abundance to gene abundance
#Then I will create a handful of quant files
#First set up all your libraries
library(stringr)
library(readr)
library(assertr)
library(tximport)
library(GenomicFeatures)
library(DESeq2)
library(ggplot2)
library(regionReport)
library(rhdf5)
library(edgeR)
library(readr)
library(tibble)
library( "genefilter" )
library(gplots)
library(pheatmap)
library(RColorBrewer)
library(UpSetR)
library(TCseq)
library(cluster)
library(EnhancedVolcano)

#Set up your working directories
indir <- "./"
outputdir <- "../analysis/DESeq2/"
##Create a transcript to gene key using the gff file for the shrew
sortxdb <- makeTxDbFromGFF("../data/0_refs/GCF_000181275.1_SorAra2.0_genomic.gff.gz")
k <- keys(sortxdb, keytype="TXNAME")       
sortx2gene <- select(sortxdb, k, "GENEID", "TXNAME")
sortx2gene <- sortx2gene[!duplicated(sortx2gene[,1]),]
sortx2gene <- na.omit(sortx2gene)
#gets rid of the XRs, which are some misc_rnas, and do not apppear to be associated with genes

cs_hyp_samples <- read.table("../data/1_ids/completecycle_hypothalamus.txt", header = T)
cs_hyp_files <-file.path("../data/4_kallisto", cs_hyp_samples$Sample_name, "abundance.tsv")
names(cs_hyp_files) <- paste0("sample_", cs_hyp_samples$Sample_name)
all(file.exists(cs_hyp_files))
#
cs_hyp.count.tsv <- tximport(cs_hyp_files, type = "kallisto", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
cs_hyp.tpm.tsv <- tximport(cs_hyp_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
write.table(cs_hyp.tpm.tsv$abundance, "cs_hyp.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)

###HYPOTHALAMUS
hyp_stages <- factor(c(cs_hyp_samples$Run))
hyp_organs <- factor(c(cs_hyp_samples$Condition))
hyp_full1 <- factor(c(cs_hyp_samples$Run))
cs_hyp_samples
hyp_stages_organ_frame <-cbind(as.data.frame(hyp_stages),as.data.frame(hyp_organs),as.data.frame(hyp_full1))
ggplot(hyp_stages_organ_frame, aes(x=hyp_full1))+
  geom_bar(stat = 'count')+
  theme_bw()
#3 RNA Quality Control
vis_hyp_samples <- read.table("~/CompleteShrew/data/samples/completecycle_hypothalamus_4vis.txt", header = T)
ggplot(vis_hyp_samples,aes(x=RIN,y=Reads_prefilt, color= Stage))+
  geom_point()+
  theme_bw()
#3B removing outlier with high reads
hyp_stages <- factor(c(cs_hyp_samples$Run))
hyp_stages <- hyp_stages[-6]
hyp_organs <- factor(c(cs_hyp_samples$Condition))
hyp_organs <- hyp_organs[-6]
hyp_full1 <- factor(c(cs_hyp_samples$Run))
hyp_full1 <- hyp_full1[-6]
hyp_stages_organ_frame <-cbind(as.data.frame(hyp_stages),as.data.frame(hyp_organs),as.data.frame(hyp_full1))
#everything looks pretty good
#Make DESeq onject
colnames(cs_hyp.count.tsv$counts) <- c("Stg1_1","Stg1_2","Stg1_3","Stg1_4","Stg1_5","Stg2_1","Stg2_2","Stg2_3","Stg2_4","Stg3_1","Stg3_2","Stg3_3","Stg3_4","Stg3_5","Stg4_1","Stg4_2","Stg4_3","Stg4_4","Stg4_5","Stg5_1","Stg5_2","Stg5_3","Stg5_4","Stg5_5")
dds_hyp_all <- DESeqDataSetFromMatrix(round(subset(cs_hyp.count.tsv$counts, select=-Stg2_4)), DataFrame(hyp_stages_organ_frame), ~ hyp_full1)
mcols(dds_hyp_all) <- cbind(mcols(dds_hyp_all), row.names(cs_hyp.count.tsv$counts))
rownames(dds_hyp_all) <- row.names(cs_hyp.count.tsv$counts)
#Run DESeq
dds_hyp_all <- DESeq(dds_hyp_all)
vst_dds_hyp_all <- vst(dds_hyp_all)
#Generate PCA
pcaData_hyp_all<- plotPCA(vst_dds_hyp_all,intgroup=c("hyp_stages","hyp_organs"), ntop=300, returnData=TRUE)
ggplot(pcaData_hyp_all, aes(x = PC1, y = PC2, color = factor(hyp_stages))) +
  geom_point(size=2)+
  theme_bw()r
#Look at distances if desired
hypsampleDists <- dist(t(assay(vst_dds_hyp_all)))
hypsampleDistMatrix <- as.matrix(hypsampleDists)
colnames(hypsampleDistMatrix) <- NULL
##make the heatmap
pheatmap(hypsampleDistMatrix, clustering_distance_rows=hypsampleDists,
         clustering_distance_cols = hypsampleDists, color = colorRampPalette(rev(brewer.pal(n = 9, name ="Reds")))(255))
hyptopVarGenes <- head( order( rowVars( assay(vst_dds_hyp_all) ), decreasing=TRUE ), 40 )
heatmap.2( assay(vst_dds_hyp_all)[ hyptopVarGenes, ], scale="row", 
           trace="none", dendrogram="column", 
           col = colorRampPalette( rev(brewer.pal(9, "RdBu")) )(255),
)

###Quick brain mass phenotype
brainmass <- c(0.2990,	0.2590,	0.2650,	0.3150,	0.2620, 0.2433,0.2305,0.2316,0.2130,0.2002,0.2482,0.2410,0.3786,0.2731,0.2727,0.2627,0.2710,0.2347,0.2614,0.2527, 0.2488)
brainstage <- c(rep("Stage1",5),rep("Stage2",3),rep("Stage3",4),rep("Stage4",4),rep("Stage5",5))
brainmass_df <- data.frame(brainmass,brainstage)
ggplot(brainmass_df,aes(x=brainstage,y=brainmass))+
  geom_boxplot()+
  scale_y_continuous(name="brainMass", limits=c(0.1,.4))+
  theme_bw()

#hypothalamus DESeq analysis
#hypothalamus 2-4
lfc <- 0
hyp_full1
hyp24res <- results(dds_hyp_all, contrast = c("hyp_full1","Stage4","Stage2"))
hyp24resSig <- subset(hyp24res,hyp24res$padj<.05)
hyp24resSigLog <- subset(hyp24resSig,abs(hyp24resSig$log2FoldChange)>=lfc)
hyp24up <- subset(hyp24resSigLog,(hyp24resSigLog$log2FoldChange)>=0)
hyp24down <- subset(hyp24resSigLog,(hyp24resSigLog$log2FoldChange)<=0)
DESeq2::plotMA(hyp24res, ylim = c (-4,4)); #drawLines()
as.data.frame(rownames(hyp24down))
rownames(hyp24down)
as.data.frame(rownames(hyp24up))
Reduce(intersect, list(row.names(hyp24resSig),row.names(hip24resSig),row.names(cor13resSig)))
length(hyp24up$log2FoldChange)
length(hyp24down$log2FoldChange)
#write.table(as.data.frame(rownames(hyp24down)), file='../analysis/DESeq/hyp24down', quote=FALSE, sep='\t')
#write.table(as.data.frame(rownames(hyp24up)), file='../analysis/DESeq/hyp24up', quote=FALSE, sep='\t')
subset(hyp24resSigLog,(rownames(hyp24resSig))=="CACNA1")
plotCounts(dds_hyp_all, gene="FOXO4", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="BCL2L1", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="FOS", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="NFKBIA", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="CTSK", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="FANCD2", intgroup="hyp_full1")

###RUN LISTS THROUGH KEGG online, can then make figures below
hyp24_kegg <- read_delim("../analysis/DESeq/KEGG_workable.tsv", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
hyp24_kegg
hyp24_kegg <- subset(hyp24_kegg, PValue <= .05)
library(dplyr)
library(ggplot2)

hyp24_kegg$Term <- factor(hyp24_kegg$Term, levels = hyp24_kegg$Term[order(-hyp24_kegg$Count)])
ggplot(hyp24_kegg , aes(Term,Count,color = Direction, size = -log(PValue)))+
  geom_point()+
  scale_size(range = c(3, 10))+
  coord_flip()+
  theme_bw()
ggplot(hyp24_kegg, aes(x=Term,y=Count, fill=Direction))+
  geom_bar(stat = 'identity', linewidth=1)+
  coord_flip()+
  scale_fill_manual(values=c("#FE4F49","#5059FF","#5059FF","#5059FF", "#5059FF"))+
  scale_y_continuous(breaks = seq(-5, 5, len = 6))+
  theme_bw()

#Apoptosis gene graph
a_hyp <- DESeq2::counts(dds_hyp_all, normalized=TRUE)
b_hyp <- c(rep("Stg1",5),rep("Stg2",3),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5))
NGF_hyp <-a_hyp[rownames(a_hyp) %in% "NGF", ]
NGF_hyp <-as.data.frame(NGF_hyp)
NGF_hyp2 <- cbind(NGF_hyp,b_hyp)
ggplot(NGF_hyp2,aes(y=NGF_hyp,x=b_hyp))+
  geom_point()+
  scale_y_continuous(trans = 'log10')+
  theme_bw()
NFKBIA_hyp <-a_hyp[rownames(a_hyp) %in% "NFKBIA", ]
NFKBIA_hyp <-as.data.frame(NFKBIA_hyp)
NFKBIA_hyp2 <- as.data.frame(cbind(NFKBIA_hyp,b_hyp))
ggplot(NFKBIA_hyp2,aes(y=NFKBIA_hyp,x=b_hyp))+
  geom_point()+
  scale_y_continuous(trans = 'log10')+
  theme_bw()
BCL2L1_hyp <-a_hyp[rownames(a_hyp) %in% "BCL2L1", ]
BCL2L1_hyp <-as.data.frame(BCL2L1_hyp)
BCL2L1_hyp2 <- as.data.frame(cbind(BCL2L1_hyp,b_hyp))
ggplot(BCL2L1_hyp2,aes(y=BCL2L1_hyp,x=b_hyp))+
  geom_point()+
  scale_y_continuous(trans = 'log10')+
  theme_bw()
FOS_hyp <-a_hyp[rownames(a_hyp) %in% "FOS", ]
FOS_hyp <-as.data.frame(FOS_hyp)
FOS_hyp2 <- as.data.frame(cbind(FOS_hyp,b_hyp))
ggplot(FOS_hyp2,aes(y=FOS_hyp,x=b_hyp))+
  geom_point()+
  scale_y_continuous(trans = 'log10')+
  theme_bw()
CTSK_hyp <-a_hyp[rownames(a_hyp) %in% "CTSK", ]
CTSK_hyp <-as.data.frame(CTSK_hyp)
CTSK_hyp2 <- as.data.frame(cbind(CTSK_hyp,b_hyp))
ggplot(CTSK_hyp2,aes(y=CTSK_hyp,x=b_hyp))+
  geom_point()+
  scale_y_continuous(trans = 'log10')+
  theme_bw()
#try it with z-score on same graph in red
a_hyp <- DESeq2::counts(dds_hyp_all, normalized=TRUE)
NGF_hyp <-a_hyp[rownames(a_hyp) %in% "NGF", ]
BCL2L1_hyp <-a_hyp[rownames(a_hyp) %in% "BCL2L1", ]
NFKBIA_hyp <-a_hyp[rownames(a_hyp) %in% "NFKBIA", ]
FOS_hyp <-a_hyp[rownames(a_hyp) %in% "FOS", ]
CTSK_hyp <-a_hyp[rownames(a_hyp) %in% "CTSK", ]
c_hyp  <- rep(c(c(rep("Stg1",5),rep("Stg2",3),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5))),5)
antiapop <- rbind(NGF_hyp,NFKBIA_hyp,BCL2L1_hyp,FOS_hyp,CTSK_hyp)
antiapop
z <- matrix(0,nrow(antiapop),ncol(antiapop))
z
for (i in 1:nrow(antiapop)) {
  for (j in 1:ncol(antiapop)) {
    z[i,j] <- (antiapop[i,j]-mean(antiapop[i,]))/sd(antiapop[i,])
  }
}
z
zz  <- matrix(0,nrow(z),5)
for (i in 1:nrow(z)) {
  zz[i,1] <- sum(z[i,1],z[i,2],z[i,3],z[i,4],z[i,5])/5
  zz[i,2] <- sum(z[i,6],z[i,7],z[i,8])/3
  zz[i,3] <- sum(z[i,9],z[i,10],z[i,11],z[i,12],z[i,13])/5
  zz[i,4] <- sum(z[i,14],z[i,15],z[i,16],z[i,17],z[i,18])/5
  zz[i,5] <- sum(z[i,19],z[i,20],z[i,21],z[i,22],z[i,23])/5
}
zz
real_antiapop <- as.matrix(rbind(as.matrix(zz[1,]),as.matrix(zz[2,]),as.matrix(zz[3,]),as.matrix(zz[4,]),as.matrix(zz[5,])))
close <-as.matrix(c(rep("NGF",5),rep("NFKBIA",5),rep("BCL2L1",5),rep("FOS",5),rep("CTSK",5)))
vclose<-as.matrix(c(rep("Hypothalamus",25)))
soclose<- rep(c("Stg1","Stg2","Stg3","Stg4","Stg5"),5)
realreal <-  cbind(as.data.frame(real_antiapop),as.data.frame(close),as.data.frame(soclose),as.data.frame(vclose))
colnames(realreal) <- c("Expr","Gene","Stage","Tissue")
realreal
#
ggplot(realreal) +
  geom_point(aes(x = Stage, colour = "#FE4F49", y =Expr,shape=Gene),size=5)+
  geom_line(aes(x = Stage, colour = "#FE4F49", y =Expr,group=Gene),size=1)+
  scale_shape_manual(values=c(15,16,17,18,3))+
  theme_bw()

#GABAnergic Figure
a_hyp <- DESeq2::counts(dds_hyp_all, normalized=TRUE)
CACNA1D_hyp <-a_hyp[rownames(a_hyp) %in% "CACNA1D", ]
GNB3_hyp <-a_hyp[rownames(a_hyp) %in% "GNB3", ]
GABRQ_hyp <-a_hyp[rownames(a_hyp) %in% "GABRQ", ]
GABRE_hyp <-a_hyp[rownames(a_hyp) %in% "GABRE", ]
c_hyp  <- rep(c(c(rep("Stg1",5),rep("Stg2",3),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5))),4)
antiapop <- rbind(CACNA1D_hyp,GNB3_hyp,GABRQ_hyp,GABRE_hyp)
antiapop
z <- matrix(0,nrow(antiapop),ncol(antiapop))
z
for (i in 1:nrow(antiapop)) {
  for (j in 1:ncol(antiapop)) {
    z[i,j] <- (antiapop[i,j]-mean(antiapop[i,]))/sd(antiapop[i,])
  }
}
z
zz  <- matrix(0,nrow(z),5)
for (i in 1:nrow(z)) {
  zz[i,1] <- sum(z[i,1],z[i,2],z[i,3],z[i,4],z[i,5])/5
  zz[i,2] <- sum(z[i,6],z[i,7],z[i,8])/3
  zz[i,3] <- sum(z[i,9],z[i,10],z[i,11],z[i,12],z[i,13])/5
  zz[i,4] <- sum(z[i,14],z[i,15],z[i,16],z[i,17],z[i,18])/5
  zz[i,5] <- sum(z[i,19],z[i,20],z[i,21],z[i,22],z[i,23])/5
}
zz
real_antiapop <- as.matrix(rbind(as.matrix(zz[1,]),as.matrix(zz[2,]),as.matrix(zz[3,]),as.matrix(zz[4,])))
close <-as.matrix(c(rep("CACNA1D",5),rep("GNB3",5),rep("GABRQ",5),rep("GABRE",5)))
vclose<-as.matrix(c(rep("Hypothalamus",20)))
soclose<- rep(c("Stg1","Stg2","Stg3","Stg4","Stg5"),4)
realreal <-  cbind(as.data.frame(real_antiapop),as.data.frame(close),as.data.frame(soclose),as.data.frame(vclose))
colnames(realreal) <- c("Expr","Gene","Stage","Tissue")
realreal
#
ggplot(realreal) +
  geom_point(aes(x = Stage, colour = "#5059FF", y =Expr,shape=Gene),size=5)+
  geom_line(aes(x = Stage, colour = "#5059FF", y =Expr,group=Gene),size=1)+
  scale_color_manual(values=c("#5059FF","#5059FF","#5059FF","#5059FF"))+
  scale_shape_manual(values=c(15,16,17,18))+
  theme_bw()

####MAKE VOLCANO PLOT FOR HYPTHALAMUS
hyp24resX <-  results(dds_hyp_all, contrast = c("hyp_full1","Stage4","Stage2"))
for (i in 1:length(hyp24resX$padj)) {
  if  (hyp24resX$padj[i]<1e-15 & !is.na (hyp24resX$padj[i])) {
    hyp24resX$padj[i] <- 1e-15
  }
  if (hyp24resX$log2FoldChange[i]>6 & !is.na (hyp24resX$log2FoldChange[i])) {
    hyp24resX$log2FoldChange[i] <- 6
  }
  if (hyp24resX$log2FoldChange[i]< -6 & !is.na (hyp24resX$log2FoldChange[i])) {
    hyp24resX$log2FoldChange[i] <- -6
  }
}

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements
keyvals.shape <- ifelse(
  abs(hyp24resX$log2FoldChange) == 6, 17,
  ifelse(hyp24resX$padj==1e-15, 17,
         16))
keyvals.shape[is.na(keyvals.shape)] <- 1
names(keyvals.shape)[keyvals.shape == 16] <- '<log1e-15'
names(keyvals.shape)[keyvals.shape == 17] <- 'Off-Graph'
###
keyvals <- ifelse(
  hyp24resX$padj > 0.05, 'grey',
  ifelse(hyp24resX$log2FoldChange <= -1.58, 'red',
         ifelse(hyp24resX$log2FoldChange >= 1.58, 'blue',
                ifelse(hyp24resX$log2FoldChange >= 0, 'lightblue',
                       'pink'))))
keyvals
keyvals[is.na(keyvals)] <- 'grey'
str(keyvals)
names(keyvals)[keyvals == 'grey'] <- 'NotSig'
names(keyvals)[keyvals == 'pink'] <- 'DownRegulated'
names(keyvals)[keyvals == 'blue'] <- 'Upregulated High Effect'
names(keyvals)[keyvals == 'lightblue'] <- 'Upregulated'
names(keyvals)[keyvals == 'red'] <- 'Downregulated High Effect'
length(keyvals)
length(keyvals.shape)
####
hyp24res[rownames(hyp24res) %in% "CTSK",]
hyp24resX[rownames(hyp24resX) %in% "CTSK",]
EnhancedVolcano(hyp24resX,
                lab = rownames(hyp24resX),
                xlim=c(-6 ,6),
                ylim=c(0,4.5),
                x = 'log2FoldChange',
                y = 'padj',
                shapeCustom = keyvals.shape,
                selectLab = c('BCL2L1','NGF','FOS','NFKBIA','CTSK','GABRQ','GABRE','CACNA1D','GNB3'),
                #selectLab = c("SIRT1",'CREBBP','FOXO','G6PC','PPARA','RXRG','FOXO1','RXRA','PER1','PPARD','MBP','PLP','FABP5','SLC27A2','FGF21'),
                #selectLab = c('SIRT1','CREBBP','FOXO3','G6PC','PPARA','RXRG','FOXO1','RXRA','PER1','PPARD','MBP','PLP','FABP5','SLC27A2','FGF21','SIRT1'),
                pCutoff = .05,
                FCcutoff = 1.58,
                colCustom = keyvals,
                pointSize = 4.0,
                legendPosition = 'none',
                drawConnectors = TRUE,
                gridlines.major  = FALSE,
                widthConnectors = 0.5)

###Timeseq Hypothalamus with ABS LFC 0.5 change)
con_hyp3 <- c(rep("Stg1",5),rep("Stg2",3),rep("Stg3",5),rep("Stg4",5),rep("Stg5",5))
hyp3timeseq <- data.frame(sampleID = 1:23, group = c(1, 1, 1,1,1,2, 2, 2, 3, 3, 3,3,3,4, 4, 4, 4,4,5, 5, 5, 5,5),
                          timepoint = con_hyp3)
gf <- data.frame(chr = c(rep('chr1', 15296), rep('chr2', 2000), rep('chr4', 2000)),
                 start = rep(100, 19296),
                 end = rep(500, 19296),
                 id = row.names(cs_hyp.tpm.tsv$abundance))
library(TCseq)
hyp3tca <- TCA(design = hyp3timeseq, counts = round(DESeq2::counts(dds_hyp_all, normalized=TRUE)), gf)
hyp3tca <- DBanalysis(hyp3tca, filter.type = "raw", filter.value = 10, samplePassfilter = 2)
hyp3tca <- timecourseTable(hyp3tca, value = "expression",  lib.norm = FALSE, filter = TRUE,abs.fold = 0.5)
hyp3_t <- tcTable(hyp3tca)
#Gene counts: 19k to 14377 to 789
#Determine how many clusters using timeclust2 with clusGap
clusGap(hyp3_t,
        FUNcluster = timeclust2,
        algo = 'cm',
        K.max = 20,
        B = 20)
##12 for Hypothalamus .5LFC
chyp3tca <- timeclust(hyp3tca, algo = "cm", k = 12, standardize = TRUE)
XXXchyp3_px <-timeclustplot(chyp3tca, value = "z-score(PRKM)", cols = 2,cl.color = "gray50",membership.color= hcl.colors(30, "Viridis",rev = TRUE))
hyp3cxxx<-clustResults(chyp3tca)
hyp_clusters<-as.data.frame(hyp3cxxx@cluster)
write.table(hyp3cxxx@membership, "../analysis/TimeSeq/MemberShip", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
#print individual cluster figures
rownames(subset(hyp_clusters,`hyp3cxxx@cluster` == 2))
hyp3cqqq <- 1
print(XXXchyp3_px[[hyp3cqqq]])

#########################
#KEGG for TCSEQ figures
library(readr)
total_workaable <- read_delim("../analysis/TimeSeq/DAVID_workable2.txt", 
                              delim = "\t", escape_double = FALSE, 
                              trim_ws = TRUE)
total_workaable <- subset(total_workaable, PValue <= .05)
library(dplyr)
library(ggplot2)

total_workaable$Term <- factor(total_workaable$Term, levels = total_workaable$Term[order(-total_workaable$PValue)])
ggplot(total_workaable , aes(Term,-log(PValue),color = Count))+
  geom_point(size=9)+
  scale_color_viridis_c(option="A", name="continuous\nvalue")+
  coord_flip()+
  theme_bw()


