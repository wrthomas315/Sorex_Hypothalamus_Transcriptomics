###Running EVE experiment using only adults in analysis
###purpose is to look for evolutionary shift in expression in shrew lineage
###Using a regrowth adult brain regions (Stage4)

### STEP1: Start by loading libraries you will need
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
library(tibble)
library(tidyr)
library(AnnotationHub)
library(pheatmap)
library(RColorBrewer)


### STEP 2  RESULTS FOR EVE BRANCH SHIFT
library(readr)
ADULTresultsSOR <- read_table2("../analysis/results_branchshift/result", 
                               col_names = FALSE)
ADULTresultsSOR_P = ecdf(ADULTresultsSOR$X2)
plot(ADULTresultsSOR_P)
pchisq(.004903, df=1, lower.tail=FALSE)
vector <- c()
for (i in 1:length(ADULTresultsSOR$X2)){
  vector[i]<- pchisq(ADULTresultsSOR$X2[i], df=1, lower.tail=FALSE)
}
vector2 <- c()
for (i in 1:length(ADULTresultsSOR$X2)){
  vector2[i]<- p.adjust(vector[i], method = "bonferroni", n = length(ADULTresultsSOR$X2))
}
ADULTresultsSOR_pval_matrx <- cbind.data.frame(ADULTresultsSOR$X1, vector,vector2)`
write.table(ADULTresultsSOR_pval_matrx, "../analysis/results_branchshift/ADULTresultsSOR_pval_matrx", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
ADULTresultsSOR_sig <- subset(ADULTresultsSOR_pval_matrx,ADULTresultsSOR_pval_matrx$vector2<.05)
write.table(ADULTresultsSOR_sig, "../analysis/results_branchshift/ADULTresultsSOR_sig", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
#make figure using 
####loading in data make a fun little graph
adult_nopap_hypoth <- read.table("~/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/no_pap/adult_nopap_hypothALV", quote="\"", comment.char="")
adult_nopap_hypot <- adult_nopap_hypoth[,-1]
rownames(adult_nopap_hypot) <- adult_nopap_hypoth[,1]
maybe <- t(adult_nopap_hypot)
maybe <- as.data.frame(maybe)
SPECIES <- c("capHir","oviAri","bosTur","susScr","sorAra","cavApe","hetGla","fukMec","perMan","micOch","musMus","ratNor","speTri","macMul","homSap","panTro")
SPECIES2 <- factor(c(rep("capHir", 3),rep("oviAri",7),rep("bosTur",6),rep("susScr",8),rep("sorAra",5),rep("cavApe",8),rep("hetGla",6),rep("fukMec",8),rep("perMan",8),rep("micOch",8),rep("musMus",4),rep("ratNor",7),rep("speTri",5),rep("macMul",7),rep("homSap",4),rep("panTro",3)))
SPECIES2 <- factor(SPECIES2, levels = SPECIES)
eve4graph <- cbind(maybe,as.data.frame(SPECIES2))
ggplot(eve4graph,aes(x = SPECIES2, y = ENSBTAG00000011392)) +
  geom_boxplot()+
  ylab("TPM") +
  theme_bw()
QUICKSUB  <-subset(newMerged_subset2,newMerged_subset2$vector2 < .05)

#make figure
expDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000011154))
expDat$V2 <- as.numeric(expDat$V2)
for (i in 1:length(tree_key$key1)) {
  expDat$V1[expDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
facet_plot(vegfa_p, panel="Gene", data=expDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=order),size=.6, outlier.shape=NA)+theme_tree2()

###Quick look at intersection between branch shift and stage 4 vs 2 changes
shrew_sigShrewBS <- read_table("project_Adult_EVE/analysis/results_2022_10_04_BS2ALV/shrew_sigShrewBS", 
                               col_names = FALSE)
Reduce(intersect, list(row.names(hyp24resSig),shrew_sigShrewBS$X3))
#CCDC22, OG0011531, ENSBTAG00000013277 ??? ENSBTAG00000011154
plotCounts(dds_hyp_all, gene="CCDC22", intgroup="hyp_full1")
ggplot(eve4graph,aes(x = SPECIES2, y = ENSBTAG00000047586)) +
  geom_boxplot()+
  ylab("TPM") +
  theme_bw()
CCDC22_hyp <-a_hyp[rownames(a_hyp) %in% "CCDC22", ]
CCDC22_hyp <-as.data.frame(CCDC22_hyp)
CCDC22_hyp2 <- cbind(CCDC22_hyp,b_hyp)
plusser_CDC <- CCDC22_hyp2 %>%
  mutate( type=ifelse(b_hyp=="Stg4","Highlighted","Normal")) %>%
  ggplot( aes(y=CCDC22_hyp,x=b_hyp,fill=type, alpha=type))+
    geom_boxplot(size =1)+
    scale_fill_manual(values=c("lightslateblue", "grey")) +
    scale_alpha_manual(values=c(1,0.1)) +
    theme_bw()+ theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank())
expDatCCDC22 <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000013277))
expDatCCDC22$V2 <- as.numeric(expDatCCDC22$V2)
for (i in 1:length(tree_key$key1)) {
  expDatCCDC22$V1[expDatCCDC22$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
expDatCCDC22
CCDC22_facet <- vegfa_p+geom_facet(panel="CCDC22", data=expDatCCDC22, geom_boxplot, mapping = aes(x=V2,group=label,colour=fore),size=.4, outlier.shape=NA)+theme_tree2(legend.position = 'none')
plot_list(CCDC22_facet, plusser_CDC, tag_levels="A", widths=c(.7, .3))

#Look at FAM57B, OG0007926 , ENSBTAG00000005033
#this
#https://www.biorxiv.org/content/10.1101/2021.03.02.433619v1.full
#The 16p11.2 homologs fam57ba and doc2a generate certain brain and body phenotypes
plotCounts(dds_hyp_all, gene="DOC2A", intgroup="hyp_full1")
ggplot(eve4graph,aes(x = SPECIES2, y = ENSBTAG00000005033)) +
  geom_boxplot()+
  ylab("TPM") +
  theme_bw()
FAM57B_hyp <-a_hyp[rownames(a_hyp) %in% "FAM57B", ]
FAM57B_hyp <-as.data.frame(FAM57B_hyp)
FAM57B_hyp2 <- cbind(FAM57B_hyp,b_hyp)
plusser_FAM57B <- FAM57B_hyp2 %>%
  mutate( type=ifelse(b_hyp=="Stg4","Highlighted","Normal")) %>%
  ggplot( aes(y=FAM57B_hyp,x=b_hyp,fill=type, alpha=type))+
  geom_boxplot(size =1)+
  scale_fill_manual(values=c("lightslateblue", "grey")) +
  scale_alpha_manual(values=c(1,0.1)) +
  theme_bw()+ theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank())
expDatFAM57B <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000005033))
expDatFAM57B$V2 <- as.numeric(expDatFAM57B$V2)
for (i in 1:length(tree_key$key1)) {
  expDatFAM57B$V1[expDatFAM57B$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
FAM57B_facet <- vegfa_p+geom_facet(panel="FAM57B", data=expDatFAM57B, geom_boxplot, mapping = aes(x=V2,group=label,colour=fore),size=.6, outlier.shape=NA)+theme_tree2(legend.position = 'none')
plot_list(FAM57B_facet, plusser_FAM57B, tag_levels="A", widths=c(.7, .3))

#GPR3, OG0008043, ENSBTAG00000002505
#https://consensus.app/details/these-results-thus-indicate-gpr3-mediator-cgps-postnatal-tanaka/8800cec514755448849f27b9bb58d27e/
#https://consensus.app/details/findings-reveal-deletion-gpr3-reduced-pathology-mouse-huang/58ef9fb1988353bbbb4a826d0903eb1a/
#https://www.nature.com/articles/srep14953
plotCounts(dds_hyp_all, gene="GPR3", intgroup="hyp_full1")
ggplot(eve4graph,aes(x = SPECIES2, y = ENSBTAG00000002505)) +
  geom_boxplot()+
  ylab("TPM") +
  theme_bw()
GPR3_hyp <-a_hyp[rownames(a_hyp) %in% "GPR3", ]
GPR3_hyp <-as.data.frame(GPR3_hyp)
GPR3_hyp2 <- cbind(GPR3_hyp,b_hyp)
plusser_GPR3 <- GPR3_hyp2 %>%
  mutate( type=ifelse(b_hyp=="Stg4","Highlighted","Normal")) %>%
  ggplot( aes(y=GPR3_hyp,x=b_hyp,fill=type, alpha=type))+
  geom_boxplot(size =1)+
  scale_fill_manual(values=c("lightslateblue", "grey")) +
  scale_alpha_manual(values=c(1,0.1)) +
  theme_bw()+ theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank())
expDatGPR3 <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000002505))
expDatGPR3$V2 <- as.numeric(expDatGPR3$V2)
for (i in 1:length(tree_key$key1)) {
  expDatGPR3$V1[expDatGPR3$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
GPR3_facet <- vegfa_p+geom_facet(panel="GPR3", data=expDatGPR3, geom_boxplot, mapping = aes(x=V2,group=label,colour=fore),size=.6, outlier.shape=NA)+theme_tree2(legend.position = 'none')
plot_list(GPR3_facet, plusser_GPR3, tag_levels="A", widths=c(.7, .3))

#LMX1A, OG0011346, ENSBTAG00000012025
#https://www.pnas.org/doi/10.1073/pnas.1520387113
plotCounts(dds_hyp_all, gene="LMX1B", intgroup="hyp_full1")
ggplot(eve4graph,aes(x = SPECIES2, y = ENSBTAG00000012025)) +
  geom_boxplot()+
  ylab("TPM") +
  theme_bw()
LMX1A_hyp <-a_hyp[rownames(a_hyp) %in% "LMX1A", ]
LMX1A_hyp <-as.data.frame(LMX1A_hyp)
LMX1A_hyp2 <- cbind(LMX1A_hyp,b_hyp)
plusser_LMX1A <- LMX1A_hyp2 %>%
  mutate( type=ifelse(b_hyp=="Stg4","Highlighted","Normal")) %>%
  ggplot( aes(y=LMX1A_hyp,x=b_hyp,fill=type, alpha=type))+
  geom_boxplot(size =1)+
  scale_fill_manual(values=c("lightslateblue", "grey")) +
  scale_alpha_manual(values=c(1,0.1)) +
  theme_bw()+ theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank())
expDatLMX1A <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000012025))
expDatLMX1A$V2 <- as.numeric(expDatLMX1A$V2)
for (i in 1:length(tree_key$key1)) {
  expDatLMX1A$V1[expDatLMX1A$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
LMX1A_facet <- vegfa_p+geom_facet(panel="LMX1A", data=expDatLMX1A, geom_boxplot, mapping = aes(x=V2,group=label,colour=fore),size=.6, outlier.shape=NA)+theme_tree2(legend.position = 'none')
plot_list(LMX1A_facet, plusser_LMX1A, tag_levels="A", widths=c(.7, .3))

#PAQR4, OG0009143, ENSBTAG00000004727
plotCounts(dds_hyp_all, gene="RXFP2", intgroup="hyp_full1")
ggplot(eve4graph,aes(x = SPECIES2, y = ENSBTAG00000004727)) +
  geom_boxplot()+
  ylab("TPM") +
  theme_bw()
PAQR4_hyp <-a_hyp[rownames(a_hyp) %in% "PAQR4", ]
PAQR4_hyp <-as.data.frame(PAQR4_hyp)
PAQR4_hyp2 <- cbind(PAQR4_hyp,b_hyp)
plusser_PAQR4 <- PAQR4_hyp2 %>%
  mutate( type=ifelse(b_hyp=="Stg4","Highlighted","Normal")) %>%
  ggplot( aes(y=PAQR4_hyp,x=b_hyp,fill=type, alpha=type))+
  geom_boxplot(size =1)+
  scale_fill_manual(values=c("lightslateblue", "grey")) +
  scale_alpha_manual(values=c(1,0.1)) +
  theme_bw()+ theme(legend.position="none",axis.title.x=element_blank(),axis.title.y=element_blank())
expDatPAQR4 <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000004727))
expDatPAQR4$V2 <- as.numeric(expDatPAQR4$V2)
for (i in 1:length(tree_key$key1)) {
  expDatPAQR4$V1[expDatPAQR4$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
PAQR4_facet <- vegfa_p+geom_facet(panel="PAQR4", data=expDatPAQR4, geom_boxplot, mapping = aes(x=V2,group=label,colour=fore),size=.6, outlier.shape=NA)+theme_tree2(legend.position = 'none')
plot_list(PAQR4_facet, plusser_PAQR4, tag_levels="A", widths=c(.7, .3))

###plot with all three
vegfa_p+geom_facet(panel="CCDC22", data=expDatCCDC22, geom_boxplot, mapping = aes(x=V2,group=label,colour=order),size=.4, outlier.shape=NA)+
  geom_facet(panel="FAM57B", data=expDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=order),size=.4, outlier.shape=NA)+
  geom_facet(panel="GPR3", data=expDatGPR3, geom_boxplot, mapping = aes(x=V2,group=label,colour=order),size=.4, outlier.shape=NA)+
  geom_facet(panel="LMX1A", data=expDatLMX1A, geom_boxplot, mapping = aes(x=V2,group=label,colour=order),size=.4, outlier.shape=NA)+
  geom_facet(panel="PAQR4", data=expDatPAQR4, geom_boxplot, mapping = aes(x=V2,group=label,colour=order),size=.4, outlier.shape=NA)+
  theme_tree2()

###or with them next to each other
vegfa_p+geom_facet(panel="CCDC22", data=expDatCCDC22, geom_boxplot, mapping = aes(x=V2,group=label,colour=order),size=.4, outlier.shape=NA)+theme_tree2(legend.position = 'none')+plusser
###Diversidication and plasticity
ADULTresultsDIV <- read_table2("../analysis/results_diversity/resultsDIV_beta", 
                               col_names = FALSE)
ADULTresultsDIV_P = ecdf(ADULTresultsDIV$X2)
length(ADULTresultsDIV$X2)
plot(ADULTresultsDIV_P)
pchisq(.004903, df=1, lower.tail=FALSE)
vectorDIV <- c()
for (i in 1:length(ADULTresultsDIV$X2)){
  vectorDIV[i]<- pchisq(ADULTresultsDIV$X2[i], df=1, lower.tail=FALSE)
}
vectorDIV2 <- c()
for (i in 1:length(ADULTresultsDIV$X2)){
  vectorDIV2[i]<- p.adjust(vectorDIV[i], method = "bonferroni", n = length(ADULTresultsDIV$X2))
}
ADULTresultsDIV$X3
ADULTresultsDIV_pval_matrx <- cbind.data.frame(ADULTresultsDIV$X1, ADULTresultsDIV$X3, vectorDIV,vectorDIV2)
ADULTresultsDIV_pval_matrx$`ADULTresultsDIV$X1`
write.table(ADULTresultsDIV_pval_matrx, "../analysis/results_diversity/ADULTresultsDIV_pval_matrx", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
ADULTresultsDIV_sig <- subset(ADULTresultsDIV_pval_matrx,ADULTresultsDIV_pval_matrx$vectorDIV2<.05)
write.table(ADULTresultsDIV_sig, "../analysis/results_diversity/ADULTresultsDIV_sig", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(subset(ADULTresultsDIV_sig,ADULTresultsDIV_sig$`ADULTresultsDIV$X3`>1), "../analysis/results_diversity/ADULTresultsDIV_sigHighBetaPlasticity", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(subset(ADULTresultsDIV_sig,ADULTresultsDIV_sig$`ADULTresultsDIV$X3`<1), "../analysis/results_diversity/ADULTresultsDIV_sigLowBetaSelection", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
ggplot(eve4graph,aes(x = SPECIES2, y = ENSBTAG00000015700)) +
  geom_boxplot()+
  ylab("TPM") +
  theme_bw()

###DROPOUT overlap
dropoutDIV <- read_table2("../analysis/results_dropout/results_dropoutput", 
                               col_names = FALSE)
dropoutDIV_P = ecdf(ADULTresultsDIV$X2)
length(dropoutDIV$X2)
plot(dropoutDIV_P)
pchisq(.004903, df=1, lower.tail=FALSE)
vectordropoutDIV <- c()
for (i in 1:length(dropoutDIV$X2)){
  vectordropoutDIV[i]<- pchisq(dropoutDIV$X2[i], df=1, lower.tail=FALSE)
}
vectordropoutDIV2 <- c()
for (i in 1:length(dropoutDIV$X2)){
  vectordropoutDIV2[i]<- p.adjust(vectordropoutDIV[i], method = "bonferroni", n = length(dropoutDIV$X2))
}
dropoutDIV_pval_matrx <- cbind.data.frame(dropoutDIV$X1, dropoutDIV$X3, vectordropoutDIV,vectordropoutDIV2)
dropoutDIV_pval_matrx$`dropoutDIV$X1`
write.table(dropoutDIV_pval_matrx, "../analysis/results_dropout/dropoutDIV_pval_matrx", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
dropoutDIV_sig <- subset(dropoutDIV_pval_matrx,dropoutDIV_pval_matrx$vectordropoutDIV2<.05)
write.table(dropoutDIV_sig, "../analysis/results_dropout/dropoutDIV_sig", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(subset(dropoutDIV_sig,dropoutDIV_sig$`dropoutDIV$X3` > 1), "../analysis/results_dropout/dropoutDIV_sigHighBetaPlasticity", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(subset(dropoutDIV_sig,dropoutDIV_sig$`dropoutDIV$X3` < 1), "../analysis/results_dropout/dropoutDIV_sigLowBetaSelection", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
subset(dropoutDIV_sig,dropoutDIV_sig$`dropoutDIV$X3` < 1)
dropoutDIV_sig$`dropoutDIV$X3` < 1
tail(dropoutDIV_sig[order(dropoutDIV_sig$`dropoutDIV$X3`),],n=length(dropoutDIV_sig$vectordropoutDIV))
ggplot(eve4graph,aes(x = SPECIES2, y = ENSBTAG00000013723)) +
  geom_boxplot()+
  ylab("TPM") +
  theme_bw()

#check for interesect between BS and DO
Reduce(intersect, list(row.names(hyp24resSig),dropoutDIV_sig$`dropoutDIV$X1`))
shrew_sigShrewDO <- read_table("../analysis/results_dropout/shrew_sigShrewDO", 
                               col_names = FALSE)
#no intersect between branchshift and DO
Reduce(intersect, list(shrew_sigShrewBS$X3,shrew_sigShrewDO$X3))
#only single gene between branchshift and ARHGAP32 OG0006830 ENSG00000134909 ENSBTAG00000015905
Reduce(intersect, list(rownames(hyp24resSig),shrew_sigShrewDO$X3))
ggplot(eve4graph,aes(x = SPECIES2, y = ENSBTAG00000015905)) +
  geom_boxplot()+
  ylab("TPM") +
  theme_bw()
eve_variances <-eve4graph %>% summarise_if(is.numeric, var)
eve_variances$Smallest_Col<-colnames(eve_variances)[apply(eve_variances,1,which.min)]
#looking to plot a housekeeping type gene with no shrew selection or plasticity ENSBTAG00000005364
plotCounts(dds_hyp_all, gene="HMBS", intgroup="hyp_full1")
expDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000005364))
expDat$V2 <- as.numeric(expDat$V2)
for (i in 1:length(tree_key$key1)) {
  expDat$V1[expDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
vegfa_p$data$order
facet_plot(vegfa_p, panel="HMBS", data=expDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=order),size=.6, outlier.shape=NA)+theme_tree2()
#example of BS (VEGFA)
expDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000005339))
expDat$V2 <- as.numeric(expDat$V2)
for (i in 1:length(tree_key$key1)) {
  expDat$V1[expDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
facet_plot(vegfa_p, panel="VEGFA", data=expDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=order),size=.6, outlier.shape=NA)+theme_tree2()
#example of DIVlowbeta EXD2 ENSBTAG00000010052 or ENSBTAG00000000497
expDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000010052))
expDat$V2 <- as.numeric(expDat$V2)
for (i in 1:length(tree_key$key1)) {
  expDat$V1[expDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
facet_plot(vegfa_p, panel="EXD2", data=expDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=order),size=.6, outlier.shape=NA)+theme_tree2()
#example of Divhighbeta (IGF2)
expDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000013066))
expDat$V2 <- as.numeric(expDat$V2)
for (i in 1:length(tree_key$key1)) {
  expDat$V1[expDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
facet_plot(vegfa_p, panel="IGF2", data=expDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=order),size=.6, outlier.shape=NA)+theme_tree2()

###Heatmap of important genes
expDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000001668,eve4graph$ENSBTAG00000013054,eve4graph$ENSBTAG00000010943,eve4graph$ENSBTAG00000005339,eve4graph$ENSBTAG00000020223,eve4graph$ENSBTAG00000044071,eve4graph$ENSBTAG00000021664,eve4graph$ENSBTAG00000016711,eve4graph$ENSBTAG00000009872,eve4graph$ENSBTAG00000047202,eve4graph$ENSBTAG00000012407,eve4graph$ENSBTAG00000011154,eve4graph$ENSBTAG00000012393,eve4graph$ENSBTAG00000004838,eve4graph$ENSBTAG00000004736,eve4graph$ENSBTAG00000004879,eve4graph$ENSBTAG00000002505,eve4graph$ENSBTAG00000002414,eve4graph$ENSBTAG00000002497,eve4graph$ENSBTAG00000012594,eve4graph$ENSBTAG00000000510,eve4graph$ENSBTAG00000010727,eve4graph$ENSBTAG00000009928,eve4graph$ENSBTAG00000007923,eve4graph$ENSBTAG00000000039,eve4graph$ENSBTAG00000020124,eve4graph$ENSBTAG00000013277,eve4graph$ENSBTAG00000008497,eve4graph$ENSBTAG00000013867,eve4graph$ENSBTAG00000016075,eve4graph$ENSBTAG00000046277,eve4graph$ENSBTAG00000006280,eve4graph$ENSBTAG00000004494))
for (i in 1:length(tree_key$key1)) {
  expDat$V1[expDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
expDat$V2 <- as.numeric(expDat$V2)
expDat$V3 <- as.numeric(expDat$V3)
expDat$V4 <- as.numeric(expDat$V4)
expDat$V5 <- as.numeric(expDat$V5)
expDat$V6 <- as.numeric(expDat$V6)
expDat$V7 <- as.numeric(expDat$V7)
expDat$V8 <- as.numeric(expDat$V8)
expDat$V9 <- as.numeric(expDat$V9)
expDat$V10 <- as.numeric(expDat$V10)
expDat$V11 <- as.numeric(expDat$V11)
expDat$V12 <- as.numeric(expDat$V12)
expDat$V13 <- as.numeric(expDat$V13)
expDat$V14 <- as.numeric(expDat$V14)
expDat$V15 <- as.numeric(expDat$V15)
expDat$V16 <- as.numeric(expDat$V16)
expDat$V17 <- as.numeric(expDat$V17)
expDat$V18 <- as.numeric(expDat$V18)
expDat$V19 <- as.numeric(expDat$V19)
expDat$V20 <- as.numeric(expDat$V20)
expDat$V21 <- as.numeric(expDat$V21)
expDat$V22 <- as.numeric(expDat$V22)
expDat$V23 <- as.numeric(expDat$V23)
expDat$V24 <- as.numeric(expDat$V24)
expDat$V25 <- as.numeric(expDat$V25)
expDat$V26 <- as.numeric(expDat$V26)
expDat$V27 <- as.numeric(expDat$V27)
expDat$V28 <- as.numeric(expDat$V28)
expDat$V29 <- as.numeric(expDat$V29)
expDat$V30 <- as.numeric(expDat$V30)
expDat$V31 <- as.numeric(expDat$V31)
expDat$V32 <- as.numeric(expDat$V32)
expDat$V33 <- as.numeric(expDat$V33)
expDat$V34 <- as.numeric(expDat$V34)

mean_expDat <-aggregate(.~V1, data=expDat, mean)
hyp_heat <- as.data.frame(mean_expDat[,2:34])
hyp_heat2 <- sapply(hyp_heat, function(hyp_heat) (hyp_heat-mean(hyp_heat))/sd(hyp_heat))
rownames(hyp_heat2)<-mean_expDat1[,1]
colnames(hyp_heat2)<-c("WNT7A","MFSD2A","SLC22A23","VEGFA","CASQ1","HRH2","TACR2","PPIF","SPHK2","GRIN1","FOXA2","SDC1","AGT","IRX5","GRB2","FOXO4","GPR3","KCNJ10","ELOVL2","MRPS6","ATG101","ATG4D","MLST8","LAMTOR3","SIRT7","HTRA2","CCDC22","RGS14","DLGAP3","EMC8","RGS4","RBFOX3","B4GALNT1")
gheatmap(vegfa_p, hyp_heat2,offset = 2,width = 6,colnames_angle=270,hjust=0,colnames_offset_y =.3)+
  scale_fill_viridis_c(option="A", name="continuous\nvalue")+ ylim(-2, 20)
gheatmap(vegfa_p, hyp_heat2,offset = 2,width = 6,colnames_angle=270,hjust=0,colnames_offset_y =.3)+
  scale_fill_viridis_c(option="A", name="continuous\nvalue")+ ylim(-2, 20)+theme_tree2(legend.position="none")


###Diversification and plasticity
library(readr)
ADULTresultsDIV <- read_table2("../analysis/results_diversity/resultsDIV_beta", 
                               col_names = FALSE)
ADULTresultsDIV
ADULTresultsDIV_P = ecdf(ADULTresultsDIV$X2)
length(ADULTresultsDIV$X2)
plot(ADULTresultsDIV_P)
pchisq(.004903, df=1, lower.tail=FALSE)
vectorDIV <- c()
for (i in 1:length(ADULTresultsDIV$X2)){
  vectorDIV[i]<- pchisq(ADULTresultsDIV$X2[i], df=1, lower.tail=FALSE)
}
vectorDIV2 <- c()
for (i in 1:length(ADULTresultsDIV$X2)){
  vectorDIV2[i]<- p.adjust(vectorDIV[i], method = "bonferroni", n = length(ADULTresultsDIV$X2))
}
ADULTresultsDIV$X3
ADULTresultsDIV_pval_matrx <- cbind.data.frame(ADULTresultsDIV$X1, ADULTresultsDIV$X3, vectorDIV,vectorDIV2)
ADULTresultsDIV_pval_matrx$`ADULTresultsDIV$X1`
write.table(ADULTresultsDIV_pval_matrx, "../analysis/results_diversity/ADULTresultsDIV_pval_matrx", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
ADULTresultsDIV_sig <- subset(ADULTresultsDIV_pval_matrx,ADULTresultsDIV_pval_matrx$vectorDIV2<.05)
write.table(ADULTresultsDIV_sig, "../analysis/results_diversity/ADULTresultsDIV_sig", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(subset(ADULTresultsDIV_sig,ADULTresultsDIV_sig$`ADULTresultsDIV$X3`>1), "../analysis/results_diversity/ADULTresultsDIV_sigHighBetaPlasticity", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
write.table(subset(ADULTresultsDIV_sig,ADULTresultsDIV_sig$`ADULTresultsDIV$X3`<1), "../analysis/results_diversity/ADULTresultsDIV_sigLowBetaSelection", na = "NA", col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
tail(ADULTresultsDIV_sig[order(ADULTresultsDIV_sig$`ADULTresultsDIV$X3`),],n=length(ADULTresultsDIV_sig$vectorDIV))
ADULTresultsDIV_sig
ggplot(eve4graph,aes(x = SPECIES2, y = ENSBTAG00000017886)) +
  geom_boxplot()+
  ylab("TPM") +
  theme_bw()
ENSBTAG00000021328
ggplot(eve4graph,aes(x = SPECIES2, y = ENSBTAG00000015024)) +
  geom_boxplot()+
  ylab("TPM") +
  theme_bw()
plotCounts(dds_hyp_all, gene="IGF2", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="KCNJ5", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="CALD1", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="DUSP27", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="FOXP3", intgroup="hyp_full1")
plotCounts(dds_hyp_all, gene="CABP7", intgroup="hyp_full1")
#justs looking at number of transcripts
length(adult_nopap_hypoth[1,])
sum(adult_nopap_hypoth[1,])
tpm_len <-  c()
for (i in 2:98){
  tpm_len[i-1]<- sum(adult_nopap_hypoth[,i])
}
tpm_len
tpmByspec <- cbind(as.data.frame(tpm_len),as.data.frame(SPECIES2))
tpmByspec

#####Quality Control
#UNIQUE HITS
#adult_nopap_hypot <- adult_nopap_hypoth[,-1]
#rownames(adult_nopap_hypot) <- adult_nopap_hypoth[,1]
#maybe <- t(adult_nopap_hypot)
#maybe <- as.data.frame(maybe)
maybe$ENSBTAG00000017701
maybe.colsums <- ifelse(maybe > 1, 1, 0) %>% as.data.frame
maybe.colsums2 <- rowSums(maybe.colsums)
maybe.colsums2
SPECIES <- c("capHir","oviAri","bosTur","susScr","sorAra","cavApe","hetGla","fukMec","perMan","micOch","musMus","ratNor","speTri","macMul","homSap","panTro")
SPECIES2 <- factor(c(rep("capHir", 3),rep("oviAri",7),rep("bosTur",6),rep("susScr",8),rep("sorAra",5),rep("cavApe",8),rep("hetGla",6),rep("fukMec",8),rep("perMan",8),rep("micOch",8),rep("musMus",4),rep("ratNor",7),rep("speTri",5),rep("macMul",7),rep("homSap",4),rep("panTro",3)))
SPECIES2 <- factor(SPECIES2, levels = SPECIES)
maybe4graph <- cbind(maybe.colsums2,as.data.frame(SPECIES2))
maybe4graph
ggplot(maybe4graph,aes(x=SPECIES2,y=maybe.colsums2))+
  geom_point()+
  theme_bw()
#distributions of data 
adult_nopap_hypoth
totalspechist <-as_tibble(adult_nopap_hypoth) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 2:98)
specRep <- rep(SPECIES2,length(adult_nopap_hypoth$V2))
length(specRep)
length(totalspechist$sample)
totalspechist$spec <- specRep
totalspechist$spec<-as.character(totalspechist$spec)
#tree_key$key1[13] <-"Ictidomys_tridecemlineatus"
for (i in 1:length(tree_key$key1)) {
  totalspechist$spec[totalspechist$spec==tree_key$key2[i]]<-tree_key$key1[i]
}
totalspechist$spec <- as.character(totalspechist$spec)
totalspechist$spec
library(ggridges)
ggrI <-ggplot(dists,aes(x=V3,y=V1,fill=V4))+
  geom_density_ridges()+
  hexpand(.2, direction = -1)+
  #stat_summary(fun.y=median, geom="point", size=2, color="red")+
  theme_bw()
dists <- as.data.frame(cbind(as.character(totalspechist$spec),totalspechist$counts))
dists$V2 <- as.numeric(dists$V2)
dists$V3 <- as.numeric(log2(dists$V2+1))
library(data.table)
setDT(dists)[, V4 := mean(V2), by = V1]
dists[1:96,]
mean(as.numeric(levels(factor(dists$V4))))
gg_abun<-vegfa_p + geom_facet(mapping = aes(x=V3,group=label,fill=V4),geom = geom_density_ridges,data = dists, panel="Abundance",color='grey80', lwd=.3)+
  scale_fill_viridis_c(direction = -1)+theme_tree2()
facet_widths(gg_abun,c(.45,.65))
facet_plot(vegfa_p, panel="Abundances(log2(TPM+1))", data=dists, geom_density_ridges, mapping = aes(x=V3,group=label,fill=V4))+theme_tree2()+scale_fill_viridis_c(direction = -1)
plot_list(vegfa_p+theme_tree2(legend.position="none")+coord_cartesian(clip="off"),ggrI)
vegfa_p+theme2_tree2(legend.position="none")
theme2
####
library(ggimage)
library(phyloseq)
library(ggtree)
library(TDbook)
library(tidyr)
library(rphylopic)
update.packages("rphylopic")

tree <- read.tree("/Users/bill/project_Adult_EVE/analysis/results_2022_10_04_BS2ALV/hypot_ALVAREZ_TreePruned4ggtree.nh")
tree2 <- read.tree("/Users/bill/project_Adult_EVE/analysis/results_2022_10_04_BS2ALV/hypot_ALVAREZ_TreePruned.nh")
#gete <- ggtree(tree,branch.length="none")#+ geom_text(aes(label=tip.label))
#ggtree(tree,branch.length="none")#+ geom_text(aes(label=tip.label))
d <- ggimage::phylopic_uid(tree$tip.label)
d$uid[11:16] <-  "f7d6d04c-73fa-4bf3-8c94-48134e6857b9"
d$uid[10] <-  "f7d6d04c-73fa-4bf3-8c94-48134e6857b9"
d$uid[1:9] <- "f7d6d04c-73fa-4bf3-8c94-48134e6857b9"
d$order <- c(rep("Artiodactyla",4),rep("Eulipotyphla",1),rep("Rodentia",8),rep("Primates",3))
d$labeling <-  d$name
options(phylopic_width=128)
#e  <- as.dataframe(c(""))
gete <-ggtree(tree,branch.length="none") %<+% d + 
  geom_tiplab(aes(image=uid,colour=order), geom="phylopic", offset=0.04, size =.04) #+
  #$geom_text(aes(label=labeling,colour=order), vjust=-.3,size = 1) 
gete
  #geom_tiplab(aes(label=label), offset = .2) + xlim(NA, 7) #+
  #scale_color_viridis_c()
vegfa_p <-flip(gete, 22, 18) %>% flip(5, 19)
vegfa_p
#change expression data to match tree
tree_key <- data.frame(key1=tree$tip.label,
                 key2=tree2$tip.label)
expDat$V1
for (i in 1:length(tree_key$key1)) {
  expDat$V1[expDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
expDat$V1
#
facet_plot(vegfa_p, panel="Expression", data=expDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=order), outlier.shape=NA)+scale_y_continuous()#+xlim(0,400)
?facet_plot
#
expDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000005339))
expDat$V2 <- as.numeric(expDat$V2)
#
facet_plot(vegfa_p, panel="Expression", data=expDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=order), outlier.shape=NA)+theme_tree2()
###
tree_y <-  function(ggtree, data){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  left_join(select(data, label), select(ggtree$data, label, y)) %>%
    pull(y)
}

###BBB WNT7A VEGFA AGT Facet plot generators
VEGFAexpDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000005339))
VEGFAexpDat$V2 <- as.numeric(VEGFAexpDat$V2)
for (i in 1:length(tree_key$key1)) {
  VEGFAexpDat$V1[VEGFAexpDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
facet_plot(vegfa_p, panel="VEGFA", data=VEGFAexpDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=fore),size=.6, outlier.shape=NA)+theme_tree2()
WNT7AexpDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000001668))
WNT7AexpDat$V2 <- as.numeric(WNT7AexpDat$V2)
for (i in 1:length(tree_key$key1)) {
  WNT7AexpDat$V1[WNT7AexpDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
facet_plot(vegfa_p, panel="WNT7A", data=WNT7AexpDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=fore),size=.6, outlier.shape=NA)+theme_tree2()
AGTexpDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000012393))
AGTexpDat$V2 <- as.numeric(AGTexpDat$V2)
for (i in 1:length(tree_key$key1)) {
  AGTexpDat$V1[AGTexpDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
facet_plot(vegfa_p, panel="AGT", data=AGTexpDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=fore),size=.6, outlier.shape=NA)+theme_tree2()
#
BBB_facet <- vegfa_p + geom_facet(geom = geom_boxplot, data = WNT7AexpDat,  panel = 'WNT7A',mapping = aes(x=V2,group=label,colour=fore)) + geom_facet(geom = geom_boxplot, data = VEGFAexpDat,  panel = 'VEGFA',mapping = aes(x=V2,group=label,colour=fore))+theme_tree2(legend.position = 'none')+ geom_facet(geom = geom_boxplot, data = AGTexpDat,  panel = 'AGT',mapping = aes(x=V2,group=label,colour=fore))+theme_tree2(legend.position = 'none')
facet_widths(BBB_facet,c(.1,.3,.3,.3))

###Leptin IRX5 FOXO4 GRB2
IRX5expDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000004838))
IRX5expDat$V2 <- as.numeric(IRX5expDat$V2)
for (i in 1:length(tree_key$key1)) {
  IRX5expDat$V1[IRX5expDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
facet_plot(vegfa_p, panel="IRX5", data=IRX5expDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=fore),size=.6, outlier.shape=NA)+theme_tree2()
FOXO4expDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000004879))
FOXO4expDat$V2 <- as.numeric(FOXO4expDat$V2)
for (i in 1:length(tree_key$key1)) {
  FOXO4expDat$V1[FOXO4expDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
facet_plot(vegfa_p, panel="FOXO4", data=FOXO4expDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=fore),size=.6, outlier.shape=NA)+theme_tree2()
GRB2expDat <- as.data.frame(cbind(as.character(eve4graph$SPECIES2),eve4graph$ENSBTAG00000004736))
GRB2expDat$V2 <- as.numeric(GRB2expDat$V2)
for (i in 1:length(tree_key$key1)) {
  GRB2expDat$V1[GRB2expDat$V1==tree_key$key2[i]]<-tree_key$key1[i]
}
facet_plot(vegfa_p, panel="GRB2", data=AGTexpDat, geom_boxplot, mapping = aes(x=V2,group=label,colour=fore),size=.6, outlier.shape=NA)+theme_tree2()
#
LEP_facet <- vegfa_p + geom_facet(geom = geom_boxplot, data = FOXO4expDat,  panel = 'FOXO4',mapping = aes(x=V2,group=label,colour=fore)) + geom_facet(geom = geom_boxplot, data = IRX5expDat,  panel = 'IRX5',mapping = aes(x=V2,group=label,colour=fore))+theme_tree2(legend.position = 'none')+ geom_facet(geom = geom_boxplot, data = GRB2expDat,  panel = 'GRB2',mapping = aes(x=V2,group=label,colour=fore))+theme_tree2(legend.position = 'none')
facet_widths(LEP_facet,c(.1,.3,.3,.3))

###Average LFC between shrew and species mean for all significant genes
SORsigPlus<- eve4graph[c(ADULTresultsSOR_sig$`ADULTresultsSOR$X1`,"SPECIES2")]
ADULTresultsSOR_sig$`ADULTresultsSOR$X1`
insig <- subset(ADULTresultsSOR_pval_matrx,ADULTresultsSOR_pval_matrx$vector2>.05)[,1]
SORinsiPLUS<- eve4graph[c(insig,"SPECIES2")]
eve_lfc_vec <- c()
xhrew <- c()
xrest <- c()
for (i in 1:1) {
  xhrew <- sum(SORsigPlus[25:29,i])/5
  xrest <- (sum(SORsigPlus[1:24,i])+sum(SORsigPlus[30:97,i]))/92
  eve_lfc_vec[i] <- xhrew/xrest
}
eve_lfc_vec
insi_lfc_vec <- c()
xhrew <- c()
xrest <- c()
for (i in 1:length(insig)) {
  xhrew <- sum(SORinsiPLUS[25:29,i])/5
  xrest <- (sum(SORinsiPLUS[1:24,i])+sum(SORinsiPLUS[30:97,i]))/92
  insi_lfc_vec[i] <- xhrew/xrest
}
mean(eve_lfc_vec)
min(eve_lfc_vec)
max(eve_lfc_vec)
insi_lfc_vec
mean(insi_lfc_vec)
min(insi_lfc_vec)
max(insi_lfc_vec)
eve_lfc_vec
#200 is mac
SORsigPlus[200]