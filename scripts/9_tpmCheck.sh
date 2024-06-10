library(tidyr)
#I want to look at the distributions of ortho tpms by species
#hypothalamus
#Lets start by looking at fuk
orthhypfuk <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.fuk.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypfuk) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:8) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Mac
orthhypmac <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.mac.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypmac) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:7) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Cav
orthhypcav <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.cav.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypcav) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:8) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Pap
orthhyppap <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.pap.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhyppap) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:8) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Pan
orthhyppan <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.pan.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhyppan) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:3) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Hom
orthhyphom <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.hom.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhyphom) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:3) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Het
orthhyphet <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.het.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhyphet) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:4) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Ict
orthhypict <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.ict.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypict) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:5) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Mic
orthhypmic <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.mic.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypmic) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:8) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Bos
orthhypbos <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.bos.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypbos) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:6) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Per
orthhypper <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.per.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypper) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:8) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Sus
orthhypsus <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.sus.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypsus) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:8) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Ovi
orthhypovi <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.ovi.hypot.tpm.tsv", header=FALSE)
orthhypovi
as_tibble(orthhypovi) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:7) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Cap
orthhypcap <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.cap.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypcap) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:3) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Mus
orthhypmus <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.mus.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypmus) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:4) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Rat
orthhyrat <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.rat.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhyprat) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:7) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Sor
orthhysor <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.sor.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypsor) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:5) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#complete species hypothalamus
orthhyptot <-cbind(rowMeans(orthhypfuk),rowMeans(orthhypmac),,rowMeans(orthhypcav),rowMeans(orthhyppap),rowMeans(orthhyppan),rowMeans(orthhyphom),rowMeans(orthhyphet),rowMeans(orthhypict),rowMeans(orthhypmic),rowMeans(orthhypbos),rowMeans(orthhypper),rowMeans(orthhypsus),rowMeans(orthhypovi),rowMeans(orthhypcap),rowMeans(orthhypmus),rowMeans(orthhyprat),rowMeans(orthhypsor))
colnames(orthhyptot) <- c("Fuk","Mac","Cav","Pap","Pan","Hom","Het","Ict","Mic","Bos","Per","Sus","Ovi","Cap","Mus","Rat","Sor")
as_tibble(orthhyptot) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:17) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
tpmhyptot <- cbind(rowMeans(fuk.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(mac.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(cav.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(pap.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(pan.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(hom.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(het.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(ict.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(mic.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(bos.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(per.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(sus.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(ovi.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(cap.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(mus.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(rat.hypot.kallisto.txi.tpm.tsv$abundance),rowMeans(sor.hypot.kallisto.txi.tpm.tsv$abundance))
colnames(tpmhyptot) <- c("Fuk","Mac","Cav","Pap","Pan","Hom","Het","Ict","Mic","Bos","Per","Sus","Ovi","Cap","Mus","Rat","Sor")
as_tibble(tpmhyptot) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:17) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  xlab("log2(tpm+1)")+
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
counthyptot <- cbind(rowMeans(fuk.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(mac.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(cav.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(pap.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(pan.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(hom.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(het.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(ict.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(mic.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(bos.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(per.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(sus.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(ovi.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(cap.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(mus.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(rat.hypot.kallisto.txi.tpm.tsv$counts),rowMeans(sor.hypot.kallisto.txi.tpm.tsv$counts))
colnames(counthyptot) <- c("Fuk","Mac","Cav","Pap","Pan","Hom","Het","Ict","Mic","Bos","Per","Sus","Ovi","Cap","Mus","Rat","Sor")
as_tibble(counthyptot) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:17) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)

#UNIQUE HITS
cs_hipp.colsums <- ifelse(cs_hip.count.tsv$counts > 1, 1, 0) %>% as.data.frame
ifelse(cs_hip.count.tsv$counts > 1, 1, 0)
library(tidyverse)
cs_hipp.colsums2 <- colSums(cs_hipp.colsums)
hipp_more_test <- cbind(vis_hip_samples, cs_hipp.colsums2)
ggplot(hipp_more_test,aes(x=RIN,y=cs_hipp.colsums2, color = Stage))+
  geom_point()+
  scale_y_continuous(name="Unique Hits", limits=c(14000,15025))+
  scale_x_continuous(name="RIN")+
  theme_bw()


#######
#######
orthhypfuk <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/tpms/hypot/fuk.hypot.kallisto.txi.tpm.tsv", header=FALSE)
as_tibble(orthhypfuk) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:8) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Mac
orthhypmac <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.mac.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypmac) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:7) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Cav
orthhypcav <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.cav.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypcav) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:8) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Pap
orthhyppap <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.pap.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhyppap) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:8) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Pan
orthhyppan <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.pan.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhyppan) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:3) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Hom
orthhyphom <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.hom.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhyphom) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:3) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Het
orthhyphet <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.het.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhyphet) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:4) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Ict
orthhypict <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.ict.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypict) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:5) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Mic
orthhypmic <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.mic.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypmic) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:8) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Bos
orthhypbos <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.bos.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypbos) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:6) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Per
orthhypper <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.per.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypper) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:8) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Sus
orthhypsus <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.sus.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypsus) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:8) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Ovi
orthhypovi <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.ovi.hypot.tpm.tsv", header=FALSE)
orthhypovi
as_tibble(orthhypovi) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:7) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Cap
orthhypcap <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.cap.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypcap) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:3) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Mus
orthhypmus <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.mus.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypmus) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:4) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Rat
orthhyrat <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.rat.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhyprat) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:7) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#Sor
orthhysor <- read.table("/Users/bill/project_Adult_EVE/data/step1_transcripts2genes/orthotpms/hypot/ortho.sor.hypot.tpm.tsv", header=FALSE)
as_tibble(orthhypsor) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:5) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)