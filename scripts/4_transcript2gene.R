###EVE experiment with only adults in analysis
###Here just taking our transcripts and converting them to gene level expression
###purpose is to look for evolutionary shift in expression in shrew lineague
###Using a regrowth adult brain regions (Stage4)

###Some of this is copied over from mini_eve.R
###Updating and cleaning it here with adult tissues

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

#STEP2: Then make your txdata bases to conver kalisto transcript outputs to gene level
#Sorex araneus(shrew) txdb 1/17
sortxdb <- makeTxDbFromGFF("../data/0_refs/GCF_000181275.1_SorAra2.0_genomic.gff.gz")
keytypes(sortxdb)
sor_k <- keys(sortxdb, keytype="TXNAME")       
sortx2gene <- select(sortxdb, sor_k, "CDSNAME", "TXNAME")
sortx2gene <- sortx2gene[!duplicated(sortx2gene[,1]),]
sortx2gene <- na.omit(sortx2gene)
#Homo sapiens (human) 2/17
homtxdb <- makeTxDbFromGFF("../data/0_refs/Homo_sapiens.GRCh38.104.gtf.gz")
hom_k <- keys(homtxdb, keytype="TXNAME")       
homtx2gene <- select(homtxdb, hom_k, "GENEID", "TXNAME")
#Bos taurus (cow) 3/17
bostxdb <- makeTxDbFromGFF("../data/0_refs/Bos_taurus.ARS-UCD1.2.104.gtf.gz")
bos_k <- keys(bostxdb, keytype="TXNAME")       
bostx2gene <- select(bostxdb, bos_k, "GENEID", "TXNAME")
#Cavia porcellus (guinea pig) 4/17
cavtxdb <- makeTxDbFromGFF("../data/0_refs/Cavia_porcellus.Cavpor3.0.104.gtf.gz")
cav_k <- keys(cavtxdb, keytype="TXNAME")       
cavtx2gene <- select(cavtxdb, cav_k, "GENEID", "TXNAME")
#Rattus norvegicus (rat) 5/17
rattxdb <- makeTxDbFromGFF("../data/0_refs/Rattus_norvegicus.Rnor_6.0.104.gtf.gz")
rat_k <- keys(rattxdb, keytype="TXNAME")       
rattx2gene <- select(rattxdb, rat_k, "GENEID", "TXNAME")
#Sus scrofa (pig) 6/17
sustxdb <- makeTxDbFromGFF("../data/0_refs/Sus_scrofa.Sscrofa11.1.104.gtf.gz")
sus_k <- keys(sustxdb, keytype="TXNAME")       
sustx2gene <- select(sustxdb, sus_k, "GENEID", "TXNAME")
#Macaca mulatta (rhesus macaque) 7/17
mactxdb <- makeTxDbFromGFF("../data/0_refs/Macaca_mulatta.Mmul_10.104.gtf.gz")
mac_k <- keys(mactxdb, keytype="TXNAME")       
mactx2gene <- select(mactxdb, mac_k, "GENEID", "TXNAME")
#Mus musculus (mouse) 8/17
mustxdb <- makeTxDbFromGFF("../data/0_refs/Mus_musculus.GRCm39.104.gtf.gz")
mus_k <- keys(mustxdb, keytype="TXNAME")       
mustx2gene <- select(mustxdb, mus_k, "GENEID", "TXNAME")
#Papio anubis (olive baboon) 9/17
paptxdb <- makeTxDbFromGFF("../data/0_refs/Papio_anubis.Panu_3.0.104.gtf.gz")
pap_k <- keys(paptxdb, keytype="TXNAME")       
paptx2gene <- select(paptxdb, pap_k, "GENEID", "TXNAME")
#Pan troglodytes (chimp) 10/17
pantxdb <- makeTxDbFromGFF("../data/0_refs/Pan_troglodytes.Pan_tro_3.0.104.gtf.gz")
pan_k <- keys(pantxdb, keytype="TXNAME")       
pantx2gene <- select(pantxdb, pan_k, "GENEID", "TXNAME")
#Fukoyms damarensis (molerat) 11/17
fuktxdb <- makeTxDbFromGFF("../data/0_refs/Fukomys_damarensis.DMR_v1.0.104.gtf.gz")
fuk_k <- keys(fuktxdb, keytype="TXNAME")       
fuktx2gene <- select(fuktxdb, fuk_k, "GENEID", "TXNAME")
#Capra hircus (goat) 12/17
captxdb <- makeTxDbFromGFF("../data/0_refs/Capra_hircus.ARS1.104.gtf.gz")
cap_k <- keys(captxdb, keytype="TXNAME")       
captx2gene <- select(captxdb, cap_k, "GENEID", "TXNAME")
#Heterocephalus glaber (male) (naked mole rat) 13/17
hettxdb <- makeTxDbFromGFF("../data/0_refs/Heterocephalus_glaber_male.HetGla_1.0.104.gtf.gz")
het_k <- keys(hettxdb, keytype="TXNAME")       
hettx2gene <- select(hettxdb, het_k, "GENEID", "TXNAME")
#Microtus ochrogaster (prairie vole) 14/17
mictxdb <- makeTxDbFromGFF("../data/0_refs/Microtus_ochrogaster.MicOch1.0.104.gtf.gz")
mic_k <- keys(mictxdb, keytype="TXNAME")       
mictx2gene <- select(mictxdb, mic_k, "GENEID", "TXNAME")
#Ovis aries (sheep) 15/17
ovitxdb <- makeTxDbFromGFF("../data/0_refs/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.104.gtf.gz")
ovi_k <- keys(ovitxdb, keytype="TXNAME")       
ovitx2gene <- select(ovitxdb, ovi_k, "GENEID", "TXNAME")
#Ictidomystridecemlineatus (squirrel) 16/17
icttxdb <- makeTxDbFromGFF("../data/0_refs/Ictidomys_tridecemlineatus.SpeTri2.0.104.gtf.gz")
ict_k <- keys(icttxdb, keytype="TXNAME")       
icttx2gene <- select(icttxdb, ict_k, "GENEID", "TXNAME")
#Peromyscus maniculatus (deer mouse) 17/17
pertxdb <- makeTxDbFromGFF("../data/0_refs/Peromyscus_maniculatus_bairdii.HU_Pman_2.1.104.gtf.gz")
per_k <- keys(pertxdb, keytype="TXNAME")       
pertx2gene <- select(pertxdb, per_k, "GENEID", "TXNAME")

###STEP3: Read in your samples and make sure files exist
###For hypothalamus
#now lets upload them samples in this script
cav_hy_samples <- read.table("../data/1_ids/specs/cav_ids.txt", header = T)
cap_hy_samples <- read.table("../data/1_ids/specs/cap_ids.txt", header = T)
adult_hyp_samples <- read.table("../data/1_ids/specs/adult_hyp_samples.txt", header = T)
ict_hy_samples <- read.table("../data/1_ids/specs/ict_ids.txt", header = T)
pan_hy_samples <- read.table("../data/1_ids/specs/pan_ids.txt", header = T)
sus_hy_samples <- read.table("../data/1_ids/specs/sus_ids.txt", header = T)
mus_hy_samples <- read.table("../data/1_ids/specs/mus_ids.txt", header = T)
bos_hy_samples <- read.table("../data/1_ids/specs/bos_ids.txt", header = T)
fuk_hy_samples <- read.table("../data/1_ids/specs/fuk_ids.txt", header = T)
het_hy_samples <- read.table("../data/1_ids/specs/het_ids.txt", header = T)
hom_hy_samples <- read.table("../data/1_ids/specs/hom_ids.txt", header = T)
mac_hy_samples <- read.table("../data/1_ids/specs/mac_ids.txt", header = T)
mic_hy_samples <- read.table("../data/1_ids/specs/mic_ids.txt", header = T)
ovi_hy_samples <- read.table("../data/1_ids/specs/ovi_ids.txt", header = T)
pap_hy_samples <- read.table("../data/1_ids/specs/pap_ids.txt", header = T)
per_hy_samples <- read.table("../data/1_ids/specs/per_ids.txt", header = T)
rat_hy_samples <- read.table("../data/1_ids/specs/rat_ids.txt", header = T)
adult_hyp_samples
#check to see if they exist
cav_hy_files <-file.path("../data/4_kallisto", cav_hy_samples$Run, "abundance.tsv")
cap_hy_files <-file.path("../data/4_kallisto", cap_hy_samples$Run, "abundance.tsv")
sor_hy_files <-file.path("../data/4_kallisto", adult_hyp_samples$Run, "abundance.tsv")
ict_hy_files <-file.path("../data/4_kallisto", ict_hy_samples$Run, "abundance.tsv")
pan_hy_files <-file.path("../data/4_kallisto", pan_hy_samples$Run, "abundance.tsv")
sus_hy_files <-file.path("../data/4_kallisto", sus_hy_samples$Run, "abundance.tsv")
mus_hy_files <-file.path("../data/4_kallisto", mus_hy_samples$Run, "abundance.tsv")
fuk_hy_files <-file.path("../data/4_kallisto", fuk_hy_samples$Run, "abundance.tsv")
bos_hy_files <-file.path("../data/4_kallisto", bos_hy_samples$Run, "abundance.tsv")
het_hy_files <-file.path("../data/4_kallisto", het_hy_samples$Run, "abundance.tsv")
hom_hy_files <-file.path("../data/4_kallisto", hom_hy_samples$Run, "abundance.tsv")
mac_hy_files <-file.path("../data/4_kallisto", mac_hy_samples$Run, "abundance.tsv")
mic_hy_files <-file.path("../data/4_kallisto", mic_hy_samples$Run, "abundance.tsv")
ovi_hy_files <-file.path("../data/4_kallisto", ovi_hy_samples$Run, "abundance.tsv")
pap_hy_files <-file.path("../data/4_kallisto", pap_hy_samples$Run, "abundance.tsv")
per_hy_files <-file.path("../data/4_kallisto", per_hy_samples$Run, "abundance.tsv")
rat_hy_files <-file.path("../data/4_kallisto", rat_hy_samples$Run, "abundance.tsv")
names(cav_hy_files) <- paste0("sample_", cav_hy_samples$Run)
names(cap_hy_files) <- paste0("sample_", cap_hy_samples$Run)
names(sor_hy_files) <- paste0("sample_", adult_hyp_samples$Run)
names(ict_hy_files) <- paste0("sample_", ict_hy_samples$Run)
names(pan_hy_files) <- paste0("sample_", pan_hy_samples$Run)
names(sus_hy_files) <- paste0("sample_", sus_hy_samples$Run)
names(mus_hy_files) <- paste0("sample_", mus_hy_samples$Run)
names(fuk_hy_files) <- paste0("sample_", fuk_hy_samples$Run)
names(bos_hy_files) <- paste0("sample_", bos_hy_samples$Run)
names(het_hy_files) <- paste0("sample_", het_hy_samples$Run)
names(hom_hy_files) <- paste0("sample_", hom_hy_samples$Run)
names(mac_hy_files) <- paste0("sample_", mac_hy_samples$Run)
names(mic_hy_files) <- paste0("sample_", mic_hy_samples$Run)
names(ovi_hy_files) <- paste0("sample_", ovi_hy_samples$Run)
names(pap_hy_files) <- paste0("sample_", pap_hy_samples$Run)
names(per_hy_files) <- paste0("sample_", per_hy_samples$Run)
names(rat_hy_files) <- paste0("sample_", rat_hy_samples$Run)
all(file.exists(cav_hy_files))
all(file.exists(cap_hy_files))
all(file.exists(sor_hy_files))
all(file.exists(ict_hy_files))
all(file.exists(pan_hy_files))
all(file.exists(sus_hy_files))
all(file.exists(mus_hy_files))
all(file.exists(fuk_hy_files))
all(file.exists(bos_hy_files))
all(file.exists(het_hy_files))
all(file.exists(hom_hy_files))
all(file.exists(mac_hy_files))
all(file.exists(ovi_hy_files))
all(file.exists(mic_hy_files))
all(file.exists(pap_hy_files))
all(file.exists(per_hy_files))
all(file.exists(rat_hy_files))
#
cav.hypot.kallisto.txi.tpm.tsv <- tximport(cav_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = cavtx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
cap.hypot.kallisto.txi.tpm.tsv <- tximport(cap_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = captx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
sor.hypot.kallisto.txi.tpm.tsv <- tximport(sor_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = sortx2gene, ignoreAfterBar=TRUE)
ict.hypot.kallisto.txi.tpm.tsv <- tximport(ict_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = icttx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
pan.hypot.kallisto.txi.tpm.tsv <- tximport(pan_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = pantx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
sus.hypot.kallisto.txi.tpm.tsv <- tximport(sus_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = sustx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
mus.hypot.kallisto.txi.tpm.tsv <- tximport(mus_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = mustx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
fuk.hypot.kallisto.txi.tpm.tsv <- tximport(fuk_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = fuktx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
bos.hypot.kallisto.txi.tpm.tsv <- tximport(bos_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = bostx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
het.hypot.kallisto.txi.tpm.tsv <- tximport(het_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = hettx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
hom.hypot.kallisto.txi.tpm.tsv <- tximport(hom_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = homtx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
mac.hypot.kallisto.txi.tpm.tsv <- tximport(mac_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = mactx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
mic.hypot.kallisto.txi.tpm.tsv <- tximport(mic_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = mictx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
pap.hypot.kallisto.txi.tpm.tsv <- tximport(pap_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = paptx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
per.hypot.kallisto.txi.tpm.tsv <- tximport(per_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = pertx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
ovi.hypot.kallisto.txi.tpm.tsv <- tximport(ovi_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = ovitx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
rat.hypot.kallisto.txi.tpm.tsv <- tximport(rat_hy_files, type = "kallisto", countsFromAbundance = "lengthScaledTPM", tx2gene = rattx2gene, ignoreAfterBar=TRUE,ignoreTxVersion=TRUE)
#write them to file
write.table(cav.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/cav.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(cap.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/cap.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(sor.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/sor.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(ict.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/ict.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(pan.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/pan.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(mus.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/mus.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(sus.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/sus.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(fuk.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/fuk.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(bos.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/bos.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(het.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/het.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(hom.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/hom.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(mac.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/mac.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(mic.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/mic.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(pap.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/pap.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(per.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/per.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(ovi.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/ovi.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
write.table(rat.hypot.kallisto.txi.tpm.tsv$abundance, "../data/5_tpms/rat.hypot.kallisto.txi.tpm.tsv", na = "NA", col.names = TRUE, row.names = TRUE, sep = "\t", quote = FALSE)
#Can look at the count structures if desired
cav.hypot.kallisto.txi.tpm.tsv
ggplot(sor.hypot.kallisto.txi.tpm.tsv, aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
#
as_tibble(cav.hypot.kallisto.txi.tpm.tsv)
as_tibble(sor.cor.kallisto.txi.tpm.tsv$abundance) %>%
  pivot_longer(names_to = "sample", values_to = "counts", cols = 1:5) %>%
  ggplot(aes(x = log2(counts + 1), fill = sample)) +
  geom_histogram(bins = 30) +
  facet_wrap(~ sample)
