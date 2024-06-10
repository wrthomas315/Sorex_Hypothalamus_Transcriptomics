#!/bin/bash

mkdir ../data/0_refs

#make list of transcriptomes
cat ../data/Adult_Table.txt | awk '{print $8}' | sort -u | grep -v "TranscriptomeRef" | tail -n 16 > ../data/0_refs/transcrList.txt
#get them, difficult to put in loop so manual
#Mus
wget https://ftp.ensembl.org/pub/release-104/gtf/mus_musculus/Mus_musculus.GRCm39.104.gtf.gz ../data/0_refs/
#Homo
wget https://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz ../data/0_refs/
#Macaca
wget https://ftp.ensembl.org/pub/release-104/gtf/macaca_mulatta/Macaca_mulatta.Mmul_10.104.gtf.gz ../data/0_refs/
#Cavia
wget https://ftp.ensembl.org/pub/release-104/gtf/cavia_porcellus/Cavia_porcellus.Cavpor3.0.104.gtf.gz ../data/0_refs/
#Papio
wget https://ftp.ensembl.org/pub/release-104/gtf/papio_anubis/Papio_anubis.Panu_3.0.104.gtf.gz ../data/0_refs/
#Pan
wget https://ftp.ensembl.org/pub/release-104/gtf/pan_troglodytes/Pan_troglodytes.Pan_tro_3.0.104.gtf.gz
#Heterocephalus
wget https://ftp.ensembl.org/pub/release-104/gtf/heterocephalus_glaber_male/Heterocephalus_glaber_male.HetGla_1.0.104.gtf.gz ../data/0_refs/
#Ictidomys
wget https://ftp.ensembl.org/pub/release-112/gtf/ictidomys_tridecemlineatus/Ictidomys_tridecemlineatus.SpeTri2.0.112.gtf.gz ../data/0_refs/
#Fukomys
wget https://ftp.ensembl.org/pub/release-104/gtf/fukomys_damarensis/Fukomys_damarensis.DMR_v1.0.104.gtf.gz ../data/0_refs/
#Microtus
wget https://ftp.ensembl.org/pub/release-104/gtf/microtus_ochrogaster/Microtus_ochrogaster.MicOch1.0.104.gtf.gz ../data/0_refs/
#Bos
wget https://ftp.ensembl.org/pub/release-104/gtf/bos_taurus/Bos_taurus.ARS-UCD1.2.104.gtf.gz  ../data/0_refs/
#Peromyscus
wget https://ftp.ensembl.org/pub/release-104/gtf/peromyscus_maniculatus_bairdii/Peromyscus_maniculatus_bairdii.HU_Pman_2.1.104.gtf.gz ../data/0_refs/
#Sus
wget https://ftp.ensembl.org/pub/release-104/gtf/sus_scrofa/Sus_scrofa.Sscrofa11.1.104.gtf.gz ../data/0_refs/
#Ovis
wget https://ftp.ensembl.org/pub/release-104/gtf/ovis_aries_rambouillet/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.104.gtf.gz../data/0_refs/
#Capra
wget https://ftp.ensembl.org/pub/release-104/gtf/capra_hircus/Capra_hircus.ARS1.104.gtf.gz ../data/0_refs/
#Rattus
wget http://ftp.ensembl.org/pub/release-104/gtf/rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.104.gtf.gz ../data/0_refs/

#make list of peptide references
cat ../data/Adult_Table.txt | awk '{print $7}' | sort -u | grep -v "TranscriptomeRef" | tail -n 16 > ../data/0_refs/peptideList.txt
#get them
#Mus
wget https://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/pep/Mus_musculus.GRCm39.pep.all.fa.gz ../data/0_refs/
#Homo
wget https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz ../data/0_refs/
#Macaca
wget https://ftp.ensembl.org/pub/release-104/fasta/macaca_mulatta/pep/Macaca_mulatta.Mmul_10.pep.all.fa.gz ../data/0_refs/
#Cavia
wget https://ftp.ensembl.org/pub/release-104/fasta/cavia_porcellus/pep/Cavia_porcellus.Cavpor3.0.pep.all.fa.gz ../data/0_refs/
#Papio
wget https://ftp.ensembl.org/pub/release-104/fasta/papio_anubis/pep/Papio_anubis.Panu_3.0.pep.all.fa.gz ../data/0_refs/
#Pan 
wget https://ftp.ensembl.org/pub/release-104/fasta/pan_troglodytes/pep/Pan_troglodytes.Pan_tro_3.0.pep.all.fa.gz ../data/0_refs/
#Heterocephalus
wget https://ftp.ensembl.org/pub/release-104/fasta/heterocephalus_glaber_male/pep/Heterocephalus_glaber_male.HetGla_1.0.pep.all.fa.gz ../data/0_refs/
#Ictidomys
wget https://ftp.ensembl.org/pub/release-104/fasta/ictidomys_tridecemlineatus/pep/Ictidomys_tridecemlineatus.SpeTri2.0.pep.all.fa.gz ../data/0_refs/
#Fukomys
wget https://ftp.ensembl.org/pub/release-104/fasta/fukomys_damarensis/pep/Fukomys_damarensis.DMR_v1.0.pep.all.fa.gz ../data/0_refs/
#Microtus
wget https://ftp.ensembl.org/pub/release-104/fasta/microtus_ochrogaster/pep/Microtus_ochrogaster.MicOch1.0.pep.all.fa.gz ../data/0_refs/
#Bos
wget https://ftp.ensembl.org/pub/release-104/fasta/bos_taurus/pep/Bos_taurus.ARS-UCD1.2.pep.all.fa.gz ../data/0_refs/
#Peromyscus
wget https://ftp.ensembl.org/pub/release-104/fasta/peromyscus_maniculatus_bairdii/pep/Peromyscus_maniculatus_bairdii.HU_Pman_2.1.pep.all.fa.gz ../data/0_refs/
#Sus
wget https://ftp.ensembl.org/pub/release-104/fasta/sus_scrofa/pep/Sus_scrofa.Sscrofa11.1.pep.all.fa.gz ../data/0_refs/
#Ovis
wget https://ftp.ensembl.org/pub/release-104/fasta/ovis_aries_rambouillet/pep/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.pep.all.fa.gz ../data/0_refs/
#Capra
wget https://ftp.ensembl.org/pub/release-104/fasta/capra_hircus/pep/Capra_hircus.ARS1.pep.all.fa.gz ../data/0_refs/
#Rattus
wget http://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/pep/Rattus_norvegicus.Rnor_6.0.pep.all.fa.gz ../data/0_refs/


#and finally get cDNA files
cat ../data/Adult_Table.txt | awk '{print $9}' | sort -u | grep -v "CDNA_file" | tail -n 16 > ../data/0_refs/cDNAList.txt
#get them
#Mus
wget https://ftp.ensembl.org/pub/release-104/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz ../data/0_refs/
#Homo
wget https://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz ../data/0_refs/
#Macaca
wget https://ftp.ensembl.org/pub/release-104/fasta/macaca_mulatta/cdna/Macaca_mulatta.Mmul_10.cdna.all.fa.gz ../data/0_refs/
#Cavia
wget https://ftp.ensembl.org/pub/release-104/fasta/cavia_porcellus/cdna/Cavia_porcellus.Cavpor3.0.cdna.all.fa.gz ../data/0_refs/
#Papio
wget https://ftp.ensembl.org/pub/release-104/fasta/papio_anubis/cdna/Papio_anubis.Panu_3.0.cdna.all.fa.gz ../data/0_refs/
#Pan 
wget https://ftp.ensembl.org/pub/release-104/fasta/pan_troglodytes/cdna/Pan_troglodytes.Pan_tro_3.0.cdna.all.fa.gz ../data/0_refs/
#Heterocephalus
wget https://ftp.ensembl.org/pub/release-104/fasta/heterocephalus_glaber_male/cdna/Heterocephalus_glaber_male.HetGla_1.0.cdna.all.fa.gz ../data/0_refs/
#Ictidomys
wget https://ftp.ensembl.org/pub/release-104/fasta/ictidomys_tridecemlineatus/cdna/Ictidomys_tridecemlineatus.SpeTri2.0.cdna.all.fa.gz ../data/0_refs/
#Fukomys
wget https://ftp.ensembl.org/pub/release-104/fasta/fukomys_damarensis/cdna/Fukomys_damarensis.DMR_v1.0.cdna.all.fa.gz ../data/0_refs/
#Microtus
wget https://ftp.ensembl.org/pub/release-104/fasta/microtus_ochrogaster/cdna/Microtus_ochrogaster.MicOch1.0.cdna.all.fa.gz ../data/0_refs/
#Bos
wget https://ftp.ensembl.org/pub/release-104/fasta/bos_taurus/cdna/Bos_taurus.ARS-UCD1.2.cdna.all.fa.gz ../data/0_refs/
#Peromyscus
wget https://ftp.ensembl.org/pub/release-104/fasta/peromyscus_maniculatus_bairdii/cdna/Peromyscus_maniculatus_bairdii.HU_Pman_2.1.cdna.all.fa.g z../data/0_refs/
#Sus
wget https://ftp.ensembl.org/pub/release-104/fasta/sus_scrofa/cdna/Sus_scrofa.Sscrofa11.1.cdna.all.fa.gz ../data/0_refs/
#Ovis
wget https://ftp.ensembl.org/pub/release-104/fasta/ovis_aries_rambouillet/cdna/Ovis_aries_rambouillet.Oar_rambouillet_v1.0.cdna.all.fa.gz ../data/0_refs/
#Capra
wget https://ftp.ensembl.org/pub/release-104/fasta/capra_hircus/cdna/Capra_hircus.ARS1.cdna.all.fa.gz ../data/0_refs/
#Rattus
wget http://ftp.ensembl.org/pub/release-104/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.Rnor_6.0.cdna.all.fa.gz ../data/0_refs/
