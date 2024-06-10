#!/bin/bash

#need to download sratoolkit.2.11.0-mac64
#set path
export PATH=$PATH:/path/to/sratoolkit.2.11.0-mac64/bin/
#Make id lists from Supplemental table
mkdir ../data/1_ids/
#Paired reads
awk '{print $1 FS $5 FS $6}' ../data/Adult_Table.txt | grep "Pair" | awk '{print $3}' > ../data/1_ids/paired_ids.txt
#Single reads
awk '{print $1 FS $5 FS $6}' ../data/Adult_Table.txt | grep "Single" | awk '{print $3}' > ../data/1_ids/single_ids.txt
#by species/genus
awk '{print $1}' ../data/Adult_Table.txt | grep -v "Spe"| sort -u | tail -n 17 | grep -v "_polionotus" > ../data/1_ids/specs.txt
../data/1_ids/specs/
while read j; do
	spec="${j:0:3}"
	echo $spec
	awk '{print $1 FS $5 FS $6}' ../data/Adult_Table.txt | grep "$spec" | awk '{print $3}' > ../data/1_ids/specs/$spec.txt
done <../data/1_ids/specs.txt

#paired end reads SRA dump
while read j; do
        echo $j
	prefetch -v $j
	fastq-dump -I --split-files --outdir "../data/2_SRA/" $j".sra"
	#use --split-files if paired
done <../data/1_ids/paired_ids.txt

#single end reads SRA dump
while read j; do
        echo $j
        prefetch -v $j
        fastq-dump -I --outdir "../data/2_SRA/" $j".sra"
        #use --split-files if paired
done <../data/1_ids/single_ids.txt
