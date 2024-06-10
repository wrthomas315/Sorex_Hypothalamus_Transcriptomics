#!/bin/bash

#mass create kallisto indexes
#while read j; do
#	gunzip "../data/0_refs/"$j".gz"
#	HEAD=$(echo $j | sed 's/fa//g')
#	REF=../data/0_refs/$j
#	IDX="../data/0_refs/"$HEAD"idx"
#	kallisto index -i $IDX $REF
#	gzip ../data/0_refs/$j
#done<../data/0_refs/cDNAList.txt

#make output directory for kallisto
#mkdir ../data/4_kallisto

#if paired
while read j; do
	echo $j
	REF=$(cat ../data/Adult_Table.txt | grep "$j" | awk '{print $9}')
	HEAD=$(echo $REF | sed 's/fa//g')
	IDX="../data/0_refs/"$HEAD"idx"
	kallisto quant -i $IDX -o ../data/4_kallisto/$j "../data/3_trimmed/"$j"_1.trimmed.fastq" "../data/3_trimmed/"$j"_2.trimmed.fastq"
done<../data/1_ids/paired_ids.txt

#if single
#note if this fails check read size, might need to change -l flag
while read j; do
	echo $j
        REF=$(cat ../data/Adult_Table.txt | grep "$j" | awk '{print $9}')
        HEAD=$(echo $REF | sed 's/fa//g')
        IDX="../data/0_refs/"$HEAD"idx"
	kallisto quant -l 100 -s .00000000001 -i $IDX -o ../data/4_kallisto/$j --single "../data/3_trimmed/"$j".trimmed.fastq"
#done <../data/1_ids/single_ids.txt
