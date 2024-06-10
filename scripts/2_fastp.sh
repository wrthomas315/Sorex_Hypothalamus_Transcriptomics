#!/bin/bash

#Paired
while read j; do
        echo $j
        fastp --in1 "../data/2_SRA/"$j"_1.fastq" --in2 "../data/2_SRA/"$j"_2.fastq" --out1 "../data/3_trimmed/"$j"_1.trimmed.fastq" --out2 "../data/3_trimmed/"$j"_2.trimmed.fastq" --cut_tail --html $j".html" --json $j".json" 2> $j".log"
done <../data/1_ids/paired_ids.txt

#Single
while read j; do
        echo $j
        fastp --in1 "../data/2_SRA/"$j"_1.fastq" --out1 "../data/3_trimmed/"$j"_1.trimmed.fastq" --cut_tail --html $j".html" --json $j".json" 2> $j".log"
done <../data/1_ids/single_ids.txt
