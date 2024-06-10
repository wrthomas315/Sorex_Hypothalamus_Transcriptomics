#!/bin/bash

#move this to results folder for branch shift then run it
while getopts f:s:g:o: flag
do
    case "${flag}" in
        f) filename=${OPTARG};;
        s) species_num=${OPTARG};;
        g) genelist=${OPTARG};;
        o) output=${OPTARG};;
    esac
done
var=1
header=$(($species_num+$var))

#how to convert the output
# take all of the spaces of LRT file and turn them into linebreaks, then save as a new file .LTR
tr ' ' '\n' < "../analysis/results_diversity/"$filename > "../analysis/results_diversity/"$filename".LRTs"
# use awk to convert it out of scientific notation
awk '{printf "%8.6f\n", $1}' "../analysis/results_diversity/"$filename".LRTs" > "../analysis/results_diversity/"$filename".LRTs.nosci"
#paste it next to a file with just the ENST gene names
paste "../data/10_genelist/"$genelist "../analysis/results_diversity/"$filename".LRTs.nosci" > "../analysis/results_diversity/"$output
