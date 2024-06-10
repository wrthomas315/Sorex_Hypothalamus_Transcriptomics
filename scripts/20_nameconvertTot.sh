#!/bin/bash

awk '{print $1}' ../analysis/results_branchshift/ADULTresultsSOR_pval_matrx > ../analysis/results_branchshift/adaptTOT

while read j; do
	p=$(cat ../data/6_keys/bos_key |  grep $j| awk '{print $1'})
	cat ../data/6_keys/final_sor_key | grep $p >> ../analysis/results_branchshift/shrew_TOTShrewBS
	cat ../data/6_keys/hom_key | grep $p >> ../analysis/results_branchshift/ENST_TOTShrewBS
done < ../analysis/results_branchshift/adaptTOT
