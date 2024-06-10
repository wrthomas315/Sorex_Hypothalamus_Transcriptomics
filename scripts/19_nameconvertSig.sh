#!/bin/bash

awk '{print $1}' ../analysis/results_branchshift/ADULTresultsSOR_sig > ../analysis/results_branchshift/adapt

while read j; do
	p=$(cat ../data/6_keys/bos_key |  grep $j| awk '{print $1'})
	#in terms of shrew genes for comparing to DESeq2 results
	cat ../data/6_keys/final_sor_key | grep $p >> ../analysis/results_branchshift/shrew_sigShrewBS
	#in ters of human ENSTs
	cat ../data/6_keys/hom_key | grep $p >> ../analysis/results_branchshift/ENST_sigShrewBS
done < ../analysis/results_branchshift/adapt
