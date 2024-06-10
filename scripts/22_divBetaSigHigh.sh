#!/bin/bash

awk '{print $1}' ../analysis/results_diversity/ADULTresultsDIV_sigHighBetaPlasticity > ../analysis/results_diversity/plast2

while read j; do
	p=$(cat ../data/6_keys/bos_key |  grep $j| awk '{print $1'})
	cat ../data/6_keys/final_sor_key | grep $p >> ../analysis/results_diversity/shrew_sigHighBetaPlasticity
	cat ../data/6_keys/hom_key | grep $p >> ../analysis/results_diversity/ENST_sigHighBetaPlasticity
done < ../analysis/results_diversity/plast2
