#!/bin/bash

awk '{print $1}' ../analysis/results_dropout/dropoutDIV_sig > ../analysis/results_dropout/adapt

while read j; do
	p=$(cat ../6_keys/bos_key |  grep $j| awk '{print $1'})
	cat ../6_keys/final_sor_key | grep $p >> ../analysis/results_dropout/shrew_sigShrewDO
	cat ../6_keys/hom_key | grep $p >> ../analysis/results_dropout/ENST_sigShrewDO
done < ../analysis/results_dropout/adapt
