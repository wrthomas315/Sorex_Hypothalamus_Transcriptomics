#!/bin/bash

awk '{print $1}' dropoutDIV_sig > adapt

while read j; do
	p=$(cat /Users/bill/project_Adult_EVE/data/step1_transcripts2genes/keys/hypot/no_pap/bos_key |  grep $j| awk '{print $1'})
	cat /Users/bill/project_Adult_EVE/data/step1_transcripts2genes/keys/hypot/no_pap/final_sor_key | grep $p >> shrew_sigShrewDO
	cat /Users/bill/project_Adult_EVE/data/step1_transcripts2genes/keys/hypot/no_pap/hom_key | grep $p >> ENST_sigShrewDO
done < adapt
#paste shrew_sigLowBetaSelection plast > Shrew_sigLowBetaSelection
#paste ENST_sigLowBetaSelection plast > eNST_sigLowBetaSelection
