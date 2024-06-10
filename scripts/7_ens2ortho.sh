#!/bin/bash

#need to run the below for each species, changing the amount of columns you keep dependent on how many individuals per species you have
#only bos keeps the var2, the rest remove var2 from awk output as they will not be pasted in the next step
while read j; do
        f=$(awk -v var=$j '$1 == var {print $2}' ../data/6_keys/bos_key)
        awk -v var=$j -v var2=$f '$1 == var2 {print var2 FS $2 FS $3 FS $4 FS $5 FS $6 FS $7}' ../data/5_tpms/bos.hypot.kallisto.txt ../data/7_orthotpms/ortho.bos.hypot.tpm.tsv
done <../data/0_refs/primary_transcripts/Orthofinder/Res_Date/singlecopyres2.txt

