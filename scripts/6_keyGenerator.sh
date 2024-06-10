#!/bin/bash

#Need to create keys to pull out desired 1:1 orthologs
#will need to change path below or run line by line within folder
cd ../data//0_refs/primary_transcripts/OrthoFinder/Results_Date/
ls Single_Copy_Orthologue_Sequences | sed 's/.fa//g' > singlecopres2.txt
bash keycreate.sh
#location
/gpfs/scratch/withomas/hippoc_primary_transcripts/OrthoFinder/Results_Sep29/
#put whereever your ressults are
#need to change id and species for each individual
#cavPor:ENSCPOG0
#susScr:ENSSSCG0
#homSap:ENSG0
#fukMec:ENSFDAG0
#macMul:ENSMMUG0
#musMus:ENSMUSG0
#panTro:ENSPTRG0
#ratNor:ENSRNOG0
#sorAra:P_
#capHir:ENSCHIG0
#hetGla:ENSHGLG0
#micOch:ENSMOCG0
#oviAri:ENSOARG0
#ictTri:ENSSTOG0
#perMan:ENSPEMG0
#bosTau:ENSBTAG0
#script
#!/bin/bash
folder="Single_Copy_Orthologue_Sequences/"
#id for sorex
id="P_"
#species id variable
spec="sorAra"
while read j; do
        echo $j
        f=$(cat $folder$j".fa" | grep "$id")
        echo $j $f >> "../data/6_keys/gc_"$spec".txt"
done <singlecopyres2.txt

#after check to make sure that there are no duplicates in each file
#test to make sure it works and that there are no duplicates, do so for everyfile
awk 'a[$2]++{print $2}' gc_$species_key.txt
awk 'a[$1]++{print $1}' gc_$species_key.txt
#trim each file of version number
cat gc_$species.txt | cut -f1 -d . > almost
sed 's/>//g' almost > $species_key
mv $species_key ../data/6_keys/
#congrats on the keys, remember to run for each species

#Shrew keys can struggle depending on cluster because of additional ".1"s at end, so sometimes need to run below
while read j; do
        echo $j
        cat ../data/0_refs/GCF_000181275.1_SorAra2.0_genomic.gtf | grep $j > "../data/6_keys/"$j
        awk -v var=$j '{print var FS $10}' "../data/6_keys/"$j | sort | uniq >> ../data/6_keys//no_pap_almost_sor_key
done</Users/bill/project_Adult_EVE/data/step1_transcripts2genes/keys/hypot/no_pap/sorXPPs.txt

cat no_pap_almost_sor_key | tr -d '"' > soralmost2
cat soralmost2 | tr -d ';' > soralmosst3
awk '{print $1}' rat_key >> justOs.txt
paste justOs.txt soralmosst3 > soralmost4
awk '{print $1 FS $2 FS $3}' soralmost4 > final_sor_key
