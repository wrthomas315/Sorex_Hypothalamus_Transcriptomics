#!/bin/bash

#need to download orthofinder from github

#change path below and run this in folder with fa files
cd ../data/0_refs/
for f in *pep.all.fa ; do python /orthofinder_tutorial/OrthoFinder/tools/primary_transcript.py $f ; done

#then run orthofinder on the primary transcript folder output
source activate /path/to/.conda/envs/orthofinderenvironment

orthofinder -f primary_transcripts/
