# Dehnels_Seasonal_RNAseq
Scripts and code to reproduce RNAseq analysis for looking at selection in gene expression for the shrew hypothalamus, and comparing these results to seasonal changes in relation to Dehnel's phenomenon.

###Goals and strategy

The objective of this project is to evaluate evolutionary changes in shrew hypothalamus expression by comparing to RNA sequencing of other mammals found on the Sequencing Read Archives (SRA). Genes found to be divergent can be inferred to be altered due to selective processes. These genes can then be compared to those that change throughout the size plasticity of Dehnel’s phenomenon in *Sorex araneus*. First, we will have to download reads from the SRA, access the quality of our RNA-seq data, filter low quality reads and trim adapters, map to the transcriptome of each species and quantify abundance. This will be followed by identifying orthologous genes and then normalization across species. Then we will 1) analyze differential expression between stages of Dehnel’s phenomenon using DESeq2, 2) characterize temporal patterns in expression using TCSeq. Overall, these analyses will show us evolutionary divergence and seasonal change in expression, as well as the intersection between these two sets of analyses.
Note: There are various spots you could start at with this analysis using this github. You can download all of the references and original reads and start at script 
0_reference_download.sh. However, as sequencing reads are not stored here, you could alternatively start at the post quantification steps, as abundance files are included 4/5 (still have to download references 0_reference_download.sh). OR! You could start right at the EVE and DESeq2 analyses (9_ and beyond)

### Data

RNA-seq analyses require alignment to a reference and quantification of reads. The genomes and original unfiltered reads can be downloaded as described below. First, download the references for each species, as well as the SRA sequencing data.

```
#get refs
bash 0_reference_download.sh
#get SRA
bash 1_SRA_download.sh
```

### Quality control, filtering, trimming
Here we will trim adapters from our reads and remove low quality reads using default settings and fastp. Will need to download fastp to your local environment (https://github.com/OpenGene/fastp).

```
bash 2_fastp.sh
```

### Mapping and quantification
Reads that have went through quality control are then mapped to the reference transcriptome and quantified using pseudoalignment. This method does not directly map reads to the genome, but can infer counts despite similarities between different coding regions (https://pachterlab.github.io/kallisto/about).

```
bash 3_kallisto.sh
```

Now need to convert kallisto abundance outputs to gene counts with transcripts to gene functions in R. (Note: Typically run this in R studio)
```
R 4_transcriptgene.R
```

### Orthology
RNAseq data is almost ready to be used. Now we need to find orthologous genes that are 1:1. Then we need to pull these out and put into a file that is useable for our EVE analyses. (Note: this will create the orthofinder result folder in ../0_refs)
Note: Probably nest to run this one line by line. Need to change path to your installation of orthofinder.
```
bash 5_orthofinder.sh
```

### Generate keys
Note: Need to run this one for each species, not automated. Something I can fix in the future.
```
6_keyGenerator.sh
```

Use these keys to change ENS count files to orthologs
Note: Need to run this one for each species, not automated. Something I can fix in the future.
```
7_ens2ortho.sh
```

And then combine into one file needed for EVE analyses
Note: Need to run this one for each species, not automated. Something I can fix in the future.
```
8_EVEexpressionGen.sh
```

And then can check to see how normalization through count distributions (Note: Here it is ran with papio, which is when we determined to remove it and go back to orthofinder; Typically run this in R studio)
```
R 9_tpmCheck.R
```

The other two inputs for EVE analyses are tree and individual files. Note, you  will need to create multiple of these for the dropout analysis to make sure that shrews were infer to be selectively upregulated in the shrew hypothalamus isn't doing this due to pervasive selection
Note: When species not available use closest relative in genus; need to download etetool kit to manipulate trees
```
python 10_treeprune.py
```

### Analyses
Note need to install EVE
### Branch shift
```
bash 11_EVE_branchshift.sh
```

### Diversity/Plasticity
```
bash 12_EVE_diversity.sh
```

### Shrew Dropout
```
bash 13_EVE_dropout.sh
```

### Convert outputs for each test so they are useable in R codes
```
bash 14_ConvertOuput_branch.sh -f BSThetaTestLRTsadult_hypothal_BS2ALV.res -s 16 -g genelist -o result
bash 15_ConvertOuput_div.sh -f betaTestLRTsadulthypot_DIVALV.res -s 16 -g genelist -o result
bash 16_ConvertOuput_dropout.sh -f  -s 15 -g genelist -o result
```

### Additionally, diversity and drop out need to have there outputs converted to get beta values, and to annotate those that are significant. 
```
bash 17_ConvertBeta_div.sh
bash 18_ConvertBeta_dropout.sh
```

### These should be run with the final R script where we determine significance and beta values, which then annotates them to make datasets comparable between differential expression and EVE
```
#branch sig
bash 19_nameconvertSig.sh
#branch all genes
bash 20_nameconvertTot.sh
#diversity Low Beta
bash 21_divBetaSigLow.sh 
#diversoty High Beta
bash 22_divBetaSigHigh.sh 
#dropout (only need diversifying, not plasticity)
bash 23_dropBetaSigDiv.sh
```

### Seasonal transcriptomics (produces both outputs for DESeq2 and TCSeq) 
(Note: run in tandem with below eve R analysis/name converters in Rstudio)
```
24_Dehnels_Hypothalamus.R
```

### Now compare results in a final R script (Note: run in tandem with above DESeq/name converters in Rstudio)
```
25_eve.R
```

### Statistical analyses and visualization for cell viability.
```
R 26_CV_t_test.R
```


### DAVID Geneset Enrichment
DAVID was done online at the below link. In a perfect world these should be scripted, however, due to conflicts in packages and Rversions they were not.
https://david.ncifcrf.gov/summary.jsp
