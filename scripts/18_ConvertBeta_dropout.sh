#!/bin/bash

awk '{print $4}' ../analysis/results_dropout/indivBetaMLsadulthypot_DIVALVdropout.res | sed '/^$/d' | awk '{printf "%8.6f\n", $1}' > ../analysis/results_dropout/p
paste ../analysis/results_dropout/results ../analysis/results_dropout/p > ../analysis/results_dropout/resultsDIV_beta
