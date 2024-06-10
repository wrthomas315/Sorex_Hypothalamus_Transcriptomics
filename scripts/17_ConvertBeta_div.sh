#!/bin/bash

awk '{print $4}' ../analysis/results_diversity/indivBetaMLparamsadulthypot_DIVALV.res | sed '/^$/d' | awk '{printf "%8.6f\n", $1}' > ../analysis/results_diversity/p
paste ../analysis/results_diversity/results ../analysis/results_diversity/p > ../analysis/results_diversity/resultsDIV_beta
