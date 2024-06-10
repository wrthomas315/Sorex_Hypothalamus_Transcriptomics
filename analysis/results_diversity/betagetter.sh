#!/bin/bash

awk '{print $4}' indivBetaMLparamsadulthypot_DIVALV.res | sed '/^$/d' | awk '{printf "%8.6f\n", $1}' > p
paste results p > resultsDIV_beta
