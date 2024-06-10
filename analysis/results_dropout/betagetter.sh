#!/bin/bash

awk '{print $4}' indivBetaMLparamsadulthypot_DIVALVdropout.res | sed '/^$/d' | awk '{printf "%8.6f\n", $1}' > p
paste dropoutput p > results_dropoutput
