#!/bin/bash

THREAD=$1

bash ../../../scripts/run_minimap2.sh . ../../../data/d4_green_algae_r94/reads.fasta ../../../data/d4_green_algae_r94/ref.fa ${THREAD}
