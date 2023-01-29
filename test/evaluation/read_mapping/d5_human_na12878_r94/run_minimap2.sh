#!/bin/bash

THREAD=$1

bash ../../../scripts/run_minimap2.sh . ../../../data/d5_human_na12878_r94/reads.fasta ../../../data/d5_human_na12878_r94/ref.fa ${THREAD}
