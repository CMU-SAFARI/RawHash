#!/bin/bash

THREAD=$1

bash ../../../scripts/run_minimap2.sh . ../../../data/d1_sars-cov-2_r94/reads.fasta ../../../data/d1_sars-cov-2_r94/ref.fa ${THREAD}
