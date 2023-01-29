#!/bin/bash

THREAD=$1

bash ../../../scripts/run_minimap2.sh . ../../../data/d2_ecoli_r94/reads.fasta ../../../data/d2_ecoli_r94/ref.fa ${THREAD}

