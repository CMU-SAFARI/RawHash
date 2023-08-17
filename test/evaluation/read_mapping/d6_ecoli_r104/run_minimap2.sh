#!/bin/bash

THREAD=$1

bash ../../../scripts/run_minimap2.sh . ../../../data/d6_ecoli_r104/reads.fasta ../../../data/d6_ecoli_r104/ref.fa ${THREAD}
