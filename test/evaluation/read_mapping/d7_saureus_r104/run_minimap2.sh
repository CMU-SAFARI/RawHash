#!/bin/bash

THREAD=$1

bash ../../../scripts/run_minimap2.sh . ../../../data/d7_saureus_r104/reads.fasta ../../../data/d7_saureus_r104/ref.fa ${THREAD}
