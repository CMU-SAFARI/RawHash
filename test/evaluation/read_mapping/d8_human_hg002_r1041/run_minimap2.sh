#!/bin/bash

THREAD=$1

bash ../../../scripts/run_minimap2.sh . ../../../data/d7_human_hg002_r1041/reads.fasta ../../../data/d7_human_hg002_r1041/ref.fa ${THREAD}
