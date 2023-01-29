#!/bin/bash

THREAD=$1

bash ../../../scripts/run_minimap2.sh . ../../../data/d3_yeast_r94/reads.fasta ../../../data/d3_yeast_r94/ref.fa ${THREAD}
