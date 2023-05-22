#!/bin/bash

THREAD=$1

#contamination
OUTDIR="./uncalled/"
FAST5="../../data/contamination/fast5_files/"
REF="../../data/d1_sars-cov-2_r94/ref.fa"
PREFIX="contamination"
mkdir -p ${OUTDIR}

#Default
bash ../../scripts/run_uncalled.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${THREAD} > "${OUTDIR}/${PREFIX}_uncalled.out" 2> "${OUTDIR}/${PREFIX}_uncalled.err"

# #contamination (negative control)
# OUTDIR="./uncalled/"
# FAST5="../../data/d4_human_na12878_r94/fast5_files/"
# REF="../../data/contamination/ref.fa"
# PREFIX="contamination_neg"
# mkdir -p ${OUTDIR}

# #Default
# bash ../../scripts/run_uncalled.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${THREAD} > "${OUTDIR}/${PREFIX}_uncalled.out" 2> "${OUTDIR}/${PREFIX}_uncalled.err"
