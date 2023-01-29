#!/bin/bash

THREAD=$1

#d1_sars-cov-2_r94
OUTDIR="./uncalled/"
FAST5="../../../data/d1_sars-cov-2_r94/fast5_files/"
REF="../../../data/d1_sars-cov-2_r94/ref.fa"
PREFIX="d1_sars-cov-2_r94"
mkdir -p ${OUTDIR}

#Default
bash ../../../scripts/run_uncalled.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${THREAD} > "${OUTDIR}/${PREFIX}_uncalled.out" 2> "${OUTDIR}/${PREFIX}_uncalled.err"
