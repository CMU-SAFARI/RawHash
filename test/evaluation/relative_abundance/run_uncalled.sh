#!/bin/bash

THREAD=$1

#relative_abundance
OUTDIR="./uncalled/"
FAST5="../../data/relative_abundance/fast5_files/"
REF="../../data/relative_abundance/ref.fa"
PREFIX="relative_abundance"
mkdir -p ${OUTDIR}

#Default
bash ../../scripts/run_uncalled.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${THREAD} > "${OUTDIR}/${PREFIX}_uncalled.out" 2> "${OUTDIR}/${PREFIX}_uncalled.err"
