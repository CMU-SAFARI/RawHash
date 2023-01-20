#!/bin/bash

THREAD=32

#d3_yeast_r94
OUTDIR="./uncalled/"
FAST5="../../../data/d3_yeast_r94/fast5_files/"
REF="../../../data/d3_yeast_r94/ref.fa"
PREFIX="d3_yeast_r94"
mkdir -p ${OUTDIR}

#Default
bash ../../../scripts/run_uncalled.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${THREAD} > "${OUTDIR}/${PREFIX}_uncalled.out" 2> "${OUTDIR}/${PREFIX}_uncalled.err"
