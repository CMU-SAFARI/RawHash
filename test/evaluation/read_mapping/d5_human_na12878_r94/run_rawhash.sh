#!/bin/bash

THREAD=$1

#d5_human_na12878_r94
OUTDIR="./rawhash/"
FAST5="../../../data/d5_human_na12878_r94/fast5_files/"
REF="../../../data/d5_human_na12878_r94/ref.fa"
PORE="../../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PREFIX="d5_human_na12878_r94"
mkdir -p ${OUTDIR}

#The following is the run using default parameters:
PRESET="fast"
bash ../../../scripts/run_rawhash.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} > "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.err"
