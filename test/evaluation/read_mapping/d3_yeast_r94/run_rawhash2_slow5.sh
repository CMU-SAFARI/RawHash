#!/bin/bash

THREAD=$1

#d3_yeast_r94
OUTDIR="./rawhash2/"
FAST5="../../../data/d3_yeast_r94/slow5_files/"
REF="../../../data/d3_yeast_r94/ref.fa"
PORE="../../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PRESET="sensitive"
mkdir -p ${OUTDIR}

#The following is the run using default parameters:
PREFIX="d3_yeast_r94_slow5"
bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"

#Minimizers
PREFIX="d3_yeast_r94_slow5_w3"
PARAMS="-w 3"
bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"
