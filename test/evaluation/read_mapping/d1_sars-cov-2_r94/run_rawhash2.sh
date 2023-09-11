#!/bin/bash

THREAD=$1

#d1_sars-cov-2_r94
OUTDIR="./rawhash2/"
FAST5="../../../data/d1_sars-cov-2_r94/fast5_files/"
REF="../../../data/d1_sars-cov-2_r94/ref.fa"
PORE="../../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PRESET="viral"
mkdir -p ${OUTDIR}

#Viral preset (default for viral genomes)
PREFIX="d1_sars-cov-2_r94"
bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"

#Minimizers
PREFIX="d1_sars-cov-2_r94_w3"
PARAMS="-w 3"
bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} ${PARAMS} > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"
