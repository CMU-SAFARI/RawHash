#!/bin/bash

THREAD=$1

#relative_abundance
OUTDIR="./rawhash2/"
FAST5="../../data/relative_abundance/fast5_files/"
# FAST5="../../data/relative_abundance/test_fast5_files/"
REF="../../data/relative_abundance/ref.fa"
PORE="../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PRESET="fast"
mkdir -p ${OUTDIR}

#The following is the run using default parameters:
# PREFIX="relative_abundance"
# PARAMS="--best-chains 5 --w-threshold 0.5 --map-model ../read_mapping/feature_set.tflite"
# bash ../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"

#Minimizers
PREFIX="relative_abundance_w3"
PARAMS="-w 3 --best-chains 5 --w-threshold 0.5 --map-model ../read_mapping/feature_set.tflite"
bash ../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"
