#!/bin/bash

THREAD=$1

#contamination
OUTDIR="./rawhash2/"
FAST5="../../data/contamination/fast5_files/"
REF="../../data/d1_sars-cov-2_r94/ref.fa"
PORE="../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PRESET="viral"
mkdir -p ${OUTDIR}

#Viral preset (default for viral genomes)
PREFIX="contamination"
PARAMS="--depletion"
bash ../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"

#Minimizers
PREFIX="contamination_w3"
PARAMS="-w 3 --depletion"
bash ../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"
