#!/bin/bash

THREAD=$1

#d2_ecoli_r94
OUTDIR="./rawhash2/"
# FAST5="../../../data/d2_ecoli_r94/fast5_files/"
FAST5="../../../data/d2_ecoli_r94/small_fast5_dir/"
REF="../../../data/d2_ecoli_r94/ref.fa"
PORE="../../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PRESET="sensitive"
mkdir -p ${OUTDIR}

PREFIX="d2_ecoli_r94"
# bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"
bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} #> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"

# #The following is the run using default parameters:
# PREFIX="d2_ecoli_r94"
# PARAMS="--best-chains 5 --w-threshold 0.5 --map-model ../w0_bc0_mc5_occ01_sensitive_logistic.tflite"
# bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" # > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"

# #The following is the run using default parameters:
# PREFIX="d2_ecoli_r94"
# PARAMS="--best-chains 5 --w-threshold 0.5 --map-model ../features_all_all_m_logistic.tflite"
# bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"

# #Minimizers
# PREFIX="d2_ecoli_r94_w3"
# PARAMS="-w 3 --best-chains 5 --w-threshold 0.5 --map-model ../features_all_all_m_logistic.tflite"
# bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"
