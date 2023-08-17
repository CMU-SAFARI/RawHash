#!/bin/bash

THREAD=$1

#d6_ecoli_r104
OUTDIR="./rawhash/"
FAST5="../../../data/d6_ecoli_r104/fast5_files/"
REF="../../../data/d6_ecoli_r104/ref.fa"
PORE="../../../../extern/kmer_models/dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt"
PRESET="sensitive"
mkdir -p ${OUTDIR}

#The following is the run using default parameters:
PREFIX="d6_ecoli_r104"
bash ../../../scripts/run_rawhash.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} > "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.err"

#Minimizers
PREFIX="d6_ecoli_r104_w3"
PARAMS="-w 3"
bash ../../../scripts/run_rawhash.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} ${PARAMS} > "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.err"

# PREFIX="d6_ecoli_r104_w5"
# PARAMS="-w 5"
# bash ../../../scripts/run_rawhash.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} ${PARAMS} > "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.err"
