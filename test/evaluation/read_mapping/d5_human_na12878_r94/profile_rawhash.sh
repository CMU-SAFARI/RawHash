#!/bin/bash

THREAD=$1

#d5_human_na12878_r94
OUTDIR="./rawhash/"
FAST5="../../../data/d5_human_na12878_r94/fast5_files/MinION3_20161013_FNFAB42260_MN20093_sequencing_run_Chip98_Genomic_R9_4_480bps_0.fast5"
REF="../../../data/d5_human_na12878_r94/ref.fa"
PORE="../../../../extern/kmer_models/r9.4_180mv_450bps_6mer/template_median68pA.model"
PREFIX="d5_human_na12878_r94_profile_"${THREAD}
mkdir -p ${OUTDIR}

#The following is the run using default parameters:
PRESET="fast"
bash ../../../scripts/run_rawhash.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} > "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.err"
