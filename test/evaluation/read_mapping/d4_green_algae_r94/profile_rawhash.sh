#!/bin/bash

#Please make sure that you compile rawhash with the profiling option enabled. In /src/Makefile, you should enable the two lines below the line "# For profiling"

THREAD=$1

#d4_green_algae_r94
OUTDIR="./rawhash/"
FAST5="../../../data/d4_green_algae_r94/fast5_files/PCT0062_20180831_0004A30B00232394_1_E5_H5_sequencing_run_PBNP18L0092_0831_A1_29025_read_111_ch_177_strand.fast5"
REF="../../../data/d4_green_algae_r94/ref.fa"
PORE="../../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PREFIX="d4_green_algae_r94_profile_"${THREAD}
mkdir -p ${OUTDIR}

#The following is the run using default parameters:
PRESET="fast"
bash ../../../scripts/run_rawhash.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} > "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.err"
