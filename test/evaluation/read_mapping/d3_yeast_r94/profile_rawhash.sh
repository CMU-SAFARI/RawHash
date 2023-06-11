#!/bin/bash

#Please make sure that you compile rawhash with the profiling option enabled. In /src/Makefile, you should enable the two lines below the line "# For profiling"

THREAD=$1

#d3_yeast_r94
OUTDIR="./rawhash/"
FAST5="../../../data/d3_yeast_r94/fast5_files/"
REF="../../../data/d3_yeast_r94/ref.fa"
PORE="../../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PREFIX="d3_yeast_r94_profile_"${THREAD}
mkdir -p ${OUTDIR}

#The following is the run using default parameters:
PRESET="sensitive"
bash ../../../scripts/run_rawhash.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} > "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.err"
