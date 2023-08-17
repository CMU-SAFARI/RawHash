#!/bin/bash

#Please make sure that you compile rawhash with the profiling option enabled. In /src/Makefile, you should enable the two lines below the line "# For profiling"

THREAD=$1

#d6_ecoli_r104
OUTDIR="./rawhash/"
FAST5="../../../data/d6_ecoli_r104/fast5_files/nomultiplex_r0b0_0.fast5"
REF="../../../data/d6_ecoli_r104/ref.fa"
PORE="../../../../extern/kmer_models/dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt"
PREFIX="d6_ecoli_r104_profile_"${THREAD}
mkdir -p ${OUTDIR}

#The following is the run using default parameters:
PRESET="sensitive"
bash ../../../scripts/run_rawhash.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} > "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash_${PRESET}.err"
