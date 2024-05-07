#!/bin/bash

#Please make sure that you compile rawhash2 with the profiling option enabled. In /src/Makefile, you should enable the two lines below the line "# For profiling"

THREAD=$1

#d7_saureus_r104
OUTDIR="./rawhash2/"
FAST5="../../../data/d7_saureus_r104/fast5_files/FAT95629_pass_barcode01_16032940_0.fast5"
REF="../../../data/d7_saureus_r104/ref.fa"
PORE="../../../../extern/kmer_models/dna_r10.4.1_e8.2_400bps/9mer_levels_v1.txt"
PRESET="sensitive"
mkdir -p ${OUTDIR}

#The following is the run using default parameters:
PREFIX="d7_saureus_r104_profile_"${THREAD}
PARAMS="--r10"
bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"
