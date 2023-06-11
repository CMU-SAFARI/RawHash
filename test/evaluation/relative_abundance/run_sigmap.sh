#!/bin/bash

THREAD=$1

#relative_abundance
OUTDIR="./sigmap/"
FAST5="../../data/relative_abundance/fast5_files/"
REF="../../data/relative_abundance/ref.fa"
PORE="../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PREFIX="relative_abundance"
mkdir -p ${OUTDIR}

#Default parameters
bash ../../scripts/run_sigmap.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${THREAD} > "${OUTDIR}/${PREFIX}_sigmap.out" 2> "${OUTDIR}/${PREFIX}_sigmap.err"
