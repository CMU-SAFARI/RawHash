#!/bin/bash

THREAD=$1

#contamination
OUTDIR="./sigmap/"
FAST5="../../data/contamination/fast5_files/"
REF="../../data/d1_sars-cov-2_r94/ref.fa"
PORE="../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PREFIX="contamination"
mkdir -p ${OUTDIR}

#Default parameters
bash ../../scripts/run_sigmap.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${THREAD} > "${OUTDIR}/${PREFIX}_sigmap.out" 2> "${OUTDIR}/${PREFIX}_sigmap.err"

# #contamination (negative control)
# OUTDIR="./sigmap/"
# FAST5="../../data/d4_human_na12878_r94/fast5_files/"
# REF="../../data/contamination/ref.fa"
# PORE="../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
# PREFIX="contamination_neg"
# mkdir -p ${OUTDIR}

# #Default parameters
# bash ../../scripts/run_sigmap.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${THREAD} > "${OUTDIR}/${PREFIX}_sigmap.out" 2> "${OUTDIR}/${PREFIX}_sigmap.err"
