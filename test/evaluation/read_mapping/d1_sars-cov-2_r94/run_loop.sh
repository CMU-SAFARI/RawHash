#!/bin/bash

THREAD=$1

#d1_sars-cov-2_r94
OUTDIR="./rawhash2/"
FAST5="../../../data/d1_sars-cov-2_r94/fast5_files/"
REF="../../../data/d1_sars-cov-2_r94/ref.fa"
PORE="../../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PRESET="viral"

for n in {1..5}
do 
    #sbatch -c 32 --wrap="bash run_rawhash2.sh 32 ${n}"
    for e in {1..6}
    do 
        sbatch -c 32 --wrap="bash run_rawhash2.sh 32 ${e} ${n} 5" 
    done 
done 

