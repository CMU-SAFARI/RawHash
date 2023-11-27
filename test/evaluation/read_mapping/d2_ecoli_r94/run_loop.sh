#!/bin/bash

FAST5="../../../data/d2_ecoli_r94/fast5_files/"
REF="../../../data/d2_ecoli_r94/ref.fa"
PORE="../../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PRESET="sensitive"
PREFIX="d2_ecoli_r94"
#d2_ecoli_r94

for e in {3..7}
do 
    #sbatch -c 32 --wrap="bash run_rawhash2.sh 32 ${n}"
    for n in {2..5}
    do 
        
        sbatch -c 32 --wrap="bash run_rawhash2.sh 32 ${n} ${e}" 
    done 
done ""

#sbatch -c 32 -e="test_5_7err.txt" --wrap="bash run_rawhash2.sh 32 5 7" 