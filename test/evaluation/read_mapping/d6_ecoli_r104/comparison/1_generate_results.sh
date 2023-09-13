#!/bin/bash

# uncalled pafstats -r ../true_mappings.paf --annotate ../uncalled/d6_ecoli_r104_uncalled.paf > d6_ecoli_r104_uncalled_ann.paf 2> d6_ecoli_r104_uncalled.throughput
# uncalled pafstats -r ../true_mappings.paf --annotate ../sigmap/d6_ecoli_r104_sigmap.paf > d6_ecoli_r104_sigmap_ann.paf 2> d6_ecoli_r104_sigmap.throughput

#UNCALLED and Sigmap cannot run on R10. So we simply use their R9 results but
# these results are not shown in the output of the second scripts as this would not be an accurate comparison.
# This is simply done just to the scripts below successfully. Make sure their R9 results are generated before running the scripts below.
ln -s ../../d2_ecoli_r94/comparison/d2_ecoli_r94_rawhash_sensitive_ann.paf d6_ecoli_r104_rawhash_sensitive_ann.paf
ln -s ../../d2_ecoli_r94/comparison/d2_ecoli_r94_rawhash_sensitive.throughput d6_ecoli_r104_rawhash_sensitive.throughput 
ln -s ../../d2_ecoli_r94/comparison/d2_ecoli_r94_sigmap_ann.paf d6_ecoli_r104_sigmap_ann.paf
ln -s ../../d2_ecoli_r94/comparison/d2_ecoli_r94_sigmap.throughput d6_ecoli_r104_sigmap.throughput 
ln -s ../../d2_ecoli_r94/comparison/d2_ecoli_r94_uncalled.throughput d6_ecoli_r104_uncalled.throughput 
ln -s ../../d2_ecoli_r94/comparison/d2_ecoli_r94_uncalled_ann.paf d6_ecoli_r104_uncalled_ann.paf 

uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d6_ecoli_r104_rawhash2_r10sensitive.paf > d6_ecoli_r104_rawhash2_r10sensitive_ann.paf 2> d6_ecoli_r104_rawhash2_r10sensitive.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d6_ecoli_r104_w3_rawhash2_r10sensitive.paf > d6_ecoli_r104_w3_rawhash2_r10sensitive_ann.paf 2> d6_ecoli_r104_w3_rawhash2_r10sensitive.throughput

python ../../../../scripts/compare_pafs.py d6_ecoli_r104_uncalled_ann.paf d6_ecoli_r104_sigmap_ann.paf d6_ecoli_r104_rawhash_sensitive_ann.paf d6_ecoli_r104_rawhash2_r10sensitive_ann.paf > d6_ecoli_r104_rawhash2_r10sensitive.comparison
python ../../../../scripts/compare_pafs.py d6_ecoli_r104_uncalled_ann.paf d6_ecoli_r104_sigmap_ann.paf d6_ecoli_r104_rawhash_sensitive_ann.paf d6_ecoli_r104_w3_rawhash2_r10sensitive_ann.paf > d6_ecoli_r104_w3_rawhash2_r10sensitive.comparison
