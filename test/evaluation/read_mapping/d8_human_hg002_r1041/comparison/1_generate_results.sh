#!/bin/bash

# uncalled pafstats -r ../true_mappings.paf --annotate ../uncalled/d7_human_hg002_r1041_uncalled.paf > d7_human_hg002_r1041_uncalled_ann.paf 2> d7_human_hg002_r1041_uncalled.throughput
# uncalled pafstats -r ../true_mappings.paf --annotate ../sigmap/d7_human_hg002_r1041_sigmap.paf > d7_human_hg002_r1041_sigmap_ann.paf 2> d7_human_hg002_r1041_sigmap.throughput

#UNCALLED and Sigmap cannot run on R10. So we simply use their R9 results but
# these results are not shown in the output of the second scripts as this would not be an accurate comparison.
# This is simply done just to the scripts below successfully. Make sure their R9 results are generated before running the scripts below.
ln -s ../../d5_human_na12878_r94/comparison/d5_human_na12878_r94_rawhash_sensitive_ann.paf d7_human_hg002_r1041_rawhash_sensitive_ann.paf
ln -s ../../d5_human_na12878_r94/comparison/d5_human_na12878_r94_rawhash_sensitive.throughput d7_human_hg002_r1041_rawhash_sensitive.throughput 
ln -s ../../d5_human_na12878_r94/comparison/d5_human_na12878_r94_sigmap_ann.paf d7_human_hg002_r1041_sigmap_ann.paf
ln -s ../../d5_human_na12878_r94/comparison/d5_human_na12878_r94_sigmap.throughput d7_human_hg002_r1041_sigmap.throughput 
ln -s ../../d5_human_na12878_r94/comparison/d5_human_na12878_r94_uncalled.throughput d7_human_hg002_r1041_uncalled.throughput 
ln -s ../../d5_human_na12878_r94/comparison/d5_human_na12878_r94_uncalled_ann.paf d7_human_hg002_r1041_uncalled_ann.paf 

uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d7_human_hg002_r1041_rawhash2_r10sensitive.paf > d7_human_hg002_r1041_rawhash2_r10sensitive_ann.paf 2> d7_human_hg002_r1041_rawhash2_r10sensitive.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d7_human_hg002_r1041_w3_rawhash2_r10sensitive.paf > d7_human_hg002_r1041_w3_rawhash2_r10sensitive_ann.paf 2> d7_human_hg002_r1041_w3_rawhash2_r10sensitive.throughput

python ../../../../scripts/compare_pafs.py d7_human_hg002_r1041_uncalled_ann.paf d7_human_hg002_r1041_sigmap_ann.paf d7_human_hg002_r1041_rawhash_sensitive_ann.paf d7_human_hg002_r1041_rawhash2_r10sensitive_ann.paf > d7_human_hg002_r1041_rawhash2_r10sensitive.comparison
python ../../../../scripts/compare_pafs.py d7_human_hg002_r1041_uncalled_ann.paf d7_human_hg002_r1041_sigmap_ann.paf d7_human_hg002_r1041_rawhash_sensitive_ann.paf d7_human_hg002_r1041_w3_rawhash2_r10sensitive_ann.paf > d7_human_hg002_r1041_w3_rawhash2_r10sensitive.comparison
