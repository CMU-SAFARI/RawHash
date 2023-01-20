#!/bin/bash

uncalled pafstats -r ../true_mappings.paf --annotate ../uncalled/d5_human_na12878_r94_uncalled.paf > d5_human_na12878_r94_uncalled_ann.paf 2> d5_human_na12878_r94_uncalled.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../sigmap/d5_human_na12878_r94_sigmap.paf > d5_human_na12878_r94_sigmap_ann.paf 2> d5_human_na12878_r94_sigmap.throughput

uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/d5_human_na12878_r94_rawhash_sensitive.paf > d5_human_na12878_r94_rawhash_sensitive_ann.paf 2> d5_human_na12878_r94_rawhash_sensitive.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/d5_human_na12878_r94_rawhash_fast.paf > d5_human_na12878_r94_rawhash_fast_ann.paf 2> d5_human_na12878_r94_rawhash_fast.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/d5_human_na12878_r94_rawhash_faster.paf > d5_human_na12878_r94_rawhash_faster_ann.paf 2> d5_human_na12878_r94_rawhash_faster.throughput

python ../../../../scripts/compare_pafs.py d5_human_na12878_r94_uncalled_ann.paf d5_human_na12878_r94_sigmap_ann.paf d5_human_na12878_r94_rawhash_sensitive_ann.paf > d5_human_na12878_r94_rawhash_sensitive.comparison
python ../../../../scripts/compare_pafs.py d5_human_na12878_r94_uncalled_ann.paf d5_human_na12878_r94_sigmap_ann.paf d5_human_na12878_r94_rawhash_fast_ann.paf > d5_human_na12878_r94_rawhash_fast.comparison
python ../../../../scripts/compare_pafs.py d5_human_na12878_r94_uncalled_ann.paf d5_human_na12878_r94_sigmap_ann.paf d5_human_na12878_r94_rawhash_faster_ann.paf > d5_human_na12878_r94_rawhash_faster.comparison
