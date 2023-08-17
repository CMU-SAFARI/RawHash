#!/bin/bash

uncalled pafstats -r ../true_mappings.paf --annotate ../uncalled/d4_green_algae_r94_uncalled.paf > d4_green_algae_r94_uncalled_ann.paf 2> d4_green_algae_r94_uncalled.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../sigmap/d4_green_algae_r94_sigmap.paf > d4_green_algae_r94_sigmap_ann.paf 2> d4_green_algae_r94_sigmap.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/d4_green_algae_r94_rawhash_fast.paf > d4_green_algae_r94_rawhash_fast_ann.paf 2> d4_green_algae_r94_rawhash_fast.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/d4_green_algae_r94_w3_rawhash_fast.paf > d4_green_algae_r94_w3_rawhash_fast_ann.paf 2> d4_green_algae_r94_w3_rawhash_fast.throughput
# uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/d4_green_algae_r94_w5_rawhash_fast.paf > d4_green_algae_r94_w5_rawhash_fast_ann.paf 2> d4_green_algae_r94_w5_rawhash_fast.throughput


python ../../../../scripts/compare_pafs.py d4_green_algae_r94_uncalled_ann.paf d4_green_algae_r94_sigmap_ann.paf d4_green_algae_r94_rawhash_fast_ann.paf > d4_green_algae_r94_rawhash_fast.comparison
python ../../../../scripts/compare_pafs.py d4_green_algae_r94_uncalled_ann.paf d4_green_algae_r94_sigmap_ann.paf d4_green_algae_r94_w3_rawhash_fast_ann.paf > d4_green_algae_r94_w3_rawhash_fast.comparison
# python ../../../../scripts/compare_pafs.py d4_green_algae_r94_uncalled_ann.paf d4_green_algae_r94_sigmap_ann.paf d4_green_algae_r94_w5_rawhash_fast_ann.paf > d4_green_algae_r94_w5_rawhash_fast.comparison
