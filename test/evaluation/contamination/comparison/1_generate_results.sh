#!/bin/bash

uncalled pafstats -r ../true_mappings.paf --annotate ../uncalled/contamination_uncalled.paf > contamination_uncalled_ann.paf 2> contamination_uncalled.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../sigmap/contamination_sigmap.paf > contamination_sigmap_ann.paf 2> contamination_sigmap.throughput

uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/contamination_rawhash_viral.paf > contamination_rawhash_viral_ann.paf 2> contamination_rawhash_viral.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/contamination_rawhash_sensitive.paf > contamination_rawhash_sensitive_ann.paf 2> contamination_rawhash_sensitive.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/contamination_rawhash_fast.paf > contamination_rawhash_fast_ann.paf 2> contamination_rawhash_fast.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/contamination_rawhash_faster.paf > contamination_rawhash_faster_ann.paf 2> contamination_rawhash_faster.throughput

python ../../../scripts/compare_pafs.py contamination_uncalled_ann.paf contamination_sigmap_ann.paf contamination_rawhash_viral_ann.paf > contamination_rawhash_viral.comparison
python ../../../scripts/compare_pafs.py contamination_uncalled_ann.paf contamination_sigmap_ann.paf contamination_rawhash_sensitive_ann.paf > contamination_rawhash_sensitive.comparison
python ../../../scripts/compare_pafs.py contamination_uncalled_ann.paf contamination_sigmap_ann.paf contamination_rawhash_fast_ann.paf > contamination_rawhash_fast.comparison
python ../../../scripts/compare_pafs.py contamination_uncalled_ann.paf contamination_sigmap_ann.paf contamination_rawhash_faster_ann.paf > contamination_rawhash_faster.comparison
