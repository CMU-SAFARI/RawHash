#!/bin/bash

uncalled pafstats -r ../true_mappings.paf --annotate ../uncalled/contamination_uncalled.paf > contamination_uncalled_ann.paf 2> contamination_uncalled.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../sigmap/contamination_sigmap.paf > contamination_sigmap_ann.paf 2> contamination_sigmap.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/contamination_rawhash_viral.paf > contamination_rawhash_viral_ann.paf 2> contamination_rawhash_viral.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/contamination_rawhash2_viral.paf > contamination_rawhash2_viral_ann.paf 2> contamination_rawhash2_viral.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/contamination_w3_rawhash2_viral.paf > contamination_w3_rawhash2_viral_ann.paf 2> contamination_w3_rawhash2_viral.throughput

python ../../../scripts/compare_pafs.py contamination_uncalled_ann.paf contamination_sigmap_ann.paf contamination_rawhash_viral_ann.paf contamination_rawhash2_viral_ann.paf > contamination_rawhash2_viral.comparison
python ../../../scripts/compare_pafs.py contamination_uncalled_ann.paf contamination_sigmap_ann.paf contamination_rawhash_viral_ann.paf contamination_w3_rawhash2_viral_ann.paf > contamination_w3_rawhash2_viral.comparison
