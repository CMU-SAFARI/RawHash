#!/bin/bash

uncalled pafstats -r ../true_mappings.paf --annotate ../uncalled/d3_yeast_r94_uncalled.paf > d3_yeast_r94_uncalled_ann.paf 2> d3_yeast_r94_uncalled.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../sigmap/d3_yeast_r94_sigmap.paf > d3_yeast_r94_sigmap_ann.paf 2> d3_yeast_r94_sigmap.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/d3_yeast_r94_rawhash_sensitive.paf > d3_yeast_r94_rawhash_sensitive_ann.paf 2> d3_yeast_r94_rawhash_sensitive.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d3_yeast_r94_rawhash2_sensitive.paf > d3_yeast_r94_rawhash2_sensitive_ann.paf 2> d3_yeast_r94_rawhash2_sensitive.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d3_yeast_r94_w3_rawhash2_sensitive.paf > d3_yeast_r94_w3_rawhash2_sensitive_ann.paf 2> d3_yeast_r94_w3_rawhash2_sensitive.throughput

python ../../../../scripts/compare_pafs.py d3_yeast_r94_uncalled_ann.paf d3_yeast_r94_sigmap_ann.paf d3_yeast_r94_rawhash_sensitive_ann.paf d3_yeast_r94_rawhash2_sensitive_ann.paf > d3_yeast_r94_rawhash2_sensitive.comparison
python ../../../../scripts/compare_pafs.py d3_yeast_r94_uncalled_ann.paf d3_yeast_r94_sigmap_ann.paf d3_yeast_r94_rawhash_sensitive_ann.paf d3_yeast_r94_w3_rawhash2_sensitive_ann.paf > d3_yeast_r94_w3_rawhash2_sensitive.comparison

#POD5
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/d3_yeast_r94_pod5_rawhash_sensitive.paf > d3_yeast_r94_pod5_rawhash_sensitive_ann.paf 2> d3_yeast_r94_pod5_rawhash_sensitive.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d3_yeast_r94_pod5_rawhash2_sensitive.paf > d3_yeast_r94_pod5_rawhash2_sensitive_ann.paf 2> d3_yeast_r94_pod5_rawhash2_sensitive.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d3_yeast_r94_pod5_w3_rawhash2_sensitive.paf > d3_yeast_r94_pod5_w3_rawhash2_sensitive_ann.paf 2> d3_yeast_r94_pod5_w3_rawhash2_sensitive.throughput

python ../../../../scripts/compare_pafs.py d3_yeast_r94_uncalled_ann.paf d3_yeast_r94_sigmap_ann.paf d3_yeast_r94_pod5_rawhash_sensitive_ann.paf d3_yeast_r94_pod5_rawhash2_sensitive_ann.paf > d3_yeast_r94_pod5_rawhash2_sensitive.comparison
python ../../../../scripts/compare_pafs.py d3_yeast_r94_uncalled_ann.paf d3_yeast_r94_sigmap_ann.paf d3_yeast_r94_pod5_rawhash_sensitive_ann.paf d3_yeast_r94_pod5_w3_rawhash2_sensitive_ann.paf > d3_yeast_r94_pod5_w3_rawhash2_sensitive.comparison
