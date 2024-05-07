#!/bin/bash

uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d6_ecoli_r104_rawhash2_sensitive.paf > d6_ecoli_r104_rawhash2_sensitive_ann.paf 2> d6_ecoli_r104_rawhash2_sensitive.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d6_ecoli_r104_w3_rawhash2_sensitive.paf > d6_ecoli_r104_w3_rawhash2_sensitive_ann.paf 2> d6_ecoli_r104_w3_rawhash2_sensitive.throughput

python ../../../../scripts/analye_rh2_paf.py d6_ecoli_r104_rawhash2_sensitive_ann.paf > d6_ecoli_r104_rawhash2_sensitive.comparison
python ../../../../scripts/analye_rh2_paf.py d6_ecoli_r104_w3_rawhash2_sensitive_ann.paf > d6_ecoli_r104_w3_rawhash2_sensitive.comparison