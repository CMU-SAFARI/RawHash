#!/bin/bash

uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d7_saureus_r104_rawhash2_sensitive.paf > d7_saureus_r104_rawhash2_sensitive_ann.paf 2> d7_saureus_r104_rawhash2_sensitive.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash2/d7_saureus_r104_w3_rawhash2_sensitive.paf > d7_saureus_r104_w3_rawhash2_sensitive_ann.paf 2> d7_saureus_r104_w3_rawhash2_sensitive.throughput

python ../../../../scripts/test_paf.py d7_saureus_r104_rawhash2_sensitive_ann.paf > d7_saureus_r104_rawhash2_sensitive.comparison
python ../../../../scripts/test_paf.py d7_saureus_r104_w3_rawhash2_sensitive_ann.paf > d7_saureus_r104_w3_rawhash2_sensitive.comparison
