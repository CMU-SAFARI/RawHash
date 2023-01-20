#!/bin/bash

/usr/bin/time -vpo relative_abundance_rawhash_fast_sequenceuntil.time rawhash --sequence-until -t 32 ../rawhash/relative_abundance_rawhash_fast.ind ../../../data/random_community/fast5_files/ > relative_abundance_rawhash_fast_sequenceuntil.paf

/usr/bin/time -vpo relative_abundance_rawhash_fast.time rawhash -t 32 ../rawhash/relative_abundance_rawhash_fast.ind ../../../data/random_community/fast5_files/ > relative_abundance_rawhash_fast.paf

/usr/bin/time -vpo relative_abundance_rawhash_faster_sequenceuntil.time rawhash --sequence-until -t 32 ../rawhash/relative_abundance_rawhash_faster.ind ../../../data/random_community/fast5_files/ > relative_abundance_rawhash_faster_sequenceuntil.paf

/usr/bin/time -vpo relative_abundance_rawhash_faster.time rawhash -t 32 ../rawhash/relative_abundance_rawhash_faster.ind ../../../data/random_community/fast5_files/ > relative_abundance_rawhash_faster.paf
