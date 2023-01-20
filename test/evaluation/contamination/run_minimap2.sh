#!/bin/bash

#contamination (negative control)
minimap2 -x map-ont -t 32 -o true_mappings_neg.paf ../../data/contamination/ref.fa ../../data/d4_human_na12878_r94/reads.fasta
