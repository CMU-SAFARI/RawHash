#!/bin/bash

minimap2 -x map-ont -t 32 -o true_mappings.paf ../../data/relative_abundance/ref.fa ../../data/d1_sars-cov-2_r94/reads.fasta ../../data/d2_ecoli_r94/reads.fasta ../../data/d3_yeast_r94/reads.fasta ../../data/d4_green_algae_r94/reads.fasta ../../data/d5_human_na12878_r94/reads.fasta
