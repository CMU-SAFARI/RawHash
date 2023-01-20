#!/bin/bash

mkdir -p fast5_files

cd fast5_files;

#The following commands will create symbolic links from already generated fast5 files to here
find ../../d1_sars-cov-2_r94/fast5_files/ -type f -name "*.fast5" | xargs -i{} ln -s {} .
find ../../d2_ecoli_r94/fast5_files/ -type f -name "*.fast5" | xargs -i{} ln -s {} .
find ../../d3_yeast_r94/fast5_files/ -type f -name "*.fast5" | xargs -i{} ln -s {} .
find ../../d4_green_algae_r94/fast5_files/ -type f -name "*.fast5" | xargs -i{} ln -s {} .
find ../../d5_human_na12878_r94/fast5_files/ -type f -name "*.fast5" | xargs -i{} ln -s {} .
