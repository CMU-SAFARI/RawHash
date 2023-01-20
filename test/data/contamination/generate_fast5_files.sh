#!/bin/bash

mkdir -p fast5_files

cd fast5_files;

#The following commands will create symbolic links from already generated fast5 files to here
find ../../d5_human_na12878_r94/fast5_files/ -type f -name "*.fast5" | xargs -i{} ln -s {} .
find ../../d1_sars-cov-2_r94/fast5_files/ -type f -name "*.fast5" | xargs -i{} ln -s {} .
