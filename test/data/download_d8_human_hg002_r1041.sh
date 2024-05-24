#!/bin/bash

mkdir -p d7_human_hg002_r1041/slow5_files/
cd d7_human_hg002_r1041

wget -qO- https://sra-pub-src-2.s3.amazonaws.com/SRR23215365/hg2_subsample_slow5.tar.1 | tar xv; mv hg2_subsample.blow5* slow5_files;

fastq-dump -Z --fasta 0 SRR23215363 | awk '{if(substr($1,1,1) == ">"){print ">"$2}else{print $0}}' > reads.fasta

# fastq-dump -Z --fasta 0 ERR11768772 | awk '{if(substr($1,1,1) == ">"){print ">"$2}else{print $0}}' > reads2.fasta

#Downloading CHM13v2 (hs1) Human reference genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz; gunzip hs1.fa.gz; mv hs1.fa ref.fa

cd ..
