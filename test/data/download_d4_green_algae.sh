#!/bin/bash

mkdir -p d4_green_algae_r94/fast5_files/
cd d4_green_algae_r94

#Download FAST5 from AWS
wget -qO- https://sra-pub-src-2.s3.amazonaws.com/ERR3237140/Chlamydomonas_0.tar.gz.1 | tar xzv;

find ./Chlamydomonas_0/reads/downloads/pass/ -type f -name '*.fast5' | xargs -i{} mv {} fast5_files/

awk 'BEGIN{line = 0}{line++; if(line %4 == 1){print ">"substr($1,2)}else if(line % 4 == 2){print $0}}' ./Chlamydomonas_0/reads/downloads/pass/*.fastq > reads.fasta

rm -rf Chlamydomonas_0

#Downloading C.reinhardtii v.5.5 reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/595/GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5/GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.fna.gz; gunzip GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.fna.gz; mv GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.fna ref.fa

cd ..
