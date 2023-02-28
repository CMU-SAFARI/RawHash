#!/bin/bash

mkdir -p d1_sars-cov-2_r94/fast5_files/
cd d1_sars-cov-2_r94

#Download FAST5 and FASTQ Caddecentre
wget -qO-  https://cadde.s3.climb.ac.uk/SP1-raw.tgz | tar -xzv; rm README

#Moving the fast5 files to fast5_files
find ./SP1-fast5-mapped -type f -name '*.fast5' | xargs -i{} mv {} ./fast5_files/; rm -rf SP1-fast5-mapped; 

#Converting FASTQ to FASTA files
awk 'BEGIN{line = 0}{line++; if(line %4 == 1){print ">"substr($1,2)}else if(line % 4 == 2){print $0}}' SP1-mapped.fastq > reads.fasta; rm SP1-mapped.fastq

#Downloading SARS-CoV-2 (covid) reference genome from RefSeq
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/858/895/GCF_009858895.2_ASM985889v3/GCF_009858895.2_ASM985889v3_genomic.fna.gz; gunzip GCF_009858895.2_ASM985889v3_genomic.fna.gz; mv GCF_009858895.2_ASM985889v3_genomic.fna ref.fa; 

cd ..
