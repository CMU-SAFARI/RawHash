#!/bin/bash

mkdir -p d2_ecoli_r94/fast5_files/
cd d2_ecoli_r94

#Download FAST5 from AWS | Unzip; Mv FAST5 files into the fast5_files directory. NCBI SRA Link: https://trace.ncbi.nlm.nih.gov/Traces/?run=ERR9127551
wget -qO- https://sra-pub-src-2.s3.amazonaws.com/ERR9127551/ecoli_r9.tar.gz.1 | tar xzv; mv r9/f5s/RefStrains210914_NK/f5s/barcode02/*.fast5 fast5_files; rm -rf r9

#Download FASTQ from SRA (Note: fastq-dump should exist in your path.) | #Processing the FASTA file so that read names contain the read ids as stored in FAST5 files
fastq-dump -Z --fasta 0 ERR9127551 | awk '{if(substr($1,1,1) == ">"){print ">"$2}else{print $0}}' > reads.fasta

#Downloading Escherichia coli CFT073, complete genome (https://www.ncbi.nlm.nih.gov/nuccore/AE014075.1/); Unzip; Change name;
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/007/445/GCA_000007445.1_ASM744v1/GCA_000007445.1_ASM744v1_genomic.fna.gz; gunzip GCA_000007445.1_ASM744v1_genomic.fna.gz; mv GCA_000007445.1_ASM744v1_genomic.fna ref.fa

cd ..
