#!/bin/bash

mkdir -p d6_ecoli_r104/fast5_files/
cd d6_ecoli_r104

#Download FAST5 from AWS | Unzip; Mv FAST5 files into the fast5_files directory. NCBI SRA Link: https://trace.ncbi.nlm.nih.gov/Traces/?run=ERR9127552
wget -qO- https://sra-pub-src-2.s3.amazonaws.com/ERR9127552/Ecoli_r10.4.tar.gz.1 | tar xzv; mv ./mnt/data/analysis/nick/q20_refstrains/data/r10.4/fast5s/211123Ecoli_Q20_112/f5s/*.fast5 fast5_files; rm -rf ./mnt

# Optional: fast5 -> pod5 conversion. Requires the pod5 tool: https://github.com/nanoporetech/pod5-file-format 
# pod5 convert fast5 -r -t 32 --output-one-to-one ./fast5_files/ ./fast5_files/ ./pod5_files/

#Download FASTQ from SRA (Note: fastq-dump should exist in your path.) | #Processing the FASTA file so that read names contain the read ids as stored in FAST5 files
fastq-dump -Z --fasta 0 ERR9127552 | awk '{if(substr($1,1,1) == ">"){print ">"$2}else{print $0}}' > reads.fasta

#Downloading Escherichia coli CFT073, complete genome (https://www.ncbi.nlm.nih.gov/nuccore/AE014075.1/); Unzip; Change name;
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/007/445/GCA_000007445.1_ASM744v1/GCA_000007445.1_ASM744v1_genomic.fna.gz; gunzip GCA_000007445.1_ASM744v1_genomic.fna.gz; mv GCA_000007445.1_ASM744v1_genomic.fna ref.fa

cd ..
