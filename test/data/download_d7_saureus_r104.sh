#!/bin/bash

mkdir -p d7_saureus_r104/fast5_files/
cd d7_saureus_r104/fast5_files

#Download FAST5 from AWS | Unzip; NCBI SRA Link: https://trace.ncbi.nlm.nih.gov/Traces/?run=SRR21386013
wget -qO- https://sra-pub-src-1.s3.amazonaws.com/SRR21386013/S_aureus_JKD6159_ONT_R10.4_fast5.tar.gz.1 | tar xzv;

cd ..;

#Basecalll the signals using dorado
dorado basecaller dna_r10.4.1_e8.2_400bps_sup@v4.1.0 fast5_files > reads.bam
samtools fasta reads.bam > reads.fasta
#TODO: Provide the direct link to download basecalled reads

#DownloadingStaphylococcus aureus subsp. aureus JKD6159, complete genome; Unzip; Change name;
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/144/955/GCF_000144955.2_ASM14495v2/GCF_000144955.2_ASM14495v2_genomic.fna.gz; gunzip GCF_000144955.2_ASM14495v2_genomic.fna.gz; mv GCF_000144955.2_ASM14495v2_genomic.fna ref.fa

cd ..
