#!/bin/bash

mkdir -p d5_human_na12878_r94/fast5_files/
cd d5_human_na12878_r94

#Download FAST5 from https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md
wget -qO- http://s3.amazonaws.com/nanopore-human-wgs/rel6/MultiFast5Tars/FAB42260-4177064552_Multi_Fast5.tar | tar xv

find ./UBC -type f -name '*.fast5' | xargs -i{} mv {} ./fast5_files/;

#Download FASTQ
wget -qO- http://s3.amazonaws.com/nanopore-human-wgs/rel6/FASTQTars/FAB42260-4177064552_Multi.tar | tar xv

zcat UBC/FAB42260-4177064552_Multi/fastq/*.fastq.gz | awk 'BEGIN{line = 0}{line++; if(line %4 == 1){print ">"substr($1,2,36)}else if(line % 4 == 2){print $0}}' > reads.fasta

rm -rf UBC;

#Downloading CHM13v2 (hs1) Human reference genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hs1/bigZips/hs1.fa.gz; gunzip hs1.fa.gz; mv hs1.fa ref.fa

cd ..
