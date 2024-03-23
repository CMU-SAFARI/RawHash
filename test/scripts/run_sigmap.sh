#!/bin/bash

#Runs the indexing and mapping step of sigmap, generates the output file as well as the corresponding .time files

OUTDIR=$1 #Path to the output directory to store all the files to be generated
PREFIX=$2 #A string prefix that you want to attach to the file names to use as an identifier for your current run (e.g., mytestrun)
SIGNALS=$3 #Path to the directory that contains the fast5 files
REF=$4 #Path to the reference genome
PORE=$5 #Path to the k-mer model file
THREAD=$6 #Number of threads to use
PARAMS=$7 #(optional -- you can keep it empty) custom parameters to set on top of the default parameters

/usr/bin/time -vpo "${OUTDIR}/${PREFIX}_sigmap_index.time" sigmap -i -p "${PORE}" -o "${OUTDIR}/${PREFIX}_sigmap.ind" -r ${REF} ${PARAMS}
/usr/bin/time -vpo "${OUTDIR}/${PREFIX}_sigmap_map.time" sigmap -m -t ${THREAD} -r ${REF} -p "${PORE}" -o "${OUTDIR}/${PREFIX}_sigmap.paf" -x "${OUTDIR}/${PREFIX}_sigmap.ind" -s ${SIGNALS} ${PARAMS}
