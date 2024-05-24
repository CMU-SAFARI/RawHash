#!/bin/bash

#Runs the indexing and mapping step of UNCALLED, generates the output file as well as the corresponding .time files

OUTDIR=$1 #Path to the output directory to store all the files to be generated
PREFIX=$2 #A string prefix that you want to attach to the file names to use as an identifier for your current run (e.g., mytestrun)
SIGNALS=$3 #Path to the directory that contains the fast5 files
REF=$4 #Path to the reference genome
THREAD=$5 #Number of threads to use
PARAMS=$6 #(optional -- you can keep it empty) custom parameters to set on top of the default parameters

/usr/bin/time -vpo "${OUTDIR}/${PREFIX}_uncalled_index.time" uncalled index ${PARAMS} -o "${OUTDIR}/${PREFIX}_uncalled.ind" ${REF}
/usr/bin/time -vpo "${OUTDIR}/${PREFIX}_uncalled_map.time" uncalled map -t ${THREAD} ${PARAMS} "${OUTDIR}/${PREFIX}_uncalled.ind" ${SIGNALS} > "${OUTDIR}/${PREFIX}_uncalled.paf"
