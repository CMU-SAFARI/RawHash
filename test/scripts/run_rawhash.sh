#!/bin/bash

#Runs the indexing and mapping step of rawhash, generates the output file as well as the corresponding .time files

OUTDIR=$1 #Path to the output directory to store all the files to be generated
PREFIX=$2 #A string prefix that you want to attach to the file names to use as an identifier for your current run (e.g., mytestrun)
SIGNALS=$3 #Path to the directory that contains the fast5 files
REF=$4 #Path to the reference genome
PORE=$5 #Path to the k-mer model file
PRESETX=$6 #Default preset of rawhash for the run (e.g., viral)
THREAD=$7 #Number of threads to use
PARAMS=$8 #(optional -- you can keep it empty) custom parameters to set on top of the default parameters

/usr/bin/time -vpo "${OUTDIR}/${PREFIX}_rawhash_index_${PRESETX}.time" rawhash -x ${PRESETX} -t ${THREAD} ${PARAMS} -p "${PORE}" -d "${OUTDIR}/${PREFIX}_rawhash_${PRESETX}.ind" ${REF}
/usr/bin/time -vpo "${OUTDIR}/${PREFIX}_rawhash_map_${PRESETX}.time" rawhash -x ${PRESETX} -t ${THREAD} ${PARAMS} -o "${OUTDIR}/${PREFIX}_rawhash_${PRESETX}.paf" "${OUTDIR}/${PREFIX}_rawhash_${PRESETX}.ind" ${SIGNALS}
