#!/bin/bash

OUTDIR=$1
READS=$2
REF=$3
THREAD=$4

minimap2 -x map-ont -t ${THREAD} -o "${OUTDIR}/true_mappings.paf" ${REF} ${READS}
