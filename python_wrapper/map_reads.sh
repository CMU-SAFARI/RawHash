#!/usr/bin/env bash
set -eux

export PATH=~/rawhash_project/rawhash2-env/bin:$PATH
export PATH=~/rawhash_project/extra_tools/ont-guppy/bin/":$PATH"
# export PATH="$PATH:/home/mmordig/rawhash_project/tools/cmake-3.28.0-rc5-linux-x86_64/bin"
export LD_LIBRARY_PATH="/home/mmordig/rawhash_project/rawhash2_new/build/src"
#:"${LD_LIBRARY_PATH+}"
export PATH=~/rawhash_project/rawhash2_new/build/bin:"$PATH"

echo $LD_LIBRARY_PATH

THREAD=32

# OUTDIR="./rawhash2/"
OUTDIR="/home/mmordig/rawhash_project/rawhash2_new/example_out"
#SIGNALS="../../../data/d2_ecoli_r94/fast5_files/"
#SIGNALS="../../../data/d2_ecoli_r94/small_fast5_dir/"
SIGNALS="/home/mmordig/rawhash_project/rawhash2_new/test/data/d2_ecoli_r94/small_slow5_files/"
# # basecalled reads
# READS=
REF="../../../data/d2_ecoli_r94/ref.fa"
# PORE="../../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PORE="/home/mmordig/rawhash_project/rawhash2_new/extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PRESETX="sensitive"

# map to positive and revcomp strand
PARAMS=
PREFIX="forwrev"

# # only map to the positive strand
# PARAMS="--no-revcomp-query --rev-query"
# PREFIX="forwonly"

mkdir -p ${OUTDIR}


# if complaining about shared lib, need to recompile target
# index
rawhash2 -x ${PRESETX} -t ${THREAD} -p "${PORE}" -d "${OUTDIR}/${PREFIX}_rawhash2_${PRESETX}.ind" ${PARAMS} ${REF}
# map
rawhash2 -x ${PRESETX} -t ${THREAD} -o "${OUTDIR}/${PREFIX}_rawhash2_${PRESETX}.paf" ${PARAMS} "${OUTDIR}/${PREFIX}_rawhash2_${PRESETX}.ind" ${SIGNALS}


# awk '$5 == "-" {print $1}' /home/mmordig/rawhash_project/rawhash2_new/test/evaluation/read_mapping/d2_ecoli_r94/true_mappings.paf > revstrand_readids.txt
# slow5tools skim --rid /home/mmordig/rawhash_project/rawhash2_new/test/data/d2_ecoli_r94/small_slow5_files/barcode02_r0barcode02b0_0.blow5 > small_read_ids.txt
# # intersect both files
# awk 'NR==FNR{a[$1];next}($1 in a)' revstrand_readids.txt small_read_ids.txt > revstrand_readids_intersect.txt

## diff -y -W 250 forwonly_rawhash2_sensitive.paf forwrev_rawhash2_sensitive.paf | less -S
## sort revstrand_readids_intersect.txt | uniq -d # check no duplicates
# slow5tools get /home/mmordig/rawhash_project/rawhash2_new/test/data/d2_ecoli_r94/small_slow5_files/barcode02_r0barcode02b0_0.blow5 -t $(nproc) --list revstrand_readids_intersect.txt -o revcomp_reads.blow5
# awk 'NR==FNR{a[$1];next}($1 in a)' revstrand_readids_intersect.txt \
#     /home/mmordig/rawhash_project/rawhash2_new/test/evaluation/read_mapping/d2_ecoli_r94/true_mappings.paf \
#     > revcomp_true_mappings.paf

# bash ../../../scripts/run_minimap2.sh . ../../../data/d2_ecoli_r94/reads.fasta ../../../data/d2_ecoli_r94/ref.fa ${THREAD}
# minimap2 -x map-ont -t ${THREAD} -o "${OUTDIR}/true_mappings.paf" ${REF} ${READS}