#!/bin/bash -eux

THREAD=$1

#d2_ecoli_r94
OUTDIR="./rawhash2/"
# FAST5="../../../data/d2_ecoli_r94/fast5_files/"
FAST5="../../../data/d2_ecoli_r94/small_fast5_dir/"
REF="../../../data/d2_ecoli_r94/ref.fa"
PORE="../../../../extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model"
PRESET="sensitive"
mkdir -p ${OUTDIR}

PREFIX="d2_ecoli_r94"
# bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"
bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} #> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"
exit 0

FAST5=""

# port=39383
port="$(cut -d' ' -f1 /home/mmordig/rawhash_project/ru_python/example_run/server_run/ont_device_server_port.txt)"
[[ $port =~ ^[0-9]+$ ]] || (echo "Port is not a number: $port"; exit 1)
echo "Connecting to port $port"

PARAMS="--ru-server-port ${port}"
# export RU_AUTH_METHOD=NO_AUTH # note: also needs to be set on the server
# export RU_AUTH_METHOD="SSL";
# export RU_SSL_DEFAULT_BASE_PATH=~/rawhash_project/readuntil_fake/python_readuntil_client/generate_certs/generated;
export RU_AUTH_METHOD="SSL"; export RU_SSL_DEFAULT_BASE_PATH=~/rawhash_project/readuntil_fake/python_readuntil_client/generate_certs/generated; export MINKNOW_API_CLIENT_KEY=~/rawhash_project/readuntil_fake/python_readuntil_client/generate_certs/generated/client_key.pem; export MINKNOW_API_CLIENT_CERTIFICATE_CHAIN=~/rawhash_project/readuntil_fake/python_readuntil_client/generate_certs/generated/client_cert.pem; export RU_SERVER_CERTIFICATE_TARGET_NAME_OVERRIDE=localhost;
export SPDLOG_LEVEL=trace
bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} "${FAST5}" ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}"

# start the server
# clear
# cd ~/rawhash_project/ru_python/example_run
# export PATH=~/rawhash_project/readuntil_fake/build/:$PATH
# # export PATH=~/rawhash_project/readuntil_fake/build_docker/:$PATH
#
# bash ~/rawhash_project/readuntil_fake/python_readuntil_client/run_fake_server.sh


# Received the following args: rawhash2 -x sensitive -t 32 -o ./rawhash2//d2_ecoli_r94_rawhash2_sensitive.paf ./rawhash2//d2_ecoli_r94_rawhash2_sensitive.ind ../../../data/d2_ecoli_r94/small_fast5_dir/

# #The following is the run using default parameters:
# PREFIX="d2_ecoli_r94"
# PARAMS="--best-chains 5 --w-threshold 0.5 --map-model ../w0_bc0_mc5_occ01_sensitive_logistic.tflite"
# bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" # > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"

# #The following is the run using default parameters:
# PREFIX="d2_ecoli_r94"
# PARAMS="--best-chains 5 --w-threshold 0.5 --map-model ../features_all_all_m_logistic.tflite"
# bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"

# #Minimizers
# PREFIX="d2_ecoli_r94_w3"
# PARAMS="-w 3 --best-chains 5 --w-threshold 0.5 --map-model ../features_all_all_m_logistic.tflite"
# bash ../../../scripts/run_rawhash2.sh ${OUTDIR} ${PREFIX} ${FAST5} ${REF} ${PORE} ${PRESET} ${THREAD} "${PARAMS}" > "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.out" 2> "${OUTDIR}/${PREFIX}_rawhash2_${PRESET}.err"
