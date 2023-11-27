OUTDIR="./rawhash2/"

for n in {1..5}
do 
    #sbatch -c 32 --wrap="bash run_rawhash2.sh 32 ${n}"
    for e in {4..8}
    do 
        #sbatch -c 32 --wrap="bash run_rawhash2.sh 32 ${n} ${e}" 
        #uncalled pafstats -r ../evaluation/read_mapping/d1_sars-cov-2_r94/true_mappings.paf -n 5000 --annotate ../evaluation/read_mapping/d1_sars-cov-2_r94/rawhash2/n_${n}_e_${e}/d1_sars-cov-2_r94_rawhash2_viral.paf > rawhash2_n_${n}_e_${e}_blqann.paf
        python test_paf.py rawhash2_n_${n}_e_${e}_blqann.paf
    done 
done 