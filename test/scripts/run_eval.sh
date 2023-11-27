#!/bin/bash

# for e in {3..7}
# do
#     for n in {2..5}
#     do 
#     #sbatch -c 32 --wrap="uncalled pafstats -r ../evaluation/read_mapping/d2_ecoli_r94/true_mappings.paf -n 5000 --annotate ../evaluation/read_mapping/d2_ecoli_r94/rawhash2/n{$n}_e{$e}_bl_reg/d2_ecoli_r94_rawhash2_sensitive.paf > rawhash2_ec_n{$n}_e{$e}_bl_reg.paf"
#     #python test_paf.py rawhash2_ec_n{$n}_e{$e}_bl_reg.paf > python_eval_n{$n}_e{$e}_bl_reg.txt 
#     python from_txt_exel.py  python_eval_n{$n}_e{$e}_bl_reg.txt $n $e
#     done
# done

for e in {3..7}
do 
    #sbatch -c 32 --wrap="uncalled pafstats -r ../evaluation/read_mapping/d2_ecoli_r94/true_mappings.paf -n 5000 --annotate ../evaluation/read_mapping/d2_ecoli_r94/rawhash2/e{$e}_reg_min_w5/d2_ecoli_r94_rawhash2_sensitive.paf > rawhash2_ec_e{$e}_reg_min_w5.paf"
    #python test_paf.py rawhash2_ec_e{$e}_reg_min_w5.paf > python_eval_e{$e}_reg_minw5.txt
    python from_txt_exel.py  python_eval_e{$e}_reg_minw5.txt $e
done

#sbatch -c 32 --wrap="uncalled pafstats -r ../evaluation/read_mapping/d2_ecoli_r94/true_mappings.paf -n 5000 --annotate ../evaluation/read_mapping/d2_ecoli_r94/rawhash2/n{2}_e{3}_M1/d2_ecoli_r94_rawhash2_sensitive.paf > rawhash2_ec_n{2}_e{3}_M1.paf"
#python test_paf.py rawhash2_ec_n{2}_e{3}_M1.paf > python_eval_n{2}_e{3}_M1.txt
#python from_txt_exel.py  python_eval_n{2}_e{3}_M1.txt 2 3