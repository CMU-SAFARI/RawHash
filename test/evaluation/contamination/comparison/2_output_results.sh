#!/bin/bash

echo "Uncalled throughput (mean/median)"
grep "BP per sec:" contamination_uncalled.throughput
echo "Sigmap throughput (mean/median)"
grep "BP per sec:" contamination_sigmap.throughput
echo "RawHash (Viral) throughput (mean/median)"
grep "BP per sec:" contamination_rawhash_viral.throughput
echo "RawHash (Sensitive) throughput (mean/median)"
grep "BP per sec:" contamination_rawhash_sensitive.throughput
echo "RawHash (Fast) throughput (mean/median)"
grep "BP per sec:" contamination_rawhash_fast.throughput
echo "RawHash (Faster) throughput (mean/median)"
grep "BP per sec:" contamination_rawhash_faster.throughput

echo;
grep "Uncalled Mean time per read" contamination_rawhash_faster.comparison
grep "Sigmap Mean time per read" contamination_rawhash_faster.comparison
echo "RawHash viral:"
grep "RawHash Mean time per read" contamination_rawhash_viral.comparison
echo "RawHash sensitive:"
grep "RawHash Mean time per read" contamination_rawhash_sensitive.comparison
echo "RawHash fast:"
grep "RawHash Mean time per read" contamination_rawhash_fast.comparison
echo "RawHash faster:"
grep "RawHash Mean time per read" contamination_rawhash_faster.comparison

echo;
echo '(Indexing) Timing and memory usage results:'
for i in `echo ../*/*index*.time`; do echo $i; 
	awk '{
		if(NR == 2){
			time = $NF
		}else if(NR == 3){
			time += $NF
			printf("CPU Time: %.2f\n", time)
		} else if(NR == 10){
			printf("Memory (GB): %.2f\n", $NF/1000000)
			print ""
		}
	}' $i; done

echo '(Mapping) Timing and memory usage results:'
for i in `echo ../*/*map.time`; do echo $i; 
	awk '{
		if(NR == 2){
			time = $NF
		}else if(NR == 3){
			time += $NF
			printf("CPU Time: %.2f\n", time)
		} else if(NR == 10){
			printf("Memory (GB): %.2f\n", $NF/1000000)
			print ""
		}
	}' $i; done
for i in `echo ../*/*_map_*.time`; do echo $i; 
	awk '{
		if(NR == 2){
			time = $NF
		}else if(NR == 3){
			time += $NF
			printf("CPU Time: %.2f\n", time)
		} else if(NR == 10){
			printf("Memory (GB): %.2f\n", $NF/1000000)
			print ""
		}
	}' $i; done


# echo;
# echo "Minimap2 mapped ratio:"
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*" && $16 == "s2:i:0"){count++;}tot++}END{print "Human mapped ratio: " count/269507}' ../true_mappings_neg.paf
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*" && $16 == "s2:i:0"){count++;}tot++}END{print "Covid mapped ratio: " count/1382016}' ../true_mappings.paf
# echo "Uncalled mapped ratio:"
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*"){count++;}tot++}END{print "Human mapped ratio: " count/269507}' ../uncalled/contamination_neg_uncalled.paf
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*"){count++;}tot++}END{print "Covid mapped ratio: " count/1382016}' ../../read_mapping/d5_sars-cov-2_r94/uncalled/d5_sars-cov-2_r94_uncalled.paf
# echo "Sigmap mapped ratio:"
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*"){count++;}tot++}END{print "Human mapped ratio: " count/269507}' ../sigmap/contamination_neg_sigmap.paf
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*"){count++;}tot++}END{print "Covid mapped ratio: " count/1382016}' ../../read_mapping/d5_sars-cov-2_r94/sigmap/d5_sars-cov-2_r94_sigmap.paf
# echo "RawHash (viral) mapped ratio:"
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*"){count++;}tot++}END{print "Human mapped ratio: " count/269507}' ../RawHash/contamination_neg_rawhash_viral.paf
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*"){count++;}tot++}END{print "Covid mapped ratio: " count/1382016}' ../../read_mapping/d5_sars-cov-2_r94/RawHash/d5_sars-cov-2_r94_rawhash_viral.paf
# echo "RawHash (sensitive) mapped ratio:"
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*"){count++;}tot++}END{print "Human mapped ratio: " count/269507}' ../RawHash/contamination_neg_rawhash_sensitive.paf
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*"){count++;}tot++}END{print "Covid mapped ratio: " count/1382016}' ../../read_mapping/d5_sars-cov-2_r94/RawHash/d5_sars-cov-2_r94_rawhash_sensitive.paf
# echo "RawHash (fast) mapped ratio:"
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*"){count++;}tot++}END{print "Human mapped ratio: " count/269507}' ../RawHash/contamination_neg_rawhash_fast.paf
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*"){count++;}tot++}END{print "Covid mapped ratio: " count/1382016}' ../../read_mapping/d5_sars-cov-2_r94/RawHash/d5_sars-cov-2_r94_rawhash_fast.paf
# echo "RawHash (faster) mapped ratio:"
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*"){count++;}tot++}END{print "Human mapped ratio: " count/269507}' ../RawHash/contamination_neg_rawhash_faster.paf
# awk 'BEGIN{count=0; tot=0;}{if($3 != "*" && $4 != "*"){count++;}tot++}END{print "Covid mapped ratio: " count/1382016}' ../../read_mapping/d5_sars-cov-2_r94/RawHash/d5_sars-cov-2_r94_rawhash_faster.paf

echo;
grep "Uncalled precision:" contamination_rawhash_faster.comparison
grep "Sigmap precision:" contamination_rawhash_faster.comparison
echo "RawHash viral:"
grep "RawHash precision:" contamination_rawhash_viral.comparison
echo "RawHash sensitive:"
grep "RawHash precision:" contamination_rawhash_sensitive.comparison
echo "RawHash fast:"
grep "RawHash precision:" contamination_rawhash_fast.comparison
echo "RawHash faster:"
grep "RawHash precision:" contamination_rawhash_faster.comparison

echo;
grep "Uncalled recall:" contamination_rawhash_faster.comparison
grep "Sigmap recall:" contamination_rawhash_faster.comparison
echo "RawHash viral:"
grep "RawHash recall:" contamination_rawhash_viral.comparison
echo "RawHash sensitive:"
grep "RawHash recall:" contamination_rawhash_sensitive.comparison
echo "RawHash fast:"
grep "RawHash recall:" contamination_rawhash_fast.comparison
echo "RawHash faster:"
grep "RawHash recall:" contamination_rawhash_faster.comparison

echo;
grep "Uncalled F-1 score:" contamination_rawhash_faster.comparison
grep "Sigmap F-1 score:" contamination_rawhash_faster.comparison
echo "RawHash viral:"
grep "RawHash F-1 score:" contamination_rawhash_viral.comparison
echo "RawHash sensitive:"
grep "RawHash F-1 score:" contamination_rawhash_sensitive.comparison
echo "RawHash fast:"
grep "RawHash F-1 score:" contamination_rawhash_fast.comparison
echo "RawHash faster:"
grep "RawHash F-1 score:" contamination_rawhash_faster.comparison
