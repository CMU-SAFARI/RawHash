#!/bin/bash

echo "Uncalled throughput (mean/median)"
grep "BP per sec:" d4_green_algae_r94_uncalled.throughput
echo "Sigmap throughput (mean/median)"
grep "BP per sec:" d4_green_algae_r94_sigmap.throughput
echo "RawHash throughput (mean/median)"
grep "BP per sec:" d4_green_algae_r94_rawhash_fast.throughput
echo "RawHash2 throughput (mean/median)"
grep "BP per sec:" d4_green_algae_r94_rawhash2_sensitive.throughput
echo "RawHash2 (Minimizer w = 3) throughput (mean/median)"
grep "BP per sec:" d4_green_algae_r94_w3_rawhash2_sensitive.throughput

echo;
grep "Uncalled Mean time per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "Sigmap Mean time per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash Mean time per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash2 Mean time per read" d4_green_algae_r94_rawhash2_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean time per read" d4_green_algae_r94_w3_rawhash2_sensitive.comparison


echo;
grep "Uncalled Mean # of sequenced bases per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash Mean # of sequenced bases per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash2 Mean # of sequenced bases per read" d4_green_algae_r94_rawhash2_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean # of sequenced bases per read" d4_green_algae_r94_w3_rawhash2_sensitive.comparison

echo;
grep "Sigmap Mean # of sequenced chunks per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash Mean # of sequenced chunks per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash2 Mean # of sequenced chunks per read" d4_green_algae_r94_rawhash2_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean # of sequenced chunks per read" d4_green_algae_r94_w3_rawhash2_sensitive.comparison

echo;
grep "Uncalled Mean (only mapped) # of sequenced bases per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash Mean (only mapped) # of sequenced bases per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash2 Mean (only mapped) # of sequenced bases per read" d4_green_algae_r94_rawhash2_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean (only mapped) # of sequenced bases per read" d4_green_algae_r94_w3_rawhash2_sensitive.comparison

echo;
grep "Sigmap Mean (only mapped) # of sequenced chunks per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash Mean (only mapped) # of sequenced chunks per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash2 Mean (only mapped) # of sequenced chunks per read" d4_green_algae_r94_rawhash2_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean (only mapped) # of sequenced chunks per read" d4_green_algae_r94_w3_rawhash2_sensitive.comparison

echo;
grep "Uncalled Mean (only unmapped) # of sequenced bases per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash Mean (only unmapped) # of sequenced bases per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash2 Mean (only unmapped) # of sequenced bases per read" d4_green_algae_r94_rawhash2_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean (only unmapped) # of sequenced bases per read" d4_green_algae_r94_w3_rawhash2_sensitive.comparison

echo;
grep "Sigmap Mean (only unmapped) # of sequenced chunks per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash Mean (only unmapped) # of sequenced chunks per read" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash2 Mean (only unmapped) # of sequenced chunks per read" d4_green_algae_r94_rawhash2_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean (only unmapped) # of sequenced chunks per read" d4_green_algae_r94_w3_rawhash2_sensitive.comparison

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

echo;
grep "Uncalled precision:" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "Sigmap precision:" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash precision:" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash2 precision:" d4_green_algae_r94_rawhash2_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 precision:" d4_green_algae_r94_w3_rawhash2_sensitive.comparison

echo;
grep "Uncalled recall:" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "Sigmap recall:" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash recall:" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash2 recall:" d4_green_algae_r94_rawhash2_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 recall:" d4_green_algae_r94_w3_rawhash2_sensitive.comparison

echo;
grep "Uncalled F-1 score:" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "Sigmap F-1 score:" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash F-1 score:" d4_green_algae_r94_rawhash2_sensitive.comparison
grep "RawHash2 F-1 score:" d4_green_algae_r94_rawhash2_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 F-1 score:" d4_green_algae_r94_w3_rawhash2_sensitive.comparison
