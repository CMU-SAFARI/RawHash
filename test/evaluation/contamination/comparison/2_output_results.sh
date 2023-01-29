#!/bin/bash

echo "Uncalled throughput (mean/median)"
grep "BP per sec:" contamination_uncalled.throughput
echo "Sigmap throughput (mean/median)"
grep "BP per sec:" contamination_sigmap.throughput
echo "RawHash throughput (mean/median)"
grep "BP per sec:" contamination_rawhash_viral.throughput

echo;
grep "Uncalled Mean time per read" contamination_rawhash_viral.comparison
grep "Sigmap Mean time per read" contamination_rawhash_viral.comparison
grep "RawHash Mean time per read" contamination_rawhash_viral.comparison

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
grep "Uncalled precision:" contamination_rawhash_viral.comparison
grep "Sigmap precision:" contamination_rawhash_viral.comparison
grep "RawHash precision:" contamination_rawhash_viral.comparison

echo;
grep "Uncalled recall:" contamination_rawhash_viral.comparison
grep "Sigmap recall:" contamination_rawhash_viral.comparison
grep "RawHash recall:" contamination_rawhash_viral.comparison

echo;
grep "Uncalled F-1 score:" contamination_rawhash_viral.comparison
grep "Sigmap F-1 score:" contamination_rawhash_viral.comparison
grep "RawHash F-1 score:" contamination_rawhash_viral.comparison
