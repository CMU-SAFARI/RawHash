#!/bin/bash

echo "Uncalled throughput (mean/median)"
grep "BP per sec:" d2_ecoli_r94_uncalled.throughput
echo "Sigmap throughput (mean/median)"
grep "BP per sec:" d2_ecoli_r94_sigmap.throughput
echo "RawHash throughput (mean/median)"
grep "BP per sec:" d2_ecoli_r94_rawhash_sensitive.throughput
# echo "RawHash (pod5) throughput (mean/median)"
# grep "BP per sec:" d2_ecoli_r94_pod5_rawhash_sensitive.throughput

echo;
grep "Uncalled Mean time per read" d2_ecoli_r94_rawhash_sensitive.comparison
grep "Sigmap Mean time per read" d2_ecoli_r94_rawhash_sensitive.comparison
grep "RawHash Mean time per read" d2_ecoli_r94_rawhash_sensitive.comparison
# grep "RawHash (pod5) Mean time per read" d2_ecoli_r94_pod5_rawhash_sensitive.comparison

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
grep "Uncalled precision:" d2_ecoli_r94_rawhash_sensitive.comparison
grep "Sigmap precision:" d2_ecoli_r94_rawhash_sensitive.comparison
grep "RawHash precision:" d2_ecoli_r94_rawhash_sensitive.comparison
# grep "RawHash precision:" d2_ecoli_r94_pod5_rawhash_sensitive.comparison

echo;
grep "Uncalled recall:" d2_ecoli_r94_rawhash_sensitive.comparison
grep "Sigmap recall:" d2_ecoli_r94_rawhash_sensitive.comparison
grep "RawHash recall:" d2_ecoli_r94_rawhash_sensitive.comparison
# grep "RawHash recall:" d2_ecoli_r94_pod5_rawhash_sensitive.comparison

echo;
grep "Uncalled F-1 score:" d2_ecoli_r94_rawhash_sensitive.comparison
grep "Sigmap F-1 score:" d2_ecoli_r94_rawhash_sensitive.comparison
grep "RawHash F-1 score:" d2_ecoli_r94_rawhash_sensitive.comparison
# grep "RawHash F-1 score:" d2_ecoli_r94_pod5_rawhash_sensitive.comparison
