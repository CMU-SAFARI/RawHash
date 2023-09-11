#!/bin/bash

echo "RawHash2 throughput (mean/median)"
grep "BP per sec:" d6_ecoli_r104_rawhash2_r10sensitive.throughput
echo "(Minimizer w = 3)"
echo "RawHash2 throughput (mean/median)"
grep "BP per sec:" d6_ecoli_r104_w3_rawhash2_r10sensitive.throughput

echo;
grep "RawHash2 Mean time per read" d6_ecoli_r104_rawhash2_r10sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean time per read" d6_ecoli_r104_w3_rawhash2_r10sensitive.comparison


echo;
grep "RawHash2 Mean # of sequenced bases per read" d6_ecoli_r104_rawhash2_r10sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean # of sequenced bases per read" d6_ecoli_r104_w3_rawhash2_r10sensitive.comparison

echo;
grep "RawHash2 Mean # of sequenced chunks per read" d6_ecoli_r104_rawhash2_r10sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean # of sequenced chunks per read" d6_ecoli_r104_w3_rawhash2_r10sensitive.comparison

echo;
grep "RawHash2 Mean (only mapped) # of sequenced bases per read" d6_ecoli_r104_rawhash2_r10sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean (only mapped) # of sequenced bases per read" d6_ecoli_r104_w3_rawhash2_r10sensitive.comparison

echo;
grep "RawHash2 Mean (only mapped) # of sequenced chunks per read" d6_ecoli_r104_rawhash2_r10sensitive.comparison
grep "RawHash2 Mean (only mapped) # of sequenced chunks per read" d6_ecoli_r104_w3_rawhash2_r10sensitive.comparison

echo;
grep "RawHash2 Mean (only unmapped) # of sequenced bases per read" d6_ecoli_r104_rawhash2_r10sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean (only unmapped) # of sequenced bases per read" d6_ecoli_r104_w3_rawhash2_r10sensitive.comparison

echo;
grep "RawHash2 Mean (only unmapped) # of sequenced chunks per read" d6_ecoli_r104_rawhash2_r10sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 Mean (only unmapped) # of sequenced chunks per read" d6_ecoli_r104_w3_rawhash2_r10sensitive.comparison

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
grep "RawHash2 precision:" d6_ecoli_r104_rawhash2_r10sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 precision:" d6_ecoli_r104_w3_rawhash2_r10sensitive.comparison

echo;
grep "RawHash2 recall:" d6_ecoli_r104_rawhash2_r10sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 recall:" d6_ecoli_r104_w3_rawhash2_r10sensitive.comparison

echo;
grep "RawHash2 F-1 score:" d6_ecoli_r104_rawhash2_r10sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash2 F-1 score:" d6_ecoli_r104_w3_rawhash2_r10sensitive.comparison
