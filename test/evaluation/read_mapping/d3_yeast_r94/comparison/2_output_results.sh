#!/bin/bash

echo "Uncalled throughput (mean/median)"
grep "BP per sec:" d3_yeast_r94_uncalled.throughput
echo "Sigmap throughput (mean/median)"
grep "BP per sec:" d3_yeast_r94_sigmap.throughput
echo "RawHash throughput (mean/median)"
grep "BP per sec:" d3_yeast_r94_rawhash_sensitive.throughput
echo "(Minimizer w = 3)"
echo "RawHash throughput (mean/median)"
grep "BP per sec:" d3_yeast_r94_w3_rawhash_sensitive.throughput
# echo "(Minimizer w = 5)"
# echo "RawHash throughput (mean/median)"
# grep "BP per sec:" d3_yeast_r94_w5_rawhash_sensitive.throughput
echo "(POD5)"
echo "RawHash throughput (mean/median)"
grep "BP per sec:" d3_yeast_r94_pod5_rawhash_sensitive.throughput
echo "(POD5, Minimizer w = 3)"
echo "RawHash throughput (mean/median)"
grep "BP per sec:" d3_yeast_r94_pod5_w3_rawhash_sensitive.throughput
# echo "(POD5, Minimizer w = 5)"
# echo "RawHash throughput (mean/median)"
# grep "BP per sec:" d3_yeast_r94_pod5_w5_rawhash_sensitive.throughput

echo;
grep "Uncalled Mean time per read" d3_yeast_r94_rawhash_sensitive.comparison
grep "Sigmap Mean time per read" d3_yeast_r94_rawhash_sensitive.comparison
grep "RawHash Mean time per read" d3_yeast_r94_rawhash_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash Mean time per read" d3_yeast_r94_w3_rawhash_sensitive.comparison
# echo "(Minimizer w = 5)"
# grep "RawHash Mean time per read" d3_yeast_r94_w5_rawhash_sensitive.comparison
echo "(POD5)"
grep "RawHash Mean time per read" d3_yeast_r94_pod5_rawhash_sensitive.comparison
echo "(POD5, Minimizer w = 3)"
grep "RawHash Mean time per read" d3_yeast_r94_pod5_w3_rawhash_sensitive.comparison
# echo "(POD5, Minimizer w = 5)"
# grep "RawHash Mean time per read" d3_yeast_r94_pod5_w5_rawhash_sensitive.comparison

# echo;
# grep "RawHash Mean gap between read anchors in the best chain" d3_yeast_r94_rawhash_sensitive.comparison
# grep "RawHash Mean gap between reference anchors in the best chain" d3_yeast_r94_rawhash_sensitive.comparison

echo;
grep "Uncalled Mean # of sequenced bases per read" d3_yeast_r94_rawhash_sensitive.comparison
grep "RawHash Mean # of sequenced bases per read" d3_yeast_r94_rawhash_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash Mean # of sequenced bases per read" d3_yeast_r94_w3_rawhash_sensitive.comparison
# echo "(Minimizer w = 5)"
# grep "RawHash Mean # of sequenced bases per read" d3_yeast_r94_w5_rawhash_sensitive.comparison
echo "(POD5)"
grep "RawHash Mean # of sequenced bases per read" d3_yeast_r94_pod5_rawhash_sensitive.comparison
echo "(POD5, Minimizer w = 3)"
grep "RawHash Mean # of sequenced bases per read" d3_yeast_r94_pod5_w3_rawhash_sensitive.comparison
# echo "(POD5, Minimizer w = 5)"
# grep "RawHash Mean # of sequenced bases per read" d3_yeast_r94_pod5_w5_rawhash_sensitive.comparison

echo;
grep "Sigmap Mean # of sequenced chunks per read" d3_yeast_r94_rawhash_sensitive.comparison
grep "RawHash Mean # of sequenced chunks per read" d3_yeast_r94_rawhash_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash Mean # of sequenced chunks per read" d3_yeast_r94_w3_rawhash_sensitive.comparison
# echo "(Minimizer w = 5)"
# grep "RawHash Mean # of sequenced chunks per read" d3_yeast_r94_w5_rawhash_sensitive.comparison
echo "(POD5)"
grep "RawHash Mean # of sequenced chunks per read" d3_yeast_r94_pod5_rawhash_sensitive.comparison
echo "(POD5, Minimizer w = 3)"
grep "RawHash Mean # of sequenced chunks per read" d3_yeast_r94_pod5_w3_rawhash_sensitive.comparison
# echo "(POD5, Minimizer w = 5)"
# grep "RawHash Mean # of sequenced chunks per read" d3_yeast_r94_pod5_w5_rawhash_sensitive.comparison

echo;
grep "Uncalled Mean (only mapped) # of sequenced bases per read" d3_yeast_r94_rawhash_sensitive.comparison
grep "RawHash Mean (only mapped) # of sequenced bases per read" d3_yeast_r94_rawhash_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash Mean (only mapped) # of sequenced bases per read" d3_yeast_r94_w3_rawhash_sensitive.comparison
# echo "(Minimizer w = 5)"
# grep "RawHash Mean (only mapped) # of sequenced bases per read" d3_yeast_r94_w5_rawhash_sensitive.comparison
echo "(POD5)"
grep "RawHash Mean (only mapped) # of sequenced bases per read" d3_yeast_r94_pod5_rawhash_sensitive.comparison
echo "(POD5, Minimizer w = 3)"
grep "RawHash Mean (only mapped) # of sequenced bases per read" d3_yeast_r94_pod5_w3_rawhash_sensitive.comparison
# echo "(POD5, Minimizer w = 5)"
# grep "RawHash Mean (only mapped) # of sequenced bases per read" d3_yeast_r94_pod5_w5_rawhash_sensitive.comparison

echo;
grep "Sigmap Mean (only mapped) # of sequenced chunks per read" d3_yeast_r94_rawhash_sensitive.comparison
grep "RawHash Mean (only mapped) # of sequenced chunks per read" d3_yeast_r94_rawhash_sensitive.comparison
grep "RawHash Mean (only mapped) # of sequenced chunks per read" d3_yeast_r94_w3_rawhash_sensitive.comparison
# echo "(Minimizer w = 5)"
# grep "RawHash Mean (only mapped) # of sequenced chunks per read" d3_yeast_r94_w5_rawhash_sensitive.comparison
echo "(POD5)"
grep "RawHash Mean (only mapped) # of sequenced chunks per read" d3_yeast_r94_pod5_rawhash_sensitive.comparison
echo "(POD5, Minimizer w = 3)"
grep "RawHash Mean (only mapped) # of sequenced chunks per read" d3_yeast_r94_pod5_w3_rawhash_sensitive.comparison
# echo "(POD5, Minimizer w = 5)"
# grep "RawHash Mean (only mapped) # of sequenced chunks per read" d3_yeast_r94_pod5_w5_rawhash_sensitive.comparison

echo;
grep "Uncalled Mean (only unmapped) # of sequenced bases per read" d3_yeast_r94_rawhash_sensitive.comparison
grep "RawHash Mean (only unmapped) # of sequenced bases per read" d3_yeast_r94_rawhash_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash Mean (only unmapped) # of sequenced bases per read" d3_yeast_r94_w3_rawhash_sensitive.comparison
# echo "(Minimizer w = 5)"
# grep "RawHash Mean (only unmapped) # of sequenced bases per read" d3_yeast_r94_w5_rawhash_sensitive.comparison
echo "(POD5)"
grep "RawHash Mean (only unmapped) # of sequenced bases per read" d3_yeast_r94_pod5_rawhash_sensitive.comparison
echo "(POD5, Minimizer w = 3)"
grep "RawHash Mean (only unmapped) # of sequenced bases per read" d3_yeast_r94_pod5_w3_rawhash_sensitive.comparison
# echo "(POD5, Minimizer w = 5)"
# grep "RawHash Mean (only unmapped) # of sequenced bases per read" d3_yeast_r94_pod5_w5_rawhash_sensitive.comparison

echo;
grep "Sigmap Mean (only unmapped) # of sequenced chunks per read" d3_yeast_r94_rawhash_sensitive.comparison
grep "RawHash Mean (only unmapped) # of sequenced chunks per read" d3_yeast_r94_rawhash_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash Mean (only unmapped) # of sequenced chunks per read" d3_yeast_r94_w3_rawhash_sensitive.comparison
# echo "(Minimizer w = 5)"
# grep "RawHash Mean (only unmapped) # of sequenced chunks per read" d3_yeast_r94_w5_rawhash_sensitive.comparison
echo "(POD5)"
grep "RawHash Mean (only unmapped) # of sequenced chunks per read" d3_yeast_r94_pod5_rawhash_sensitive.comparison
echo "(POD5, Minimizer w = 3)"
grep "RawHash Mean (only unmapped) # of sequenced chunks per read" d3_yeast_r94_pod5_w3_rawhash_sensitive.comparison
# echo "(POD5, Minimizer w = 5)"
# grep "RawHash Mean (only unmapped) # of sequenced chunks per read" d3_yeast_r94_pod5_w5_rawhash_sensitive.comparison

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
grep "Uncalled precision:" d3_yeast_r94_rawhash_sensitive.comparison
grep "Sigmap precision:" d3_yeast_r94_rawhash_sensitive.comparison
grep "RawHash precision:" d3_yeast_r94_rawhash_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash precision:" d3_yeast_r94_w3_rawhash_sensitive.comparison
# echo "(Minimizer w = 5)"
# grep "RawHash precision:" d3_yeast_r94_w5_rawhash_sensitive.comparison
echo "(POD5)"
grep "RawHash precision:" d3_yeast_r94_pod5_rawhash_sensitive.comparison
echo "(POD5, Minimizer w = 3)"
grep "RawHash precision:" d3_yeast_r94_pod5_w3_rawhash_sensitive.comparison
# echo "(POD5, Minimizer w = 5)"
# grep "RawHash precision:" d3_yeast_r94_pod5_w5_rawhash_sensitive.comparison

echo;
grep "Uncalled recall:" d3_yeast_r94_rawhash_sensitive.comparison
grep "Sigmap recall:" d3_yeast_r94_rawhash_sensitive.comparison
grep "RawHash recall:" d3_yeast_r94_rawhash_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash recall:" d3_yeast_r94_w3_rawhash_sensitive.comparison
# echo "(Minimizer w = 5)"
# grep "RawHash recall:" d3_yeast_r94_w5_rawhash_sensitive.comparison
echo "(POD5)"
grep "RawHash recall:" d3_yeast_r94_pod5_rawhash_sensitive.comparison
echo "(POD5, Minimizer w = 3)"
grep "RawHash recall:" d3_yeast_r94_pod5_w3_rawhash_sensitive.comparison
# echo "(POD5, Minimizer w = 5)"
# grep "RawHash recall:" d3_yeast_r94_pod5_w5_rawhash_sensitive.comparison

echo;
grep "Uncalled F-1 score:" d3_yeast_r94_rawhash_sensitive.comparison
grep "Sigmap F-1 score:" d3_yeast_r94_rawhash_sensitive.comparison
grep "RawHash F-1 score:" d3_yeast_r94_rawhash_sensitive.comparison
echo "(Minimizer w = 3)"
grep "RawHash F-1 score:" d3_yeast_r94_w3_rawhash_sensitive.comparison
# echo "(Minimizer w = 5)"
# grep "RawHash F-1 score:" d3_yeast_r94_w5_rawhash_sensitive.comparison
echo "(POD5)"
grep "RawHash F-1 score:" d3_yeast_r94_pod5_rawhash_sensitive.comparison
echo "(POD5, Minimizer w = 3)"
grep "RawHash F-1 score:" d3_yeast_r94_pod5_w3_rawhash_sensitive.comparison
# echo "(POD5, Minimizer w = 5)"
# grep "RawHash F-1 score:" d3_yeast_r94_pod5_w5_rawhash_sensitive.comparison
