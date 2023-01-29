#!/bin/bash

uncalled pafstats -r ../true_mappings.paf --annotate ../uncalled/relative_abundance_uncalled.paf > relative_abundance_uncalled_ann.paf 2> relative_abundance_uncalled.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../sigmap/relative_abundance_sigmap.paf > relative_abundance_sigmap_ann.paf 2> relative_abundance_sigmap.throughput
uncalled pafstats -r ../true_mappings.paf --annotate ../rawhash/relative_abundance_rawhash_fast.paf > relative_abundance_rawhash_fast_ann.paf 2> relative_abundance_rawhash_fast.throughput

python ../../../scripts/compare_pafs.py relative_abundance_uncalled_ann.paf relative_abundance_sigmap_ann.paf relative_abundance_rawhash_fast_ann.paf > relative_abundance_rawhash_fast.comparison

for i in '../true_mappings.paf' '../uncalled/relative_abundance_uncalled.paf' '../sigmap/relative_abundance_sigmap.paf' '../rawhash/relative_abundance_rawhash_fast.paf' ; do
if test -f $i; then
fname=`basename $i | sed 's/.paf/.abundance/g'`;
awk 'BEGIN{tot=0; covid=0; human=0; ecoli=0; yeast=0; algae=0; totsum=0; covidsum=0; humansum=0; ecolisum=0; yeastsum=0; algaesum=0;}{
	if(substr($6,1,5) == "ecoli" && $3 != "*" && $4 != "*" && $4 >= $3){ecoli++; count++; ecolisum+=$4-$3}
	else if(substr($6,1,5) == "yeast" && $3 != "*" && $4 != "*" && $4 >= $3){yeast++; count++; yeastsum+=$4-$3}
	else if(substr($6,1,11) == "green_algae" && $3 != "*" && $4 != "*" && $4 >= $3){algae++; count++; algaesum+=$4-$3}
	else if(substr($6,1,5) == "human" && $3 != "*" && $4 != "*" && $4 >= $3){human++; count++; humansum+=$4-$3}
	else if(substr($6,1,5) == "covid" && $3 != "*" && $4 != "*" && $4 >= $3){covid++; count++; covidsum+=$4-$3}
	}END{
		totsum = ecolisum + yeastsum + algaesum + humansum + covidsum;
		print "Ratio of bases: covid:" covidsum/totsum " ecoli: " ecolisum/totsum " yeast: " yeastsum/totsum " green_algae: " algaesum/totsum " human: " humansum/totsum;
	}' $i > $fname
fi done;
