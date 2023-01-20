#!/bin/bash

awk '{if(substr($1,1,1) == ">"){print ">ecoli_"substr($1,2)}else{print $0}}' ../d2_ecoli_r94/ref.fa > ref.fa
awk '{if(substr($1,1,1) == ">"){print ">yeast_"substr($1,2)}else{print $0}}' ../d3_yeast_r94/ref.fa >> ref.fa
awk '{if(substr($1,1,1) == ">"){print ">green_algae_"substr($1,2)}else{print $0}}' ../d4_green_algae_r94/ref.fa >> ref.fa
awk '{if(substr($1,1,1) == ">"){print ">human_"substr($1,2)}else{print $0}}' ../d5_human_na12878_r94/ref.fa >> ref.fa
awk '{if(substr($1,1,1) == ">"){print ">covid_"substr($1,2)}else{print $0}}' ../d1_sars-cov-2_r94/ref.fa >> ref.fa
