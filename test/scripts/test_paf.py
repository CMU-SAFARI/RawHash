import os
import sys
import argparse
import subprocess
import fileinput

from statistics import median
from statistics import mean

if (len(sys.argv) < 2):
  print("usage: test_paf.py rawhash_ann.paf")
  sys.exit(1)

fps = []
for tool in [1]:
  fps.append(open(sys.argv[tool]))

rawhash2_tp = 0
rawhash2_fp = 0
rawhash2_fn = 0
rawhash2_tn = 0

rawhash2_time_per_chunk = []
rawhash2_time_per_read = []
rawhash2_maplast_pos = []
rawhash2_maplast_chunk = []
rawhash2_umaplast_pos = []
rawhash2_umaplast_chunk = []
rawhash2_refgap = []
rawhash2_readgap = []
for line in fps[0]:
  cols = line.rstrip().split()
  if (len(cols) == 21):
    mt = float(cols[12].split(":")[2])
    lastpos = int(cols[1])
    if (cols[20].split(":")[2] != 'na'):
      rawhash2_time_per_read.append(mt)
      if(cols[2] != '*'):
        rawhash2_maplast_pos.append(lastpos)
      else:
        rawhash2_umaplast_pos.append(lastpos)
    chunk = int(cols[13].split(":")[2])
    if(cols[2] != '*'):
      rawhash2_maplast_chunk.append(chunk)
    else:
      rawhash2_umaplast_chunk.append(chunk)
    cm = int(cols[15].split(":")[2])
    nc = int(cols[16].split(":")[2])
    s1 = float(cols[17].split(":")[2])
    s2 = float(cols[18].split(":")[2])
    sm = float(cols[19].split(":")[2])
    if (cols[20].split(":")[2] == 'tp'):
      rawhash2_tp += 1
      rawhash2_time_per_chunk.append(mt / chunk)
    if (cols[20].split(":")[2] == 'fp' or cols[20].split(":")[2] == 'na'):
      rawhash2_fp += 1
      rawhash2_time_per_chunk.append(mt / chunk)
    if (cols[20].split(":")[2] == 'fn'):
      rawhash2_fn += 1
      rawhash2_time_per_chunk.append(mt / chunk)
    if (cols[20].split(":")[2] == 'tn'):
      rawhash2_tn += 1
      rawhash2_time_per_chunk.append(mt / chunk)
  if (len(cols) == 15):
    mt = float(cols[12].split(":")[2])
    if (cols[14].split(":")[2] != 'na'):
      rawhash2_time_per_read.append(mt)
    if (cols[14].split(":")[2] == 'fn'):
      rawhash2_fn += 1
    if (cols[14].split(":")[2] == 'tn'):
      rawhash2_tn += 1

fps[0].close()

print("RawHash2 TP: " + str(rawhash2_tp))
print("RawHash2 FP: " + str(rawhash2_fp))
print("RawHash2 FN: " + str(rawhash2_fn))
print("RawHash2 TN: " + str(rawhash2_tn))
rawhash2_precision = rawhash2_tp / (rawhash2_tp + rawhash2_fp)
print("RawHash2 precision: " + str(rawhash2_precision))
rawhash2_recall = rawhash2_tp / (rawhash2_tp + rawhash2_fn)
print("RawHash2 recall: " + str(rawhash2_recall))
print("RawHash2 F-1 score: " + str(2 * rawhash2_precision * rawhash2_recall / (rawhash2_precision + rawhash2_recall)))
print("RawHash2 Mean time per mapped read : " + str(mean(rawhash2_time_per_chunk)))
print("RawHash2 Median time per mapped read : " + str(median(rawhash2_time_per_chunk)))
print("RawHash2 Mean time per unmapped read : " + str(mean(rawhash2_time_per_read)))
print("RawHash2 Median time per unmapped read : " + str(median(rawhash2_time_per_read)))
print("RawHash2 Mean time per read : " + str(mean(rawhash2_time_per_read + rawhash2_time_per_chunk)))
print("RawHash2 Median time per read : " + str(median(rawhash2_time_per_read + rawhash2_time_per_chunk)))
print("RawHash2 Mean # of sequenced bases per read : " + str(mean(rawhash2_maplast_pos + rawhash2_umaplast_pos)))
print("RawHash2 Mean # of sequenced chunks per read : " + str(mean(rawhash2_maplast_chunk + rawhash2_umaplast_chunk)))

print("RawHash2 Mean (only mapped) # of sequenced bases per read : " + str(mean(rawhash2_maplast_pos)))
print("RawHash2 Mean (only mapped) # of sequenced chunks per read : " + str(mean(rawhash2_maplast_chunk)))

print("RawHash2 Mean (only unmapped) # of sequenced bases per read : " + str(mean(rawhash2_umaplast_pos)))
print("RawHash2 Mean (only unmapped) # of sequenced chunks per read : " + str(mean(rawhash2_umaplast_chunk)))

print("#Done with RawHash2\n")