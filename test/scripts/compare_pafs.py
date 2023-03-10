import os
import sys
import argparse
import subprocess
import fileinput

from statistics import median
from statistics import mean

if (len(sys.argv) < 4):
  print("usage: eval_rawhash.py uncalled_ann.paf sigmap_ann.paf rawhash_ann.paf")
  sys.exit(1)

fps = []
for tool in [1, 2, 3]:
  fps.append(open(sys.argv[tool]))

uncalled_tp = 0
uncalled_fp = 0
uncalled_fn = 0
uncalled_tn = 0

uncalled_chunk = []
uncalled_time_per_read = []
for line in fps[0]:
  cols = line.rstrip().split()
  mt = float(cols[14].split(":")[2])
  if (cols[15].split(":")[2] != 'na'):
    uncalled_time_per_read.append(mt)
  if (cols[15].split(":")[2] == 'tp'):
    uncalled_tp += 1
  if (cols[15].split(":")[2] == 'fp' or cols[15].split(":")[2] == 'na'):
    uncalled_fp += 1
  if (cols[15].split(":")[2] == 'fn'):
    uncalled_fn += 1
  if (cols[15].split(":")[2] == 'tn'):
    uncalled_tn += 1
print("Uncalled TP: " + str(uncalled_tp))
print("Uncalled FP: " + str(uncalled_fp))
print("Uncalled FN: " + str(uncalled_fn))
print("Uncalled TN: " + str(uncalled_tn))
uncalled_precision = uncalled_tp / (uncalled_tp + uncalled_fp)
print("Uncalled precision: " + str(uncalled_precision))
uncalled_recall = uncalled_tp / (uncalled_tp + uncalled_fn)
print("Uncalled recall: " + str(uncalled_recall))
print("Uncalled F-1 score: " + str(2 * uncalled_precision * uncalled_recall / (uncalled_precision + uncalled_recall)))
print("Uncalled Mean time per read : " + str(mean(uncalled_time_per_read)))
print("Uncalled Median time per read : " + str(median(uncalled_time_per_read)))
print("#Done with uncalled\n")

sigmap_tp = 0
sigmap_fp = 0
sigmap_fn = 0
sigmap_tn = 0

sigmap_time_per_chunk = []
sigmap_time_per_read = []
for line in fps[1]:
  cols = line.rstrip().split()
  if (len(cols) == 24):
    mt = float(cols[12].split(":")[2])
    if (cols[23].split(":")[2] != 'na'):
      sigmap_time_per_read.append(mt)
    chunk = int(cols[13].split(":")[2])
    cm = int(cols[15].split(":")[2])
    nc = int(cols[16].split(":")[2])
    s1 = float(cols[17].split(":")[2])
    s2 = float(cols[18].split(":")[2])
    sm = float(cols[19].split(":")[2])
    ad = float(cols[20].split(":")[2])
    at = float(cols[21].split(":")[2])
    aq = float(cols[22].split(":")[2])
    if (cols[23].split(":")[2] == 'tp'):
      sigmap_tp += 1
      sigmap_time_per_chunk.append(mt / chunk)
    if (cols[23].split(":")[2] == 'fp' or cols[23].split(":")[2] == 'na'):
      sigmap_fp += 1
      sigmap_time_per_chunk.append(mt / chunk)
    if (cols[23].split(":")[2] == 'fn'):
      sigmap_fn += 1
      sigmap_time_per_chunk.append(mt / chunk)
    if (cols[23].split(":")[2] == 'tn'):
      sigmap_tn += 1
      sigmap_time_per_chunk.append(mt / chunk)
  if (len(cols) == 15):
    mt = float(cols[12].split(":")[2])
    if (cols[14].split(":")[2] != 'na'):
      sigmap_time_per_read.append(mt)
    if (cols[14].split(":")[2] == 'fn'):
      sigmap_fn += 1
    if (cols[14].split(":")[2] == 'tn'):
      sigmap_tn += 1
print("Sigmap TP: " + str(sigmap_tp))
print("Sigmap FP: " + str(sigmap_fp))
print("Sigmap FN: " + str(sigmap_fn))
print("Sigmap TN: " + str(sigmap_tn))
sigmap_precision = sigmap_tp / (sigmap_tp + sigmap_fp)
print("Sigmap precision: " + str(sigmap_precision))
sigmap_recall = sigmap_tp / (sigmap_tp + sigmap_fn)
print("Sigmap recall: " + str(sigmap_recall))
print("Sigmap F-1 score: " + str(2 * sigmap_precision * sigmap_recall / (sigmap_precision + sigmap_recall)))
print("Sigmap Mean time per mapped read : " + str(mean(sigmap_time_per_chunk)))
print("Sigmap Median time per mapped read : " + str(median(sigmap_time_per_chunk)))
print("Sigmap Mean time per unmapped read : " + str(mean(sigmap_time_per_read)))
print("Sigmap Median time per unmapped read : " + str(median(sigmap_time_per_read)))
print("Sigmap Mean time per read : " + str(mean(sigmap_time_per_read + sigmap_time_per_chunk)))
print("Sigmap Median time per read : " + str(median(sigmap_time_per_read + sigmap_time_per_chunk)))
print("#Done with sigmap\n")

rawhash_tp = 0
rawhash_fp = 0
rawhash_fn = 0
rawhash_tn = 0

rawhash_time_per_chunk = []
rawhash_time_per_read = []
for line in fps[2]:
  cols = line.rstrip().split()
  if (len(cols) == 23):
    mt = float(cols[12].split(":")[2])
    if (cols[22].split(":")[2] != 'na'):
      rawhash_time_per_read.append(mt)
    chunk = int(cols[13].split(":")[2])
    cm = int(cols[15].split(":")[2])
    nc = int(cols[16].split(":")[2])
    s1 = float(cols[17].split(":")[2])
    s2 = float(cols[18].split(":")[2])
    sm = float(cols[19].split(":")[2])
    at = float(cols[20].split(":")[2])
    aq = float(cols[21].split(":")[2])
    if (cols[22].split(":")[2] == 'tp'):
      rawhash_tp += 1
      rawhash_time_per_chunk.append(mt / chunk)
    if (cols[22].split(":")[2] == 'fp' or cols[22].split(":")[2] == 'na'):
      rawhash_fp += 1
      rawhash_time_per_chunk.append(mt / chunk)
    if (cols[22].split(":")[2] == 'fn'):
      rawhash_fn += 1
      rawhash_time_per_chunk.append(mt / chunk)
    if (cols[22].split(":")[2] == 'tn'):
      rawhash_tn += 1
      rawhash_time_per_chunk.append(mt / chunk)
  if (len(cols) == 15):
    mt = float(cols[12].split(":")[2])
    if (cols[14].split(":")[2] != 'na'):
      rawhash_time_per_read.append(mt)
    if (cols[14].split(":")[2] == 'fn'):
      rawhash_fn += 1
    if (cols[14].split(":")[2] == 'tn'):
      rawhash_tn += 1
print("RawHash TP: " + str(rawhash_tp))
print("RawHash FP: " + str(rawhash_fp))
print("RawHash FN: " + str(rawhash_fn))
print("RawHash TN: " + str(rawhash_tn))
rawhash_precision = rawhash_tp / (rawhash_tp + rawhash_fp)
print("RawHash precision: " + str(rawhash_precision))
rawhash_recall = rawhash_tp / (rawhash_tp + rawhash_fn)
print("RawHash recall: " + str(rawhash_recall))
print("RawHash F-1 score: " + str(2 * rawhash_precision * rawhash_recall / (rawhash_precision + rawhash_recall)))
print("RawHash Mean time per mapped read : " + str(mean(rawhash_time_per_chunk)))
print("RawHash Median time per mapped read : " + str(median(rawhash_time_per_chunk)))
print("RawHash Mean time per unmapped read : " + str(mean(rawhash_time_per_read)))
print("RawHash Median time per unmapped read : " + str(median(rawhash_time_per_read)))
print("RawHash Mean time per read : " + str(mean(rawhash_time_per_read + rawhash_time_per_chunk)))
print("RawHash Median time per read : " + str(median(rawhash_time_per_read + rawhash_time_per_chunk)))
print("#Done with RawHash\n")

for fp in fps:
  fp.close()
