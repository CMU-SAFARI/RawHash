import os
import sys
import argparse
import subprocess
import fileinput

from statistics import median
from statistics import mean

if (len(sys.argv) < 5):
  print("usage: compare_pafs.py uncalled_ann.paf sigmap_ann.paf rawhash_ann.paf rawhash2_ann.paf")
  sys.exit(1)

fps = []
for tool in [1, 2, 3, 4]:
  fps.append(open(sys.argv[tool]))

uncalled_tp = 0
uncalled_fp = 0
uncalled_fn = 0
uncalled_tn = 0

uncalled_chunk = []
uncalled_time_per_read = []
uncalled_maplast_pos = []
uncalled_umaplast_pos = []
for line in fps[0]:
  cols = line.rstrip().split()
  mt = float(cols[14].split(":")[2])
  lastpos = int(cols[1])
  if (cols[15].split(":")[2] != 'na'):
    uncalled_time_per_read.append(mt)
    if(cols[2] != '*'):
      uncalled_maplast_pos.append(lastpos)
    else:
      uncalled_umaplast_pos.append(lastpos)
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
if (uncalled_tp + uncalled_fp == 0):
   uncalled_precision = 0
else:
  uncalled_precision = uncalled_tp / (uncalled_tp + uncalled_fp)
print("Uncalled precision: " + str(uncalled_precision))
if (uncalled_tp + uncalled_fn == 0):
  uncalled_recall = 0
else:
   uncalled_recall = uncalled_tp / (uncalled_tp + uncalled_fn)
print("Uncalled recall: " + str(uncalled_recall))
if (uncalled_precision + uncalled_recall == 0):
  print("Uncalled F-1 score: 0")
else:
  print("Uncalled F-1 score: " + str(2 * uncalled_precision * uncalled_recall / (uncalled_precision + uncalled_recall)))
print("Uncalled Mean time per read : " + str(mean(uncalled_time_per_read)))
print("Uncalled Median time per read : " + str(median(uncalled_time_per_read)))
print("Uncalled Mean # of sequenced bases per read : " + str(mean(uncalled_maplast_pos + uncalled_umaplast_pos)))
print("Uncalled Mean (only mapped) # of sequenced bases per read : " + str(mean(uncalled_maplast_pos)))
print("Uncalled Mean (only unmapped) # of sequenced bases per read : " + str(mean(uncalled_umaplast_pos)))
print("#Done with uncalled\n")

fps[0].close()

sigmap_tp = 0
sigmap_fp = 0
sigmap_fn = 0
sigmap_tn = 0

sigmap_time_per_chunk = []
sigmap_time_per_read = []
sigmap_maplast_chunk = []
sigmap_umaplast_chunk = []
for line in fps[1]:
  cols = line.rstrip().split()
  if (len(cols) == 24):
    mt = float(cols[12].split(":")[2])
    if (cols[23].split(":")[2] != 'na'):
      sigmap_time_per_read.append(mt)
    chunk = int(cols[13].split(":")[2])
    if(cols[2] != '*'):
      sigmap_maplast_chunk.append(chunk)
    else:
      sigmap_umaplast_chunk.append(chunk)
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
if (sigmap_tp + sigmap_fp == 0):
    sigmap_precision = 0
else:
  sigmap_precision = sigmap_tp / (sigmap_tp + sigmap_fp)
print("Sigmap precision: " + str(sigmap_precision))

if (sigmap_tp + sigmap_fn == 0):
  sigmap_recall = 0
else:
  sigmap_recall = sigmap_tp / (sigmap_tp + sigmap_fn)
print("Sigmap recall: " + str(sigmap_recall))
if (sigmap_precision + sigmap_recall == 0):
  print("Sigmap F-1 score: 0")
else:
  print("Sigmap F-1 score: " + str(2 * sigmap_precision * sigmap_recall / (sigmap_precision + sigmap_recall)))
print("Sigmap Mean time per mapped read : " + str(mean(sigmap_time_per_chunk)))
print("Sigmap Median time per mapped read : " + str(median(sigmap_time_per_chunk)))
print("Sigmap Mean time per unmapped read : " + str(mean(sigmap_time_per_read)))
print("Sigmap Median time per unmapped read : " + str(median(sigmap_time_per_read)))
print("Sigmap Mean time per read : " + str(mean(sigmap_time_per_read + sigmap_time_per_chunk)))
print("Sigmap Median time per read : " + str(median(sigmap_time_per_read + sigmap_time_per_chunk)))
print("Sigmap Mean # of sequenced chunks per read : " + str(mean(sigmap_maplast_chunk + sigmap_umaplast_chunk)))
print("Sigmap Mean (only mapped) # of sequenced chunks per read : " + str(mean(sigmap_maplast_chunk)))
print("Sigmap Mean (only unmapped) # of sequenced chunks per read : " + str(mean(sigmap_umaplast_chunk)))
print("#Done with sigmap\n")

fps[1].close()

rawhash_tp = 0
rawhash_fp = 0
rawhash_fn = 0
rawhash_tn = 0

rawhash_time_per_chunk = []
rawhash_time_per_read = []
rawhash_maplast_pos = []
rawhash_maplast_chunk = []
rawhash_umaplast_pos = []
rawhash_umaplast_chunk = []
rawhash_refgap = []
rawhash_readgap = []
for line in fps[2]:
  cols = line.rstrip().split()
  if (len(cols) == 23):
    mt = float(cols[12].split(":")[2])
    lastpos = int(cols[1])
    if (cols[22].split(":")[2] != 'na'):
      rawhash_time_per_read.append(mt)
      if(cols[2] != '*'):
        rawhash_maplast_pos.append(lastpos)
      else:
        rawhash_umaplast_pos.append(lastpos)
    chunk = int(cols[13].split(":")[2])
    if(cols[2] != '*'):
      rawhash_maplast_chunk.append(chunk)
    else:
      rawhash_umaplast_chunk.append(chunk)
    cm = int(cols[15].split(":")[2])
    nc = int(cols[16].split(":")[2])
    s1 = float(cols[17].split(":")[2])
    s2 = float(cols[18].split(":")[2])
    sm = float(cols[19].split(":")[2])
    at = float(cols[20].split(":")[2])
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
if (rawhash_tp + rawhash_fp == 0):
    rawhash_precision = 0
else:
  rawhash_precision = rawhash_tp / (rawhash_tp + rawhash_fp)
print("RawHash precision: " + str(rawhash_precision))
if (rawhash_tp + rawhash_fn == 0):
  rawhash_recall = 0
else:
  rawhash_recall = rawhash_tp / (rawhash_tp + rawhash_fn)
print("RawHash recall: " + str(rawhash_recall))
if (rawhash_precision + rawhash_recall == 0):
  print("RawHash F-1 score: 0")
else:
  print("RawHash F-1 score: " + str(2 * rawhash_precision * rawhash_recall / (rawhash_precision + rawhash_recall)))
print("RawHash Mean time per mapped read : " + str(mean(rawhash_time_per_chunk)))
print("RawHash Median time per mapped read : " + str(median(rawhash_time_per_chunk)))
print("RawHash Mean time per unmapped read : " + str(mean(rawhash_time_per_read)))
print("RawHash Median time per unmapped read : " + str(median(rawhash_time_per_read)))
print("RawHash Mean time per read : " + str(mean(rawhash_time_per_read + rawhash_time_per_chunk)))
print("RawHash Median time per read : " + str(median(rawhash_time_per_read + rawhash_time_per_chunk)))
print("RawHash Mean # of sequenced bases per read : " + str(mean(rawhash_maplast_pos + rawhash_umaplast_pos)))
print("RawHash Mean # of sequenced chunks per read : " + str(mean(rawhash_maplast_chunk + rawhash_umaplast_chunk)))

print("RawHash Mean (only mapped) # of sequenced bases per read : " + str(mean(rawhash_maplast_pos)))
print("RawHash Mean (only mapped) # of sequenced chunks per read : " + str(mean(rawhash_maplast_chunk)))

print("RawHash Mean (only unmapped) # of sequenced bases per read : " + str(mean(rawhash_umaplast_pos)))
print("RawHash Mean (only unmapped) # of sequenced chunks per read : " + str(mean(rawhash_umaplast_chunk)))

print("#Done with RawHash\n")

fps[2].close()

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
for line in fps[3]:
  cols = line.rstrip().split()
  if (len(cols) == 20):
    mt = float(cols[12].split(":")[2])
    lastpos = int(cols[1])
    if (cols[19].split(":")[2] != 'na'):
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
    # s2 = float(cols[18].split(":")[2])
    sm = float(cols[18].split(":")[2])
    if (cols[19].split(":")[2] == 'tp'):
      rawhash2_tp += 1
      rawhash2_time_per_chunk.append(mt / chunk)
    if (cols[19].split(":")[2] == 'fp' or cols[19].split(":")[2] == 'na'):
      rawhash2_fp += 1
      rawhash2_time_per_chunk.append(mt / chunk)
    if (cols[19].split(":")[2] == 'fn'):
      rawhash2_fn += 1
      rawhash2_time_per_chunk.append(mt / chunk)
    if (cols[19].split(":")[2] == 'tn'):
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
print("RawHash2 TP: " + str(rawhash2_tp))
print("RawHash2 FP: " + str(rawhash2_fp))
print("RawHash2 FN: " + str(rawhash2_fn))
print("RawHash2 TN: " + str(rawhash2_tn))
if (rawhash2_tp + rawhash2_fp == 0):
    rawhash2_precision = 0
else:
  rawhash2_precision = rawhash2_tp / (rawhash2_tp + rawhash2_fp)
print("RawHash2 precision: " + str(rawhash2_precision))
if (rawhash2_tp + rawhash2_fn == 0):
  rawhash2_recall = 0
else:
  rawhash2_recall = rawhash2_tp / (rawhash2_tp + rawhash2_fn)
print("RawHash2 recall: " + str(rawhash2_recall))
if (rawhash2_precision + rawhash2_recall == 0):
  print("RawHash2 F-1 score: 0")
else:
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

fps[3].close()