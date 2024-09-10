#!/usr/bin/env python

import pandas as pd
import numpy as np
import argparse
import sys

import pandas as pd

import pandas as pd

def read_paf(file_path, is_input_paf):
    previous_query_id = None  # To track the previous read ID
    mapped_pairs = set()
    unmapped_pairs = set()
    bases = []
    signals = []
    mt_ms = []

    mt_col_index = None  # To store the index of the mt tag column
    sl_col_index = None  # To store the index of the sl tag column

    with open(file_path, 'r') as file:
        for line_number, line in enumerate(file):
            parts = line.strip().split('\t')
            query_id = parts[0]
            reference_id = parts[5]

            # Dynamically find the mt column on the first line
            if line_number == 0 and is_input_paf:
                for index, part in enumerate(parts):
                    if part.startswith('mt:f:'):
                        mt_col_index = index
                    elif part.startswith('sl:i:'):
                        sl_col_index = index
                if mt_col_index is None:
                    raise ValueError("mt tag not found in the input PAF file.")
                continue  # Skip the first line after finding mt column index

            # Accuracy metrics: count all mappings
            if reference_id != '*':
                mapped_pairs.add((query_id, reference_id))
            else:
                unmapped_pairs.add((query_id, reference_id))

            # Throughput calculations: check if current read ID is different from the previous one
            if is_input_paf and query_id != previous_query_id:
                if mt_col_index is not None and len(parts) > mt_col_index:
                    mt_value = float(parts[mt_col_index].split(':')[-1])
                    mt_ms.append(mt_value)
                    bases.append(int(parts[1]))
                    if sl_col_index is not None and len(parts) > sl_col_index:
                        signals.append(int(parts[sl_col_index].split(':')[-1]))

            previous_query_id = query_id  # Update previous_query_id to the current read ID

    return mapped_pairs, unmapped_pairs, bases, mt_ms, signals

def classify_pairs(input_mapped, input_unmapped, minimap2_mapped, minimap2_unmapped):
    results = []
    counts = {'tp': 0, 'fp': 0, 'fn': 0, 'tn': 0}
    all_pairs = input_mapped.union(input_unmapped, minimap2_mapped, minimap2_unmapped)

    for pair in all_pairs:
        if pair in input_mapped:
            if pair in minimap2_mapped:
                results.append((*pair, 'rf:Z:tp'))
                counts['tp'] += 1
            else:
                results.append((*pair, 'rf:Z:fp'))
                counts['fp'] += 1
        elif pair in minimap2_mapped:
            results.append((*pair, 'rf:Z:fn'))
            counts['fn'] += 1
        else:
            results.append((*pair, 'rf:Z:tn'))
            counts['tn'] += 1

    return results, counts

import numpy as np
import argparse
import sys

def compute_throughput(bases, signals, mt_ms):
    if not bases:
        return 0, 0, 0, 0  # Avoid division by zero if no data
    bps = [1000 * b / mt for b, mt in zip(bases, mt_ms)]  # Calculate bases per second
    if not signals:
        return np.mean(bps), np.median(bps), np.mean(mt_ms), np.median(mt_ms), 0, 0
    sps = [1000 * s / mt for s, mt in zip(signals, mt_ms)]  # Calculate signals per second
    mean_bps = np.mean(bps)
    median_bps = np.median(bps)
    mean_mt = np.mean(mt_ms)
    median_mt = np.median(mt_ms)
    mean_sps = np.mean(sps)
    median_sps = np.median(sps)
    return mean_bps, median_bps, mean_mt, median_mt, mean_sps, median_sps

def process_files(input_path, minimap2_path):
    input_mapped, input_unmapped, bases, mt_ms, signals = read_paf(input_path, True)
    minimap2_mapped, minimap2_unmapped, _, _, _ = read_paf(minimap2_path, False)
    
    classified_results, counts = classify_pairs(input_mapped, input_unmapped, minimap2_mapped, minimap2_unmapped)

    total = sum(counts.values())
    ratios = {k: v / total for k, v in counts.items()}
    sys.stderr.write(f"Counts: {counts}\n")
    sys.stderr.write(f"Ratios: {ratios}\n")
    sys.stderr.write(f"TP: {counts['tp']}, FP: {counts['fp']}, FN: {counts['fn']}, TN: {counts['tn']}\n")
    sys.stderr.write(f"Precision: {counts['tp'] / (counts['tp'] + counts['fp']) if counts['tp'] + counts['fp'] > 0 else 0}\n")
    sys.stderr.write(f"Recall: {counts['tp'] / (counts['tp'] + counts['fn']) if counts['tp'] + counts['fn'] > 0 else 0}\n")
    sys.stderr.write(f"F1 Score: {2 * counts['tp'] / (2 * counts['tp'] + counts['fp'] + counts['fn']) if 2 * counts['tp'] + counts['fp'] + counts['fn'] > 0 else 0}\n")

    mean_bps, median_bps, mean_mt, median_mt, mean_sps, median_sps = compute_throughput(bases, signals, mt_ms)
    sys.stderr.write(f"Speed            Mean    Median\n")
    sys.stderr.write(f"BP per sec: %9.2f %9.2f\n" % (mean_bps, median_bps))
    sys.stderr.write(f"Signals per sec: %9.2f %9.2f\n" % (mean_sps, median_sps))
    sys.stderr.write(f"MS to map:  %9.2f %9.2f\n" % (mean_mt, median_mt))

    # for result in classified_results:
    #     print('\t'.join(map(str, result)))

def parse_arguments():
    parser = argparse.ArgumentParser(description='Compare PAF files to classify read pairs.')
    parser.add_argument('input_path', type=str, help='Path to input.paf file.')
    parser.add_argument('minimap2_path', type=str, help='Path to minimap2.paf file.')
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    process_files(args.input_path, args.minimap2_path)
