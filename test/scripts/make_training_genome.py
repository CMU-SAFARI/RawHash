import sys

def load_truth_set(filename):
    truth_mapping = {}
    with open(filename, 'r') as file:
        for line in file:
            if line.strip():  # Skip empty lines
                parts = line.strip().split('\t')
                read_name = parts[0]
                qs = int(parts[2])
                qe = int(parts[3])
                strand = parts[4]
                chromosome = parts[5]
                rs = int(parts[7])
                re = int(parts[8]) 
                # position_start = int(parts[7])  # Starting position of the mapping
                truth_mapping[read_name] = (qs, qe, rs, re, strand, chromosome)
    return truth_mapping

def replace_nan_with_zero(s):
    return s.replace('-nan', '0')

def check_mapping(truth_info, feature_info):
    tqs, tqe, trs, tre, tstrand, tchr = truth_info
    qs, qe, rs, re, strand, chr = feature_info

    if tstrand == strand and tchr == chr:
        return 1
    return 0

def label_feature_set(truth_set_filename, feature_set_filename, output_filename):
    truth_mapping = load_truth_set(truth_set_filename)
    with open(feature_set_filename, 'r') as feature_file:
        lines = feature_file.readlines()  # Read all lines into memory

    # Skip the first 8 lines and the last 6 lines
    lines = lines[8:-3]

    with open(output_filename, 'w') as output_file:
        for line in lines:
            if line.strip():  # Skip empty lines
                parts = line.strip().split('\t')
                # Apply replacement for '-nan' values in the parts
                # parts = [replace_nan_with_zero(part) for part in parts]
                read_name = parts[2]
                qs = int(parts[3])
                qe = int(parts[4])
                chromosome = parts[5]
                rs = int(parts[6])
                re = int(parts[7])
                strand = parts[8]
                # position = int(parts[4])
                if read_name in truth_mapping:
                    tqs, tqe, trs, tre, tstrand, tchr = truth_mapping[read_name]
                    is_correct = check_mapping((tqs, tqe, trs, tre, tstrand, tchr), (qs, qe, rs, re, strand, chromosome))
                    output_line = '\t'.join(parts) + '\t' + str(is_correct) + '\n'
                    output_file.write(output_line)
                else:
                    output_file.write('\t'.join(parts) + '\t0\n')  # Default to incorrect if not found in truth set

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python make_training.py <truth_set.paf> <feature_set.txt> <output_file.txt>")
        sys.exit(1)
    label_feature_set(sys.argv[1], sys.argv[2], sys.argv[3])
