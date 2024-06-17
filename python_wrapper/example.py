#%%

# Make sure to compile rawhash properly, both the CLI and the shared lib!!! Otherwise, the kernel may crash

from pathlib import Path
import sys
sys.path.extend(["/home/mmordig/rawhash_project/rawhash2_new/python_wrapper"])
from rawhash_wrapper import get_read_signal, load_rawhash_wrapper_lib, get_rawhash_mapper, prepare_signal_for_rawhash, get_rawhash_alignments, parse_paf_line

#%%
rawhash_args = [
    "-x", "sensitive", "-t", "8", "-p", 
    "/home/mmordig/rawhash_project/rawhash2_new/extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model", 
]
# rawhash_args += [
#     "/home/mmordig/rawhash_project/rawhash2_new/test/data/d2_ecoli_r94/ref.fa", 
#     "no-revcomp-query",
# ]
rawhash_args += [
    "/home/mmordig/rawhash_project/rawhash2_new/example_out/forwrev_rawhash2_sensitive.ind", 
]
# forwrev_rawhash2_sensitive.paf
# 0006a106-77eb-4752-b405-38d1635718fa
# /home/mmordig/rawhash_project/rawhash2_new/example_out/forwonly_rawhash2_sensitive.ind
# /home/mmordig/rawhash_project/rawhash2_new/example_out/forwrev_rawhash2_sensitive.ind


#%%

slow5_filename = "/home/mmordig/rawhash_project/rawhash2_new/test/data/d2_ecoli_r94/small_slow5_files/barcode02_r0barcode02b0_0.blow5"
# we know rawhash finds an alignment for this read from the file /home/mmordig/rawhash_project/rawhash2_new/example_out/forwrev_rawhash2_sensitive.paf
read_id = "0006a106-77eb-4752-b405-38d1635718fa"

raw_signal = get_read_signal(slow5_filename, read_id)
print(f"Raw signal shape: {raw_signal.shape}")

rawhash_lib_dir = Path("/home/mmordig/rawhash_project/rawhash2_new")
load_rawhash_wrapper_lib(rawhash_lib_dir)
mapper = get_rawhash_mapper(rawhash_args)
mapper.print_idx_info()
raw_signal = prepare_signal_for_rawhash(raw_signal)
# raw_signal = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)
print(f"Raw signal shape: {raw_signal.shape}")
alignments = get_rawhash_alignments(mapper, raw_signal)
#%%

for (i, alignment) in enumerate(alignments):
    print(f"Alignment {i}: {alignment}")
    
gt_paf_filename = "/home/mmordig/rawhash_project/rawhash2_new/example_out/forwrev_rawhash2_sensitive.paf"
with open(gt_paf_filename) as f:
    for line in f:
        rid, gt_alignment = parse_paf_line(line)
        if rid == read_id:
            print(f"GT Alignment for {rid}:\n             {gt_alignment}")
            break
#%%
