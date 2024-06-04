#%%
from pathlib import Path
import cppyy
import os

# os.chdir("/home/mmordig/rawhash_project/rawhash2_new/build/extern")

header_file = Path("/home/mmordig/rawhash_project/rawhash2_new/src/rawhash_wrapper.hpp")
library_file = Path("/home/mmordig/rawhash_project/rawhash2_new/build/src/librawhash2_wrapper.so")
cppyy.include(str(header_file))
# cppyy.add_library_path("/home/mmordig/rawhash_project/rawhash2_new/build/extern/hdf5/lib")
cppyy.add_library_path("/home/mmordig/rawhash_project/rawhash2_new/build/src") # for shared libs
cppyy.load_library("librawhash2_wrapper")

# cppyy.add_library_path(str(library_file.parent))
# cppyy.load_library(str(library_file.name))
# cppyy.load_library(str(library_file))

# cppyy.gbl.RawHashMapper(5, [b"hello"])
args = [
    "-x", "sensitive", "-t", "8", "-p", 
    "/home/mmordig/rawhash_project/rawhash2_new/extern/kmer_models/legacy/legacy_r9.4_180mv_450bps_6mer/template_median68pA.model", 
]
# args += [
#     "/home/mmordig/rawhash_project/rawhash2_new/test/data/d2_ecoli_r94/ref.fa", 
#     "no-revcomp-query",
# ]
args += [
    "/home/mmordig/rawhash_project/rawhash2_new/example_out/forwrev_rawhash2_sensitive.ind", 
]
# forwrev_rawhash2_sensitive.paf
# 0006a106-77eb-4752-b405-38d1635718fa
# /home/mmordig/rawhash_project/rawhash2_new/example_out/forwonly_rawhash2_sensitive.ind
# /home/mmordig/rawhash_project/rawhash2_new/example_out/forwrev_rawhash2_sensitive.ind
args = ["my_dummy_program"] + args
mapper = cppyy.gbl.RawHashMapper(len(args), args)
mapper.idx_info()

#%%
slow5_filename = "/home/mmordig/rawhash_project/rawhash2_new/test/data/d2_ecoli_r94/small_slow5_files/barcode02_r0barcode02b0_0.blow5"
# we know rawhash finds an alignment since from the file /home/mmordig/rawhash_project/rawhash2_new/example_out/forwrev_rawhash2_sensitive.paf
read_id = "0006a106-77eb-4752-b405-38d1635718fa"
import contextlib
import pyslow5
s5 = pyslow5.Open(slow5_filename, "r")
raw_signal = s5.get_read(read_id, pA=True)["signal"]
raw_signal = raw_signal[(raw_signal > 30.) & (raw_signal < 200.)]
# with contextlib.autoclosing(pyslow5.Open(slow5_filename, "r")) as s5:
raw_signal.shape
#%%
import numpy as np

# Create a numpy array
# arr = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)
# # Ensure the numpy array is contiguous
# arr = np.ascontiguousarray(arr, dtype=np.float32)

# arr = raw_signal.astype(np.float32)
arr = np.ascontiguousarray(raw_signal, dtype=np.float32)

import cppyy.ll
# cppyy.set_debug()
# Call the C++ function with the numpy array
alignments = mapper.map(cppyy.ll.cast["float*"](arr.ctypes.data), len(arr))

class Alignment:
    def __init__(self, contig, ref_start, ref_end, is_pos_strand):
        self.contig = contig
        self.ref_start = ref_start
        self.ref_end = ref_end
        self.is_pos_strand = is_pos_strand
    @staticmethod
    def from_cppyy(alignment):
        return Alignment(
            contig=alignment.contig,
            ref_start=alignment.ref_start,
            ref_end=alignment.ref_end,
            is_pos_strand=alignment.is_pos_strand,
        )
    def __repr__(self):
        return f"Alignment(contig: {self.contig}, start: {self.ref_start}, end: {self.ref_end}, pos_strand: {self.is_pos_strand})"
    
def parse_paf_line(line):
    """returns read_id, alignment"""
    fields = line.strip().split("\t")
    return fields[0], Alignment(
        contig=fields[5],
        ref_start=int(fields[7]),
        ref_end=int(fields[8]),
        is_pos_strand=fields[4] == "+",
    )
# def format_alignment(alignment):
#     # return alignment
#     return f"Alignment(contig: {alignment.contig}, start: {alignment.ref_start}, end: {alignment.ref_end}, pos_strand: {alignment.is_pos_strand})"

for (i, alignment) in enumerate(alignments):
    print(f"Alignment {i}: {Alignment.from_cppyy(alignment)}")
    
gt_paf_filename = "/home/mmordig/rawhash_project/rawhash2_new/example_out/forwrev_rawhash2_sensitive.paf"
with open(gt_paf_filename) as f:
    for line in f:
        rid, gt_alignment = parse_paf_line(line)
        if rid == read_id:
            print(f"GT Alignment for {rid}: {gt_alignment}")
            break
#%%
1

# --rev-query