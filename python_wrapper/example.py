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
    # "-d", "example_out/ref.ind", 
    "/home/mmordig/rawhash_project/rawhash2_new/test/data/d2_ecoli_r94/ref.fa", 
    "no-revcomp-query",
]
# forwrev_rawhash2_sensitive.paf
# 0006a106-77eb-4752-b405-38d1635718fa
args = ["my_dummy_program"] + args
mapper = cppyy.gbl.RawHashMapper(len(args), args)
mapper.idx_info()
#%%
import numpy as np

# Create a numpy array
arr = np.array([1.0, 2.0, 3.0, 4.0], dtype=np.float32)

# Ensure the numpy array is contiguous
arr = np.ascontiguousarray(arr, dtype=np.float32)

import cppyy.ll
cppyy.set_debug()
# Call the C++ function with the numpy array
alignments = mapper.map(cppyy.ll.cast["float*"](arr.ctypes.data), len(arr))

def format_alignment(alignment):
    # return alignment
    return f"Alignment(contig: {alignment.contig}, start: {alignment.ref_start}, end: {alignment.ref_end}, strand: {alignment.is_pos_strand})"

for (i, alignment) in enumerate(alignments):
    print(f"Alignment {i}: {format_alignment(alignment)}")
#%%
1

# --rev-query