"""
Import this file if you want to use rawhash2 in Python.
"""

from pathlib import Path
from typing import List, TypeVar
import cppyy
import cppyy.ll
import contextlib
import pyslow5
import numpy as np

def get_read_signal(slow5_filename, read_id: str):
    with contextlib.closing(pyslow5.Open(slow5_filename, "r")) as s5:
        return s5.get_read(read_id, pA=True)["signal"]

def prepare_signal_for_rawhash(raw_signal: np.ndarray):
    """
    Remove outliers, convert to numpy
    raw_signal should be floats (with offset, range, digitisation scaling)
    
    raw_signal: np.ndarray[float]
    """
    raw_signal = raw_signal[(raw_signal > 30.) & (raw_signal < 200.)]
    return raw_signal

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
    
def load_rawhash_wrapper_lib(rawhash_lib_dir: Path):
    rawhash_lib_dir = Path(rawhash_lib_dir)
    header_file = rawhash_lib_dir / "src/rawhash_wrapper.hpp"
    # library_file = rawhash_lib_dir / "build/src/librawhash2_wrapper.so"
    cppyy.include(str(header_file))
    cppyy.add_library_path("/home/mmordig/rawhash_project/rawhash2_new/build/src") # for shared libs
    cppyy.load_library("librawhash2_wrapper")

# define type RawHashMapper

RawHashMapper = TypeVar("RawHashMapper")
def get_rawhash_mapper(rawhash_args: List[str]) -> RawHashMapper:
    args = ["my_dummy_program"] + rawhash_args
    mapper = cppyy.gbl.RawHashMapper(len(args), args)
    return mapper

def get_rawhash_alignments(mapper: RawHashMapper, raw_signal: np.ndarray) -> List[Alignment]:
    """
    raw_signal: np.ndarray[float]
    """
    raw_signal = np.ascontiguousarray(raw_signal, dtype=np.float32)
    alignments = mapper.map(cppyy.ll.cast["float*"](raw_signal.ctypes.data), len(raw_signal))
    # return alignments
    return [Alignment.from_cppyy(alignment) for alignment in alignments]

def parse_paf_line(line):
    """returns read_id, alignment"""
    fields = line.strip().split("\t")
    return fields[0], Alignment(
        contig=fields[5],
        ref_start=int(fields[7]),
        ref_end=int(fields[8]),
        is_pos_strand=fields[4] == "+",
    )