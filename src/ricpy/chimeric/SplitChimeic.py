#!/usr/bin/env python3
import argparse
import time
import pysam
from collections import namedtuple
from typing import Literal, Optional
from stream import Pipe, mapwith, apply, filter
from pybiotk.utils import reverse_seq, infer_fragment_strand, logging
from pybiotk.io import Bam, FastaFile, OpenFqGzip



class ChimericSegment:
    def __init__(self, read: pysam.AlignedSegment,
                 rule: Literal["1+-,1-+,2++,2--", "1++,1--,2+-,2-+", "+-,-+", "++,--"] = "+-,-+"):
        self.read = read
        self.name = read.query_name
        self.sequence = read.get_forward_sequence()
        self.quality = pysam.qualities_to_qualitystring(
            read.get_forward_qualities())
        self.strand = '-' if read.is_reverse else '+'
        self.length = len(self.sequence)
        self.cigar = read.cigartuples
        self.ref_name = read.reference_name
        self.aligned_sequence = read.query_alignment_sequence
        self.aligned_length = read.reference_length
        self.fragment_strand = infer_fragment_strand(self.strand, rule, read.is_read2)
        if not self.fragment_strand == self.strand:
            self.sequence = reverse_seq(self.sequence)
            self.quality = "".join(reversed(self.quality))
            
    def __repr__(self):
        return str(self.read)
    
    __str__ = __repr__        


@Pipe
def to_ChimericSegment(read: pysam.AlignedSegment, rule: Literal["1+-,1-+,2++,2--", "1++,1--,2+-,2-+", "+-,-+", "++,--"] = "+-,-+") -> ChimericSegment:
    return ChimericSegment(read, rule)
        
        
def filter_ChimericSegment(chims: ChimericSegment):
    if len(chims.cigar) > 3:
        return False
    if len(chims.cigar) == 1:
        return False
    if chims.fragment_strand == '-':
        return False
    return True


def split_ChimericSegment(chims: ChimericSegment,
                          ref_fasta: Optional[FastaFile] = None,
                          more_C: Optional[dict] = None):
    assert len(chims.cigar) == 2
    if chims.cigar[0][0] == 0:
        ref_loc = 'left'
    else:
        ref_loc = 'right'
    cigar_dict = dict(chims.cigar)
    ref_len = cigar_dict[0]
    other_len = cigar_dict[4]
    ref_name = chims.ref_name
    mapped_ref_seq = chims.aligned_sequence
    assert ref_len == chims.aligned_length
    assert ref_len + other_len == chims.length
    if ref_loc == 'left':
        pCp = 'withoutpCp'
        other_seq = chims.sequence[-other_len:]
        other_qua = chims.quality[-other_len:]
        ref_seq = chims.sequence[:ref_len]
        if other_seq[0] == 'C':
            pCp = 'withpCp'
            other_seq = other_seq[1:]
            other_qua = other_qua[1:]

    if ref_loc == 'right':
        pCp = 'withoutpCp'
        other_seq = chims.sequence[:other_len]
        other_qua = chims.quality[:other_len]
        ref_seq = chims.sequence[-ref_len:]
        if other_seq[-1] == 'C':
            pCp = 'withpCp'
            other_seq = other_seq[:-1]
            other_qua = other_qua[:-1]

    if not mapped_ref_seq == ref_seq:
        return None
    if len(other_seq) < 25:
        return None
    
    if ref_fasta is not None:
        true_ref_seq = ref_fasta.fetch(chims.ref_name)
        if not ref_seq == true_ref_seq:
            return None
        if more_C is not None:          
            if pCp == "withoutpCp":
                if ref_loc == 'left':
                    if ref_name in more_C.RC:
                        subseq = more_C.RC[ref_name]
                        ref_name = subseq.name
                        ref_seq = subseq.seq
                elif ref_loc == 'right':
                    if ref_name in more_C.LC:
                        subseq = more_C.LC[ref_name]
                        ref_name = subseq.name
                        ref_seq = subseq.seq
                        pCp = "withpCp"
    seq_name = f"{chims.name}|{ref_name}|{ref_loc}|{pCp}"
    split_info = namedtuple("split_info", ["seq_name", "ref_name", "ref_seq", "other_seq", "other_qua"])
    
    return split_info(seq_name, ref_name, ref_seq, other_seq, other_qua)


def split_ChimericSegment3(chims: ChimericSegment,
                           ref_fasta: Optional[FastaFile] = None,
                           more_C: Optional[dict] = None):
    assert len(chims.cigar) == 3
    cigar_type = [i[0] for i in chims.cigar]
    if not cigar_type == [4, 0, 4]:
        return None
    ref_loc = "mid"
    left_len = chims.cigar[0][1]
    right_len = chims.cigar[2][1]
    mapped_ref_seq = chims.aligned_sequence
    ref_name = chims.ref_name
    left = chims.sequence[:left_len]
    left_qua = chims.quality[:left_len]
    right = chims.sequence[-right_len:]
    right_qua = chims.quality[-right_len:]
    ref_seq = chims.sequence[left_len:-right_len]

    if not mapped_ref_seq == ref_seq:
        return None

    LpCp = "withoutLpCp"
    RpCp = "withoutRpCp"

    if left[-1] == 'C':
        LpCp = "withLpCp"
        left = left[:-1]
        left_qua = left_qua[:-1]

    if right[0] == 'C':
        RpCp = "withRpCp"
        right = right[1:]
        right_qua = right_qua[1:]

    if ref_fasta is not None:
        true_ref_seq = ref_fasta.fetch(chims.ref_name)
        if not ref_seq == true_ref_seq:
            return None
        if more_C is not None:
            if LpCp == "withoutLpCp" and RpCp == "withRpCp":
                if ref_name in more_C.LC:
                    subseq = more_C.LC[ref_name]
                    ref_name = subseq.name
                    ref_seq = subseq.seq
                    LpCp = "withLpCp"
            elif LpCp == "withLpCp" and RpCp == "withoutRpCp":
                if ref_name in more_C.RC:
                    subseq = more_C.RC[ref_name]
                    ref_name = subseq.name
                    ref_seq = subseq.seq
                    RpCp = "withRpCp"
            elif LpCp == "withoutLpCp" and RpCp == "withoutRpCp":
                if ref_name in more_C.LR:
                    subseq = more_C.LR[ref_name]
                    ref_name = subseq.name
                    ref_seq = subseq.seq
                    LpCp = "withLpCp"
                    RpCp = "withRpCp"
                elif ref_name in more_C.LC and ref_name not in more_C.RC:
                    subseq = more_C.LC[ref_name]
                    ref_name = subseq.name
                    ref_seq = subseq.seq
                    LpCp = "withLpCp"
                elif ref_name in more_C.RC and ref_name not in more_C.LC:
                    subseq = more_C.RC[ref_name]
                    ref_name = subseq.name
                    ref_seq = subseq.seq
                    RpCp = "withRpCp"
    if not (left or right):
        return None
    if len(left) < 5 or len(right) < 5:
        return None    
    if len(left) + len(right) < 25:
        return None
    seq_name = f"{chims.name}|{ref_name}|{ref_loc}|{LpCp}|{RpCp}"
    split_info3 = namedtuple("split_info3", ["seq_name", "ref_name", "ref_seq", "left_seq", "left_qua", "right_seq", "right_qua"])
    return split_info3(seq_name, ref_name, ref_seq, left, left_qua, right, right_qua)


def split_chimeic(bamfile:str,
                  out_fa: str,
                  other_fq: str,
                  unmap_fq: str,
                  mid_left_fq: str,
                  mid_right_fq: str,
                  mapped_info: str,
                  ref_fasta:Optional[str] = None,
                  more_C_file: Optional[str] = None,
                  rule: str = "+-,-+"):
    logging.info("start to split chimeric reads...")
    start = time.perf_counter()
    ref_fa = FastaFile(ref_fasta) if ref_fasta is not None else None
    more_C = more_C_file if more_C_file is not None else None
    mapped_reads = {}
    mapped_reads3 = {}
    out_ref = open(out_fa, "w")
    unmap = OpenFqGzip(unmap_fq)
    other = OpenFqGzip(other_fq)
    mid_left = OpenFqGzip(mid_left_fq)
    mid_right = OpenFqGzip(mid_right_fq)
    info = open(mapped_info, "w")
    total_reads_num = 0
    with Bam(bamfile) as bam:
        for chims in bam.iter_mapped(secondary=True) | to_ChimericSegment(rule) | filter(filter_ChimericSegment):
            if chims.read.is_secondary:
                total_reads_num += 1
            if len(chims.cigar) == 2:
                split_info = split_ChimericSegment(chims, ref_fa, more_C)
                if split_info is None:
                    if not chims.read.is_secondary:
                        unmap.write_entry(name=chims.name, sequence=chims.sequence, quality=chims.quality)
                        logging.info(f"{chims.name} is unmapped.")
                    continue
                if split_info.seq_name not in mapped_reads:
                    mapped_reads[split_info.seq_name] = split_info
                else:
                    if len(split_info.ref_seq) > len(mapped_reads[split_info.seq_name].ref_seq):
                        mapped_reads[split_info.seq_name] = split_info
            elif len(chims.cigar) == 3:
                split_info3 = split_ChimericSegment3(chims, ref_fa, more_C)
                if split_info3 is None:
                    if not chims.read.is_secondary:
                        unmap.write_entry(name=chims.name, sequence=chims.sequence, quality=chims.quality)
                        logging.info(f"{chims.name} is unmapped")
                    continue
                if split_info3.seq_name not in mapped_reads3:
                    mapped_reads3[split_info3.seq_name] = split_info3
                else:
                    if len(split_info3.ref_seq) > len(mapped_reads3[split_info3.seq_name].ref_seq):
                        mapped_reads3[split_info3.seq_name] = split_info3
    mapped_reads_num = len(mapped_reads) + len(mapped_reads3)
    for read in mapped_reads:
        split_info = mapped_reads[read]
        out_ref.write(f">{split_info.seq_name}\n{split_info.ref_seq}\n")
        other.write_entry(name=split_info.seq_name, sequence=split_info.other_seq, quality=split_info.other_qua)
        info.write(f"{split_info.seq_name}\tcigar_len=2\t{split_info.ref_name}\t{len(split_info.ref_seq)}\n")
    for read in mapped_reads3:
        split_info3 = mapped_reads3[read]
        out_ref.write(f">{split_info3.name}\n{split_info3.ref_seq}\n")
        mid_left.write_entry(name=split_info3.seq_name, sequence=split_info3.left_seq, quality=split_info3.left_qua)
        mid_right.write_entry(name=split_info3.seq_name, sequence=split_info3.right_seq, quality=split_info3.right_qua)
        info.write(f"{split_info3.seq_name}\tcigar_len=3\t{split_info3.ref_name}\t{len(split_info3.ref_seq)}\n")
    
    out_ref.close()
    unmap.close()
    other.close()
    mid_left.close()
    mid_right.close()
    info.close()
    end = time.perf_counter()
    logging.info(f"{mapped_reads_num} mapped reads in total {total_reads_num} reads, {mapped_reads_num*100/total_reads_num:.2f}%")
    logging.info(f"task finished in {end-start:.2f}s.")
