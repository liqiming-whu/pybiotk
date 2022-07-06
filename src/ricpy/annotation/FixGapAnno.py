#!/usr/bin/env python3
import os
import pysam
import argparse
from itertools import groupby
from typing import Tuple, List, NamedTuple, Optional
from pybiotk.utils import is_overlap
from pybiotk.intervals import merge_intervals


class Mergeblocks(NamedTuple):
    mergeblocks: Optional[List[Tuple[int, int]]] = None
    gap: Optional[Tuple[int, int]] = None
    gap_length: Optional[int] = None
    

def merge_blocks(blocks1: List[Tuple[int, int]], blocks2: List[Tuple[int, int]]) -> Mergeblocks:
    if blocks1[0][0] > blocks2[0][0]:
        blocks1, blocks2 = blocks2, blocks1
    blocks1_start = blocks1[0][0]
    blocks1_end = blocks1[-1][1]
    blocks2_start = blocks2[0][0]
    blocks2_end = blocks2[-1][1]
    mergeblocks = None
    gap = None,
    gap_length = None
    if is_overlap((blocks1_start, blocks1_end), (blocks2_start, blocks2_end)):
        gap = 0
        mergeblocks = merge_intervals(blocks1 + blocks2)
    else:
        mergeblocks = blocks1[:-1] + [(blocks1[-1][0], blocks2[0][1])] + blocks2[1:]
        gap = [(blocks1_end, blocks2_start)]
        gap_length = blocks2_start - blocks1_end
    
    return Mergeblocks(mergeblocks, gap, gap_length)
