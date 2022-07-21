# -*- coding: utf-8 -*-
from typing import List, Dict, Iterable, Iterator

from pybiotk.bx.bitset import BinnedBitSet, MAX


def binned_bitsets_from_list(lst: Iterable[List]) -> Dict:
    """Read a list into a dictionary of bitsets"""
    last_chrom = None
    last_bitset = None
    bitsets = dict()
    for bed in lst:
        chrom = bed[0]
        if chrom != last_chrom:
            if chrom not in bitsets:
                bitsets[chrom] = BinnedBitSet(MAX)
            last_chrom = chrom
            last_bitset = bitsets[chrom]
        start, end = int(bed[1]), int(bed[2])
        last_bitset.set_range(start, end - start)
    return bitsets


def unionBed3(lst: Iterable[List]) -> Iterator[List]:
    """Take the union of 3 column bed files. return a new list"""
    bitsets = binned_bitsets_from_list(lst)
    for chrom in bitsets:
        bits = bitsets[chrom]
        end = 0
        while 1:
            start = bits.next_set(end)
            if start == bits.size:
                break
            end = bits.next_clear(start)
            yield [chrom, start, end]
    del bitsets


def intersectBed3(lst1, lst2) -> Iterator[List]:
    """Take the intersection of two bed files (3 column bed files)"""
    bits1 = binned_bitsets_from_list(lst1)
    bits2 = binned_bitsets_from_list(lst2)

    bitsets = dict()
    for key in bits1:
        if key in bits2:
            bits1[key].iand(bits2[key])
            bitsets[key] = bits1[key]

    for chrom in bitsets:
        bits = bitsets[chrom]
        end = 0
        while 1:
            start = bits.next_set(end)
            if start == bits.size:
                break
            end = bits.next_clear(start)
            yield [chrom, start, end]
    bits1.clear()
    bits2.clear()
    bitsets.clear()


def subtractBed3(lst1, lst2) -> Iterator[List]:
    """subtract lst2 from lst1"""
    bitsets1 = binned_bitsets_from_list(lst1)
    bitsets2 = binned_bitsets_from_list(lst2)

    for chrom in bitsets1:
        if chrom not in bitsets1:
            continue
        bits1 = bitsets1[chrom]
        if chrom in bitsets2:
            bits2 = bitsets2[chrom]
            bits2.invert()
            bits1.iand(bits2)
        end = 0
        while 1:
            start = bits1.next_set(end)
            if start == bits1.size:
                break
            end = bits1.next_clear(start)

            yield [chrom, start, end]

    bitsets1.clear()
    bitsets2.clear()
