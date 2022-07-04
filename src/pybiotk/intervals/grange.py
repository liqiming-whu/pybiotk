# -*- coding: utf-8 -*-
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, List, Dict, Literal, Tuple, Union, Optional
from pybiotk.bx.intersection import IntervalTree
from pybiotk.utils import is_overlap
from pybiotk.intervals.merge_intervals import MergedIntervals


@dataclass
class GRange:
    chrom: str
    start: int
    end: int
    strand: Optional[Literal['+', '-']] = field(default=None)
    value: List[Any] = field(default_factory=list, init=False, repr=True, compare=False)

    def __len__(self):
        return self.end - self.start

    def __str__(self) -> str:
        return f"{self.chrom}:{self.start}-{self.end}({self.strand})"

    def is_overlap(self, other: GRange) -> bool:
        if self.chrom != other.chrom:
            return False
        if self.strand != other.strand:
            return False
        if not is_overlap((self.start, self.end), (other.start, other.end)):
            return False
        return True


class GRangeTree:
    def __init__(self, strand: bool = False):
        self.strand: bool = strand
        self._trees: Dict[Union[str, Tuple[str, ...]], IntervalTree] = {}

    def __iter__(self):
        for _, tree in sorted(self._trees.items()):
            nodes = list()
            fn = nodes.append
            tree.traverse(fn)
            for node in nodes:
                yield node.interval

    def _add_record(self, key: Union[str, Tuple[str, ...]], start: int, end: int, record: Any):
        if key not in self._trees.keys():
            self._trees[key] = IntervalTree()
        self._trees[key].insert(int(start), int(end), record)

    def add(self, obj: Any, chrom: str, start: int, end: int, strand: Optional[Literal['+', '-']] = None):
        if self.strand:
            key = (chrom, strand)
        else:
            key = chrom
        self._add_record(key, start, end, obj)

    def add_range(self, record: GRange):
        if self.strand:
            key = (record.chrom, record.strand)
        else:
            key = record.chrom
        start = record.start
        end = record.end
        self._add_record(key, start, end, record)

    def find(self, chrom: str, start: int, end: int, strand: Optional[Literal['+', '-']] = None):
        if self.strand:
            key = (chrom, strand)
        else:
            key = chrom
        overlaped = []
        if key in self._trees.keys():
            overlaped = self._trees[key].find(int(start), int(end))
        return overlaped

    def iter_chrom(self, chrom: str, strand: Optional[Literal['+', '-']] = None):
        if self.strand:
            key = (chrom, strand)
        else:
            key = chrom
        nodes = []
        tree = self._trees[key]
        fn = nodes.append
        tree.traverse(fn)
        for node in nodes:
            yield node.interval


class MergedGRange:
    def __init__(self, strand: bool = False):
        self.strand: bool = strand
        self._dict: Dict[Union[str, Tuple[str, ...]], MergedIntervals] = {}

    def __iter__(self):
        for key in self._dict:
            merge_intervals = self._dict[key]
            if self.strand:
                chrom, strand = key
            else:
                chrom = key
                strand = None
            for start, end in merge_intervals:
                yield GRange(chrom, start, end, strand)

    def _add_record(self, key: Union[str, Tuple[str, ...]], start: int, end: int):
        if key not in self._trees.keys():
            self._trees[key] = MergedIntervals()
        self._trees[key].add_interval(int(start), int(end))

    def add(self, obj: Any, chrom: str, start: int, end: int, strand: Optional[Literal['+', '-']] = None):
        if self.strand:
            key = (chrom, strand)
        else:
            key = chrom
        self._add_record(key, start, end, obj)

    def add_range(self, record: GRange):
        if self.strand:
            key = (record.chrom, record.strand)
        else:
            key = record.chrom
        start = record.start
        end = record.end
        self._add_record(key, start, end, record)

    def intersect(self, other: MergedGRange):
        assert self.strand == other.strand
        for key in self._dict:
            if key in other._dict:
                self._dict[key].intersect(other._dict[key])

    def subtract(self, other: MergedGRange):
        assert self.strand == other.strand
        for key in self._dict:
            if key in other._dict:
                self._dict[key].subtract(other._dict[key])
