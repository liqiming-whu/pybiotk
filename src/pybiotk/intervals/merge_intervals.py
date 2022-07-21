# -*- coding: utf-8 -*-
from __future__ import annotations

from typing import List, Iterable, Tuple, Iterator

from pybiotk.bx.bitset import BinnedBitSet, MAX


class MergedIntervals:
    def __init__(self, intervals: Iterable[Tuple[int, int]] = []):
        self.bitset = None
        self.intervals = intervals
        self.filled = False

    def _fill(self):
        self.bitset = BinnedBitSet(MAX)
        for interval in self.intervals:
            start, end = int(interval[0]), int(interval[1])
            self.bitset.set_range(start, end - start)
        self.filled = True

    def __iter__(self) -> Iterator[Tuple[int, int]]:
        if not self.filled:
            self._fill()
        end = 0
        while True:
            start = self.bitset.next_set(end)
            if start == self.bitset.size:
                break
            end = self.bitset.next_clear(start)
            yield start, end

    def clear(self):
        del self.bitset
        self.bitset = None

    def extend_intervals(self, intervals: Iterable[Tuple[int, int]]):
        self.intervals.extend(intervals)

    def add_interval(self, interval: Tuple[int, int]):
        self.intervals.append(interval)

    def intersect(self, other: MergedIntervals):
        if not self.filled:
            self._fill()
        if not other.filled:
            other._fill()
        self.bitset.iand(other.bitset)

    def subtract(self, other: MergedIntervals):
        if not self.filled:
            self._fill()
        if not other.filled:
            other._fill()
        other.bitset.invert()
        self.bitset.iand(other.bitset)


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    intervals_iter = MergedIntervals(intervals)
    interval_list = list(intervals_iter)
    intervals_iter.clear()
    return interval_list


def intersect_intervals(intervals1: Iterable[Tuple[int, int]], intervals2: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged_intervals1 = MergedIntervals(intervals1)
    merged_intervals2 = MergedIntervals(intervals2)
    merged_intervals1.intersect(merged_intervals2)
    merged_intervals2.clear()
    interval_list = list(merged_intervals1)
    merged_intervals1.clear()
    return interval_list


def subtract_intervals(intervals1: Iterable[Tuple[int, int]], intervals2: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    merged_intervals1 = MergedIntervals(intervals1)
    merged_intervals2 = MergedIntervals(intervals2)
    merged_intervals1.subtract(merged_intervals2)
    interval_list = list(merged_intervals1)
    merged_intervals1.clear()
    return interval_list
