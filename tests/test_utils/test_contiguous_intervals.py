import pytest
import pysam

from pybiotk.utils import get_blocks_from_read, merge_contiguous_intervals


@pytest.mark.parametrize("intervals, expected", [
    ([], []),
    ([(1, 2)], [(1, 2)]),
    ([(1, 2), (2, 3), (4, 5)], [(1, 3), (4, 5)]),
    ([(1, 2), (2, 3), (3, 4)], [(1, 4)]),
    ([(1, 3), (2, 4)], [(1, 4)]),
])
def test_merge_contiguous_intervals(intervals, expected):
    assert merge_contiguous_intervals(iter(intervals)) == expected


@pytest.mark.parametrize("intervals", [[(3, 2)], [(2, 3), (1, 2)]])
def test_invalid_intervals(intervals):
    with pytest.raises(ValueError):
        merge_contiguous_intervals(intervals)


def test_read_blocks_include_deletions_merge_insertions_and_preserve_introns():
    read = pysam.AlignedSegment()
    read.reference_start = 100
    read.cigarstring = "10M2I5M3D4M20N6M"
    assert get_blocks_from_read(read) == [(100, 122), (142, 148)]
