#!/usr/bin/env python3
"""Convert one or more BED files to a bedGraph coverage stream.

The built-in implementation uses endpoint events and accepts unsorted BED
input. Its memory use is proportional to the number of distinct interval
endpoints rather than the number of covered bases.

For very large, chromosome-grouped BED files, bedtools provides a compiled
alternative that also requires a chromosome sizes file::

    bedtools genomecov -i input.bed -bg -g chrom.sizes
"""
import argparse
import gzip
import os
import sys
from collections import defaultdict
from contextlib import nullcontext
from typing import Dict, Iterable, Iterator, TextIO, Tuple

from pybiotk.utils import ignore


Events = Dict[str, Dict[int, int]]


def _open_bed(filename: str):
    if filename == "-":
        return nullcontext(sys.stdin)
    if filename.endswith(".gz"):
        return gzip.open(filename, "rt")
    return open(filename, encoding="utf-8")


def _read_intervals(bedfile: str, stream: TextIO) -> Iterator[Tuple[str, int, int]]:
    for line_number, line in enumerate(stream, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "track", "browser")):
            continue
        fields = stripped.split("\t", 3)
        if len(fields) < 3:
            raise ValueError(f"{bedfile}:{line_number}: expected at least three BED columns")
        chrom = fields[0]
        try:
            start = int(fields[1])
            end = int(fields[2])
        except ValueError as exc:
            raise ValueError(f"{bedfile}:{line_number}: BED coordinates must be integers") from exc
        if start < 0 or end <= start:
            raise ValueError(
                f"{bedfile}:{line_number}: expected 0 <= start < end, got {start}, {end}"
            )
        yield chrom, start, end


def _collect_events(bed_list: Iterable[str]) -> Events:
    events: Events = defaultdict(lambda: defaultdict(int))
    for bedfile in bed_list:
        with _open_bed(bedfile) as stream:
            for chrom, start, end in _read_intervals(bedfile, stream):
                events[chrom][start] += 1
                events[chrom][end] -= 1
    return events


def _coverage_segments(chrom_events: Dict[int, int]) -> Iterator[Tuple[int, int, int]]:
    coverage = 0
    previous = None
    for position in sorted(position for position, delta in chrom_events.items() if delta):
        if previous is not None and previous < position and coverage > 0:
            yield previous, position, coverage
        coverage += chrom_events[position]
        previous = position
    if coverage != 0:
        raise RuntimeError("unbalanced BED interval endpoints")


def main(bed_list, header):
    if header:
        basename = os.path.basename(bed_list[0])
        name = os.path.splitext(basename)[0]
        if name.endswith(".bed"):
            name = os.path.splitext(name)[0]
        print(
            f'track type=bedGraph name="{name}" visibility=full color=0,51,51 '
            "AutoScale=on alwaysZero=on maxHeightPixels=50:50:50"
        )

    events = _collect_events(bed_list)
    output = sys.stdout
    for chrom in sorted(events):
        for start, end, coverage in _coverage_segments(events[chrom]):
            output.write(f"{chrom}\t{start}\t{end}\t{coverage}\n")


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", type=str, nargs="+", help="input BED files")
    parser.add_argument("--header", dest="header", action="store_true",
                        help="add a UCSC bedGraph track header")
    args = parser.parse_args()
    try:
        main(args.input, args.header)
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    run()
