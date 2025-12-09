#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Path: src/pybiotk/utils/bigwigfetcher.py
"""
Fetch signal value from bigwig file.
"""
import argparse
import re
import sys
import time
from typing import Sequence, Optional, TextIO, Literal

from pybiotk.io import Openbwn, Openbed
from pybiotk.utils import logging, ignore


class LocationFormatError(Exception):
    pass


def main(
    bwfiles: Sequence[str],
    location: Optional[str] = None,
    bedfile: Optional[str] = None,
    type: Literal["stats", "values"] = "stats",
    stats: Literal["mean", "max", "min", "sum", "coverage", "std"] = "mean",
    output: Optional[TextIO] = None
):
    """Fetch signal value from bigwig file.

    Args:
        bwfile (str): bigwig files.
        location (str, optional): location file. Defaults to None.
        bedfile (str, optional): bed file. Defaults to None.
        output (TextIO, optional): output file. Defaults to None.

    Raises:
        LocationFormatError: location format error.
    """
    time_start = time.perf_counter()
    with Openbwn(bwfiles) as bw, output or sys.stdout as out:
        logging.info("read bigwig file ...")
        if bedfile is not None:
            with Openbed(bedfile) as bed:
                for i in bed:
                    if type == "stats":
                        res = bw.stats(i.chrom, i.start, i.end, stat=stats)
                    elif type == "values":
                        res = bw.values(i.chrom, i.start, i.end)
                    res_str = "\t".join([str(x) for x in res])
                    out.write(f"{i.chrom}\t{i.start}\t{i.end}\t{res_str}\n")
        elif location is not None:
            location_pattern = re.compile(r"(.+):(\d+)-(\d+)\(([+-])\)")
            loc = location_pattern.match(location)
            if loc is None:
                raise LocationFormatError(f"location format error: {location}")
            else:
                chrom, start, end, strand = loc.groups()
                if type == "stats":
                    res = bw.stats(chrom, start, end, stat=stats)
                elif type == "values":
                    res = bw.values(chrom, start, end)
                res_str = "\t".join([str(x) for x in res])
                out.write(f"{chrom}\t{start}\t{end}\t{strand}\t{res_str}\n")
    time_end = time.perf_counter()
    logging.info(f"task finished in {time_end-time_start:.2f}s.")


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--bwfiles', dest="bwfiles", type=str, nargs="+", help="bigwig files.")
    parser.add_argument('-l', '--location', dest="location", type=str, default=None, help='"chromosome:start-end(+/-)". will be ignored if -b is set.')
    parser.add_argument('-b', '--bedfile', dest="bedfile", type=str, default=None, help="bed file.")
    parser.add_argument('-t', '--type', dest="type", type=str, default="stats", choices=("stats", "values"), help="output type. [stats|values].")
    parser.add_argument('-s', '--stats', dest="stats", type=str, default="mean", choices=("mean", "max", "min", "sum", "coverage", "std"), help="stats type. [mean|max|min|sum|coverage|std].")
    parser.add_argument('-o', '--output', dest='output', type=argparse.FileType('w'), default=sys.stdout, help="output file name.")

    args = parser.parse_args()
    if args.bedfile is None and args.location is None:
        parser.error("please set -l or -b option.")
    if args.bedfile is not None and args.location is not None:
        logging.warning("both -l and -b are set, -l will be ignored.")
    main(**vars(args))


if __name__ == "__main__":
    run()
