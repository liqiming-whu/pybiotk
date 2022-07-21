#!/usr/bin/env python3
import argparse
import time
from typing import Sequence

from pybiotk.io import FastaFile
from pybiotk.utils import logging, ignore


def main(filename, chromList: Sequence[str], wrap: bool = True):
    logging.info("reading fasta ....")
    start = time.perf_counter()
    with FastaFile(filename) as fa:
        fa.stdout(referenceList=chromList, wrap=wrap)
    end = time.perf_counter()
    logging.info(f"task finished in {end-start:.2f}s.")


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(dest="fasta", type=str, help="Genome fasta file.")
    parser.add_argument("-c", dest="chromlist", nargs="+", default=None, help="chrom list for split genome.")
    parser.add_argument("--wrap", dest="wrap", action="store_true", help="line-wrapped display.")
    args = parser.parse_args()
    main(args.fasta, args.chromlist, args.wrap)


if __name__ == "__main__":
    run()
