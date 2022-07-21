#!/usr/bin/env python3
import argparse

from pybiotk.io import FastqFile
from pybiotk.utils import ignore
from stream import stdout


def main(fq_list):
    for fq in fq_list:
        with FastqFile(fq) as f:
            f.to_fasta() | stdout


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str, nargs="+",
                        help="Input fastq files.")
    args = parser.parse_args()
    main(args.input)


if __name__ == "__main__":
    run()
