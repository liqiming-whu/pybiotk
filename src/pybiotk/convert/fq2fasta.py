#!/usr/bin/env python3
import argparse
import sys

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
    parser.add_argument(dest="input", type=str, nargs="*", default=(None if sys.stdin.isatty() else "-"), help="input fastq files. [stdin]")
    args = parser.parse_args()
    if not args.input:
        parser.print_help()
        sys.exit(1)
    main(args.input)


if __name__ == "__main__":
    run()
