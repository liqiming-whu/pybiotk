#!/usr/bin/env python3
import argparse
import sys

from pybiotk.io import FastxFile
from pybiotk.utils import ignore

def main(fa_list):
    for fa in fa_list:
        with FastxFile(fa) as f:
            f.to_fastq()


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(dest="input", type=str, nargs="*", default=(None if sys.stdin.isatty() else "-"), help="input fasta files. [stdin]")
    args = parser.parse_args()
    if not args.input:
        parser.print_help()
        sys.exit(1)
    main(args.input)


if __name__ == "__main__":
    run()