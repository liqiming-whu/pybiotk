#!/usr/bin/env python3
import argparse

from pybiotk.io import BamPE, BamType, BamTypeError, check_bam_type
from pybiotk.utils import logging


def main(filename: str, output: str, ordered_by_name: bool = False):
    bamtype = check_bam_type(filename)

    if bamtype is not BamType.PE:
        raise BamTypeError(f"{filename} is not a pair end bam.")

    with BamPE(filename) as bam:
        bam.ordered_by_name = ordered_by_name
        bam.to_bam_ordered_by_name(output)

    logging.info(f"output bam have been saved in {output}.")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str, help="input pair end bam.")
    parser.add_argument("-o", dest="output", type=str, required=True,
                        help="output file.")
    parser.add_argument("--ordered_by_name", dest="ordered_by_name", action="store_true",
                        help="if input bam is ordered by name.")

    args = parser.parse_args()
    main(args.input, args.output, args.ordered_by_name)


if __name__ == "__main__":
    run()
