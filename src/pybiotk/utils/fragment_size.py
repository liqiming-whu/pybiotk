#!/usr/bin/env python3
import argparse
import os
import sys
import time
from collections import Counter

from pybiotk.intervals import merge_intervals
from pybiotk.io import Bam, BamPE, BamType, check_bam_type
from pybiotk.utils import logging, blocks_len, intervals_is_overlap, ignore


def main(filename: str, output: str, ordered_by_name: bool = False):
    start = time.perf_counter()
    bamtype = check_bam_type(filename)
    length = []
    sys.stdout.write("read_name\treference_name\tblocks\tlength\n")
    if bamtype is BamType.SE:
        logging.info("SingleEnd mode ...")
        with Bam(filename) as bam:
            for read in bam.iter_mapped():
                blocks = read.get_blocks()
                read_len = blocks_len(blocks)
                length.append(read_len)
                sys.stdout.write(f"{read.query_name}\t{read.reference_name}\t{blocks}\t{read_len}\n")
    elif bamtype is BamType.PE:
        logging.info("PairEnd mode ...")
        with BamPE(filename) as bampe:
            bampe.ordered_by_name = ordered_by_name
            for read1, read2 in bampe.iter_pair(properly_paired=False):
                if read1 is not None and read2 is not None:
                    if read1.reference_name == read2.reference_name:
                        read1_blocks = read1.get_blocks()
                        read2_blocks = read2.get_blocks()
                        if not intervals_is_overlap(read1_blocks, read2_blocks):
                            logging.warning(f"{read1.query_name} read1 and read2 not overlap, use read1_len + read2_len")
                        merge_blocks = merge_intervals(read1_blocks+read2_blocks)
                        fragment_len = blocks_len(merge_blocks)
                        length.append(fragment_len)
                        sys.stdout.write(f"{read1.query_name}\t{read1.reference_name}\t{merge_blocks}\t{fragment_len}")
                    elif not read1.reference_name == read2.reference_name:
                        read1_blocks = read1.get_blocks()
                        read1_len = blocks_len(read1_blocks)
                        length.append(read1_len)
                        sys.stdout.write(f"{read1.query_name}_R1\t{read1.reference_name}\t{read1_blocks}\t{read1_len}\n")
                        read2_blocks = read2.get_blocks()
                        read2_len = blocks_len(read2_blocks)
                        length.append(read2_len)
                        sys.stdout.write(f"{read2.query_name}_R2\t{read2.reference_name}\t{read2_blocks}\t{read2_len}\n")
                if read1 is not None and read2 is None:
                    read1_blocks = read1.get_blocks()
                    read1_len = blocks_len(read1_blocks)
                    length.append(read1_len)
                    sys.stdout.write(f"{read1.query_name}_R1\t{read1.reference_name}\t{read1_blocks}\t{read1_len}\n")
                if read2 is not None and read1 is None:
                    read2_blocks = read2.get_blocks()
                    read2_len = blocks_len(read2_blocks)
                    length.append(read2_len)
                    sys.stdout.write(f"{read2.query_name}_R2\t{read2.reference_name}\t{read2_blocks}\t{read2_len}\n")
    counter = Counter(length)
    counter = list(counter.items())
    counter.sort(key=lambda x: x[0])
    with open(output, "w") as f:
        f.write("Length\tCount\n")
        for length, count in counter:
            f.write(f"{length}\t{count}\n")

    end = time.perf_counter()
    logging.info(f"task finished in {end-start:.2f}s.")


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str, help="input bam.")
    parser.add_argument("-s", dest="summary", type=str, default=os.devnull,
                        help="output summary file name.")
    parser.add_argument("--ordered_by_name", dest="ordered_by_name", action="store_true",
                        help="if input bam is ordered by name.")

    args = parser.parse_args()
    main(args.input, args.summary, args.ordered_by_name)


if __name__ == "__main__":
    run()
