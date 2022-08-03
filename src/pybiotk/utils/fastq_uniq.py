#!/usr/bin/env python3
import argparse
import time
from typing import Sequence, Literal

from pybiotk.io import FastqFile, FastqPair, OpenFqGzip
from pybiotk.utils import logging


def single_end(input_fq: str, output: str, min_len: int = 15, by: Literal["seq", "id", "name"] = "seq"):
    logging.info(f"Single end mode, by {by} ...")
    too_short_reads = 0
    output_reads = 0
    with FastqFile(input_fq) as fqi, OpenFqGzip(output) as fqo:
        for fq in fqi.uniq(by):
            if len(fq.sequence) < min_len:
                too_short_reads += 1
                continue
            fqo.write_fastx_record(fq)
            output_reads += 1
        input_reads = fqi.ptr
        fqi.ptr = None
    duplicate_ratio = (input_reads - output_reads) * 100 / input_reads
    logging.info((f"result summary:\nread pairs too short (<{min_len}nt): {too_short_reads}\ninput read pairs: {input_reads}\n"
                  f"output read pairs: {output_reads}\nduplicate ratio: {duplicate_ratio:.2f}%"))


def pair_end(input_r1: str, input_r2: str, output_r1: str, output_r2: str, min_len: int = 15, by: Literal["seq", "id", "name"] = "seq"):
    logging.info(f"Pair end mode, by {by} ...")
    too_short_reads = 0
    output_reads = 0

    with FastqPair(input_r1, input_r2) as pairfq, OpenFqGzip(output_r1) as r1, OpenFqGzip(output_r2) as r2:
        for fq1, fq2 in pairfq.uniq(by):
            if min(len(fq1.sequence), len(fq2.sequence)) < min_len:
                too_short_reads += 1
                continue
            r1.write_fastx_record(fq1)
            r2.write_fastx_record(fq2)
            output_reads += 1
        input_reads = pairfq.ptr
        pairfq.ptr = None
    duplicate_ratio = (input_reads - output_reads) * 100 / input_reads
    logging.info((f"result summary:\nread pairs too short (<{min_len}nt): {too_short_reads}\ninput read pairs: {input_reads}\n"
                  f"output read pairs: {output_reads}\nduplicate ratio: {duplicate_ratio:.2f}%"))


def main(input_files: Sequence[str], output_files: Sequence[str], min_len: int = 15, by: Literal["seq", "id", "name"] = "seq"):
    input_str = " ".join(input_files)
    output_str = " ".join(output_files)
    assert len(input_files) == len(output_files)
    logging.info(f"Processing {input_str} ...")

    start = time.perf_counter()
    if len(input_files) == 1:
        single_end(input_files[0], output_files[0], min_len, by)
    else:
        pair_end(*input_files[:2], *output_files[:2], min_len=min_len, by=by)
    end = time.perf_counter()
    logging.info(f"task finished in {end-start:.2f}s, output saved in {output_str}.")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", dest="input", type=str, nargs="+", required=True, help="input *.fastq.gz files.")
    parser.add_argument("-o", dest="output", type=str, nargs="+", required=True, help="output *.fastq.gz files.")
    parser.add_argument("-m", dest="min_len", type=int, default=15, help="min length.")
    parser.add_argument("-by", dest="by", type=str, default="seq", choices=["seq", "id", "name"], help="by seq, id or full name.")
    args = parser.parse_args()
    main(args.input, args.output, args.min_len, args.by)


if __name__ == "__main__":
    run()
