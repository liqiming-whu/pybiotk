#!/usr/bin/env python3
import argparse
import sys
import time
from typing import Sequence, Literal, Union

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


def main(input_files: Union[Sequence[str], str] = "-", output_files: Union[Sequence[str], str] = "-", min_len: int = 15, by: Literal["seq", "id", "name"] = "seq"):
    if isinstance(input_files, str):
        input_files = [input_files]
    if isinstance(output_files, str):
        output_files = [output_files]

    assert len(input_files) == len(output_files)
    input_str = " ".join(input_files)
    output_str = " ".join(output_files)
    if input_str == "-":
        input_str = "stdin"
    logging.info(f"Processing {input_str} ...")

    start = time.perf_counter()
    if len(input_files) == 1:
        single_end(input_files[0], output_files[0], min_len, by)
    else:
        pair_end(*input_files[:2], *output_files[:2], min_len=min_len, by=by)
    end = time.perf_counter()
    if output_str == "-":
        logging.info(f"task finished in {end-start:.2f}s")
    else:
        logging.info(f"task finished in {end-start:.2f}s, output saved in {output_str}.")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(dest="input", type=str, nargs="*", default=(None if sys.stdin.isatty() else "-"), help="input fastq files. [stdin]")
    parser.add_argument("-o", "--output", dest="output", type=str, nargs="+", default="-", help="output fastq files. [stdout]")
    parser.add_argument("-m", "--min-len", dest="min_len", type=int, default=15, help="min length.")
    parser.add_argument("-i", "--by-id", dest="id", action="store_true", help="by id instead of seq.")
    parser.add_argument("-n", "--by-name", dest="name", action="store_true", help="by full name instead of just id.")
    args = parser.parse_args()
    if not args.input:
        parser.print_help()
        sys.exit(1)
    by = "seq"
    if args.id:
        by = "id"
    if args.name:
        by = "name"
    main(args.input, args.output, args.min_len, by)


if __name__ == "__main__":
    run()
