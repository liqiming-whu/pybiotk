#!/usr/bin/env python3
import argparse
import sys
import time
from typing import Literal

from pybiotk.io import FastxFile, OpenFqGzip
from pybiotk.utils import logging, ignore


def fastx_rename(input_fq: str, output: str, outfmt: Literal["fastq", "fasta"] = "fastq") -> None:
    input_str = "stdin" if input_fq == "-" else input_fq
    logging.info(f"Processing {input_str} ...")
    start = time.perf_counter()
    if outfmt == "fastq":
        ostream = OpenFqGzip(output)
    else:
        ostream = open(output, "w") if output != "-" else sys.stdout
    with FastxFile(input_fq) as fqi, ostream as fqo:
        for fq in fqi.rename():
            if outfmt == "fastq":
                fqo.write_fastx_record(fq)
            else:
                fqo.write(f">{fq.name}\n{fq.sequence}\n")
        input_reads = fqi.ptr
        fqi.ptr = None
    end = time.perf_counter()
    logging.info(f"Processed {input_reads} reads in {end - start:.2f} seconds.")


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(dest="input", type=str, nargs="?", default=(None if sys.stdin.isatty() else "-"), help="input fastq/a. [stdin]")
    parser.add_argument("-o", "--output", dest="output", type=str, default="-", help="output fastq/a. [stdout]")
    parser.add_argument("--outfmt", dest="outfmt", type=str, choices=("fastq", "fasta"), default="fastq", help="output fastq or fasta. [fastq]")
    args = parser.parse_args()
    if args.input is None:
        parser.print_help()
        sys.exit(1)
    fastx_rename(args.input, args.output, args.outfmt)


if __name__ == "__main__":
    run()


