#!/usr/bin/env python3
import argparse
import sys
from typing import Literal

import pysam

from pybiotk.io import OpenFqGzip
from pybiotk.utils import reverse_seq, ignore


def reverse_fastx(filename: str, output: str, outfmt: Literal['fastq', 'fasta'] = 'fastq'):
    if outfmt == 'fastq':
        ostream = OpenFqGzip(output)
    else:
        ostream = open(output, "w") if output != "-" else sys.stdout
        
    for entry in pysam.FastxFile(filename):
        name = entry.name
        sequence = reverse_seq(entry.sequence)
        comment = entry.comment
        quality = entry.quality if entry.quality is None else "".join(reversed(str(entry.quality)))
        if outfmt == "fastq":
            ostream.write_entry(name, sequence, comment, quality)
        else:
            ostream.write(f">{name}\n{sequence}\n")
    ostream.close()
    

@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(dest="input", type=str, nargs="?", default=(None if sys.stdin.isatty() else "-"), help="input fastq/a. [stdin]")
    parser.add_argument("-o", "--output", dest="output", type=str, default="-", help="output fastq/a. [stdout]")
    parser.add_argument("--outfmt", dest="outfmt", type=str, choices=("fasta", 'fastq'), default="fastq",  help="output file format fasta or fastq, default:fastq.")
    args = parser.parse_args()
    if args.input is None:
        parser.print_help()
        sys.exit(1)
    reverse_fastx(args.input, args.output, args.outfmt)


if __name__ == "__main__":
    run()
