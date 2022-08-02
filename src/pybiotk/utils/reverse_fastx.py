#!/usr/bin/env python3
import argparse
from typing import Literal

import pysam

from pybiotk.io import OpenFqGzip
from pybiotk.utils import reverse_seq


def reverse_fastx(filename: str, output: str, outformat: Literal['fastq', 'fasta'] = 'fastq'):
    if outformat == 'fastq':
        ostream = OpenFqGzip(output)
    else:
        ostream = open(output, "w")
        
    for entry in pysam.FastxFile(filename):
        name = entry.name
        sequence = reverse_seq(entry.sequence)
        comment = entry.comment
        quality = entry.quality if entry.quality is None else "".join(reversed(str(entry.quality)))
        if outformat == "fastq":
            ostream.write_entry(name, sequence, comment, quality)
        else:
            ostream.write(f">{name}\n{sequence}\n")
    ostream.close()
    

def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", dest="input", type=str, required=True, help="input *.fastx file.")
    parser.add_argument("-o", dest="output", type=str, required=True, help="output *.fastx file.")
    parser.add_argument("-f", dest="outformat", type=str, choices=["fasta", 'fastq'], default="fastq",  help="output file format fasta or fastq, default:fastq.")
    args = parser.parse_args()
    reverse_fastx(args.input, args.output, args.outformat)


if __name__ == "__main__":
    run()
