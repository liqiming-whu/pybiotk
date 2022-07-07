#!/usr/bin/env python3
import sys
import argparse
import pysam
import numpy as np
from typing import Optional
from pybiotk.utils import logging, ignore
from pybiotk.io import count_bam_size


def main(input: str = "-", output: str = "-", bam: bool = False, bamsize: Optional[int] = None, count: int = 10000):
    mode = "wb" if bam else "w"
    alignment = []
    header = None
    if bamsize is None:
        logging.info("bamsize is unknown, it will take some time and memory to calculate ...")
        if input == "-":
            with pysam.AlignmentFile(input) as alignmentfile:
                header = alignmentfile.header
                alignment = list(alignmentfile)
                bamsize = len(alignment)
        else:
            bamsize = count_bam_size(input)
    
    index_arr = np.arange(bamsize)
    np.random.shuffle(index_arr)
    select = set(index_arr[:count])
    
    if not alignment:
        alignment = pysam.AlignmentFile(input)
        header = alignment.header
        
    with pysam.AlignmentFile(output, mode, header=header) as outf:
        for index, segment in enumerate(alignment):
            if index in select:
                outf.write(segment)


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", dest="filename", type=str, default=(None if sys.stdin.isatty() else "-"), help="the input BAM file, must with header [stdin]")
    parser.add_argument("-o", "--output", dest="file", type=str, default="-", help="Write output to FILE [stdout]")
    parser.add_argument("-b", "--bam", dest="bam", action="store_true", help="Output BAM")
    parser.add_argument("--bamsize", dest="bamsize", type=int, default=None, help="bamsize. If unknown, it will take some time and memory to calculate")
    parser.add_argument("-n", dest="count", type=int, default=10000, help="number of alignments to grab")
    args = parser.parse_args()
    if args.filename is None:
        args = parser.parse_args(['-h'])
    try:
        main(args.filename, args.file, args.bam, args.bamsize, args.count)
    except ValueError:
        logging.error("An error occurs, input BAM file must have a header, use '-h' option when using samtools view.")
        raise


if __name__ == "__main__":
    run()
