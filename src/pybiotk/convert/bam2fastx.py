#!/usr/bin/env python3
import time
import argparse
import pysam
from typing import Literal, Union, TextIO, Optional
from pybiotk.io import OpenFqGzip, Bam, BamPE, BamType, check_bam_type
from pybiotk.utils import logging


def write_recode(fileobj: Union[TextIO, OpenFqGzip], record: pysam.libcfaidx.FastxRecord, outformat: Literal['fastq', 'fasta']):
    if outformat == 'fastq':
        fileobj.write_fastx_recode(record)
    else:
        fileobj.write(f">{record.name}\n{record.sequence}\n")


def main(filename: str, prefix: str, bamtype: Optional[BamType] = None, outformat: Literal['fastq', 'fasta'] = "fasta", ordered_by_name: bool = False):
    logging.info(f"start convert {filename} to {outformat} ...")
    start = time.perf_counter()
    if bamtype is None:
        bamtype = check_bam_type(filename)

    if bamtype is BamType.SE:
        if outformat == 'fasta':
            out = open(prefix + ".fasta", "w")
        else:
            out = OpenFqGzip(prefix + ".fastq.gz")
        with Bam(filename) as bam:
            logging.info(f"writing reads to {out.name} ...")
            for record in bam.to_fastx_recorde():
                write_recode(out, record, outformat)
        out.close()
    elif bamtype is BamType.PE:
        if outformat == 'fasta':
            out_r1 = open(prefix + ".R1.fa", "w")
            out_r2 = open(prefix + ".R2.fa", "w")
            unpaired_r1 = open(prefix + ".unpaired.R1.fa", "w")
            unpaired_r2 = open(prefix + ".unpaired.R2.fa", "w")
        else:
            out_r1 = OpenFqGzip(prefix + ".R1.fq.gz")
            out_r2 = OpenFqGzip(prefix + ".R2.fq.gz")
            unpaired_r1 = OpenFqGzip(prefix + ".unpaired.R1.fq.gz")
            unpaired_r2 = OpenFqGzip(prefix + ".unpaired.R2.fq.gz")

        with BamPE(filename) as bampe:
            bampe.ordered_by_name = ordered_by_name
            logging.info(f"writing property paired reads to {out_r1.name} {out_r2.name} ...")
            for record1, record2 in bampe.to_fastx_recorde_pair():
                write_recode(out_r1, record1, outformat)
                write_recode(out_r2, record2, outformat)
            logging.info(f"writing unpaired reads to {unpaired_r1.name} ...")
            for record in bampe.to_fastx_recorde_unpaired(terminal="read1"):
                write_recode(unpaired_r1, record, outformat)
            logging.info(f"writing unpaired reads to {unpaired_r2.name} ...")
            for record in bampe.to_fastx_recorde_unpaired(terminal="read2"):
                write_recode(unpaired_r2, record, outformat)
        out_r1.close()
        out_r2.close()
        unpaired_r1.close()
        unpaired_r2.close()
        end = time.perf_counter()
        logging.info(f"task finished in {end-start:.2f}s.")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str, help="input bam.")
    parser.add_argument("--bamtype", dest="bamtype", type=str, choices=["SE", "PE"],
                        default=None, help="input bam type single end or pair end, default:autodectect.")
    parser.add_argument("-o", dest="outprefix", type=str, required=True,
                        help="output file name prefix.")
    parser.add_argument("--outformat", dest="outformat", type=str, choices=["fasta", 'fastq'],
                        default="fastq",  help="output file format fasta or fastq, default:fastq.")
    parser.add_argument("--ordered_by_name", dest="ordered_by_name", action="store_true",
                        help="if input bam is ordered by name.")

    args = parser.parse_args()
    bamtype = None if args.bamtype is None else BamType(args.bamtype)
    main(args.input, args.outprefix, bamtype, args.outformat, args.ordered_by_name)


if __name__ == "__main__":
    run()
