#!/usr/bin/env python3
import argparse
import os
import sys
import time

from pybiotk.io import GenomeFile
from pybiotk.utils import logging, ignore


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(dest="genome", type=str, help="Genome fasta file.")
    parser.add_argument("--size", dest="size", type=int, default=100000, help="dataset size.")
    parser.add_argument("--prefix", dest="prefix", type=str, default="test", help="dataset name prefix.")
    parser.add_argument("--length", dest="length", type=int, default=50, help="random sequence length.")
    parser.add_argument("--random_on", dest="random_on", nargs="+", default=None, help="random on specific chroms. default:All chroms.")
    parser.add_argument("--load2dict", dest="load", action="store_true", help="load genome into dict.")

    return parser


def main():
    parser = parse_args()
    if len(sys.argv[1:]) < 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    faindex = args.genome + ".fai"
    regenerate = False
    if not os.path.exists(faindex):
        regenerate = True
        logging.info("Fasta index file is not exist, generating ...")
    else:
        logging.info(f"Fasta index file {faindex} found.")
    index_start = time.perf_counter()
    genome = GenomeFile(args.genome)
    index_end = time.perf_counter()
    if regenerate:
        logging.info(f"File {faindex} generated in {index_end-index_start:.2f}s.")
    genome.random_on_chroms = args.random_on

    if args.load:
        load_start = time.perf_counter()
        logging.info("Load genome into dict ...")
        genome.load_into_dict()
        load_end = time.perf_counter()
        logging.info(f"Load completed in {load_end-load_start:.2f}s.")

    fastafile = args.prefix + ".fa"
    bed_file = args.prefix + ".bed"
    size = args.size
    logging.info("Generating random sequence ...")
    start = time.perf_counter()
    with open(fastafile, "w") as fa, open(bed_file, "w", encoding="utf-8") as bed:
        for i in range(size):
            region, seq = genome.random(args.length, args.random_on)
            fa.write(f">{region}\n{seq}\n")
            bed.write(f"{region.chrom}\t{region.start}\t{region.end}\t{i+1}\t.\t{region.strand}\n")            
            width = len(str(size))
            sys.stderr.write(f"\r{i+1:{width}}/{size}\t\t\t\t{(i+1)*100/size:.2f}%")
    sys.stderr.write("\n")
    genome.close()
    end = time.perf_counter()
    logging.info(f"{size} random reads generated in {end-start:.2f}s.")
    logging.info(f"Datasets saved in {fastafile}, {bed_file}.")


@ignore
def run():
    main()


if __name__ == "__main__":
    run()
