#!/usr/bin/env python3
import argparse
import re
import sys
import time
from typing import Sequence, Optional, TextIO

from numpy import block

from pybiotk.io import FastaFile, GtfFile
from pybiotk.utils import logging, ignore


class LocationFormatError(Exception):
    pass


def main(
    fasta: str,
    location:Optional[str] = None,
    gtf: Optional[str] = None, 
    regions: Optional[str] = None,
    separate: bool = False,
    gene_types: Optional[Sequence[str]] = None,
    transcript_types: Optional[Sequence[str]] = None,
    transcript_ids: Optional[Sequence[str]] = None,
    transcript_names: Optional[Sequence[str]] = None,
    gene_ids: Optional[Sequence[str]] = None,
    gene_names: Optional[Sequence[str]] = None,
    output: Optional[TextIO] = None):
    """Fetch sequences from genome fasta file.

    Args:
        fasta (str): genome fasta file.
        location (str, optional): location file. Defaults to None.
        gtf (str, optional): gtf file. Defaults to None.
        regions (str, optional): regions file. Defaults to None.
        separate (bool, optional): separate
        gene_types (Sequence[str], optional): gene types. Defaults to None.
        transcript_types (Sequence[str], optional): transcript types. Defaults to None.
        transcript_ids (Sequence[str], optional): transcript ids. Defaults to None.
        transcript_names (Sequence[str], optional): transcript names. Defaults to None.
        gene_ids (Sequence[str], optional): gene ids. Defaults to None.
        gene_names (Sequence[str], optional): gene names. Defaults to None.
        output (TextIO, optional): output file. Defaults to None.

    Returns:
        None
    """
    start = time.perf_counter()
    logging.info("loading genome fasta file ...")
    with FastaFile(fasta) as genome, output or sys.stdout as out:
        logging.info("genome fasta file loaded.")
        if gtf is not None:
            logging.info("loading gtf file ...")
            with GtfFile(gtf) as gtf:
                logging.info("gtf file loaded.")
                for transcript in gtf.iter_transcript(
                    gene_types=gene_types,
                    transcript_types=transcript_types,
                    transcript_ids=transcript_ids,
                    transcript_names=transcript_names,
                    gene_ids=gene_ids,
                    gene_names=gene_names
                    ):
                    if regions == "all":
                        blocks = [(transcript.start, transcript.end)]
                    elif regions == "exon":
                        blocks = transcript.exons()
                    elif regions == "intron":
                        blocks = transcript.introns()
                    elif regions == "fist_exon":
                        blocks = transcript.exons()[:1] if transcript.strand == "+" else transcript.exons()[-1:]
                    elif regions == "last_exon":
                        blocks = transcript.exons()[-1:] if transcript.strand == "+" else transcript.exons()[:1]

        elif location is not None:
            location_pattern = re.compile(r"(.+):(\d+)-(\d+)\(([+-])\)")
            loc = location_pattern.match(location)
            if loc is None:
                raise LocationFormatError(f"location format error: {location}")
            else:
                reference, start, end, strand = loc.groups()
                sequence = genome.fetchs(reference, start, end, strand)
                out.write(f">{location}\n{sequence}\n")
        else:
            for reference, sequence in genome:
                out.write(f">{reference}\n{sequence}\n")
    end = time.perf_counter()
    logging.info(f"task finished in {end - start:.2f} seconds.")



@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--fasta', dest="fasta", type=str, required=True, help="Genome fasta file.")
    parser.add_argument('-l', '--location', dest="location", type=str, default=None, help="chromosome:start-end(+/-). will be ignored if -r is set.")
    parser.add_argument('-g', '--gtf', dest="gtf", type=str, default=None, help="GTF file.")
    parser.add_argument('-r', '--regions', dest='regions', default='exons',
                        choices=['all', 'exons', 'introns', 'first_exon', 'last_exon', 'first_intron', 'last_intron', '5utr', '3utr', 'cds'],
                        type=str, help="annotation regions. need -g option.")
    parser.add_argument('--separate', dest='separate', action='store_true', help="separate exons or introns. when -r is set to exons or introns.")
    parser.add_argument('--gene_types', dest='gene_types', nargs="+",
                        default=None, help="choose gene types to filter gtf. need -g option.")
    parser.add_argument('--transcript_types', dest='transcript_types', nargs="+",
                        default=None, help="choose transcript types to filter gtf. need -g option.")
    parser.add_argument('--transcript_ids', dest='transcript_ids', nargs="+",
                        default=None, help="choose transcript ids to filter gtf. need -g option.")
    parser.add_argument('--transcript_names', dest='transcript_names', nargs="+",
                        default=None, help="choose transcript names to filter gtf. need -g option.")
    parser.add_argument('--gene_ids', dest='gene_ids', nargs="+",
                        default=None, help="choose gene ids to filter gtf. need -g option.")
    parser.add_argument('--gene_names', dest='gene_names', nargs="+",
                        default=None, help="choose gene names to filter gtf. need -g option.")
    parser.add_argument('-o', '--output', dest='output', type=argparse.FileType('w'),
                        default=sys.stdout, help="output file name.")
    
    args = parser.parse_args()
    if args.gtf is None and args.location is None:
        parser.error("please set -l or -g option.")
    if args.gtf is not None and args.location is not None:
        logging.warning("both -l and -g are set, -l will be ignored.")
    main(**vars(args))


if __name__ == "__main__":
    run()
