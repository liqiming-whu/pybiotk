#!/usr/bin/env python3
import argparse
import sys
from typing import Sequence, Optional, TextIO

from pybiotk.io import GtfFile
from pybiotk.utils import ignore
from stream import write


def main(istream: TextIO,
         features: Optional[Sequence[str]] = None,
         gene_types: Optional[Sequence[str]] = None,
         transcript_types: Optional[Sequence[str]] = None,
         output: Optional[TextIO] = None,
         transcript_ids: Optional[Sequence[str]] = None,
         transcript_names: Optional[Sequence[str]] = None,
         gene_ids: Optional[Sequence[str]] = None,
         gene_names: Optional[Sequence[str]] = None
         ):
    with GtfFile(istream) as gtf:
        gtf.iter(features=features,
                 gene_types=gene_types,
                 transcript_types=transcript_types,
                 transcript_ids=transcript_ids,
                 transcript_names=transcript_names,
                 gene_ids=gene_ids,
                 gene_names=gene_names
                 ) | write(output)


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtf', nargs='?', type=argparse.FileType('r'),
                        default=(None if sys.stdin.isatty() else sys.stdin),
                        help="gtf file download from Genecode, or a sorted gtf file.")
    parser.add_argument('-f', '--features', dest='features', default=['exon'], nargs="+", help="annotation features.")
    parser.add_argument('-g', '--gene_types', dest='gene_types', nargs="+",
                        default=None, help="choose gene types to filter gtf.")
    parser.add_argument('-t', '--transcript_types', dest='transcript_types', nargs="+",
                        default=None, help="choose transcript types to filter gtf.")
    parser.add_argument('-o', dest='output', type=argparse.FileType('w'),
                        default=sys.stdout, help="output file name.")
    parser.add_argument('--transcript_ids', dest='transcript_ids', nargs="+",
                        default=None, help="choose transcript ids to filter gtf.")
    parser.add_argument('--transcript_names', dest='transcript_names', nargs="+",
                        default=None, help="choose transcript names to filter gtf.")
    parser.add_argument('--gene_ids', dest='gene_ids', nargs="+",
                        default=None, help="choose gene ids to filter gtf.")
    parser.add_argument('--gene_names', dest='gene_names', nargs="+",
                        default=None, help="choose gene names to filter gtf.")
    args = parser.parse_args()
    if args.gtf is None:
        parser.print_help()
        sys.exit(1)

    main(args.gtf, args.features, args.gene_types, args.transcript_types, args.output,
         args.transcript_ids, args.transcript_names, args.gene_ids, args.gene_names)


if __name__ == "__main__":
    run()
