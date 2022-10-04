#!/usr/bin/env python3
import argparse
import sys
from typing import Literal, Sequence, Optional, TextIO

from pybiotk.io import GtfFile
from pybiotk.utils import ignore
from stream import write


def main(
        istream: TextIO,
        name: Literal["gene_id", "gene_name", "transcript_id",
                      "transcript_name"] = "transcript_id",
        feature: Literal["exon", "gene"] = "exon",
        gene_types: Optional[Sequence[str]] = None,
        transcript_types: Optional[Sequence[str]] = None,
        transcript_ids: Optional[Sequence[str]] = None,
        transcript_names: Optional[Sequence[str]] = None,
        gene_ids: Optional[Sequence[str]] = None,
        gene_names: Optional[Sequence[str]] = None,
        output: TextIO = sys.stdout,
        outfmt: Literal["bed12", "bed6", "intron", "gene_info", "trans_info"] = "bed12"):
    with GtfFile(istream) as gtf:
        if outfmt == 'bed12':
            gtf.to_bed12(name, gene_types, transcript_types, transcript_ids, transcript_names, gene_ids, gene_names) | write(output)
        elif outfmt == 'bed6':
            gtf.to_bed6(feature, name, gene_types, transcript_types, transcript_ids, transcript_names, gene_ids, gene_names) | write(output)
        elif outfmt == 'intron':
            gtf.to_intron(name, gene_types, transcript_types, transcript_ids, transcript_names, gene_ids, gene_names) | write(output)
        elif outfmt == 'gene_info':
            gtf.to_gene_info(gene_types, transcript_types, transcript_ids, transcript_names, gene_ids, gene_names) | write(output)
        elif outfmt == 'trans_info':
            gtf.to_trans_info(gene_types, transcript_types, transcript_ids, transcript_names, gene_ids, gene_names) | write(output)
        else:
            raise RuntimeError(f"{outfmt} not support.")


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('gtf', nargs='?', type=argparse.FileType('r'),
                        default=(None if sys.stdin.isatty() else sys.stdin),
                        help="gtf file download from Genecode, or a sorted gtf file.")
    parser.add_argument('-n', '--name', dest='name', default='transcript_id',
                        choices=("gene_id", "gene_name", "transcript_id", "transcript_name"),
                        help="""attribute ID by which to group bed entries.
                                'gene_id' or 'gene_name' for gene.bed6, gene_info,
                                'transcript_id' for genome.bed12, exon, intron, transcrpt_info"""
                        )
    parser.add_argument('-f', '--feature', dest='feature', default='exon',
                        choices=("exon", "gene"),
                        help="""
                        annotation feature to join, all others are filtered out:
                        'exon' genome.bed12, exon, intron, 'gene' for gene.bed6, gene_info."""
                        )
    parser.add_argument('-o', dest='output', type=argparse.FileType('w'),
                        default=sys.stdout, help="output file name.")
    parser.add_argument('--outfmt', dest='outfmt', default='bed12',
                        choices=("bed12", "bed6", "intron", "gene_info", "trans_info"),
                        help="choose output file format:bed12, bed6, intron, info")
    parser.add_argument('--gene_types', dest='gene_types', nargs="+",
                        default=None, help="choose gene types to filter bed.")
    parser.add_argument('--transcript_types', dest='transcript_types', nargs="+",
                        default=None, help="choose transcript types to filter bed.")
    parser.add_argument('--transcript_ids', dest='transcript_ids', nargs="+",
                        default=None, help="choose transcript ids to filter bed.")
    parser.add_argument('--transcript_names', dest='transcript_names', nargs="+",
                        default=None, help="choose transcript names to filter bed.")
    parser.add_argument('--gene_ids', dest='gene_ids', nargs="+",
                        default=None, help="choose gene ids to filter bed.")
    parser.add_argument('--gene_names', dest='gene_names', nargs="+",
                        default=None, help="choose gene names to filter bed.")

    args = parser.parse_args()
    if args.gtf is None:
        parser.print_help()
        sys.exit(1)

    main(args.gtf, args.name, args.feature,
         args.gene_types, args.transcript_types,
         args.transcript_ids, args.transcript_names,
         args.gene_ids, args.gene_names,
         args.output, args.outfmt)


if __name__ == "__main__":
    run()
