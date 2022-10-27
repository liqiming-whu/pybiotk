#!/usr/bin/env python3
import argparse
import os
import sys
import time
from typing import Literal, Tuple, TextIO

from pybiotk.annodb import AnnoSet, Gene
from pybiotk.intervals import GRangeTree, merge_intervals
from pybiotk.io import (
    GtfFile, Openbed, BamType, Bam, BamPE, check_bam_type
)
from pybiotk.utils import logging, infer_fragment_strand, intervals_is_overlap
from stream import groupby


def load_grangetree(
    filename: str,
    level: Literal["transcript", "gene"] = "transcript",
    tss_region: Tuple[int, int] = (-1000, 1000),
    downstream: int = 3000,
    strand: bool = True
) -> GRangeTree:

    start = time.perf_counter()
    logging.info(f"start to load {level} grangetree ...")
    grangetree = GRangeTree(strand)
    with GtfFile(filename) as gtf:
        i = 0
        if level == "transcript":
            for transcript in gtf.to_transcript():
                grange = [*transcript.tss_region(tss_region), *transcript.downstream(downstream)]
                grange.sort()
                grangetree.add(transcript, transcript.chrom, grange[0], grange[-1], transcript.strand)
                del grange
                i += 1
                if i % 1000 == 0:
                    sys.stderr.write(f"\rload {i:6} {level}s.")
        else:
            for transcripts in gtf.to_transcript() | groupby(lambda x: x.gene_id):
                gene = Gene.init_by_transcripts(transcripts)
                grange = [*gene.tss_region(tss_region), *gene.downstream(downstream)]
                grange.sort()
                grangetree.add(gene, gene.chrom, grange[0], grange[-1], gene.strand)
                del grange
                i += 1
                if i % 1000 == 0:
                    sys.stderr.write(f"\rload {i:6} {level}s.")
        sys.stderr.write(f"\rload {i:6} {level}s.\n")
    end = time.perf_counter()
    logging.info(f"grangetree load completed in {end-start}s.")
    return grangetree


def annobam(filename: str,
            file_obj: TextIO,
            grangetree: GRangeTree,
            anno_fragments: bool = False,
            tss_region: Tuple[int, int] = (-1000, 1000),
            downstream: int = 3000,
            rule: str = "1+-,1-+,2++,2--",
            ordered_by_name: bool = False
            ):
    bamtype = check_bam_type(filename)

    def anno_read(_blocks, _strand, _read):
        start = int(_blocks[0][0])
        end = int(_blocks[-1][1])
        fragment_strand = infer_fragment_strand(_strand, rule, _read.is_read2)
        genes = grangetree.find(_read.reference_name, start, end, fragment_strand)
        if not genes:
            file_obj.write(f"{_read.query_name}\t{_read.reference_name}\t{start}\t{end}\t{_blocks}\t{fragment_strand}\tIntergenic\t*\t*\t*\t*\t*\n")
        else:
            annoset = AnnoSet(gene.annotation(_blocks, tss_region, downstream) for gene in genes)
            gene_types = set(gene.gene_type for gene in genes if gene.id in set(annoset.id))
            if not ("protein_coding" in gene_types and "protein_coding" not in set(annoset.type)):
                annoset.type = gene_types
            file_obj.write(f"{_read.query_name}\t{_read.reference_name}\t{start}\t{end}\t{_blocks}\t{fragment_strand}\t{annoset}\n")

    if bamtype is BamType.PE and anno_fragments:
        with BamPE(filename) as bam:
            bam.ordered_by_name = ordered_by_name
            for read1, read2 in bam.iter_pair(properly_paired=False):
                if read1 is not None and read2 is not None:
                    strand1 = "-" if read1.is_reverse else "+"
                    strand2 = "-" if read2.is_reverse else "+"
                    blocks1 = read1.get_blocks()
                    blocks2 = read2.get_blocks()
                    if read1.reference_name == read2.reference_name and strand1 == strand2 and intervals_is_overlap(blocks1, blocks2):
                        merge_blocks = merge_intervals(blocks1 + blocks2)
                        anno_read(merge_blocks, strand1, read1)
                    else:
                        anno_read(blocks1, strand1, read1)
                        anno_read(blocks2, strand2, read2)
                if read1 is not None and read2 is None:
                    blocks = read1.get_blocks()
                    strand = '-' if read1.is_reverse else '+'
                    anno_read(blocks, strand, read1)
                if read2 is not None and read1 is None:
                    blocks = read2.get_blocks()
                    strand = '-' if read2.is_reverse else '+'
                    anno_read(blocks, strand, read2)
    else:
        with Bam(filename) as bam:
            i = 0
            for read in bam.iter_mapped():
                blocks = read.get_blocks()
                strand = '-' if read.is_reverse else '+'
                anno_read(blocks, strand, read)
                i += 1
            sys.stderr.write(f"annotate {i} reads.\n")


def annobed(filename: str,
            file_obj: TextIO,
            grangetree: GRangeTree,
            tss_region: Tuple[int, int] = (-1000, 1000),
            downstream: int = 3000
            ):
    with Openbed(filename) as bedfile:
        i = 0
        for bed in bedfile:
            genes = grangetree.find(bed.chrom, bed.start, bed.end, bed.strand)
            if not genes:
                file_obj.write(f"{bed.name}\t{bed.chrom}\t{bed.start}\t{bed.end}\t{bed.strand}\tIntergenic\t*\t*\t*\t*\t*\n")
            else:
                annoset = AnnoSet(gene.annotation([(bed.start, bed.end)], tss_region, downstream) for gene in genes)
                gene_types = set(gene.gene_type for gene in genes if gene.id in set(annoset.id))
                if not ("protein_coding" in gene_types and "protein_coding" not in set(annoset.type)):
                    annoset.type = gene_types
                file_obj.write(f"{bed.name}\t{bed.chrom}\t{bed.start}\t{bed.end}\t{bed.strand}\t{annoset}\n")
            i += 1
        sys.stderr.write(f"annotate {i} beds.\n")


def main(
    filename: str,
    outfilename: str,
    gtf_file: str,
    level: Literal["transcript", "gene"] = "transcript",
    tss_region: Tuple[int, int] = (-3000, 0),
    downstream: int = 3000,
    strand: bool = True,
    rule: str = "1+-,1-+,2++,2--",
    annofragments: bool = False,
    ordered_by_name: bool = False
):
    start = time.perf_counter()
    filetype = os.path.splitext(filename)[1]
    grangetree = load_grangetree(gtf_file, level, tss_region, downstream, strand)
    if outfilename == "-":
        ostream = sys.stdout
    else:
        ostream = open(outfilename, "w", encoding="utf-8")
    with ostream as annofile:
        if filetype == ".bam":
            annofile.write("seqname\tchrom\tstart\tend\tblocks\tstrand\tannotation\tgeneStart\tgeneEnd\tgeneName\tid\tgeneType\n")
            logging.info("start annotating, use bam mode ...")
            annobam(filename, annofile, grangetree, annofragments, tss_region, downstream, rule, ordered_by_name)
        elif filetype.startswith(".bed"):
            annofile.write("seqname\tchrom\tstart\tend\tstrand\tannotation\tgeneStart\tgeneEnd\tgeneName\tid\tgeneType\n")
            logging.info("start annotating, use bed mode ...")
            annobed(filename, annofile, grangetree, tss_region, downstream)
        else:
            raise RuntimeError(f"Unable to infer file type as bam or bed from filename: {filename}.")
    end = time.perf_counter()
    logging.info(f"task completed in {end-start:.2f}s, annofile saved in {outfilename}")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", dest="input", type=str, required=True,
                        help="input file, bam or bed. The file type will be inferred from the filename suffix ['*.bam', '*.bed'].")
    parser.add_argument("-o", "--output", dest="output", type=str, default="-", help="output file name. [stdout]")
    parser.add_argument("-g", "--gtf", dest='gtf', required=True,
                        help="gtf file download from Genecode, or a sorted gtf file.")
    parser.add_argument("-l", "--level", dest="level", type=str, default="transcript", choices=("transcript", "gene"),
                        help="annotation level, transcript or gene.")
    parser.add_argument("--tss_region", dest="tss_region", type=int, nargs="+", default=[-3000, 0],
                        help="choose region from tss.")
    parser.add_argument("--downstream", dest="downstream", type=int, default=3000, help="downstream length from tes.")
    parser.add_argument("-s", "--strand", dest="strand", action="store_true", help="require same strandedness.")
    parser.add_argument("--rule", dest="rule", type=str, default="1+-,1-+,2++,2--",
                        choices=("1+-,1-+,2++,2--", "1++,1--,2+-,2-+", "+-,-+", "++,--"),
                        help="how read(s) were stranded during sequencing. only for bam.")
    parser.add_argument("-p", "--pair", dest="pair", action="store_true",
                        help="annotate fragments instead of reads.")
    parser.add_argument("--ordered_by_name", dest="ordered_by_name", action="store_true",
                        help="if input bam is ordered by name, only for pair-end bam.")

    args = parser.parse_args()

    if not len(args.tss_region) == 2:
        parser.error("--tss_region must be a tuple of 2 elements.")
    main(args.input, args.output, args.gtf, args.level, args.tss_region, args.downstream, args.strand, args.rule, args.pair, args.ordered_by_name)


if __name__ == "__main__":
    run()
