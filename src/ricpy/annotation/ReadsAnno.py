#!/usr/bin/env python3
import os
import sys
import time
import argparse
from typing import Literal, Tuple, TextIO, Optional
from pybiotk.utils import logging, infer_fragment_strand, intervals_is_overlap, reverse_seq, blocks_len
from pybiotk.io import BamType, Bam, BamPE, GenomeFile, check_bam_type
from pybiotk.annodb import AnnoSet
from pybiotk.intervals import GRangeTree, merge_intervals
from pybiotk.utils.pyanno import load_grangetree


def annobam(filename: str,
            file_obj: TextIO,
            grangetree: GRangeTree,
            anno_fragments: bool = False,
            tss_region: Tuple[int, int] = (-1000, 1000),
            downstream: int = 3000,
            rule: str = "++,--",
            min_len: int = 0,
            ordered_by_name: bool = False,
            genomefile: Optional[str] = None,
            referencefile: Optional[str] = None,
            check_up_and_down: Optional[int] = 1000,
            filter_log: str = os.devnull,
            mapped_sequence: Optional[str] = None,
            unmapped_sequence: Optional[str] = None
            ):
    if unmapped_sequence is not None:
        unmapped = open(unmapped_sequence, "w")
    bamtype = check_bam_type(filename)
    filter_mode = False if not (genomefile is not None and referencefile is not None) else True
    
    if genomefile is not None:
        genome = GenomeFile(genomefile)
    
    if mapped_sequence is not None:
        mapped = open(mapped_sequence, "w")
    
    if filter_mode:
        reference = GenomeFile(referencefile)
        filter_log_obj = open(filter_log, "w")

    def anno_read(blocks, strand, read):
        start = int(blocks[0][0])
        end = int(blocks[-1][1])
        name = read.query_name
        chrom = read.reference_name
        fragment_strand = infer_fragment_strand(strand, rule, read.is_read2)
        genes = grangetree.find(read.reference_name, start, end, fragment_strand)
        
        if not genes:
            file_obj.write(f"{name}\t{chrom}\t{start}\t{end}\t{blocks}\t{fragment_strand}\tIntergenic\t*\t*\t*\t*\t*\n")
            if mapped_sequence is not None:
                mapped_seq = genome.fetch_blocks(chrom, blocks, fragment_strand)
                mapped.write(f"{name}\n{mapped_seq}\n")
        else:
            annoset = AnnoSet(gene.annotation(blocks, tss_region, downstream) for gene in genes)
            if filter_mode:
                ref_sequence = reference.fetch(name)
                anno_start = min(annoset.start)
                anno_end = max(annoset.end)
                gene_blocks = [(anno_start, anno_end)]
                gene_sequence = genome.fetch_blocks(chrom, gene_blocks, fragment_strand)
                if gene_sequence.find(reverse_seq(ref_sequence)) > -1:
                    filter_log_obj.write(f"{name} has reversed seqence of reference.\n")
                    file_obj.write(f"{name}\t{chrom}\t{start}\t{end}\t{blocks}\t{fragment_strand}\t{annoset}\n")
                    if mapped_sequence is not None:
                        mapped_seq = genome.fetch_blocks(chrom, blocks, fragment_strand)
                        mapped.write(f"{name}\n{mapped_seq}\n")
                    return
                
                if gene_sequence.find(ref_sequence) > -1:
                    filter_log_obj.write(f"{name} has reference seqence in gene body.\n")
                    return
                
                if check_up_and_down is not None:
                    up_blocks = [(anno_start, anno_start - check_up_and_down)]
                    down_blocks = [(anno_end, anno_end + check_up_and_down)]
                    if fragment_strand == "-":
                        up_blocks, down_blocks = down_blocks, up_blocks
                        
                    up_sequence = genome.fetch_blocks(chrom, up_blocks, fragment_strand)
                    down_sequence = genome.fetch_blocks(chrom, down_blocks, fragment_strand)
                    if up_sequence.find(ref_sequence) > -1:
                        filter_log_obj.write(f"{name} has reference in upStream {check_up_and_down}.\n")
                        return
                    if down_sequence.find(ref_sequence) > -1:
                        filter_log_obj.write(f"{name} has reference in downStream {check_up_and_down}.\n")
                        return
            file_obj.write(f"{name}\t{chrom}\t{start}\t{end}\t{blocks}\t{fragment_strand}\t{annoset}\n")
            if mapped_sequence is not None:
                mapped_seq = genome.fetch_blocks(chrom, blocks, fragment_strand)
                mapped.write(f"{name}\n{mapped_seq}\n")
                

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
                        merge_blocks = merge_intervals(blocks1, blocks2)
                        if min_len and blocks_len(merge_blocks) < min_len:
                            logging.info(f"{read1.query_name} < {min_len}nt, considered an unmapped read pair.")
                            continue
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
                    strand = '-' if read1.is_reverse else '+'
                    anno_read(blocks, strand, read2)
    else:
        with Bam(filename) as bam:
            i = 0
            for read in bam.iter(secondary=False, supplementary=False):
                if read.is_qcfail:
                    continue
                if read.is_unmapped :
                    if unmapped_sequence is not None:
                        unmapped.write(f">{read.query_name}\n{read.get_forward_sequence()}\n")
                    continue
                blocks = read.get_blocks()
                if min_len and blocks_len(blocks) < min_len:
                    logging.info(f"{read.query_name} < {min_len}nt, considered an unmapped read.")
                    if unmapped_sequence is not None:
                        unmapped.write(f">{read.query_name}\n{read.get_forward_sequence()}\n")
                    continue
                strand = '-' if read.is_reverse else '+'
                anno_read(blocks, strand, read)
                i += 1
            sys.stderr.write(f"annotate {i} reads.\n")
    if filter_mode:
        filter_log_obj.close()
    if mapped_sequence is not None:
        mapped.close()
    if unmapped_sequence is not None:
        unmapped.close()


def main(
    filename: str,
    outfilename: str,
    gtf_file: str,
    level: Literal["transcript", "gene"] = "transcript",
    tss_region: Tuple[int, int] = (-1000, 1000),
    downstream: int = 3000,
    strand: bool = True,
    rule: str = "1+-,1-+,2++,2--",
    annofragments: bool = False,
    min_len: int = 0,
    ordered_by_name: bool = False,
    genomefile: Optional[str] = None,
    referencefile: Optional[str] = None,
    check_up_and_down: Optional[int] = 1000,
    filter_log: str = os.devnull,
    mapped_sequence: Optional[str] = None,
    unmapped_sequence: Optional[str] = None,
):
    start = time.perf_counter()
    grangetree = load_grangetree(gtf_file, level, tss_region, downstream, strand)
    with open(outfilename, "w", encoding="utf-8") as annofile:
        annofile.write("seqname\tchrom\tstart\tend\tblocks\tstrand\tannotation\tgeneStart\tgeneEnd\tgeneName\tid\ttype\n")
        logging.info("start annotating, use bam mode ...")
        annobam(filename, annofile, grangetree, annofragments,
                tss_region, downstream, rule, min_len,
                ordered_by_name, genomefile, referencefile,
                check_up_and_down, filter_log, mapped_sequence, unmapped_sequence)
    end = time.perf_counter()
    logging.info(f"task completed in {end-start:.2f}s, annofile saved in {outfilename}")
    

def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", dest="input", type=str, required=True,
                        help="input bam file.")
    parser.add_argument("-o", "--output", dest="output", type=str, required=True, help="output file name.")
    parser.add_argument("-g", "--gtf", dest='gtf', required=True,
                        help="gtf file download from Genecode, or a sorted gtf file.")
    parser.add_argument("-l", "--level", dest="level", type=str, default="transcript", choices=["transcript", "gene"],
                        help="annotation level, transcript or gene.")
    parser.add_argument("--tss_region", dest="tss_region", type=int, nargs="+", default=[-1000, 1000],
                        help="choose region from tss.")
    parser.add_argument("--downstream", dest="downstream", type=int, default=3000, help="downstream length from tes.")
    parser.add_argument("-s", "--strand", dest="strand", action="store_true", help="require same strandedness.")
    parser.add_argument("--rule", dest="rule", type=str, default="++,--",
                        choices=["1+-,1-+,2++,2--", "1++,1--,2+-,2-+", "+-,-+", "++,--"],
                        help="how read(s) were stranded during sequencing. only for bam.")
    parser.add_argument("-p", "--pair", dest="pair", action="store_true",
                        help="annotate fragments instead of reads.")
    parser.add_argument("--min_len", dest="len", type=int, default=0,
                        help="length < len is considered an unmapped read.")
    parser.add_argument("--ordered_by_name", dest="ordered_by_name", action="store_true",
                        help="if input bam is ordered by name, only for pair-end bam.")
    parser.add_argument("--filter_mode", dest="filter_mode", action="store_true", help="enable filter mode, --genome and --refernce needed.")
    parser.add_argument("--filter_log", dest="filter_log", type=str, default=os.devnull, help="filter log.")
    parser.add_argument('--genome', dest="genome", type=str, default=None, help="genome sequence file")
    parser.add_argument("--reference", dest="reference", type=str, default=None, help="refernce sequence")
    parser.add_argument("--check_up_and_down", dest="check_up_and_down", type=int, default=None, help="check_upstream and downstream")
    parser.add_argument("--mapped_sequence", dest="mapped_sequence", type=str, default=None, help="mapped sequence, --genome needed")
    parser.add_argument("--unmapped_sequence", dest="unmapped_sequence", type=str, default=None, help="unmapped sequence. not support --pair mode")

    args = parser.parse_args()
    assert len(args.tss_region) == 2, "tss_region must be a tuple of 2 elements."
    
    if args.filter_mode:
        if not (args.genome is not None and args.reference is not None):
            args = parser.parse_args(['-h'])
    
    if args.mapped_sequence is not None:
        if args.genome is None:
            args = parser.parse_args(['-h'])
    
    if args.unmapped_sequence is not None:
        if args.pair is not None:
            args = parser.parse_args(['-h'])
    
    main(args.input, args.output, args.gtf, args.level, args.tss_region,
         args.downstream, args.strand, args.rule, args.pair, args.len,
         args.ordered_by_name, args.genome, args.reference, args.check_up_and_down,
         args.mapped_sequence, args.unmapped_sequence)


if __name__ == "__main__":
    run()
