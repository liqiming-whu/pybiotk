#!/usr/bin/env python3
import os
import re
import sys
import time
import argparse
import pandas as pd
from typing import Tuple, List, NamedTuple, Optional
from pybiotk.io import GenomeFile, FastaFile
from pybiotk.utils import is_overlap, reverse_seq, logging
from pybiotk.intervals import merge_intervals
from ricpy.annotation.FixGapAnno import ReadAnno



class ReadAnno2(ReadAnno):
    def __init__(self, annofile: str, unmapped_reads: str):
        self.anno = pd.read_table(annofile, header=0, index_col=0)
        self.unmapped_reads = FastaFile(unmapped_reads)
        self.fix_gap_rna_anno = None
        self.fix_circle_gap = None
        
    
    def match_reads(self, genome: GenomeFile, refseq: Optional[FastaFile] = None, check_up_and_down: Optional[int] = 1000):
        fix_gap_rna_anno = []
        read_mate_name_set = set()
        fix_circle_gap = []
        for index, anno in self.anno.iterrows():
            read_name = index
            chrom = anno["chrom"]
            start = int(anno["start"])
            end = int(anno["end"])
            blocks = eval(anno["blocks"])
            strand = anno["strand"]
            if read_name.endswith("L"):
                new_read_name = read_name.rstrip("|L")
                read_mate_name = new_read_name + "|R"
            elif read_name.endswith("R"):
                new_read_name = read_name.rstrip("|R")
                read_mate_name = new_read_name + "|L"
            read_mate_name_set.add(read_mate_name)
            
            if read_mate_name not in self.unmapped_reads:
                continue
            
            unmap_seq = self.unmapped_reads.fetch(read_mate_name)

            reference_length = genome.get_chrom_length(chrom)
            right_region_end = end + 150 if end + 150 < reference_length else reference_length
            left_region_start = start - 150 if start > 150 else 0

            left_region = genome.fetch(chrom, [(left_region_start, start)], strand)
            right_region = genome.fetch(chrom, [(end, right_region_end)], strand)
            left_pos = [i.start() for i in re.finditer(unmap_seq, left_region)]
            right_pos = [i.start() for i in re.finditer(unmap_seq, right_region)]
            if not (left_pos or right_pos):
                continue
            elif right_pos and not left_pos:
                find_pos = right_pos[0]
                new_start = start
                if strand == '+':
                    new_end = end + find_pos + len(unmap_seq)
                else:
                    new_end = end + 150 - find_pos
            elif left_pos and not right_pos:
                find_pos = left_pos[-1]
                if strand == '+':
                    new_start = start - 150 + find_pos
                else:
                    new_start = start - find_pos - len(unmap_seq)
                new_end = end
            elif left_pos and right_pos:
                if right_pos[0] <= 150 - left_pos[-1] - len(unmap_seq):
                    find_pos = right_pos[0]
                    new_start = start
                    if strand == '+':
                        new_end = end + find_pos + len(unmap_seq)
                    else:
                        new_end = end + 150 - find_pos
                else:
                    find_pos = left_pos[-1]
                    if strand == '+':
                        new_start = start - 150 + find_pos
                    else:
                        new_start = start - find_pos
                    new_end = end

            if len(blocks) > 1:
                new_blocks_p = [(new_start, blocks[0][1])] + blocks[1:-1] + [(blocks[-1][0], new_end)]
                new_blocks = []
                for block in new_blocks_p:
                    if not block[0] == block[1]:
                        new_blocks.append(block)
                del new_blocks_p
            else:
                new_blocks = [(new_start, new_end)]
            if new_start == start:
                gap = [(end, new_end-len(unmap_seq))]
                gap_length = new_end - end - len(unmap_seq)
            else:
                gap = [(new_start + len(unmap_seq), start)]
                gap_length = start - new_start - len(unmap_seq)

            gap_sequence = None
            
            ref_sequence = refseq.fetch(new_read_name) if refseq is not None else None
            
            geneStart = set(anno["geneStart"].split(","))
            geneEnd = set(anno["geneEnd"].split(","))
            if refseq is not None:
                if gap_length:
                    gap_sequence = genome.fetch(chrom, gap, strand)
                    if gap_sequence.find(ref_sequence) > -1:
                        sys.stderr.write(f"{new_read_name} ref_sequence in gap: {gap}, gap_length: {gap_length}.")
                        continue
                blocks_sequence = genome.fetch_blocks(chrom, new_blocks, strand)
                if blocks_sequence.find(reverse_seq(ref_sequence)) > -1:
                    sys.stderr.write(f"{read_name} has reversed sequence of reference in blocks.\n")
                else:
                    if "*" in geneStart or "*" in geneEnd:
                        up_blocks = [(blocks[0][0]-1000, blocks[0][0])]
                        down_blocks = [(blocks[-1][1], blocks[-1][1]+1000)]
                        if strand == "-":
                            up_blocks, down_blocks = down_blocks, up_blocks
                        up_sequence = genome.fetch_blocks(chrom, up_blocks, strand)
                        down_sequence = genome.fetch_blocks(chrom, down_blocks, strand)
                        if up_sequence.find(ref_sequence) > -1:
                            sys.stderr.write(f"{read_name}(Intergenic) has reference in blocks upStream 1k.\n")
                            continue
                        if down_sequence.find(ref_sequence) > -1:
                            sys.stderr.write(f"{read_name}(Intergenic) has reference in blocks downStream 1k.\n")
                            continue
                    else:
                        gene_start = min(list(map(int, geneStart)))
                        gene_end = max(list(map(int, geneEnd)))
                        gene_sequence = genome.fetch_blocks(chrom, [(gene_start, gene_end)], strand)
                        if gene_sequence.find(reverse_seq(ref_sequence)) > -1:
                            sys.stderr.write(f"{read_name} has reversed sequence of reference in gene body.\n")
                        else:
                            if gene_sequence.find(ref_sequence) > -1:
                                sys.stderr.write(f"{read_name} has ref_sequence in gene body.\n")
                                continue
                            if check_up_and_down is not None:
                                up_blocks = [(gene_start-check_up_and_down, gene_start)]
                                down_blocks = [(gene_end, gene_end+ check_up_and_down)]
                                if strand == "-":
                                    up_blocks, down_blocks = down_blocks, up_blocks
                                up_sequence = genome.fetch_blocks(chrom, up_blocks, strand)
                                down_sequence = genome.fetch_blocks(chrom, down_blocks, strand)
                                if up_sequence.find(ref_sequence) > -1:
                                    sys.stderr.write(f"{read_name} has reference in gene upStream {check_up_and_down}.\n")
                                    continue
                                if down_sequence.find(ref_sequence) > -1:
                                    sys.stderr.write(f"{read_name} has reference in gene downStream {check_up_and_down}.\n")
                                    continue
                
            fix_gap_rna_anno.append([new_read_name, chrom, start, end, str(new_blocks), strand,
                                 anno["annotation"], ",".join(list(geneStart), ",".join(list(geneEnd))), anno["geneName"], anno["id"], anno["type"]])
            fix_circle_gap.append([new_read_name, chrom, gap, gap_length, strand, gap_sequence])
        
        self.fix_gap_rna_anno= fix_gap_rna_anno
        self.fix_circle_gap = fix_circle_gap
        
        return read_mate_name_set
    
    def unmap_reads_pair(self, unmap_reads_file: str, genome: GenomeFile, refseq: Optional[FastaFile] = None, check_up_and_down: Optional[int] = 1000):
        reads_mate_name_set = self.match_reads(genome, refseq, check_up_and_down)
        unmap_reads_dict = {}
        for name, seq in self.unmapped_reads:
            if name not in reads_mate_name_set:
                unmap_reads_dict[read_name] = seq
        
        with open(unmap_reads_file, "w") as f:
            for read_name in unmap_reads_dict:
                if read_name.endswith("R"):
                    continue
                new_read_name = read_name.rstrip("|L")
                read_mate_name = new_read_name + "|R"

                left_seq = unmap_reads_dict[read_name]
                if read_mate_name not in unmap_reads_dict:
                    continue
                right_seq = unmap_reads_dict[read_mate_name]
                new_seq = left_seq + right_seq
                f.write(f">{new_read_name}\n{new_seq}\n")


def main(genomefa: str, anno:str, unmap:str, merge: str, fix_gap: str, fix_gap_fa: str, circle_gap: str,
         ref_fasta: Optional[str] = None, check_up_and_down: Optional[int] = 1000):
    logging.info("start fix gap annotaion ...")
    start = time.perf_counter()
    genome = GenomeFile(genomefa)
    refseq = FastaFile(ref_fasta) if ref_fasta is not None else None
    read_anno = ReadAnno2(anno, unmap)
    logging.info("start to fix gap ...")
    read_anno.unmap_reads_pair(merge, genome, refseq, check_up_and_down)
    logging.info("start to write gap info ...")
    read_anno.writeout_gap_info(
        read_anno.fix_circle_gap, circle_gap)
    logging.info("start to write reads anno ...")
    read_anno.writeout_reads_anno(
        read_anno.fix_gap_rna_anno,
        genome, fix_gap, fix_gap_fa)
    end = time.perf_counter()
    logging.info(f"task finished in {end-start:.2f}s.")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    
    parser.add_argument("-g", dest="genome",
                        type=str, required=True, help="genome fasta")
    parser.add_argument("-a", dest="anno",
                        type=str, required=True, help="rna anno file")
    parser.add_argument("-u", dest="unmap",
                        type=str, required=True, help="unmap fasta")
    parser.add_argument("-o", dest="outdir", type=str,
                        default=os.getcwd(), help="output dir")
    parser.add_argument("--reference", dest="reference", type=str, default=None, help="refernce sequence")
    parser.add_argument("--check_up_and_down", dest="check_up_and_down", type=int, default=None, help="check_upstream and downstream")
    
    args = parser.parse_args()
    outdir = args.outdir
    main(args.genome, args.anno, args.unmap, os.path.join(outdir, "merged_unmapped.fa"), os.path.join(outdir, "fix_gap_anno.tsv"),
         os.path.join(outdir, "fix_gap.fa"), os.path.join(outdir, "circle_gap_info.tsv"), args.reference, args.check_up_and_down)


if __name__ == "__main__":
    run()
