#!/usr/bin/env python3
import os
import sys
import time
import argparse
import pandas as pd
from typing import Tuple, List, NamedTuple, Optional
from pybiotk.io import GenomeFile, FastaFile
from pybiotk.utils import is_overlap, reverse_seq, logging
from pybiotk.intervals import merge_intervals


class Mergeblocks(NamedTuple):
    blocks: Optional[List[Tuple[int, int]]] = None
    gap: Optional[Tuple[int, int]] = None
    gap_length: Optional[int] = None
    

def merge_blocks(blocks1: List[Tuple[int, int]], blocks2: List[Tuple[int, int]]) -> Mergeblocks:
    if blocks1[0][0] > blocks2[0][0]:
        blocks1, blocks2 = blocks2, blocks1
    blocks1_start = blocks1[0][0]
    blocks1_end = blocks1[-1][1]
    blocks2_start = blocks2[0][0]
    blocks2_end = blocks2[-1][1]
    mergeblocks = None
    gap = None,
    gap_length = None
    if is_overlap((blocks1_start, blocks1_end), (blocks2_start, blocks2_end)):
        gap_length = 0
        mergeblocks = merge_intervals(blocks1 + blocks2)
    else:
        mergeblocks = blocks1[:-1] + [(blocks1[-1][0], blocks2[0][1])] + blocks2[1:]
        gap = [(blocks1_end, blocks2_start)]
        gap_length = blocks2_start - blocks1_end
    
    return Mergeblocks(mergeblocks, gap, gap_length)


class ReadAnno:
    def __init__(self, annofile):
        self.anno = pd.read_table(annofile, header=0, index_col=0)
        self.one_side_rna_anno = None
        self.two_rna_anno = None
        self.fix_gap_rna_anno = None
        self.fix_circle_gap = None
        
    @staticmethod
    def pd_series_to_list(series: pd.Series):
        return [series.name] + series.to_list()

    def classify_reads(self) -> List[pd.Series]:
        read_pair_names = set()
        for read_name in self.anno.index:
            read_pair_names.add(read_name.rstrip("L").rstrip('R'))
        read_pair_names = list(read_pair_names)

        one_side_rna_anno = []
        two_rna_anno = []
        for read_pair_name in read_pair_names:
            left_anno = read_pair_name + 'L'
            right_anno = read_pair_name + 'R'

            if left_anno in self.anno.index and right_anno in self.anno.index:
                two_rna_anno.append(
                    [self.anno.loc[left_anno], self.anno.loc[right_anno]])
            elif left_anno in self.anno.index:
                one_side_rna_anno.append(self.pd_series_to_list(self.anno.loc[left_anno]))
            elif right_anno in self.anno.index:
                one_side_rna_anno.append(self.pd_series_to_list(self.anno.loc[right_anno]))

        self.one_side_rna_anno = one_side_rna_anno
        return two_rna_anno


    def fix_gap(self, genome: GenomeFile, refseq: Optional[FastaFile] = None, check_up_and_down: Optional[int] = 1000):
        two_rna_anno_set = self.classify_reads()
        two_rna_anno = []
        fix_gap_rna_anno = []
        fix_circle_gap = []
        for read_L, read_R in two_rna_anno_set:
            read_L_chrom = read_L["chrom"]
            read_L_blocks = eval(read_L["blocks"])
            read_L_strand = read_L["strand"]
            read_L_gene_names = set(read_L["geneName"].split(","))
            read_R_chrom = read_R["chrom"]
            read_R_blocks = eval(read_R["blocks"])
            read_R_strand = read_R["strand"]
            read_R_gene_names = set(read_R["geneName"].split(","))
            if read_L_chrom == read_R_chrom and read_L_strand == read_R_strand and (read_L_gene_names & read_R_gene_names):
                mergeblocks = merge_blocks(read_L_blocks, read_R_blocks)
                blocks = mergeblocks.blocks
                gap = mergeblocks.gap
                gap_length = mergeblocks.gap_length
                read_name = read_L.name.rstrip("|L")
                
                if refseq is not None:
                    ref_sequence = refseq.fetch(read_name)
                
                gap_sequence = None
                if gap_length and refseq is not None:
                    gap_sequence = genome.fetch_blocks(read_L_chrom, gap, read_L_strand)
                    if gap_sequence.find(ref_sequence) > -1:
                        sys.stderr.write(f"{read_name} ref_sequence in gap: {gap}, gap_length: {gap_length}.\n")
                        continue

                geneStart = set(read_L["geneStart"].split(",")) | set(read_R["geneStart"].split(","))
                geneEnd = set(read_L["geneEnd"].split(",")) | set(read_R["geneEnd"].split(","))
                
                if refseq is not None:
                    blocks_sequence = genome.fetch_blocks(read_L_chrom, blocks, read_L_strand)
                    if blocks_sequence.find(reverse_seq(ref_sequence)) > -1:
                        sys.stderr.write(f"{read_name} has reversed sequence of reference in blocks.\n")
                    else:
                        if "*" in geneStart or "*" in geneEnd:
                            up_blocks = [(blocks[0][0]-1000, blocks[0][0])]
                            down_blocks = [(blocks[-1][1], blocks[-1][1]+1000)]
                            if read_L_strand == "-":
                                up_blocks, down_blocks = down_blocks, up_blocks
                            up_sequence = genome.fetch_blocks(read_L_chrom, up_blocks, read_L_strand)
                            down_sequence = genome.fetch_blocks(read_L_chrom, down_blocks, read_L_strand)
                            if up_sequence.find(ref_sequence) > -1:
                                sys.stderr.write(f"{read_name}(Intergenic) has reference in blocks upStream 1k.\n")
                                continue
                            if down_sequence.find(ref_sequence) > -1:
                                sys.stderr.write(f"{read_name}(Intergenic) has reference in blocks downStream 1k.\n")
                                continue
                        else:
                            start = min(list(map(int, geneStart)))
                            end = max(list(map(int, geneEnd)))
                            gene_sequence = genome.fetch_blocks(read_L_chrom, [(start, end)], read_L_strand)
                            if gene_sequence.find(reverse_seq(ref_sequence)) > -1:
                                sys.stderr.write(f"{read_name} has reversed sequence of reference in gene body.\n")
                            else:
                                if gene_sequence.find(ref_sequence) > -1:
                                    sys.stderr.write(f"{read_name} has ref_sequence in gene body.\n")
                                    continue
                                if check_up_and_down is not None:
                                    up_blocks = [(start-check_up_and_down, start)]
                                    down_blocks = [(end, end+ check_up_and_down)]
                                    if read_L_strand == "-":
                                        up_blocks, down_blocks = down_blocks, up_blocks
                                    up_sequence = genome.fetch_blocks(read_L_chrom, up_blocks, read_L_strand)
                                    down_sequence = genome.fetch_blocks(read_L_chrom, down_blocks, read_L_strand)
                                    if up_sequence.find(ref_sequence) > -1:
                                        sys.stderr.write(f"{read_name} has reference in gene upStream {check_up_and_down}.\n")
                                        continue
                                    if down_sequence.find(ref_sequence) > -1:
                                        sys.stderr.write(f"{read_name} has reference in gene downStream {check_up_and_down}.\n")
                                        continue

                if gap_length > 150:
                    two_rna_anno.append([self.pd_series_to_list(read_L), self.pd_series_to_list(read_R)])
                    sys.stderr.write(f"{read_name} gap: {gap}, gap_length: {gap_length} > 150, considered a two_rna_anno.\n")
                    continue
                read_L_anno = read_L["annotation"]
                read_R_anno = read_R["annotation"]
                priority = ("Promoter", "5UTR", "3UTR", "CDS", "Exon", "Intron", "Downstream", "Intergenic")
                if priority.index(read_L_anno) <= priority.index(read_R_anno):
                    annotation = read_L_anno
                else:
                    annotation = read_R_anno    
                gene_names = read_L_gene_names | read_R_gene_names
                ids = ",".join(list(set(read_L["id"].split(",")) | set(read_R["id"].split(","))))
                types = ",".join(list(set(read_L["type"].split(",")) | set(read_R["type"].split(","))))
                
                fix_gap_rna_anno.append([read_name, read_L_chrom, blocks[0][0], blocks[-1][1], str(blocks), read_L_strand,
                                         annotation, ",".join(list(geneStart), ",".join(list(geneEnd))), ",".join(list(gene_names)), ids, types])
                fix_circle_gap.append([read_name, read_L_chrom, gap, gap_length, read_L_strand, gap_sequence])
            
            else:
                if refseq is not None:
                    ref_sequence = refseq.fetch(read_name)
                    read_name = read_L.name.rstrip("|L")
                    left_start = read_L_blocks[0][0]
                    left_end = read_L_blocks[-1][1]
                    left_block_down = [(left_end, left_end + 1000)]
                    left_block_up = [(left_start - 1000, left_start)]
                    
                    if genome.fetch_blocks(read_L_chrom, left_block_down, read_L_strand).find(
                        ref_sequence) > -1 or genome.fetch_blocks(
                            read_L_chrom, left_block_up, read_L_strand).find(ref_sequence) > -1:
                        sys.stderr.write(f"{read_name}(two rna) ref_sequence in left part up or down 1k.\n")
                        continue
                    
                    right_start = read_R_blocks[0][0]
                    right_end = read_R_blocks[-1][1]
                    right_block_down = [(right_end, right_end + 1000)]
                    right_block_up = [(right_start - 1000, right_start)]
                    
                    if genome.fetch_blocks(read_R_chrom, right_block_up, read_R_strand).find(
                        ref_sequence) > -1 or genome.fetch_blocks(
                            read_R_chrom, right_block_down, read_R_strand).find(ref_sequence) > -1:
                        sys.stderr.write(f"{read_name}(two rna) ref_sequence in right part up or down 1k.\n")
                        continue
                    
                two_rna_anno.append([self.pd_series_to_list(read_L), self.pd_series_to_list(read_R)])
                
            self.two_rna_anno = two_rna_anno
            self.fix_gap_rna_anno = fix_gap_rna_anno
            self.fix_circle_gap = fix_circle_gap
            
    @staticmethod
    def get_fasta_from_anno(genome: GenomeFile, anno: list):
        name = anno[0]
        chrom = anno[1]
        blocks = eval(anno[4])
        strand = anno[5]
        sequence = genome.fetch_blocks(chrom, blocks, strand)    
        return f">{name}\n{sequence}\n"
    
    def writeout_reads_anno(self, anno_list: list, genome: GenomeFile, annofile: str, fasta: str, two_rna_anno=False):
        header = "seqname\tchrom\tstart\tend\tblocks\tstrand\tannotation\tgeneStart\tgeneEnd\tgeneName\tid\ttype\n"
        
        with open(annofile, "w") as annof, open(fasta, "w") as fa:
            annof.write(header)
            if two_rna_anno:
                for read_L, read_R in anno_list:
                    for anno in (read_L, read_R):
                        annof.write("\t".join(map(str, anno))+"\n")
                        fa.write(self.get_fasta_from_anno(genome, anno))
            else:
                for anno in anno_list:
                    annof.write("\t".join(map(str, anno))+"\n")
                    fa.write(self.get_fasta_from_anno(genome, anno))


    @staticmethod
    def writeout_gap_info(gap_info: list, gap_file: str):
        header = "seqname\tchrom\tcircle_gap\tcircle_gap_length\tstrand\tgap_sequence\n"
        with open(gap_file, "w") as gap_fileobj:
            gap_fileobj.write(header)
            for gap in gap_info:
                gap_fileobj.write("\t".join(map(str, gap))+"\n")



def main(genomefa: str, anno:str, fix_gap:str, fix_gap_fa: str, circle_gap: str,
         two: str, two_fa:str, one_side:str, one_side_fa:str,
         ref_fasta: Optional[str] = None, check_up_and_down: Optional[int] = 1000):
    logging.info("start fix gap annotaion ...")
    start = time.perf_counter()
    genome = GenomeFile(genomefa)
    refseq = FastaFile(ref_fasta) if ref_fasta is not None else None
    read_anno = ReadAnno(anno)
    logging.info("start to fix gap ...")
    read_anno.fix_gap(genome, refseq, check_up_and_down)
    logging.info("start to write gap info ...")
    read_anno.writeout_gap_info(
        read_anno.fix_circle_gap, circle_gap)
    logging.info("start to write reads anno ...")
    read_anno.writeout_reads_anno(
        read_anno.fix_gap_rna_anno,
        genome, fix_gap, fix_gap_fa)
    read_anno.writeout_reads_anno(
        read_anno.one_side_rna_anno,
        genome, one_side, one_side_fa)
    read_anno.writeout_reads_anno(
        read_anno.two_rna_anno,
        genome, two, two_fa, two_rna_anno=True)
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
    parser.add_argument("-o", dest="outdir", type=str,
                        default=os.getcwd(), help="output dir")
    parser.add_argument("--reference", dest="reference", type=str, default=None, help="refernce sequence")
    parser.add_argument("--check_up_and_down", dest="check_up_and_down", type=int, default=None, help="check_upstream and downstream")
    
    args = parser.parse_args()
    outdir = args.outdir
    main(args.genome, args.anno, os.path.join(outdir, "fix_gap_anno.tsv"), os.path.join(outdir, "fix_gap.fa"),
         os.path.join(outdir, "circle_gap_info.tsv"), os.path.join(outdir, "two_anno.tsv"), os.path.join(outdir, "two.fa"),
         os.path.join(outdir, "one_side_anno.tsv"), os.path.join(outdir, "one_side.fa"), args.reference, args.check_up_and_down)


if __name__ == "__main__":
    run()
