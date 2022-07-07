#!/usr/bin/env python3
import os
import pysam
import argparse
import pandas as pd
from itertools import groupby
from typing import Tuple, List, NamedTuple, Optional
from pybiotk.utils import is_overlap
from pybiotk.intervals import merge_intervals


class Mergeblocks(NamedTuple):
    mergeblocks: Optional[List[Tuple[int, int]]] = None
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

    def classify_reads(self):
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
                one_side_rna_anno.append(self.anno.loc[left_anno])
            elif right_anno in self.anno.index:
                one_side_rna_anno.append(self.anno.loc[right_anno])

        self.one_side_rna_anno = one_side_rna_anno
        return two_rna_anno


    def fix_gap(self, bed_dict, genome, pirna):
        two_rna_anno_set = self.classify_reads()
        two_rna_anno = []
        fix_gap_rna_anno = []
        fix_circle_gap = []
        for read_L, read_R in two_rna_anno_set:
            read_L_chrom = read_L[0][1]
            read_L_blocks = eval(read_L[0][4])
            read_L_strand = read_L[0][5]
            read_L_gene_name = set(i[9] for i in read_L)
            read_R_chrom = read_R[0][1]
            read_R_blocks = eval(read_R[0][4])
            read_R_strand = read_R[0][5]
            read_R_gene_name = set(i[9] for i in read_R)
            if read_L_chrom == read_R_chrom and read_L_strand == read_R_strand and read_L_gene_name == read_R_gene_name:
                blocks, gap, gap_length = merge_neighbor(read_L_blocks, read_R_blocks)

                read_name = read_L[0][0].rstrip("|L")
                if not blocks:
                    # two_rna_anno.append([read_L, read_R])
                    print(f"{read_name} gap > 100nt")
                    continue

                trans_info_set = []
                trans_name_set = set()
                for read in read_L + read_R:
                    if read[6] in trans_name_set:
                        continue
                    trans_name_set.add(read[7])
                    trans_info_set.append([read[7], read[8], read[9], read[10]])

                if gap == "overlap":
                    print(f"{read_name} overlap gap")
                    continue

                gap_sequence = None
                pirna_sequence = pirna.fetch(read_name)
                if gap_length and pirna:
                    gap_sequence = genome.fetch(read_L_chrom, gap, read_L_strand)
                    if gap_sequence.find(pirna_sequence) > -1:
                        print(f"{read_name} pirna in gap\t{gap_sequence}\t{pirna_sequence}")
                        continue
                    
                pirna_in_trans = False
                for trans_info in trans_info_set:
                    transcript_id, transcript_type, gene_name, gene_type = trans_info
                    if transcript_type == "Unknown":
                        continue
                    bed = bed_dict[transcript_id]
                    transcript = Transcript(
                        transcript_id, None, transcript_type, None,
                        gene_name, gene_type, bed[5], bed[0], bed[1], bed[2],
                        bed[6], bed[7], bed[10], bed[11]
                    )
                    transcript_sequence_with_intron = transcript.get_sequence_with_intron(genome)
                    if transcript_sequence_with_intron.find(pirna_sequence) > -1:
                        pirna_in_trans = True
                        break
                    transcript_sequence = transcript.get_sequence(genome)
                    if transcript_sequence.find(pirna_sequence) > -1:
                        pirna_in_trans = True
                        break

                if pirna_in_trans:
                    print(f"{read_name} pirna in gene")
                    continue
                    
                if len(trans_info_set) == 1 and list(trans_info_set)[0][1] == "Unknown":
                    up_1k_blocks = [(blocks[0][0]-1000, blocks[0][0])]
                    down_1k_blocks = [(blocks[-1][1], blocks[-1][1]+1000)]
                    up_sequence = genome.fetch(read_L_chrom, up_1k_blocks, read_L_strand)
                    down_sequence = genome.fetch(read_L_chrom, down_1k_blocks, read_L_strand)
                    if up_sequence.find(pirna_sequence) > -1 or down_sequence.find(pirna_sequence) > -1:
                        print(f"{read_name} pirna in upstream or downstream 1k")
                        continue
                
                read_anno_list = []
                for trans_info in trans_info_set:
                    position = self.anno_read(trans_info, blocks, bed_dict)
                    read_anno_list.append(
                        [read_name, read_L_chrom, blocks[0][0], blocks[-1][1], str(blocks), read_L_strand, position, *trans_info])
                if read_anno_list:
                    fix_gap_rna_anno.append(read_anno_list)
                    fix_circle_gap.append([read_name, read_L_chrom, gap, gap_length, read_L_strand, gap_sequence])
            else:
                if genome and pirna:
                    read_name = read_L[0][0].rstrip("|L")
                    pirna_sequence = pirna.fetch(read_name)
                    left_end = read_L_blocks[-1][1]
                    left_block = [(left_end, left_end + len(pirna_sequence) + 1)]
                    
                    if genome.fetch(read_L_chrom, left_block, read_L_strand).find(pirna_sequence):
                        print(f"{read_name} pirna in left")
                        continue
                    
                    right_start = read_R_blocks[0][0]
                    right_block = [(right_start - len(pirna_sequence) - 1, right_start)]
                    
                    if genome.fetch(read_R_chrom, right_block, read_R_strand).find(pirna_sequence):
                        print(f"{read_name} pirna in right")
                        continue
                    
                two_rna_anno.append([read_L, read_R])

            self.two_rna_anno = two_rna_anno
            self.fix_gap_rna_anno = fix_gap_rna_anno
            self.fix_circle_gap = fix_circle_gap