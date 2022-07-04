#!/usr/bin/env python3
import os
import argparse
import pickle
import time
import numpy as np
import pandas as pd
from typing import List, Optional
from stream import for_each, filter, apply, to_list
from pybiotk.annodb import MergedTranscript
from pybiotk.intervals import GRangeTree
from pybiotk.io import Openbed
from pybiotk.utils import logging, bedtools_sort
from groan.merge_transcript import load_gene


class PeaksTrees:
    def __init__(self, filenamelist: List[str]):
        self.filenamelist = filenamelist
        self.trees = [self.load_peak_tree(filename) for filename in filenamelist]

    def __len__(self):
        return len(self.filenamelist)

    @staticmethod
    def load_peak_tree(filename: str):
        tree = GRangeTree()
        with Openbed(filename) as bedfile:
            for bed in bedfile:
                tree.add(bed, bed.chrom, bed.start, bed.end)
        logging.info(f"{filename} loaded.")
        return tree

    def tes_peak(self, gene: MergedTranscript, distance: int = 1000):
        dislist = []
        for tree in self.trees:
            if gene.strand == '+':
                peaks = tree.find(gene.chrom, gene.end-distance, gene.end+distance)
                if peaks:
                    dis = max(i.end for i in peaks)-gene.end+distance
                    assert dis > 0
                    dislist.append(dis)
                else:
                    dislist.append(0)
            else:
                peaks = tree.find(gene.chrom, gene.start-distance, gene.start+distance)
                if peaks:
                    dis = gene.start+distance-min(i.start for i in peaks)
                    assert dis > 0
                    dislist.append(dis)
                else:
                    dislist.append(0)
        return np.array(dislist)

    def tss_peak(self, gene: MergedTranscript, distance: int = 1000):
        dislist = []
        for tree in self.trees:
            if gene.strand == '+':
                peaks = tree.find(gene.chrom, gene.start-distance, gene.start+distance)
                if peaks:
                    dis = gene.start+distance-min(i.start for i in peaks)
                    assert dis > 0
                    dislist.append(dis)
                else:
                    dislist.append(0)
            else:
                peaks = tree.find(gene.chrom, gene.end-distance, gene.end+distance)
                if peaks:
                    dis = max(i.end for i in peaks)-gene.end+distance
                    dislist.append(0)
                else:
                    dislist.append(0)
        return np.array(dislist)

    def downStream_peak(self, gene: MergedTranscript, downStream: int = 10000):
        if gene.strand == '+':
            start = gene.end
            if downStream > gene.after:
                downStream = gene.after*0.75
            end = start + downStream
        else:
            end = gene.start
            if downStream > gene.before:
                downStream = int(gene.before*0.75)
            start = end - downStream
            if start < 0:
                logging.warning(f"{gene.gene_name} {gene.chrom}:{start}-{end} start < 0; use 0 instead")
                start = 0
        if start > end:
            logging.warning(f"{gene.gene_name} {gene.chrom}:{start}-{end} start > end, sum(dislist)=0.")
            dislist = np.array([0]*len(self))
        else:
            dislist = []
            for tree in self.trees:
                peaks = tree.find(gene.chrom, start, end)
                if peaks:
                    if gene.strand == '+':
                        dis = max(i.end for i in peaks)-start
                        assert dis > 0
                        dislist.append(dis)
                    else:
                        dis = end-min(i.start for i in peaks)
                        assert dis > 0
                        dislist.append(dis)
                else:
                    dislist.append(0)
        return dislist


def peak_extend(gene: MergedTranscript,
                filenamelist: List[str],
                ped_tes_table: str,
                peak_in_gene_table: str = os.devnull,
                tes_distance: int = 1000,
                tss_distance: int = 1000):
    peakstree = PeaksTrees(filenamelist)
    tes_dislist = []
    gene_info = []
    logging.info("start to calculate peak in genes ...")
    start = time.perf_counter()

    gene_with_peak = []

    def cal_dislist(gene):
        tes = peakstree.tes_peak(gene, tes_distance)
        tss = peakstree.tss_peak(gene, tss_distance)
        withtss = sum(tss) > 0
        withtes = sum(tes) > 0

        gene_with_peak.append([f"{gene.gene_name}", f"{gene.chrom}:{gene.start}-{gene.end}({gene.strand})", f"{withtss}", f"{withtes}"])
        tes_dislist.append(tes)
        gene_info.append((f"{gene.gene_name}", f"{gene.chrom}:{gene.start}-{gene.end}({gene.strand})",
                          f"{gene.length()}", f"{gene.downStream()}",
                          f"{gene.gene_type}", f"{withtss}"))

    gene['+'] | for_each(cal_dislist)
    gene['-'] | for_each(cal_dislist)

    tes_df = pd.DataFrame(tes_dislist)

    g = list(zip(*gene_info))
    cols = [os.path.basename(x).split(".")[0] for x in filenamelist]
    tes_df.columns = cols
    tes_df['gene_name'] = g[0]
    tes_df['gene_loci'] = g[1]
    tes_df['gene_length'] = g[2]
    tes_df['downStream'] = g[3]
    tes_df['gene_type'] = g[4]
    tes_df['peak_in_tss'] = g[5]
    newcols = ['gene_name', 'gene_loci', 'gene_length', 'downStream', 'gene_type', 'peak_in_tss']
    newcols.extend(cols)
    tes_df.to_csv(ped_tes_table, sep='\t', index=False, columns=newcols)
    gene_df = pd.DataFrame(gene_with_peak, columns=["gene_name", "gene_loci", "peak_in_tss", "peak_in_tes"])
    gene_df.to_csv(peak_in_gene_table, sep="\t", index=False)
    end = time.perf_counter()
    logging.info(f"task finished in {end-start:.2f}s.")


def filter_gene_by_peak(gene: MergedTranscript,
                        gene_file: str,
                        table_file: str,
                        filtered_bed: Optional[str] = None,
                        filtered_df: Optional[str] = None,
                        dump: Optional[str] = None,
                        withtss: Optional[str] = None,
                        withtes: Optional[str] = None):
    logging.info("start to filter gene by peak ...")
    df = pd.read_table(gene_file, header=0)
    if withtes is not None:
        df = df.loc[df["peak_in_tes"] == eval(withtes)]
    if withtss is not None:
        df = df.loc[df["peak_in_tss"] == eval(withtss)]

    bedpath = os.devnull if filtered_bed is None else filtered_bed
    select = set(df["gene_loci"].to_list())
    with open(bedpath, 'w') as f:
        a = gene['+'] | filter(lambda x: f"{x.chrom}:{x.start}-{x.end}({x.strand})" in select) | apply(lambda x: f.write(f"{x.to_bed6()}\n")) | to_list
        b = gene['-'] | filter(lambda x: f"{x.chrom}:{x.start}-{x.end}({x.strand})" in select) | apply(lambda x: f.write(f"{x.to_bed6()}\n")) | to_list
        filtered_gene = {'+': a, '-': b}
    if filtered_bed is not None:
        bedtools_sort(filtered_bed)
    if dump is not None:
        with open(dump, "wb") as f:
            pickle.dump(filtered_gene, f)
    table_df = pd.read_table(table_file, header=0)
    table_df = table_df.loc[table_df["gene_loci"].isin(select)]
    if filtered_df is not None:
        table_df.to_csv(filtered_df, sep="\t", index=False)
    logging.info(f"save peak_in_tes: {withtes} peak_in_tss: {withtss} {len(a) + len(b)} genes.")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--load", dest="load_pickle", type=str, required=True, help="load pickle file.")
    parser.add_argument("-f", "--files", dest="files", nargs="+", required=True, type=str,
                        help="input peaks file.")
    parser.add_argument("-t", dest="table_file", type=str, required=True, help="output peak extend table.")
    parser.add_argument("-g", dest="gene_with_peak", type=str, required=True, help="peak in gene table.")
    parser.add_argument("--tes_distance", dest="tes_distance", type=int, default=1000, help="distance form tes.")
    parser.add_argument("--tss_distance", dest="tss_distance", type=int, default=1000, help="distance form tss.")
    parser.add_argument("-b", dest="bed", type=str, default=None, help="output filtered bed file."),
    parser.add_argument("-p", dest="pickle", type=str, default=None, help="output pickle dump file.")
    parser.add_argument("--filtered_table", dest="filtered_table", type=str, default=os.devnull, help="output filtered table.")
    parser.add_argument("--peak_in_tss", dest="peak_in_tss", type=str, choices=["True", "False"], default=None, help="peak in tss")
    parser.add_argument("--peak_in_tes", dest="peak_in_tes", type=str, choices=["True", "False"], default=None, help="peak in tes")    
    args = parser.parse_args()
    if not (os.path.exists(args.table_file) and os.path.exists(args.gene_with_peak)):
        peak_extend(load_gene(args.load_pickle), args.files, args.table_file, args.gene_with_peak, args.tes_distance, args.tss_distance)
    filter_gene_by_peak(load_gene(args.load_pickle), args.gene_with_peak, args.table_file, args.bed, args.filtered_table, args.pickle, args.peak_in_tss, args.peak_in_tes)


if __name__ == '__main__':
    run()
