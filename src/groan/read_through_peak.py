#!/usr/bin/env python3
import os
import argparse
import time
import pandas as pd
from typing import List
from stream import for_each
from pybiotk.annodb import MergedTranscript
from pybiotk.utils import logging
from groan.merge_transcript import load_gene
from groan.peak_qc import PeaksTrees


def peak_in_downStream(gene: MergedTranscript, filenamelist: List[str], ped_tes_table: str, downStream: int = 1000):
    peakstree = PeaksTrees(filenamelist)
    tes_dislist = []
    gene_info = []
    logging.info(f"start to calculate peak in downStream {downStream} ...")
    start = time.perf_counter()

    def cal_dislist(gene):
        tes = peakstree.downStream_peak(gene, downStream)
        tes_dislist.append(tes)
        gene_info.append((f"{gene.gene_name}", f"{gene.chrom}:{gene.start}-{gene.end}({gene.strand})",
                          f"{gene.length()}", f"{gene.downStream()}",
                          f"{gene.gene_type}"))

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
    newcols = ['gene_name', 'gene_loci', 'gene_length', 'downStream', 'gene_type']
    newcols.extend(cols)
    tes_df.to_csv(ped_tes_table, sep='\t', index=False, columns=newcols)
    end = time.perf_counter()
    logging.info(f"task finished in {end-start:.2f}s.")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--load", dest="load_pickle", type=str, required=True, help="load pickle file.")
    parser.add_argument("-f", "--files", dest="files", nargs="+", required=True, type=str,
                        help="input peaks file.")
    parser.add_argument("-t", dest="table_file", type=str, required=True, help="output peak extend table.")
    parser.add_argument("--downStream", dest="downStream", type=int, default=10000, help="downStream length.") 
    args = parser.parse_args()
    peak_in_downStream(load_gene(args.load_pickle), args.files, args.table_file, args.downStream)


if __name__ == '__main__':
    run()
