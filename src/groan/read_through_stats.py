#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Literal
from matplotlib import pyplot as plt
from pybiotk.annodb import MergedTranscript
from pybiotk.io import Openbw
from pybiotk.utils import logging
from groan.merge_transcript import load_gene
from groan.task import bw_task


def downStream_stats(bw: Openbw, t: MergedTranscript, downStream: int = 10000, shift: int = 0,
                     method: Literal['coverage', 'reads'] = 'coverage') -> Tuple[float, str]:
    if t.strand == '+':
        start = t.end
        if downStream > t.after:
            downStream = t.after*0.75
        end = start + downStream
        start = start + shift
    else:
        end = t.start
        if downStream > t.before:
            downStream = int(t.before*0.75)
        start = end - downStream
        if start < 0:
            logging.warning(f"{t.gene_name} {t.chrom}:{start}-{end} start < 0; use 0 instead")
            start = 0
        end = end - shift
    if start > end:
        logging.warning(f"{t.gene_name} {t.chrom}:{start}-{end} start > end, value=0.")
        values = 0.0
    else:
        if method == "coverage":
            values = bw.coverage(t.chrom, start, end)
        else:
            values = bw.stats(t.chrom, start, end)
    return values, f"{t.chrom}:{start}-{end}({t.strand})"


def readthrough_stats(gene: Dict[str, List[MergedTranscript]],
                      out_table: str, bw_fwd: List[str], bw_rev: List[str],
                      downStream: int = 10000, shift: int = 0,
                      method: Literal['coverage', 'reads'] = 'coverage'):
    bw_task(downStream_stats, gene, out_table, bw_fwd, bw_rev, downStream=downStream, shift=shift, method=method)


def plot_readthrough_stats(table_file: str, plot_out: str, ylabel='', xlabel='', rot=45, figsize=(10, 6)):
    df = pd.read_table(table_file)
    _, ax = plt.subplots(figsize=figsize)
    ax.boxplot(df.iloc[:, 2:], labels=df.iloc[:, 2:].columns.tolist(),
               showfliers=False, patch_artist=True, flierprops={'color': 'black'},
               boxprops={'color': 'black', 'facecolor': '#68ad81'},
               medianprops={'color': 'black'})
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"n={df.shape[0]}    ")
    plt.xticks(rotation=rot)
    plt.tight_layout()
    plt.savefig(plot_out)


def interpolation_occ(values: np.array):
    last = 1
    interp_occ = 0
    gap_len = 0
    for occ in values:
        if occ == 0 and last == 1:
            gap_len = 1
            last = 0
        elif occ == 0 and last == 0:
            gap_len += 1
        elif occ > 0:
            interp_occ += 1
            if gap_len:
                if gap_len <= 10:
                    interp_occ += gap_len
                else:
                    interp_occ += pow(gap_len, 0.7)
                last = 1
                gap_len = 0
    return interp_occ


def downStream_interpolation_coverage(bw: Openbw, t: MergedTranscript, downStream: int = 10000, shift: int = 0) -> Tuple[float, str]:
    if t.strand == '+':
        start = t.end
        if downStream > t.after:
            downStream = t.after*0.75
        end = start + downStream
        start = start + shift
    else:
        end = t.start
        if downStream > t.before:
            downStream = int(t.before*0.75)
        start = end - downStream
        if start < 0:
            logging.warning(f"{t.gene_name} {t.chrom}:{start}-{end} start < 0; use 0 instead")
            start = 0
        end = end - shift    
    values = bw.values(t.chrom, start, end)
    occ = interpolation_occ(values)
    coverage = occ / values.shape[-1]
    return coverage, f"{t.chrom}:{start}-{end}({t.strand})"


def readthrough_interpolation_coverage(
    gene: Dict[str, List[MergedTranscript]],
    out_table: str, bw_fwd: List[str], bw_rev: List[str],
    downStream: int = 10000, shift: int = 0):
    bw_task(downStream_interpolation_coverage, gene, out_table, bw_fwd, bw_rev, downStream=downStream, shift=shift)


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--load", dest="load_pickle", type=str, required=True, help="load pickle file.")
    parser.add_argument("-f", "--fwd", dest="fwd", nargs="+", required=True, type=str,
                        help="input forward bigwigs or collapsed bigwigs.")
    parser.add_argument("-r", "--rev", dest="rev", nargs="+", default=None, type=str,
                        help="input reverse bigwigs or None.")
    parser.add_argument("-t", dest="table_file", type=str, required=True, help="output expression table.")
    parser.add_argument("-d", "--downStream", dest="downStream", type=int, default=10000, help="downStream distance.")
    parser.add_argument("--downStart", dest="downStart", type=int, default=0, help="downStream distance from TES.")
    parser.add_argument("-m", "--method", dest="method", type=str, default="reads", choices=["reads", "coverage"],
                        help="calculte method.")
    parser.add_argument("--plot", dest="plot_file", type=str, default=None, help="output plot file.")
    args = parser.parse_args()
    readthrough_stats(load_gene(args.load_pickle), args.table_file, args.fwd, args.rev, args.downStream, args.downStart, args.method)
    if args.plot_file:
        plot_readthrough_stats(args.table_file, args.plot_file)


if __name__ == '__main__':
    run()
