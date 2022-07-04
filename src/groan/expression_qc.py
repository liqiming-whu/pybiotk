#!/usr/bin/env python3
import os
import argparse
import pickle
import time
import numpy as np
import pandas as pd
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from typing import Optional, Dict, List, Literal, Union
from matplotlib import pyplot as plt
from stream import to_list, filter, apply, mapwith, concat
from pybiotk.annodb import MergedTranscript
from pybiotk.io import Openbw
from pybiotk.utils import logging, bedtools_sort
from groan.merge_transcript import load_gene


def gene_stats(bw: Openbw, t: MergedTranscript,
               loci: Literal["body", "tss", "tes"],
               method: Literal["coverage", "value"],
               distance: int = 3000) -> float:
    if loci == "tss":
        if t.strand == '+':
            start = t.start - distance
            end = t.start + distance
        else:
            start = t.end - distance
            end = t.end + distance
    elif loci == "tes":
        if t.strand == '-':
            start = t.start - distance
            end = t.start + distance
        else:
            start = t.end - distance
            end = t.end + distance
    else:
        start = t.start
        end = t.end
    if method == "coverage":
        a = bw.coverage(t.chrom, start, end)
    else:
        a = bw.stats(t.chrom, start, end)
    return a


def gene_task(gene: Dict[str, List[MergedTranscript]], bwfile: str,
              strand: Literal['+', '-'] = '+',
              loci: Literal["body", "tss", "tes"] = "body",
              method: Literal['coverage', 'reads'] = 'coverage',
              distance: int = 3000):
    fwd = gene[strand]
    with Openbw(bwfile) as bw:
        gene_coverage_or_reads = partial(gene_stats, bw, loci=loci, method=method, distance=distance)
        a = fwd | mapwith(gene_coverage_or_reads) | to_list
    return a


def gene_expression(gene: Dict[str, List[MergedTranscript]],
                    table_file: str, bw_fwd: List[str], bw_rev: Optional[List[str]] = None,
                    loci: Literal["body", "tss", "tes"] = "body",
                    method: Literal['coverage', 'reads'] = 'coverage',
                    distance: int = 3000, njobs=20):
    logging.info(f"start to calculate gene {method} in {loci} +- {distance}...")
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=njobs) as pool:
        bw_rev = bw_rev if bw_rev is not None else bw_fwd
        fwd_task = [pool.submit(gene_task, gene, fwd, '+', loci, method, distance) for fwd in bw_fwd]
        rev_task = [pool.submit(gene_task, gene, rev, '-', loci, method, distance) for rev in bw_rev]
        fwd_task_out: List[List[float]] = [x.result() for x in fwd_task]
        rev_task_out: List[List[float]] = [x.result() for x in rev_task]
        for x in range(len(fwd_task_out)):
            fwd_task_out[x].extend(rev_task_out[x])

        a = np.array(fwd_task_out).transpose()
        df = pd.DataFrame(a)

        b = gene['+'] | concat(gene['-']) | mapwith(
            lambda x: (f"{x.gene_name}", f"{x.chrom}:{x.start}-{x.end}({x.strand})",
                       f"{x.length()}", f"{x.downStream()}",
                       f"{x.gene_type}")) | to_list
        b = list(zip(*b))
        cols = [os.path.basename(x).split(".")[0] for x in bw_fwd]
        df.columns = cols
        df['gene_name'] = b[0]
        df['gene_loci'] = b[1]
        df['gene_length'] = b[2]
        df['downStream'] = b[3]
        df['gene_type'] = b[4]
        newcols = ['gene_name', 'gene_loci', 'gene_length', 'downStream', 'gene_type']
        newcols.extend(cols)
        df.to_csv(table_file, sep='\t', float_format='%.4f', index=False, columns=newcols)
        end = time.perf_counter()
        logging.info(f"task finished in {end-start:.2f}s.")


def filter_gene_expression(gene: Dict[str, List[MergedTranscript]],
                           table_file: str, cols: Optional[List[int]] = None,
                           filtered_bed: Optional[str] = None,
                           dump: Optional[str] = None, value: Union[int, float] = 0.25,
                           method: Literal['coverage', 'reads'] = 'coverage',
                           filter_method: Literal["value", "top", "mean", "quantile"] = "mean"):
    logging.info(f"start to filter gene by {method}...")
    start = time.perf_counter()
    df = pd.read_table(table_file)
    cols = [x+5 for x in cols] if cols is not None else range(df.shape[1])[5:]
    if filter_method == "top":
        df['sum'] = df.iloc[:, cols].sum(axis=1)
        select = set(df.sort_values(by="sum", ascending=False)['gene_loci'][:int(value)])
    elif filter_method == "mean":
        mean_value = df.iloc[:, cols].sum(axis=1).mean()
        select = set(df['gene_loci'][(df.iloc[:, cols].sum(axis=1) >= mean_value)])
    elif filter_method == "quantile":
        assert value <= 1
        qcut = df.iloc[:, cols].sum(axis=1).quantile(value)
        select = set(df['gene_loci'][(df.iloc[:, cols].sum(axis=1) >= qcut)])
    else:
        select = set(df['gene_loci'][(df.iloc[:, cols] >= value).all(axis=1)])

    bedpath = os.devnull if filtered_bed is None else filtered_bed
    with open(bedpath, 'w') as f:
        a = gene['+'] | filter(lambda x: f"{x.chrom}:{x.start}-{x.end}({x.strand})" in select) | apply(lambda x: f.write(f"{x.to_bed6()}\n")) | to_list
        b = gene['-'] | filter(lambda x: f"{x.chrom}:{x.start}-{x.end}({x.strand})" in select) | apply(lambda x: f.write(f"{x.to_bed6()}\n")) | to_list
        filtered_gene = {'+': a, '-': b}
    if filtered_bed is not None:
        bedtools_sort(filtered_bed)
    if dump is not None:
        with open(dump, "wb") as f:
            pickle.dump(filtered_gene, f)
    end = time.perf_counter()
    if filter_method == "top":
        logging.info(f"save {method} top {len(a) + len(b)} genes.")
    elif filter_method == "mean":
        logging.info(f"save {len(a) + len(b)} > mean {method}({mean_value}) genes.")
    elif filter_method == "quantile":
        logging.info(f"save {len(a) + len(b)} {method} > {qcut}({value*100:.2f}%) genes.")
    else:
        logging.info(f"save {len(a) + len(b)} {method} > {value} genes.")
    logging.info(f"task finished in {end-start:.2f}s.")
    return filtered_gene


def plot_gene_expression(table_file: str, plot_file: str, figsize=(10, 6), xlabel: str = 'coverage', title: str = ''):
    logging.info(f"start to plot gene {xlabel}...")
    start = time.perf_counter()
    df = pd.read_table(table_file)
    fig, ax = plt.subplots(figsize=figsize)
    df.iloc[:, 5:].plot.kde(ax=ax)
    ax.set_xlabel(xlabel)
    ax.set_title(f"n={df.shape[0]}    "+title)
    fig.tight_layout()
    plt.savefig(plot_file)
    end = time.perf_counter()
    logging.info(f"task finished in {end-start:.2f}s.")


def expression_qc(gene: Dict[str, List[MergedTranscript]],
                  table_file: str,
                  bw_fwd: List[str],
                  bw_rev: Optional[List[str]] = None,
                  dump: Optional[str] = None,
                  filtered_bed: Optional[str] = None,
                  cols: Optional[List[int]] = None,
                  loci: Literal["body", "tss", "tes"] = "body",
                  method: Literal['coverage', 'reads'] = 'coverage',
                  distance: int = 3000,
                  value: Union[float, int] = 0.7,
                  filter_method: Literal["value", "top", "mean", "quantile"] = "mean",
                  plot_file: Optional[str] = None,
                  title: str = ''):
    if not os.path.exists(table_file):
        gene_expression(gene, table_file, bw_fwd, bw_rev, loci, method, distance)
    else:
        logging.info(f"using exist expression table: {table_file}.")
    filtered_gene = filter_gene_expression(gene, table_file, cols, filtered_bed, dump, value, method, filter_method)
    if plot_file is not None:
        plot_gene_expression(table_file, plot_file, xlabel=method, title=title)
    return filtered_gene


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--load", dest="load_pickle", type=str, required=True, help="load pickle file.")
    parser.add_argument("-f", "--fwd", dest="fwd", nargs="+", required=True, type=str,
                        help="input forward bigwigs or collapsed bigwigs.")
    parser.add_argument("-r", "--rev", dest="rev", nargs="+", default=None, type=str,
                        help="input reverse bigwigs or None.")
    parser.add_argument("-b", dest="bed", type=str, default=None, help="output filtered bed file.")
    parser.add_argument("-p", dest="pickle", type=str, default=None, help="output pickle dump file.")
    parser.add_argument("-t", dest="table_file", type=str, required=True, help="output expression table.")
    parser.add_argument("-c", dest="cols", type=int, nargs="+", default=None, help="choose which samples to filter, 0 based.")
    parser.add_argument("-l", "--loci", dest="loci", type=str, default="body", choices=["body", "tss", "tes"], help="calc in which loci")
    parser.add_argument("-d", "--distance", dest="distance", type=int, default=3000, help="tss/tes +- distance.")
    parser.add_argument("-m", dest="method", type=str, default="coverage", choices=["coverage", "reads"], help="choose calculate method.")
    parser.add_argument("-v", dest="value", type=float, default=0.25, help="value cutoff.")
    parser.add_argument("--filter_gene_by", dest="filter_method", type=str, default="mean", choices=["value", "top", "mean", "quantile"], help="filter method")
    parser.add_argument("--plot", dest="plot_file", type=str, default=None, help="output plot file.")
    parser.add_argument("--title", dest="title", type=str, default='', help="output plot file title.")
    args = parser.parse_args()
    expression_qc(load_gene(args.load_pickle), args.table_file,
                  args.fwd, args.rev, args.pickle,
                  args.bed, args.cols, args.loci,
                  args.method, args.distance,
                  args.value, args.filter_method,
                  args.plot_file, args.title
                  )


if __name__ == '__main__':
    run()
