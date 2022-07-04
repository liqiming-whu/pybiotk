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
from stream import to_list, filter, apply, mapwith, concat
from pybiotk.annodb import MergedTranscript
from pybiotk.io import Openbw
from pybiotk.utils import logging, bedtools_sort
from groan.merge_transcript import load_gene


def gap_score(values: np.array):
    last = 1
    score = 0
    gap_len = 0
    for occ in values:
        if occ == 0 and last == 1:
            gap_len = 1
            last = 0
        elif occ == 0 and last == 0:
            gap_len += 1
        elif occ > 0:
            if gap_len:
                score += 2 + np.log(gap_len)
                last = 1
                gap_len = 0
    return score

def downStream_signal_gap_score(bw: Openbw, t: MergedTranscript, downStream: int = 10000, shift: int = 0) -> float:
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
    return gap_score(values)
    

def downStream_signal_gap_score_task(gene: Dict[str, List[MergedTranscript]], bwfile: str,
              strand: Literal['+', '-'] = '+', downStream: int = 10000, shift: int = 0):
    fwd = gene[strand]
    with Openbw(bwfile) as bw:
        gap_score_func = partial(downStream_signal_gap_score, bw, downStream=downStream, shift=shift)
        a = fwd | mapwith(gap_score_func) | to_list
    return a


def gene_downStream_gap_score(gene: Dict[str, List[MergedTranscript]],
                    table_file: str, bw_fwd: List[str], bw_rev: Optional[List[str]] = None,
                    downStream: int = 10000, shift: int = 0, njobs=20):
    logging.info("start to calculate gene downSteam signal gap score")
    start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=njobs) as pool:
        bw_rev = bw_rev if bw_rev is not None else bw_fwd
        fwd_task = [pool.submit(downStream_signal_gap_score_task, gene, fwd, '+', downStream, shift) for fwd in bw_fwd]
        rev_task = [pool.submit(downStream_signal_gap_score_task, gene, rev, '-', downStream, shift) for rev in bw_rev]
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


def filter_gene_gap_score(gene: Dict[str, List[MergedTranscript]],
                          table_file: str, cols: Optional[List[int]] = None,
                          filtered_bed: Optional[str] = None,
                          dump: Optional[str] = None, value: Union[int, float] = 0.75,
                          filter_method: Literal["value", "top", "mean", "quantile"] = "quantile"):
    logging.info("start to filter gene by downStream gap score")
    start = time.perf_counter()
    df = pd.read_table(table_file)
    cols = [x+5 for x in cols] if cols is not None else range(df.shape[1])[5:]
    if filter_method == "top":
        df['sum'] = df.iloc[:, cols].sum(axis=1)
        select = set(df.sort_values(by="sum", ascending=True)['gene_loci'][:int(value)])
    elif filter_method == "mean":
        mean_value = df.iloc[:, cols].sum(axis=1).mean()
        select = set(df['gene_loci'][(df.iloc[:, cols].sum(axis=1) <= mean_value)])
    elif filter_method == "quantile":
        assert value <= 1
        qcut = df.iloc[:, cols].sum(axis=1).quantile(value)
        select = set(df['gene_loci'][(df.iloc[:, cols].sum(axis=1) <= qcut)])
    else:
        select = set(df['gene_loci'][(df.iloc[:, cols] <= value).all(axis=1)])

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
        logging.info(f"save gap score top {len(a) + len(b)} genes.")
    elif filter_method == "mean":
        logging.info(f"save {len(a) + len(b)} < mean gap score({mean_value}) genes.")
    elif filter_method == "quantile":
        logging.info(f"save {len(a) + len(b)} gap score < {qcut}({value*100:.2f}%) genes.")
    else:
        logging.info(f"save {len(a) + len(b)} gap score < {value} genes.")
    logging.info(f"task finished in {end-start:.2f}s.")
    return filtered_gene


def signal_continuity_qc(gene: Dict[str, List[MergedTranscript]],
                         table_file: str,
                         bw_fwd: List[str],
                         bw_rev: Optional[List[str]] = None,
                         dump: Optional[str] = None,
                         filtered_bed: Optional[str] = None,
                         cols: Optional[List[int]] = None,
                         downStream: int = 10000,
                         shift: int = 0, 
                         value: Union[float, int] = 0.7,
                         filter_method: Literal["value", "top", "mean", "quantile"] = "quantile"):
    if not os.path.exists(table_file):
        gene_downStream_gap_score(gene, table_file, bw_fwd, bw_rev, downStream, shift)
    else:
        logging.info(f"using exist downStream gap score table: {table_file}.")
    filtered_gene = filter_gene_gap_score(gene, table_file, cols, filtered_bed, dump, value, filter_method)
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
    parser.add_argument("--downStream", dest="downStream", type=int, default=10000, help="downStream cutoff.")
    parser.add_argument("--downStart", dest="downStart", type=int, default=0, help="downStream distance from TES.")
    parser.add_argument("-v", dest="value", type=float, default=0.75, help="value cutoff.")
    parser.add_argument("--filter_gene_by", dest="filter_method", type=str, default="quantile", choices=["value", "top", "mean", "quantile"], help="filter method")
    args = parser.parse_args()
    signal_continuity_qc(load_gene(args.load_pickle), args.table_file,
                  args.fwd, args.rev, args.pickle,
                  args.bed, args.cols, args.loci,
                  args.downStream, args.downStart,
                  args.value, args.filter_method,
                  )


if __name__ == '__main__':
    run()
