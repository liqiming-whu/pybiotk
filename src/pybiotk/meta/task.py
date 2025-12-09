#!/usr/bin/env python3
import os
import warnings
import numpy as np
import pandas as pd
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from typing import Optional, Dict, List, Literal, Callable, Tuple, Protocol
from scipy.stats import ttest_ind, levene
from stream import to_list, mapwith, concat, for_each, count
from pybiotk.annodb import MergedTranscript
from pybiotk.io import Openbw


class FuncProc(Protocol):
    def __call__(self, bw: Openbw, t: MergedTranscript, **kwargs) -> Tuple[float, str]:...


def single_task(func: FuncProc, gene: Dict[str, List[MergedTranscript]], bwfile: str, strand: Literal['+', '-'] = '+', **kwargs):
    fwd = gene[strand]
    with Openbw(bwfile) as bw:
        p_func = partial(func, bw, **kwargs)
        a = fwd | mapwith(p_func) | to_list
    return a


def bw_task(func: Callable, gene: Dict[str, List[MergedTranscript]],
            table_file: str, bw_fwd: List[str],
            bw_rev: Optional[List[str]] = None,
            njobs=20, **kwargs):
    with ProcessPoolExecutor(max_workers=njobs) as pool:
        bw_rev = bw_rev if bw_rev is not None else bw_fwd
        fwd_task = [pool.submit(single_task, func, gene, fwd, '+', **kwargs) for fwd in bw_fwd]
        rev_task = [pool.submit(single_task, func, gene, rev, '-', **kwargs) for rev in bw_rev]
        fwd_task_out = [x.result() for x in fwd_task]
        rev_task_out = [x.result() for x in rev_task]
        fwd_array = np.array([j[0] for i in fwd_task_out for j in i]).reshape((len(bw_fwd), -1))
        rev_array = np.array([j[0] for i in rev_task_out for j in i]).reshape((len(bw_rev), -1))
        fwd_loci = [x[1] for x in fwd_task_out[0]]
        fwd_loci.extend(x[1] for x in rev_task_out[0])
        df = pd.DataFrame(np.concatenate([fwd_array, rev_array], axis=1).transpose())
        names = gene['+'] | concat(gene['-']) | mapwith(lambda x: x.gene_name) | to_list
        cols = [os.path.basename(x).split(".")[0] for x in bw_fwd]
        df.columns = cols
        df['gene_name'] = names
        df['gene_loci'] = fwd_loci
        newcols = ['gene_name', 'gene_loci']
        newcols.extend(cols)
        df = df[newcols]
        df = df.dropna(axis=0)
        df = df.loc[df.iloc[:, 2:].sum(axis=1) > 0]
        if len(bw_fwd) >= 4 and len(bw_fwd) % 2 == 0:
            def t_test(x):
                a = [i+0.01*sum(x[2 + len(bw_fwd)//2:]) for i in x[2: 2 + len(bw_fwd)//2]]
                b = [i+0.01*sum(x[2: 2 + len(bw_fwd)//2]) for i in x[2 + len(bw_fwd)//2:]]
                return ttest_ind(a, b, equal_var=levene(a,b).pvalue > 0.05).pvalue

            def fc(x):
                a = sum(x[2: 2 + len(bw_fwd)//2]) + 0.01*sum(x[2 + len(bw_fwd)//2:])
                b = sum(x[2 + len(bw_fwd)//2:]) + 0.01*sum(x[2: 2 + len(bw_fwd)//2])
                return b / a

            def diff(x):
                a = sum(x[2: 2 + len(bw_fwd)//2]) / (len(bw_fwd)//2)
                b = sum(x[2 + len(bw_fwd)//2:]) / (len(bw_fwd)//2)
                return b - a
            df['difference'] = df.apply(diff, axis=1)
            df['foldchange'] = df.apply(fc, axis=1)
            df['pvalue'] = df.apply(t_test, axis=1)
            df = df.sort_values(by=["difference", "pvalue"], ascending=[False, True])
        elif len(bw_fwd) == 2:
            def fc(x):
                a = x[2] + 0.01*x[3]
                b = x[3] + 0.01*x[2]
                return b / a

            def diff(x):
                return x[3] - x[2]
            df['difference'] = df.apply(diff, axis=1)
            df['foldchange'] = df.apply(fc, axis=1)
            df = df.sort_values(by="difference", ascending=False)
        df.to_csv(table_file, sep='\t', float_format='%.4f', index=False)


def scale_regions(bw: Openbw, metaplot_values: np.ndarray, t: MergedTranscript,
                  upStream: int = 2000, downStream: int = 2000,
                  length: int = 5000, bins: int = 10):
    # TODO use index to calculate, not to use concatenate
    if t.strand == '-':
        upStream, downStream = downStream, upStream

    start = t.start - upStream
    if start < 0:
        start = 0
    end = t.end + downStream
    try:
        values = bw.values(t.chrom, start, end)
    except RuntimeError:
        warnings.warn(f'Invalid interval bounds! in {t.chrom}:{start}-{end}, skip it')
        values = np.zeros(metaplot_values.shape[-1])
    medion_values = bw.scale_region_values_np_values(values[upStream: -downStream], length, bins)
    start_values = bw.scale_region_values_np_values(values[:upStream], upStream, bins)
    end_values = bw.scale_region_values_np_values(values[-downStream:], downStream, bins)
    z = np.concatenate([start_values, medion_values, end_values])
    metaplot_values += z


def scale_regions_single_task(gene: Dict[str, List[MergedTranscript]],
                              bw_file: str, strand: Literal['+', '-'],
                              upStream: int = 2000, downStream: int = 2000,
                              length: int = 5000, bins: int = 10):
    assert length % bins == 0
    assert upStream % bins == 0
    assert downStream % bins == 0
    values = np.zeros((upStream+downStream+length)//bins)
    with Openbw(bw_file) as bw:
        scale_regions_func = partial(scale_regions, bw, values, upStream=upStream, downStream=downStream, length=length, bins=bins)
        gene[strand] | for_each(scale_regions_func)
    return values


def scale_regions_task(numpy_file: str, gene: Dict[str, List[MergedTranscript]],
                       bw_fwd: List[str], bw_rev: Optional[List[str]] = None,
                       upStream: int = 2000, downStream: int = 2000,
                       length: int = 5000, bins: int = 10, njobs: int = 20):

    with ProcessPoolExecutor(max_workers=njobs) as pool:
        bw_rev = bw_rev if bw_rev is not None else bw_fwd
        fwd_task_list = [pool.submit(
            scale_regions_single_task, gene, bw, '+', upStream=upStream,
            downStream=downStream, length=length, bins=bins) for bw in bw_fwd]
        rev_task_list = [pool.submit(
            scale_regions_single_task, gene, bw, '-', upStream=upStream,
            downStream=downStream, length=length, bins=bins) for bw in bw_rev]

        fwd_task_out = np.array([x.result() for x in fwd_task_list])
        rev_task_out = np.flip(np.array([x.result() for x in rev_task_list]), axis=1)
        merged_strand_out = fwd_task_out + rev_task_out
        n = (gene['+'] | count) + (gene['-'] | count)
        with open(numpy_file, 'wb') as f:
            np.save(f, merged_strand_out/n)


def reference_point(bw: Openbw, a: np.ndarray, t: MergedTranscript,
                    loci: Literal['TES', 'TSS'] = 'TES',
                    method: Literal['reads', 'coverage'] = 'reads',
                    upStream: int = 2000, downStream: int = 10000, bins: int = 10):
    if t.strand == '-':
        upStream, downStream = downStream, upStream
    point = None
    if loci == 'TES':
        point = t.end if t.strand == '+' else t.start
    else:
        point = t.start if t.strand == '+' else t.end
    start = point - upStream
    end = point + downStream
    if start < 0:
        start = 0
    values = bw.scale_region_values(t.chrom, start, end, upStream+downStream, bins) if method == 'reads' else bw.coverage_sliding_window(t.chrom, start, end, nbins=bins)
    if values.shape != a.shape:
        values = np.concatenate([np.zeros(a.shape[-1]-values.shape[-1]), values])
    a += values


def reference_point_single_task(gene: Dict[str, List[MergedTranscript]],
                                bw_file: str, strand: Literal['+', '-'],
                                loci: Literal['TES', 'TSS'] = 'TES',
                                method: Literal['reads', 'coverage'] = 'reads',
                                upStream: int = 2000,
                                downStream: int = 10000,
                                bins: int = 10):
    values = np.zeros((upStream+downStream)//bins) if method == 'reads' else np.zeros(upStream+downStream-bins+1)
    with Openbw(bw_file) as bw:
        reference_point_func = partial(reference_point, bw, values, loci=loci, method=method, upStream=upStream, downStream=downStream, bins=bins)
        gene[strand] | for_each(reference_point_func)
        return values


def reference_point_task(numpy_file: str, gene: Dict[str, List[MergedTranscript]],
                         bw_fwd: List[str], bw_rev: Optional[List[str]] = None,
                         loci: Literal['TES', 'TSS'] = 'TES',
                         method: Literal['reads', 'coverage'] = 'reads',
                         upStream: int = 2000, downStream: int = 10000,
                         bins: int = 10, njobs: int = 20):
    with ProcessPoolExecutor(max_workers=njobs) as pool:
        bw_rev = bw_rev if bw_rev is not None else bw_fwd
        fwd_task_list = [pool.submit(
            reference_point_single_task, gene, bw, '+', loci=loci,
            method=method, upStream=upStream, downStream=downStream,
            bins=bins) for bw in bw_fwd]
        rev_task_list = [pool.submit(
            reference_point_single_task, gene, bw, '-', loci=loci,
            method=method, upStream=upStream, downStream=downStream,
            bins=bins) for bw in bw_rev]
        fwd_task_out = np.array([x.result() for x in fwd_task_list])
        rev_task_out = np.flip(np.array([x.result() for x in rev_task_list]), axis=1)
        merged_strand_out = fwd_task_out + rev_task_out
        n = (gene['+'] | count) + (gene['-'] | count)
        with open(numpy_file, 'wb') as f:
            np.save(f, merged_strand_out/n)


def signal_point(bw: Openbw, a: np.ndarray, t: MergedTranscript,
                 loci: Literal['TES', 'TSS'] = 'TES',
                 upStream: int = 2000, downStream: int = 2000):
    if t.strand == '-':
        upStream, downStream = downStream, upStream
    point = None
    if loci == 'TES':
        point = t.end if t.strand == '+' else t.start
    else:
        point = t.start if t.strand == '+' else t.end
    start = point - upStream
    end = point + downStream
    if start < 0:
        start = 0
    values = bw.values(t.chrom, start, end)
    if values.shape != a.shape:
        values = np.concatenate([np.zeros(a.shape[-1]-values.shape[-1]), values])
    a += values


def signal_point_single_task(gene: Dict[str, List[MergedTranscript]],
                             bw_file: str, strand: Literal['+', '-'],
                             loci: Literal['TES', 'TSS'] = 'TSS',
                             upStream: int = 2000, downStream: int = 2000):
    values = np.zeros(upStream+downStream)
    with Openbw(bw_file) as bw:
        reference_point_func = partial(signal_point, bw, values, loci=loci, upStream=upStream, downStream=downStream)
        gene[strand] | for_each(reference_point_func)
        return values


def signal_point_task(numpy_file: str, gene: Dict[str, List[MergedTranscript]],
                      bw_fwd: str, bw_rev: str = None,
                      loci: Literal['TES', 'TSS'] = 'TSS',
                      upStream: int = 2000, downStream: int = 2000):
    with ProcessPoolExecutor(max_workers=4) as pool:
        fwd_task_list = [pool.submit(signal_point_single_task, gene, bw, '+', loci=loci, upStream=upStream, downStream=downStream) for bw in (bw_fwd, bw_rev)]
        rev_task_list = [pool.submit(signal_point_single_task, gene, bw, '-', loci=loci, upStream=upStream, downStream=downStream) for bw in (bw_rev, bw_fwd)]
    fwd_task_out = [x.result() for x in fwd_task_list]
    rev_task_out = [x.result() for x in rev_task_list]
    fwd, fwd_opposite = fwd_task_out
    rev, rev_opposite = rev_task_out
    signal = fwd + np.flip(rev)
    signal_opposite = fwd_opposite + np.flip(rev_opposite)    
    n = (gene['+'] | count) + (gene['-'] | count)
    merged_out = np.array([signal, signal_opposite])
    with open(numpy_file, 'wb') as f:
        np.save(f, merged_out/n)
