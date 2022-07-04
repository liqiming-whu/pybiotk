#!/usr/bin/env python3
import os
import warnings
import numpy as np
import pandas as pd
from functools import partial
from concurrent.futures import ProcessPoolExecutor
from typing import Optional, Dict, List, Literal, Callable, Tuple, Protocol
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
            bw_rev: Optional[List[str]] = None, njobs=20, **kwargs):
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
        df['loci'] = fwd_loci
        newcols = ['gene_name', 'loci']
        newcols.extend(cols)
        df.to_csv(table_file, sep='\t', float_format='%.4f', index=False, columns=newcols)


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
    if t.strand == '+':
        if downStream > t.after:
            downStream = t.after
    else:
        if upStream > t.before:
            upStream = t.before

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
    values = bw.values(t.chrom, start, end) if method == 'reads' else bw.coverage_sliding_windw(t.chrom, start, end, nbins=bins)
    if values.shape != a.shape:
        if t.strand == '+':
            values = np.concatenate([values, np.zeros(a.shape[-1]-values.shape[-1])])
        else:
            values = np.concatenate([np.zeros(a.shape[-1]-values.shape[-1]), values])
    a += values


def reference_point_single_task(gene: Dict[str, List[MergedTranscript]],
                                bw_file: str, strand: Literal['+', '-'],
                                a: np.ndarray, **kwargs):
    with Openbw(bw_file) as bw:
        reference_point_func = partial(reference_point, bw, a, **kwargs)
        gene[strand] | for_each(reference_point_func)
        return a


def reference_point_task(numpy_file: str, gene: Dict[str, List[MergedTranscript]],
                         bw_fwd: List[str], bw_rev: Optional[List[str]] = None,
                         loci: Literal['TES', 'TSS'] = 'TES',
                         method: Literal['reads', 'coverage'] = 'reads',
                         upStream: int = 2000, downStream: int = 10000,
                         bins: int = 10, njobs: int = 20):
    values = np.zeros(upStream+downStream) if method == 'reads' else np.zeros(upStream+downStream-bins+1)
    with ProcessPoolExecutor(max_workers=njobs) as pool:
        bw_rev = bw_rev if bw_rev is not None else bw_fwd
        fwd_task_list = [pool.submit(
            reference_point_single_task, gene, bw, '+', values, loci=loci,
            method=method, upStream=upStream, downStream=downStream,
            bins=bins) for bw in bw_fwd]
        rev_task_list = [pool.submit(
            reference_point_single_task, gene, bw, '-', values, loci=loci,
            method=method, upStream=upStream, downStream=downStream,
            bins=bins) for bw in bw_rev]
        fwd_task_out = np.array([x.result() for x in fwd_task_list])
        rev_task_out = np.flip(np.array([x.result() for x in rev_task_list]), axis=1)
        merged_strand_out = fwd_task_out + rev_task_out
        n = (gene['+'] | count) + (gene['-'] | count)
        with open(numpy_file, 'wb') as f:
            np.save(f, merged_strand_out/n)
