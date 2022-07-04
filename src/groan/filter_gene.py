#!/usr/bin/env python3
import os
import argparse
import pickle
import time
from typing import Optional, Dict, List
from stream import to_list, filter, apply
from pybiotk.annodb import MergedTranscript
from pybiotk.utils import logging, bedtools_sort
from groan.merge_transcript import load_gene


def filter_gene(merged_transcripts: Dict[str, List[MergedTranscript]],
                filtered_bed: Optional[str] = None,
                dump: Optional[str] = None,
                gene_length: int = 3000,
                downStream: int = 10000,
                gene_types: Optional[List[str]] = None,
                transcript_types: Optional[List[str]] = None):
    bedpath = os.devnull if filtered_bed is None else filtered_bed
    logging.info("start to filter gene ...")
    start = time.perf_counter()
    with open(bedpath, "w") as bed:
        iter_a = merged_transcripts['+'] | filter(lambda x: x.length() > gene_length) | filter(
            lambda x: x.downStream() > downStream)
        iter_b = merged_transcripts['-'] | filter(lambda x: x.length() > gene_length) | filter(
            lambda x: x.downStream() > downStream)
        if gene_types is not None:
            gene_types = set(gene_types)
            iter_a = iter_a | filter(lambda x: set(x.gene_type.split(",")) & gene_types)
            iter_b = iter_b | filter(lambda x: set(x.gene_type.split(",")) & gene_types)
        if transcript_types is not None:
            transcript_types = set(transcript_types)
            iter_a = iter_a | filter(lambda x: set(x.transcript_type.split(",")) & transcript_types)
            iter_b = iter_b | filter(lambda x: set(x.transcript_type.split(",")) & transcript_types)
        a = iter_a | apply(lambda x: bed.write(f"{x.to_bed6()}\n")) | to_list
        b = iter_b | apply(lambda x: bed.write(f"{x.to_bed6()}\n")) | to_list
        filtered_gene = {'+': a, '-': b}

    if filtered_bed is not None:
        bedtools_sort(filtered_bed)

    if dump is not None:
        with open(dump, "wb") as f:
            pickle.dump(filtered_gene, f)
    end = time.perf_counter()
    logging.info(f"save {len(a) + len(b)} filtered genes.")
    logging.info(f"task finished in {end-start:.2f}s.")
    return filtered_gene


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--load", dest="load_pickle", type=str, required=True, help="load pickle file.")
    parser.add_argument("-b", dest="bed", type=str, default=None, help="output filtered bed file.")
    parser.add_argument("-p", dest="pickle", type=str, default=None, help="output pickle dump file.")
    parser.add_argument("--gene_length", dest="gene_length", type=int, default=5000, help="gene_length cutoff.")
    parser.add_argument("--downStream", dest="downStream", type=int, default=10000, help="downStream cutoff.")
    parser.add_argument('--gene_types', dest='gene_types', nargs="+", default=None, help="choose gene types to filter gtf.")
    parser.add_argument('--transcript_types', dest='transcript_types', nargs="+", default=None, help="choose transcript types to filter gtf.")
    args = parser.parse_args()
    filter_gene(load_gene(args.load_pickle), args.bed, args.pickle, args.gene_length, args.downStream, args.gene_types, args.transcript_types)


if __name__ == "__main__":
    run()
