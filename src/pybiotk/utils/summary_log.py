#!/usr/bin/env python3
import argparse
import os
import re
import sys
from collections import namedtuple
from typing import Optional, Sequence, TextIO, Union

import pandas as pd

from pybiotk.utils import write_table, ignore
from stream import cat, head


def parse_cutadapt(filename: str):
    total_line = ""
    reads1_line = ""
    reads2_line = ""
    for line in cat(filename) | head(15):
        line = line.strip()
        if line.startswith("Total read pairs processed:"):
            total_line = line.replace(",", "")
        if line.startswith("Read 1 with adapter:"):
            reads1_line = line.replace(",", "")
        if line.startswith("Read 2 with adapter:"):
            reads2_line = line.replace(",", "")
    assert all([total_line, reads1_line, reads2_line])
    total = re.compile('Total read pairs processed:\\s*(\\d*)').findall(total_line)[0]
    read1 = re.compile('Read 1 with adapter:\\s*(\\d*.*)').findall(reads1_line)[0]
    read1_num = read1.split(" ")[0]
    read1_percent = read1.split(" ")[1].lstrip("(").rstrip(")")
    read2 = re.compile('Read 2 with adapter:\\s*(\\d*.*)').findall(reads2_line)[0]
    read2_num = read2.split(" ")[0]
    read2_percent = read2.split(" ")[1].lstrip("(").rstrip(")")

    cutadapt_summary = namedtuple("cutadapt_summary", [
        "total_read_pairs",
        "read1_with_adapter",
        "read2_with_adapter",
        "read1_with_adapter_percent",
        "read2_with_adapter_percent"
        ])
    return cutadapt_summary(total, read1_num, read2_num, read1_percent, read2_percent)


def parse_bowtie2(filename: str):
    lines = [line.strip() for line in cat(filename)]
    input_read_pairs = lines[0].split()[0]
    mapped_read_pairs = str(int(lines[3].split()[0]) + int(lines[4].split()[0]))
    alignment_rate = lines[5].split()[0]
    bowtie2_summary = namedtuple("bowtie2_summary", [
        "input_read_pairs",
        "mapped_read_pairs",
        "alignment_rate"
    ])

    return bowtie2_summary(input_read_pairs, mapped_read_pairs, alignment_rate)


def parse_picards_rmdup(filename: str):
    lines = ""
    for line in cat(filename) | head(10):
        line = line.strip()
        if line.startswith("Unknown Library"):
            lines = line
    assert lines
    line_toks = lines.split()
    estimated_librarysize = line_toks[-1]
    duplication_rate = f"{float(line_toks[-2]) * 100:.2f}%"
    picards_rmdup_summary = namedtuple("picards_rmdup_summary", ["duplication_rate", "estimated_librarysize"])
    return picards_rmdup_summary(duplication_rate, estimated_librarysize)


def parse_flagstat(filename: str):
    lines = [line.strip() for line in cat(filename)]
    filtered_reads = lines[0].split()[0]
    return filtered_reads


def parse_sample_name(sample_name: str):
    rep_finder = re.compile('.*(\\d)$').findall(sample_name)
    rep = f"rep{rep_finder[0]}" if rep_finder else "rep0"
    group = sample_name.rstrip(rep).rstrip("_").rstrip("-")
    sample = namedtuple("sample", ["name", "replicate", "group"])
    return sample(sample_name, rep, group)


def parse_path(filename: str):
    return os.path.splitext(os.path.basename(filename))[0].split(".")[0]


def main(filepath_or_buffer: Union[str, TextIO],
         cutadapt: Optional[Sequence[str]],
         bowtie2: Optional[Sequence[str]] = None,
         picards_rmdup: Optional[Sequence[str]] = None,
         flagstat: Optional[Sequence[str]] = None,
         reads_in_peaks: Optional[Sequence[str]] = None,
         sample_names: Optional[Sequence[str]] = None,
         replicates: Optional[Sequence[str]] = None,
         group: Optional[Sequence[str]] = None,
         ):
    sample_names = sample_names if sample_names is not None else [parse_path(i) for i in cutadapt]
    group = group if group is not None else [parse_sample_name(i).group for i in sample_names]
    replicates = replicates if replicates is not None else [parse_sample_name(i).replicate for i in sample_names]
    assert len(sample_names) == len(cutadapt)
    data = dict()
    data["samples"] = sample_names
    data["group"] = group
    data["replicates"] = replicates
    total_read_pairs = []
    read1_with_adapter = []
    read2_with_adapter = []
    read1_with_adapter_percent = []
    read2_with_adapter_percent = []
    for log in cutadapt:
        cutadapt_summary = parse_cutadapt(log)
        total_read_pairs.append(cutadapt_summary.total_read_pairs)
        read1_with_adapter.append(cutadapt_summary.read1_with_adapter)
        read2_with_adapter.append(cutadapt_summary.read2_with_adapter)
        read1_with_adapter_percent.append(cutadapt_summary.read1_with_adapter_percent)
        read2_with_adapter_percent.append(cutadapt_summary.read2_with_adapter_percent)
    data["total_read_pairs"] = total_read_pairs
    data["read1_with_adapter"] = read1_with_adapter
    data["read2_with_adapter"] = read2_with_adapter
    data["read1_with_adapter_percent"] = read1_with_adapter_percent
    data["read2_with_adapter_percent"] = read2_with_adapter_percent
    if bowtie2 is not None:
        assert len(bowtie2) == len(sample_names)
        input_read_pairs = []
        mapped_read_pairs = []
        alignment_rate = []
        for log in bowtie2:
            bowtie2_log = parse_bowtie2(log)
            input_read_pairs.append(bowtie2_log.input_read_pairs)
            mapped_read_pairs.append(bowtie2_log.mapped_read_pairs)
            alignment_rate.append(bowtie2_log.alignment_rate)
        data["input_read_pairs"] = input_read_pairs
        data["mapped_read_pairs"] = mapped_read_pairs
        data["alignment_rate"] = alignment_rate
    if picards_rmdup is not None:
        assert len(picards_rmdup) == len(sample_names)
        duplication_rate = []
        estimated_librarysize = []
        for log in picards_rmdup:
            picards_rmdup_summary = parse_picards_rmdup(log)
            duplication_rate.append(picards_rmdup_summary.duplication_rate)
            estimated_librarysize.append(picards_rmdup_summary.estimated_librarysize)
        data["duplication_rate"] = duplication_rate
        data["estimated_librarysize"] = estimated_librarysize
    if flagstat is not None:
        assert len(flagstat) == len(sample_names)
        filtered_reads = []
        for log in flagstat:
            filtered_reads.append(parse_flagstat(log))
        data["filtered_reads"] = filtered_reads
    if reads_in_peaks is not None:
        assert len(reads_in_peaks) == len(sample_names)
        peaks = []
        in_peaks = []
        for log in reads_in_peaks:
            with open(log) as l:
                toks = l.read().split()
            peaks.append(toks[0])
            in_peaks.append(toks[1])
        data["peaks_count"] = peaks
        data["reads_in_peaks"] = in_peaks
    df = pd.DataFrame(data)
    write_table(df, filepath_or_buffer)


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-o', dest='output', type=argparse.FileType('w'),
                        default=sys.stdout, help="output file name.")
    parser.add_argument('-s', "--sample_names", dest="sample_names", nargs="+",
                        type=str, default=None, help="sample_names.")
    parser.add_argument('-r', "--replicates", dest="replicates", nargs="+",
                        type=str, default=None, help="replicates.")
    parser.add_argument('-g', "--group", dest="group", nargs="+",
                        type=str, default=None, help="group.")
    parser.add_argument('--cutadapt', dest="cutadapt", nargs="+", type=str,
                        required=True, help="cutadapt log.")
    parser.add_argument('--bowtie2', dest="bowtie2", nargs="+", type=str,
                        default=None, help="bowtie2 log.")
    parser.add_argument('--picards_rmdup', dest="picards_rmdup", nargs="+", type=str,
                        default=None, help="picards_rmdup.")
    parser.add_argument('--flagstat', dest="flagstat", nargs="+", type=str,
                        default=None, help="flagstat.")
    parser.add_argument('--reads_in_peaks', dest="reads_in_peaks", nargs="+", type=str,
                        default=None, help="reads_in_peaks.")
    args = parser.parse_args()
    main(args.output, args.cutadapt, args.bowtie2, args.picards_rmdup, args.flagstat, args.reads_in_peaks, args.sample_names, args.replicates, args.group)


if __name__ == "__main__":
    run()
