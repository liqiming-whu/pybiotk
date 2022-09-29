# -*- coding: utf-8 -*-
import importlib
import importlib.resources
import logging
import os
import subprocess
import sys
import warnings
from functools import wraps
from io import TextIOWrapper
from types import ModuleType
from typing import List, Dict, Sequence, Tuple, Literal, Iterator, Iterable, Optional, Callable, Union, TextIO

import pandas as pd

logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s", datefmt='%Y-%m-%d %A %H:%M:%S', level=logging.INFO)


def reverse_seq(seq: str) -> str:
    seq = seq.upper()
    complement = {
        "A": "T", "C": "G", "G": "C", "T": "A",
        "N": "N", "W": "W", "S": "S", "K": "M",
        "M": "K", "Y": "R", "R": "Y", "H": "D",
        "B": "V", "V": "B", "D": "H", "&": "&",
    }
    return "".join(reversed([complement[i] for i in seq]))


def is_overlap(interval1: Tuple[int, int], interval2: Tuple[int, int]) -> bool:
    x1, y1 = interval1
    x2, y2 = interval2
    if y1 <= x2 or y2 <= x1:
        return False
    else:
        return True


def intervals_is_overlap(intervals1: List[Tuple[int, int]], intervals2: List[Tuple[int, int]]) -> bool:
    if not is_overlap((intervals1[0][0], intervals1[-1][1]), (intervals2[0][0], intervals2[-1][1])):
        return False
    for interval1 in intervals1:
        for interval2 in intervals2:
            if is_overlap(interval1, interval2):
                return True
    return False


def cigar_tuples2blocks(start: int, cigartuples: Sequence[Tuple[int, int]], shift: int = 0) -> List[Tuple[int, int]]:
    tmp_abs_start = shift + start
    blocks = []
    for cigar, length in cigartuples:
        if cigar in [0, 2]:
            blocks.append((tmp_abs_start, tmp_abs_start + length))
            tmp_abs_start += length
        elif cigar in [3, 6, 10]:
            tmp_abs_start += length
    merged_blocks = []
    tmp_start = None
    block_num = len(blocks)
    for i in range(block_num):
        if tmp_start is None:
            tmp_start = blocks[i][0]
        if not ((i < block_num - 1) and (blocks[i + 1][0] <= blocks[i][1])):
            merged_blocks.append((tmp_start, blocks[i][1]))
            tmp_start = None
    return merged_blocks


def merge_blocks(blocks1: List[Tuple[int,int]], blocks2: List[Tuple[int,int]]) -> Tuple[bool, List[Tuple[int,int]]]:
    if not intervals_is_overlap(blocks1, blocks2):
        return False, sorted(blocks1 + blocks2, key=lambda x: x[0])
    block1_num = len(blocks1)
    block2_num = len(blocks2)
    js1_li = list()
    js2_li = list()
    if block1_num > 1:
        js1_li = [(blocks1[i][1], blocks1[i + 1][0]) for i in range(len(blocks1) - 1)]
    if block2_num > 1:
        js2_li = [(blocks2[i][1], blocks2[i + 1][0]) for i in range(len(blocks2) - 1)]
    for js1_s, js1_e in js1_li:
        for js2_s, js2_e in js2_li:
            if (js1_s == js2_s) and (js1_e == js2_e):
                continue
            if not ((js1_s >= js1_e) or (js2_s >= js2_e)):
                conflict = True
                return conflict, None
    conflict = False
    js_li = sorted(set(js1_li + js2_li))
    start = min([x[0] for x in blocks1 + blocks2])
    end = max([x[1] for x in blocks1 + blocks2])
    if js_li:
        new_blocks = list()
        for i in range(len(js_li)):
            if i == 0:
                new_blocks.append((start, js_li[i][0]))
            else:
                new_blocks.append((js_li[i - 1][1], js_li[i][0]))
            if i == (len(js_li) - 1):
                new_blocks.append((js_li[i][1], end))
        return conflict, new_blocks
    else:
        return conflict, [(start, end)]


def split_discontinuous(discontinuous_list: Iterable[int]) -> Iterator[List[int]]:
    continuous_list = []
    for p in discontinuous_list:
        if not continuous_list:
            continuous_list.append(p)
        elif p - continuous_list[-1] == 1:
            continuous_list.append(p)
        else:
            yield continuous_list
            continuous_list = [p]
    if continuous_list:
        yield continuous_list


def read_table(filepath_or_buffer: Union[str, TextIO], header: Optional[int] = 0, sep: Optional[str] = None, **kwargs) -> pd.DataFrame:
    if filepath_or_buffer == "-":
        filepath_or_buffer = sys.stdin
    if isinstance(filepath_or_buffer, TextIOWrapper):
        filename = filepath_or_buffer.name
    else:
        filename = filepath_or_buffer
    ext = os.path.splitext(filename)[1]
    try:
        if ext == ".csv":
            deli = sep if sep else ","
            df = pd.read_csv(filepath_or_buffer, delimiter=deli, header=header, **kwargs)
        elif ext == ".xlsx":
            df = pd.read_excel(filepath_or_buffer, header=header, **kwargs)
        else:
            deli = sep if sep else "\t"
            df = pd.read_table(filepath_or_buffer, delimiter=deli, header=header, **kwargs)
    except pd.errors.EmptyDataError:
        df = pd.DataFrame()
    if isinstance(filepath_or_buffer, TextIOWrapper):
        filepath_or_buffer.close()
    return df


def write_table(out_df: pd.DataFrame, path_or_buf: Union[str, TextIO], index: bool = False, header: bool = True, sep: str = "\t", **kwargs):
    if path_or_buf == "-":
        path_or_buf = sys.stdout
    if isinstance(path_or_buf, TextIOWrapper):
        filename = path_or_buf.name
    else:
        filename = path_or_buf
    ext = os.path.splitext(filename)[1]
    if ext == ".csv":
        out_df.to_csv(path_or_buf, index=index, header=header, sep=",", **kwargs)
    elif ext == ".xlsx":
        out_df.to_excel(path_or_buf, index=index, header=header, **kwargs)
    else:
        out_df.to_csv(path_or_buf, index=index, header=header, sep=sep, **kwargs)

    if isinstance(path_or_buf, TextIOWrapper):
        path_or_buf.close()


def blocks_len(blocks: Iterable[Tuple[int, int]]) -> int:
    return sum(y-x for (x, y) in blocks)


def infer_fragment_strand(strand: Literal["+", "-"], rule: str = "1+-,1-+,2++,2--", is_read2: bool = False):
    fragment_strand = None
    rules = ["1+-,1-+,2++,2--", "1++,1--,2+-,2-+", "+-,-+", "++,--"]
    assert rule in rules
    if rule == rules[0] or rule == rules[2]:
        if strand == "-":
            fragment_strand = "+"
        elif strand == "+":
            fragment_strand = "-"
    else:
        if strand == "+":
            fragment_strand = "+"
        elif strand == "-":
            fragment_strand = "-"

    if is_read2:
        if fragment_strand == "+":
            fragment_strand = "-"
        else:
            fragment_strand = "+"
    return fragment_strand


def ignore(func: Callable):
    @wraps(func)
    def wrapper(*args, **kargs):
        try:
            traceback = func(*args, **kargs)
        except BrokenPipeError:
            sys.stdout = None
            sys.exit(0)
        except KeyboardInterrupt:
            sys.stdout = None
            sys.exit(0)
        return traceback
    return wrapper


def load_chrom_length_dict(filepath_or_buffer: Union[str, TextIOWrapper]) -> Dict[str, int]:
    if isinstance(filepath_or_buffer, TextIOWrapper):
        buffer = filepath_or_buffer
    else:
        buffer = open(filepath_or_buffer)
    chrom_length_dict = {}
    with buffer as f:
        for line in f:
            fileds = line.split()
            chrom_length_dict[fileds[0]] = int(fileds[1])

    return chrom_length_dict


def package_resource(package: Union[str, ModuleType], resource: Union[str, os.PathLike]):
    with importlib.resources.path(package, resource) as path:
        pathstr = str(path)
    return pathstr


def default_hg38_chrom_length_dict() -> Dict[str, int]:
    return load_chrom_length_dict(importlib.resources.open_text("pybiotk.data", "hg38.chrom.sizes"))


def default_mm10_chrom_length_dict() -> Dict[str, int]:
    return load_chrom_length_dict(importlib.resources.open_text("pybiotk.data", "mm10.chrom.sizes"))


def bedtools_sort(path: str):
    logging.info(f"start to sort {path}...")
    tmp_path = path + ".tmp"
    os.rename(path, tmp_path)
    cmd = f"bedtools sort -i {tmp_path} > {path}"
    try:
        subprocess.check_call(cmd, shell=True)
        os.remove(tmp_path)
    except subprocess.CalledProcessError:
        os.rename(tmp_path, path)
        warnings.warn(f"An error occurred while executing bedtools, bedtools may not be installed. {path} is unsorted")
