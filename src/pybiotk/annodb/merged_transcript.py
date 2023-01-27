 # -*- coding: utf-8 -*-
from __future__ import annotations

import itertools
import os
from collections import deque
from functools import partial
from typing import List, Dict, Tuple, Literal, Iterable, Iterator, Sequence, Optional, TYPE_CHECKING

from pybiotk.intervals import merge_intervals
from pybiotk.io.bed import Bed6
from pybiotk.utils import bedtools_sort
from stream import window, to_list, filter, transform, mapwith, apply, uniq, flatten

if TYPE_CHECKING:
    from pybiotk.io import GtfFile
    from pybiotk.annodb import Transcript


class MergedTranscript:
    def __init__(
        self,
        transcript_id: Optional[str] = None,
        transcript_name: Optional[str] = None,
        transcript_type: Optional[str] = None,
        gene_id: Optional[str] = None,
        gene_name: Optional[str] = None,
        gene_type: Optional[str] = None,
        chrom: Optional[str] = None,
        start: int = 0,
        end: int = 0,
        strand: Optional[Literal['+', '-']] = None,
        cds_start: Optional[int] = None,
        cds_end: Optional[int] = None,
        exons: Iterable[Tuple[int, int]] = (),
        count: int = 0,
        before: Optional[int] = None,
        after: Optional[int] = None,
    ):
        self.transcript_id = transcript_id
        self.transcript_name = transcript_name
        self.transcript_type = transcript_type
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_type = gene_type
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.cds_start = int(cds_start) if cds_start is not None else cds_start
        self.cds_end = int(cds_end) if cds_end is not None else cds_end
        self.exons = exons
        self.count = count
        self.before = int(before) if before is not None else before
        self.after = int(after) if after is not None else after

    def __repr__(self) -> str:
        return f"{self.gene_name}:{self.chrom}:{self.start}-{self.end}({self.strand})"

    __str__ = __repr__

    @classmethod
    def init_by_transcripts(cls, transcripts: Iterable[Transcript]):
        min_start = float("inf")
        max_end = 0
        min_cds_start = float("inf")
        max_cds_end = 0
        transcript_ids = []
        transcript_names = []
        transcript_types = set()
        gene_ids = set()
        gene_names = set()
        gene_types = set()
        chroms = set()
        strands = set()
        exons = []
        for transcript in transcripts:
            transcript_ids.append(transcript.transcript_id)
            transcript_names.append(transcript.transcript_name)
            transcript_types.add(transcript.transcript_type)
            gene_ids.add(transcript.gene_id)
            gene_names.add(transcript.gene_name)
            gene_types.add(transcript.gene_type)
            chroms.add(transcript.chrom)
            strands.add(transcript.strand)
            exons.extend(transcript.exons())
            if min_start > transcript.start:
                min_start = transcript.start
            if max_end < transcript.end:
                max_end = transcript.end
            if transcript.cds_start is None and transcript.cds_end is None:
                continue
            if min_cds_start > transcript.cds_start:
                min_cds_start = transcript.cds_start
            if max_cds_end < transcript.cds_end:
                max_cds_end = transcript.cds_end

        exons = merge_intervals(set(exons))
        start = min_start
        end = max_end
        if min_cds_start == float("inf"):
            min_cds_start = None
        if max_cds_end == 0:
            max_cds_end = None
        cds_start = min_cds_start
        cds_end = max_cds_end
        count = len(transcript_ids)
        transcript_id = ",".join(transcript_ids)
        transcript_name = ",".join(transcript_names)
        transcript_type = ",".join(list(transcript_types))
        gene_id = ",".join(list(gene_ids))
        gene_name = ",".join(list(gene_names))
        gene_type = ",".join(list(gene_types))
        chrom = ",".join(list(chroms))
        strand = ",".join(list(strands))
        return cls(transcript_id, transcript_name, transcript_type,
                   gene_id, gene_name, gene_type,
                   chrom, start, end, strand,
                   cds_start, cds_end, exons, count)

    def to_bed6(self) -> Bed6:
        return Bed6(self.chrom, self.start, self.end, self.gene_name, str(self.count), self.strand)

    def upStream(self):
        if self.strand == '-':
            return self.after
        else:
            return self.before

    def downStream(self):
        if self.strand == '-':
            return self.before
        else:
            return self.after

    def length(self):
        return self.end - self.start


def group_overlap_transcripts(iterable: Iterable[Transcript]) -> Iterator[Tuple[Transcript, ...]]:
    a = deque(itertools.islice(iterable, 1))
    max_end = 0
    for i in iterable:
        before, later = a[-1], i
        max_end = max(max_end, before.end)
        if (before.chrom == later.chrom) and (later.start <= max_end):
            a.append(i)
        else:
            yield tuple(a)
            if not before.chrom == later.chrom:
                max_end = 0
            a.clear()
            a.append(i)
    if a:
        yield tuple(a)
        a.clear()


def add_before_and_after(x: Tuple[MergedTranscript, ...]):
    x[0].after = x[1].start - x[0].end
    x[1].before = x[0].after


def add_chrom_ends_before_and_after(x: MergedTranscript, chrom_length_dict: Optional[Dict[str, int]] = None):
    if x.before is None:
        x.before = x.start
    if x.after is None:
        if chrom_length_dict is not None:
            try:
                x.after = chrom_length_dict[x.chrom] - x.end
            except KeyError:
                x.after = 0
        else:
            x.after = 0


def merge_transcripts(gtf: GtfFile, strand: Optional[Literal["+", "-"]] = "+",
                      escape_gene_types: Sequence[str] = (),
                      escape_gene_name_startswith: Tuple[str] = (),
                      chrom_length_dict: Optional[Dict[str, int]] = None
                      ) -> List[MergedTranscript]:
    escape = set(escape_gene_types)
    add_chrom_ends_before_and_after_partial = partial(add_chrom_ends_before_and_after, chrom_length_dict=chrom_length_dict)
    merged_transcripts = gtf.to_transcript() | filter(lambda x: not x.gene_name.startswith(
        escape_gene_name_startswith)) | filter(lambda x: x.gene_type not in escape) | filter(
            lambda x: x.strand == strand if strand else True) | transform(
                group_overlap_transcripts) | mapwith(MergedTranscript.init_by_transcripts) | filter(
                    lambda x: x.strand in {'+', '-'}) | window | filter(
                    lambda x: x[0].chrom == x[-1].chrom) | apply(add_before_and_after) | flatten | uniq(
                        lambda x: x.gene_name) | apply(add_chrom_ends_before_and_after_partial)
    return merged_transcripts


def merge_transcripts_groupby_strand(
    gtf: GtfFile,
    escape_gene_types: Sequence[str] = (),
    escape_gene_name_startswith: Tuple[str] = (),
    chrom_length_dict: Optional[Dict[str, int]] = None,
    savebed: Optional[str] = None,
    remove_strand_overlap: bool = False
) -> Dict[str, List[MergedTranscript]]:
    bedpath = savebed if savebed is not None else os.devnull
    with open(bedpath, "w") as bed:
        if remove_strand_overlap:
            fwd = []
            rev = []
            for x in merge_transcripts(gtf, None, escape_gene_types, escape_gene_name_startswith, chrom_length_dict):
                bed.write(f"{x.to_bed6()}\n")
                if x.strand == '+':
                    fwd.append(x)
                elif x.strand == '-':
                    rev.append(x)
                else:
                    raise RuntimeError(f"Wrong strand for {x.gene_name}.")
            merged_transcripts_dict = {'+': fwd, '-': rev}
        else:
            merged_transcripts_dict = {
                '+': merge_transcripts(gtf, '+', escape_gene_types, escape_gene_name_startswith, chrom_length_dict) | apply(
                    lambda x: bed.write(f"{x.to_bed6()}\n")) | to_list,
                '-': merge_transcripts(gtf, '-', escape_gene_types, escape_gene_name_startswith, chrom_length_dict) | apply(
                    lambda x: bed.write(f"{x.to_bed6()}\n")) | to_list,
            }
    if savebed is not None:
        bedtools_sort(savebed)
    return merged_transcripts_dict
