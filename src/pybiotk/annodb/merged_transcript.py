 # -*- coding: utf-8 -*-
from __future__ import annotations

import itertools
import os
import warnings
from collections import deque
from functools import partial
from typing import List, Dict, Tuple, Literal, Iterable, Iterator, Sequence, Optional, TYPE_CHECKING

from pybiotk.annodb.anno import GFeature
from pybiotk.intervals import merge_intervals
from pybiotk.io.bed import Bed6
from pybiotk.utils import bedtools_sort
from stream import window, to_list, filter, transform, mapwith, apply, uniq, flatten, skip_while

if TYPE_CHECKING:
    from pybiotk.io import GtfFile
    from pybiotk.annodb import Transcript


class MergedTranscript(GFeature):
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
        starts: Iterable[int] = (),
        ends: Iterable[int] = (),
        cds_starts: Iterable[int] = (),
        cds_ends: Iterable[int] = (),
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
        self.starts = list(starts)
        self.ends = list(ends)
        self.cds_starts = list(cds_starts)
        self.cds_ends = list(cds_ends)
        self._exons = list(exons)
        self._introns = None
        self._cds_exons = None
        self._utr5_exons = None
        self._utr3_exons = None
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
        starts = []
        ends = []
        cds_starts = []
        cds_ends = []
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
            starts.append(transcript.start)
            ends.append(transcript.end)
            if min_start > transcript.start:
                min_start = transcript.start
            if max_end < transcript.end:
                max_end = transcript.end
            if transcript.cds_start is None and transcript.cds_end is None:
                continue
            cds_starts.append(transcript.cds_start)
            cds_ends.append(transcript.cds_end)
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
                   cds_start, cds_end, starts, ends, 
                   cds_starts, cds_ends, exons, count)

    def _classify_exons(self):
        if self.cds_start is not None and self.cds_end is not None:
            utr5_exons = []
            utr3_exons = []
            cds_exons = []
            for exon in self.exons():
                if exon[1] <= self.cds_start:
                    utr5_exons.append(exon)
                elif exon[0] >= self.cds_end:
                    utr3_exons.append(exon)
                elif self.cds_start <= exon[0] and exon[1] <= self.cds_end:
                    cds_exons.append(exon)
                elif exon[0] < self.cds_start < exon[1] <= self.cds_end:
                    utr5_exons.append((exon[0], self.cds_start))
                    cds_exons.append((self.cds_start, exon[1]))
                elif self.cds_start <= exon[0] < self.cds_end < exon[1]:
                    cds_exons.append((exon[0], self.cds_end))
                    utr3_exons.append((self.cds_end, exon[1]))
                elif exon[0] < self.cds_start and self.cds_end < exon[1]:
                    utr5_exons.append((exon[0], self.cds_start))
                    cds_exons.append((self.cds_start, self.cds_end))
                    utr3_exons.append((self.cds_end, exon[1]))

            if self.strand == '-':
                utr5_exons, utr3_exons = utr3_exons, utr5_exons
                self._utr5_exons = utr5_exons
                self._cds_exons = cds_exons
                self._utr3_exons = utr3_exons
            loginfo = f"merged exons: {self.exons()}\ncds_start: {self.cds_start}, cds_end: {self.cds_end}"
            if not cds_exons:
                warnings.warn(f"cds_exons of {self.gene_id} does not exist, please check:\n{loginfo}")
        else:
            self._utr5_exons = []
            self._cds_exons = []
            self._utr3_exons = []

    def is_protein_coding(self) -> bool:
        if self.gene_type == 'protein_coding':
            return True
        else:
            return False

    def exons(self) -> List[Tuple[int, int]]:
        return self._exons

    def introns(self) -> List[Tuple[int, int]]:
        if self._introns is None:
            self._introns = [self.start, *(self.exons() | flatten), self.end] | window(2, 2) | skip_while(lambda x: x[0] == x[1]) | to_list
        return self._introns

    def tss_region(self, region: Tuple[int, int] = (-1000, 1000)) -> Tuple[int, ...]:
        if self.strand == '+':
            region = tuple(i+self.start for i in region)
        else:
            region = tuple(reversed(tuple(self.end-i for i in region)))
        return region

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

    downstream = downStream

    def cds_exons(self) -> List[Tuple[int, int]]:
        if self._cds_exons is None:
            self._classify_exons()
        return self._cds_exons

    def utr5_exons(self) -> List[Tuple[int, int]]:
        if self._utr5_exons is None:
            self._classify_exons()
        return self._utr5_exons

    def utr3_exons(self) -> List[Tuple[int, int]]:
        if self._utr3_exons is None:
            self._classify_exons()
        return self._utr3_exons

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
