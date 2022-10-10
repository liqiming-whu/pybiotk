# -*- coding: utf-8 -*-
import sys
import time
from collections import deque, defaultdict
from enum import Enum
from typing import (
    AbstractSet,
    Set,
    DefaultDict,
    Literal,
    Tuple,
    Optional,
    Iterator
)

import pysam

from pybiotk.utils import logging


class BamType(Enum):
    SE = "SE"
    PE = "PE"


def check_bam_type(filename: str) -> BamType:
    logging.info(f"checking bam type: {filename} ...")
    with pysam.AlignmentFile(filename) as bam:
        try:
            read = next(bam)
        except StopIteration:
            logging.warning(f"empty bam file: {filename}")
            return None
        if read.is_read1 or read.is_read2:
            bamtype = BamType.PE
        else:
            bamtype = BamType.SE
        sys.stderr.write(f"bamtype is {bamtype.value}.\n")
    return bamtype


def count_bam_size(filename: str) -> int:
    with pysam.AlignmentFile(filename) as bam:
        count = bam.count()
    return count


class BamTypeError(RuntimeError):
    pass


class Bam(pysam.AlignmentFile):
    def iter(self, flags: Optional[AbstractSet[int]] = None,
             remove_flags: Optional[AbstractSet[int]] = None,
             secondary=False,
             supplementary=False,
             onlymapped=False) -> Iterator[pysam.AlignedSegment]:
        for read in self:
            if onlymapped:
                if read.is_unmapped:
                    continue
            if not secondary:
                if read.is_secondary:
                    continue
            if not supplementary:
                if read.is_supplementary:
                    continue
            if flags is not None:
                if read.flag not in flags:
                    continue
            if remove_flags is not None:
                if read.flag in remove_flags:
                    continue
            yield read

    def iter_mapped(self, secondary=False, supplementary=False, qcut=None) -> Iterator[pysam.AlignedSegment]:
        for read in self:
            if read.is_unmapped:
                continue
            if read.is_qcfail:
                continue
            if read.is_duplicate:
                continue
            if not secondary:
                if read.is_secondary:
                    continue
            if not supplementary:
                if read.is_supplementary:
                    continue
            if qcut:
                if read.mapping_quality < qcut:
                    continue
            yield read

    def iter_unmapped(self, secondary=False, supplementary=False) -> Iterator[pysam.AlignedSegment]:
        for read in self:
            if not secondary:
                if read.is_secondary:
                    continue
            if not supplementary:
                if read.is_supplementary:
                    continue
            if read.is_qcfail:
                continue
            if read.is_unmapped:
                yield read

    @staticmethod
    def read_to_fastx_record(read) -> pysam.libcfaidx.FastxRecord:
        return pysam.libcfaidx.FastxRecord(
                name=read.query_name,
                comment="",
                sequence=read.get_forward_sequence(),
                quality=pysam.qualities_to_qualitystring(read.get_forward_qualities())
                )

    def to_fastx_record(self) -> Iterator[pysam.libcfaidx.FastxRecord]:
        for read in self.iter_mapped():
            yield self.read_to_fastx_record(read)

    def to_bam(self, filename: str, flags: Optional[AbstractSet[int]] = None,
               remove_flags: Optional[AbstractSet[int]] = None, header=True,
               secondary=False, supplementary=False):
        template = self if header else None
        with pysam.AlignmentFile(filename, mode="wb", template=template) as bam:
            for read in self.iter(flags=flags, remove_flags=remove_flags, secondary=secondary, supplementary=supplementary):
                bam.write(read)


class BamPE(Bam):
    def __init__(self, /, *args, **kwargs):
        super().__init__()
        self.ordered_by_name = False
        self.query_names: Set[str] = set()
        self.read1_set: DefaultDict[str, pysam.AlignedSegment] = defaultdict(type(None))
        self.read2_set: DefaultDict[str, pysam.AlignedSegment] = defaultdict(type(None))
        self.read1_unpaired_querynames: Set[str] = set()
        self.read2_unpaired_querynames: Set[str] = set()
        self.ptr: Optional[int] = None

    def is_ordered_by_name(self):
        self.ordered_by_name = True

    def to_dict(self):
        self.ptr = 0
        logging.warning("saving bam to dict, make sure you have enough memory...")
        start = time.perf_counter()
        for read in self.iter(secondary=False, supplementary=False):
            self.ptr += 1
            if self.ptr % 10000 == 0:
                sys.stderr.write(f"processed {self.ptr} reads.\n")
            self.query_names.add(read.query_name)
            if read.is_read1:
                self.read1_set[read.query_name] = read
            elif read.is_read2:
                self.read2_set[read.query_name] = read
            else:
                raise BamTypeError(f"{self.filename.decode('utf-8')} seems like a single end bam.")
        end = time.perf_counter()
        logging.info(f"{len(self.query_names)} reads have been saved in {end-start:.2f}s.")

    def iter_pair(self, properly_paired=False, onlymapped=True, secondary=False, supplementary=False) -> Iterator[Tuple[pysam.AlignedSegment, ...]]:
        if self.ordered_by_name:
            logging.info(f"{self.filename} is ordered by queryname, use io iter mode ...")
            d = deque(maxlen=2)
            for read in self.iter(secondary=secondary, supplementary=supplementary, onlymapped=onlymapped):
                d.append(read)
                if len(d) == 2:
                    if d[0].query_name == d[1].query_name:
                        if d[0].is_read1 and d[1].is_read2:
                            yield tuple(d)
                        elif d[0].is_read2 and d[1].is_read1:
                            yield d[1], d[0]
                        else:
                            if not properly_paired:
                                if d[0].is_read1 and d[1].is_read1:
                                    yield d[0], None
                                    yield d[1], None
                                elif d[0].is_read2 and d[1].is_read2:
                                    yield None, d[0]
                                    yield None, d[1]
                                else:
                                    raise BamTypeError(f"{d[0]} or {d[1]} in {self.filename} is neither read1 nor read2.")
                        d.clear()
                    else:
                        if not properly_paired:
                            if d[0].is_read1:
                                yield d[0], None
                            elif d[0].is_read2:
                                yield None, d[0]
                            else:
                                raise BamTypeError(f"{d[0]} in {self.filename} is neither read1 nor read2.")
            if d:
                if not properly_paired:
                    if d[-1].is_read1:
                        yield d[-1], None
                    elif d[-1].is_read2:
                        yield None, d[-1]
                    else:
                        raise BamTypeError(f"{d[-1]} in {self.filename} is neither read1 nor read2.")
        else:
            logging.info(f"{self.filename} is not ordered by queryname, use dict iter mode ...")
            if not self.query_names:
                self.to_dict()
            for query_name in list(self.query_names):
                read1 = self.read1_set[query_name]
                read2 = self.read2_set[query_name]

                if read1 is not None and read2 is None:
                    self.read1_unpaired_querynames.add(query_name)

                if read1 is None and read2 is not None:
                    self.read2_unpaired_querynames.add(query_name)

                if properly_paired:
                    if read1 is not None and read2 is not None:
                        yield read1, read2
                else:
                    yield read1, read2

    def iter_unpaired(self, terminal: Literal["read1", "read2"] = "read1") -> Iterator[pysam.AlignedSegment]:
        if self.read1_unpaired_querynames or self.read2_unpaired_querynames:
            if terminal == "read1":
                for query_name in list(self.read1_unpaired_querynames):
                    yield self.read1_set[query_name]
            else:
                for query_name in list(self.read2_unpaired_querynames):
                    yield self.read2_set[query_name]
        else:
            for read1, read2 in self.iter_pair(properly_paired=False):
                if terminal == "read1":
                    if read2 is None:
                        yield read1
                else:
                    if read1 is None:
                        yield read2

    def to_bam_ordered_by_name(self, filename: str):
        header = self.header.to_dict()
        header['HD']['SO'] = "queryname"
        if self.ordered_by_name:
            read1_list = []
            read2_list = []
            with pysam.AlignmentFile(filename, mode="wb", header=header) as bam:
                for read1, read2 in self.iter_pair(properly_paired=False):
                    if read1 is None:
                        read2_list.append(read2)
                    elif read2 is None:
                        read1_list.append(read1)
                    elif read1 is not None and read2 is not None:
                        bam.write(read1)
                        bam.write(read2)
                for read in read1_list:
                    bam.write(read)
                for read in read2_list:
                    bam.write(read)
        else:
            with pysam.AlignmentFile(filename, mode="wb", header=header) as bam:
                for read1, read2 in self.iter_pair(properly_paired=True):
                    bam.write(read1)
                    bam.write(read2)
                for read in self.iter_unpaired(terminal="read1"):
                    bam.write(read)
                for read in self.iter_unpaired(terminal="read2"):
                    bam.write(read)

    def to_fastx_record_pair(self) -> Iterator[Tuple[pysam.libcfaidx.FastxRecord, ...]]:
        for read1, read2 in self.iter_pair(properly_paired=True):
            yield self.read_to_fastx_record(read1), self.read_to_fastx_record(read2)

    def to_fastx_record_unpaired(self, terminal: Literal["read1", "read2"] = "read1") -> Iterator[pysam.libcfaidx.FastxRecord]:
        for read in self.iter_unpaired(terminal):
            yield self.read_to_fastx_record(read)
