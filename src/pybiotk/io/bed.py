# -*- coding: utf-8 -*-
from __future__ import annotations

from dataclasses import dataclass, field
from io import TextIOWrapper
from typing import List, Tuple, Literal, Iterable, Union, Optional, TextIO, TYPE_CHECKING

if TYPE_CHECKING:
    from pybiotk.io.gtf import GTF


@dataclass
class TransInfo:
    transcript_id: str = field(default=None)
    transcript_name: str = field(default=None)
    transcript_type: str = field(default=None)
    gene_id: str = field(default=None)
    gene_name: str = field(default=None)
    gene_type: str = field(default=None)
    strand: Literal['+', '-'] = field(default="+")

    @classmethod
    def init_by_gtf(cls, gtf: GTF):
        return cls(gtf.transcript_id(), gtf.transcript_name(), gtf.transcript_type(),
                   gtf.gene_id(), gtf.gene_name(), gtf.gene_type(), gtf.strand)

    def __str__(self):
        return "\t".join(str(s) for s in list(self.__dict__.values()))


@dataclass
class GeneInfo:
    gene_id: str = field(default=None)
    gene_name: str = field(default=None)
    gene_type: str = field(default=None)
    strand: Literal['+', '-'] = field(default='+')

    @classmethod
    def init_by_gtf(cls, gtf: GTF):
        return cls(gtf.gene_id(), gtf.gene_name(), gtf.gene_type(), gtf.strand)

    def __str__(self):
        return "\t".join(str(s) for s in list(self.__dict__.values()))


@dataclass
class Bed6:
    chrom: str = field(default=None)
    start: int = field(default=0)
    end: int = field(default=0)
    name: str = field(default=None)
    score: str = field(default=None)
    strand: Literal['+', '-'] = field(default='+')

    def __post_init__(self):
        self.start = int(self.start)
        self.end = int(self.end)

    @classmethod
    def init_by_gtf(cls, gtf: GTF, attr: str):
        return cls(gtf.chrom, gtf.start-1, gtf.end, attr, gtf.score, gtf.strand)

    def __str__(self):
        return "\t".join(str(s) for s in list(self.__dict__.values()))


@dataclass
class Bed12:
    chrom: str = field(default=None)
    start: int = field(default=0)
    end: int = field(default=0)
    name: str = field(default=None)
    score: int = field(default=0, repr=False)
    strand: Literal['+', '-'] = field(default='+')
    thickStart: int = field(default=0, repr=False)
    thickEnd: int = field(default=0, repr=False)
    itemRgb: int = field(default=0, repr=False)
    blockCount: int = field(default=1, repr=False)
    blockSizes: str = field(default=None, repr=False)
    blockStarts: str = field(default=None, repr=False)

    def __post_init__(self):
        self.start = int(self.start)
        self.end = int(self.end)
        self.thickStart = int(self.thickStart)
        self.thickEnd = int(self.thickEnd)

    @classmethod
    def init_by_gtf(cls, gtf: GTF, attr: str, thickStart: Optional[int] = 0, thickEnd: Optional[int] = 0):
        chrom = gtf.chrom
        start = gtf.start - 1
        end = gtf.end
        name = attr
        score = 0 if gtf.score == '.' else int(gtf.score)
        strand = gtf.strand
        thickStart = thickStart
        thickEnd = thickEnd
        itemRgb = 0
        blockCount = 1
        blockSizes = "{size},".format(size=gtf.end - gtf.start + 1)
        blockStarts = "{start},".format(start=0)
        return cls(chrom, start, end, name, score, strand, thickStart,
                   thickEnd, itemRgb, blockCount, blockSizes, blockStarts)

    def __str__(self):
        return "\t".join(str(s) for s in list(self.__dict__.values()))

    def update(self, gtf: GTF):
        self.end = gtf.end
        self.blockSizes += "{size},".format(size=gtf.end - gtf.start + 1)
        self.blockStarts += "{start},".format(start=gtf.start - 1 - self.start)
        self.blockCount += 1

    def exons(self) -> List[Tuple[int, int]]:
        exons = sorted((self.start+int(y), self.start+int(y)+int(x)) for x, y in zip(
            self.blockSizes.rstrip(',').split(','), self.blockStarts.rstrip(',').split(',')))
        return exons


@dataclass
class Intron:
    chrom: str = field(default=None)
    start: int = field(default=0)
    end: int = field(default=0)
    name: str = field(default=None)
    score: str = field(default='.')
    strand: Literal['+', '-'] = field(default='+')

    def __post_init__(self):
        self.start = int(self.start)
        self.end = int(self.end)

    @classmethod
    def init_by_gtf(cls, exon_1: GTF, exon_2: GTF, attr: str):
        return cls(exon_1.chrom, exon_1.end, exon_2.start-1, attr, '.', exon_1.strand)

    def __str__(self):
        return "\t".join(str(s) for s in list(self.__dict__.values()))


class Openbed:
    def __init__(self, filepath_or_buffer: Union[str, TextIO], mode: Literal["bed6", "bed12"] = "bed6"):
        if isinstance(filepath_or_buffer, TextIOWrapper):
            self.filename: str = filepath_or_buffer.name
            self.bed: TextIO = filepath_or_buffer
        else:
            self.filename: str = filepath_or_buffer
            self.bed: TextIO = open(filepath_or_buffer)
        self.mode: str = mode

    def __repr__(self):
        return f"Openbed(filename='{self.filename}', mode='{self.mode}')"

    __str__ = __repr__

    @staticmethod
    def str_to_bed(line: str, mode: Literal["bed6", "bed12"]) -> Union[Bed6, Bed12]:
        args = line.rstrip("\r\n").split("\t")
        args_len = len(args)
        if mode == "bed6":
            if args_len < 6:
                args.extend((6-len(args))*[''])
            if args_len > 6:
                args = args[:6]
            bed = Bed6(*args)
        elif mode == "bed12":
            if args_len < 12:
                args.extend((12-len(args))*[''])
            if args_len > 12:
                args = args[:12]
            bed = Bed12(*args)
        else:
            raise RuntimeError(f"{mode} not support.")
        return bed

    def close(self):
        self.bed.close()

    def __iter__(self) -> Iterable[Union[Bed6, Bed12]]:
        for line in self.bed:
            yield self.str_to_bed(line, self.mode)

    def __next__(self) -> Union[Bed6, Bed12]:
        return self.str_to_bed(next(self.bed), self.mode)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, traceback):
        self.close()
