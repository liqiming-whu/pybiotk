# -*- coding: utf-8 -*-
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Tuple, Set, Iterable, AbstractSet, Optional

from pybiotk.utils import blocks_len, intervals_is_overlap


@dataclass
class GenomicAnnotation:
    id: Optional[str] = field(default=None)
    name: Optional[str] = field(default=None)
    start: Optional[int] = field(default=None)
    end: Optional[int] = field(default=None)
    type: Optional[str] = field(default=None)
    detail: Set[str] = field(default_factory=set)

    def update(self, anno: str):
        self.detail.add(anno)

    @staticmethod
    def select_anno(
        annoset: AbstractSet[str],
        priority: Tuple[str, ...] = ("Promoter", "5UTR", "3UTR", "CDS", "Exon", "Intron", "Downstream", "Intergenic")
    ):
        for anno in priority:
            if anno in annoset:
                return anno

    def primary_anno(
        self,
        priority: Tuple[str, ...] = ("Promoter", "5UTR", "3UTR", "CDS", "Exon", "Intron", "Downstream", "Intergenic")
    ):
        return self.select_anno(self.detail, priority)


@dataclass
class AnnoSet:
    annoset: Iterable[GenomicAnnotation] = field(default_factory=list, repr=False)
    id: List[str] = field(init=False, default_factory=list)
    name: List[str] = field(init=False, default_factory=list)
    start: List[int] = field(init=False, default_factory=list)
    end: List[int] = field(init=False, default_factory=list)
    type: List[str] = field(init=False, default_factory=list)
    anno: List[str] = field(init=False, default_factory=list)

    def __post_init__(self):
        for anno in self.annoset:
            self.id.append(anno.id)
            self.name.append(anno.name)
            self.start.append(anno.start)
            self.end.append(anno.end)
            self.type.append(anno.type)
            self.anno.append(anno.primary_anno())

    def primary_anno(
        self,
        priority: Tuple[str, ...] = ("Promoter", "5UTR", "3UTR", "CDS", "Exon", "Intron", "Downstream", "Intergenic")
    ) -> str:
        return GenomicAnnotation.select_anno(set(self.anno), priority)

    def __str__(self) -> str:
        _id = ",".join(self.id)
        _name = ",".join(set(self.name))
        _type = ",".join(set(self.type))
        _start = ",".join(str(i) for i in set(self.start))
        _end = ",".join(str(i) for i in set(self.end))
        _anno = self.primary_anno()
        return f"{_anno}\t{_start}\t{_end}\t{_name}\t{_id}\t{_type}"


class GFeature(ABC):
    @abstractmethod
    def is_protein_coding(self) -> bool: ...

    @abstractmethod
    def exons(self) -> List[Tuple[int, int]]: ...

    @abstractmethod
    def introns(self) -> List[Tuple[int, int]]: ...

    @abstractmethod
    def tss_region(self, region: Tuple[int, int] = (-1000, 1000)) -> Tuple[int, int]: ...

    @abstractmethod
    def downstream(self, down: int = 3000) -> Tuple[int, int]: ...

    @abstractmethod
    def cds_exons(self) -> List[Tuple[int, int]]: ...

    @abstractmethod
    def utr5_exons(self) -> List[Tuple[int, int]]: ...

    @abstractmethod
    def utr3_exons(self) -> List[Tuple[int, int]]: ...

    def exons_len(self) -> int:
        return blocks_len(self.exons())

    def introns_len(self) -> int:
        return blocks_len(self.introns())

    def cds_len(self) -> int:
        return blocks_len(self.cds_exons())

    def utr5_len(self) -> int:
        return blocks_len(self.utr5_exons())

    def utr3_len(self) -> int:
        return blocks_len(self.utr3_exons())

    def anno(self, blocks: List[Tuple[int, int]], region: Tuple[int, int] = (-1000, 1000), down: int = 3000) -> Set[str]:
        anno = []
        tss = intervals_is_overlap(blocks, [self.tss_region(region=region)])
        introns = self.introns()
        if introns:
            introns = intervals_is_overlap(blocks, introns)
        downstream = intervals_is_overlap(blocks, [self.downstream(down=down)])
        if self.is_protein_coding():
            utr5_exons = self.utr5_exons()
            if utr5_exons:
                utr5_exons = intervals_is_overlap(blocks, utr5_exons)
            utr3_exons = self.utr3_exons()
            if utr3_exons:
                utr3_exons = intervals_is_overlap(blocks, utr3_exons)
            cds_exons = self.cds_exons()
            if cds_exons:
                cds_exons = intervals_is_overlap(blocks, cds_exons)
            exons = False
        else:
            exons = intervals_is_overlap(blocks, self.exons())
            utr5_exons = False
            utr3_exons = False
            cds_exons = False

        if tss:
            anno.append("Promoter")
        if self.is_protein_coding():
            if utr5_exons:
                anno.append("5UTR")
            if utr3_exons:
                anno.append("3UTR")
            if cds_exons:
                anno.append("CDS")
        elif exons:
            anno.append("Exon")
        if introns:
            anno.append("Intron")
        if downstream:
            anno.append("Downstream")
        if not anno:
            anno.append("Intergenic")
        return set(anno)
