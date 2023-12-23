# -*- coding: utf-8 -*-
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Tuple, Set, Iterable, AbstractSet, Optional, Literal

from pybiotk.utils import blocks_len, intervals_is_overlap


@dataclass
class GenomicAnnotation:
    id: Optional[str] = field(default=None)
    name: Optional[str] = field(default=None)
    start: Optional[int] = field(default=None)
    end: Optional[int] = field(default=None)
    strand: Optional[Literal['+', '-']] = field(default=None)
    type: Optional[str] = field(default=None)
    detail: Set[str] = field(default_factory=set)

    def update(self, anno: str):
        self.detail.add(anno)

    @staticmethod
    def select_anno(
        annoset: AbstractSet[str],
        priority: Tuple[str, ...] = ("TSS", "TES", "StartCondon", "StopCondon", "5UTR", "3UTR", "CDS", "Exon", "Intron", "Upstream", "Downstream", "Intergenic")
    ):
        for anno in priority:
            if anno in annoset:
                return anno

    def primary_anno(
        self,
        priority: Tuple[str, ...] = ("TSS", "TES", "StartCondon", "StopCondon", "5UTR", "3UTR", "CDS", "Exon", "Intron", "Upstream", "Downstream", "Intergenic")
    ):
        return self.select_anno(self.detail, priority)


@dataclass
class AnnoSet:
    annoset: Iterable[GenomicAnnotation] = field(default_factory=list, repr=False)
    id: List[str] = field(init=False, default_factory=list)
    name: List[str] = field(init=False, default_factory=list)
    start: List[int] = field(init=False, default_factory=list)
    end: List[int] = field(init=False, default_factory=list)
    strand: List[str] = field(init=False, default_factory=list)
    type: List[str] = field(init=False, default_factory=list)
    anno: List[str] = field(init=False, default_factory=list)

    def __post_init__(self):
        first = []
        second = []
        third = []
        fourth = []
        other = []
        for anno in self.annoset:
            if not {"Upstream", "Downstream", "Intergenic"} & anno.detail:
                if not anno.primary_anno() == "Intron":
                    first.append(anno)
                else:
                    third.append(anno)
            else:
                if {"TSS", "TES", "StartCondon", "StopCondon", "5UTR", "3UTR", "CDS", "Exon"} & anno.detail:
                    second.append(anno)
                elif "Intron" in anno.detail:
                    fourth.append(anno)
                else:
                    other.append(anno)
        if first:
            for anno in first:
                self.id.append(anno.id)
                self.name.append(anno.name)
                self.start.append(anno.start)
                self.end.append(anno.end)
                self.strand.append(anno.strand)
                self.type.append(anno.type)
                self.anno.append(anno.primary_anno())
        elif second:
            for anno in second:
                self.id.append(anno.id)
                self.name.append(anno.name)
                self.start.append(anno.start)
                self.end.append(anno.end)
                self.strand.append(anno.strand)
                self.type.append(anno.type)
                self.anno.append(anno.primary_anno())
        elif third:
            for anno in third:
                self.id.append(anno.id)
                self.name.append(anno.name)
                self.start.append(anno.start)
                self.end.append(anno.end)
                self.strand.append(anno.strand)
                self.type.append(anno.type)
                self.anno.append(anno.primary_anno())
        elif fourth:
            for anno in fourth:
                self.id.append(anno.id)
                self.name.append(anno.name)
                self.start.append(anno.start)
                self.end.append(anno.end)
                self.strand.append(anno.strand)
                self.type.append(anno.type)
                self.anno.append(anno.primary_anno())
        else:
            for anno in other:
                self.id.append(anno.id)
                self.name.append(anno.name)
                self.start.append(anno.start)
                self.end.append(anno.end)
                self.strand.append(anno.strand)
                self.type.append(anno.type)
                self.anno.append(anno.primary_anno())

    def primary_anno(
        self,
        priority: Tuple[str, ...] = ("TSS", "TES", "StartCondon", "StopCondon", "5UTR", "3UTR", "CDS", "Exon", "Intron", "Upstream", "Downstream", "Intergenic")
    ) -> str:
        anno = GenomicAnnotation.select_anno(set(self.anno), priority)
        return anno

    def __str__(self) -> str:
        _id = ",".join(self.id)
        _name = ",".join(set(self.name))
        _type = ",".join(set(self.type))
        _start = ",".join(str(i) for i in set(self.start))
        _end = ",".join(str(i) for i in set(self.end))
        _strand = ",".join(str(i) for i in set(self.strand))
        _anno = self.primary_anno()
        return f"{_anno}\t{_start}\t{_end}\t{_strand}\t{_name}\t{_id}\t{_type}"


class GFeature(ABC):
    @abstractmethod
    def is_protein_coding(self) -> bool: ...

    @abstractmethod
    def exons(self) -> List[Tuple[int, int]]: ...

    @abstractmethod
    def introns(self) -> List[Tuple[int, int]]: ...

    @abstractmethod
    def tss(self) -> int: ...

    @abstractmethod
    def tes(self) -> int: ...
    
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
    
    @abstractmethod
    def length(self) -> int: ...

    def exons_len(self) -> int:
        return blocks_len(self.exons())
    
    def exons_count(self) -> int:
        return len(self.exons())

    def introns_len(self) -> int:
        return blocks_len(self.introns())
    
    def introns_count(self) -> int:
        return len(self.introns())

    def cds_len(self) -> int:
        return blocks_len(self.cds_exons())

    def utr5_len(self) -> int:
        return blocks_len(self.utr5_exons())

    def utr3_len(self) -> int:
        return blocks_len(self.utr3_exons())

    def anno(self, blocks: List[Tuple[int, int]], region: Tuple[int, int] = (-3000, 0), down: int = 3000,
             anno_tss: bool = False, anno_tes: bool = False,
             anno_start_condon: bool = False, anno_stop_condon: bool = False) -> Set[str]:
        anno = []
        st = self.tss_region(region=region)
        end = self.downstream(down=down)
        pos = sorted([*st, *end])
        if blocks[-1][1] < pos[0] or blocks[0][0] > pos[3]:
            anno.append("Intergenic")
        upstream = intervals_is_overlap(blocks, [st])
        downstream = intervals_is_overlap(blocks, [end])
        
        if anno_tss:
            tss_block = (self.tss(), self.tss()+1)
            tss = intervals_is_overlap(blocks, [tss_block])
        if anno_tes:
            tes_block = (self.tes(), self.tes()+1)
            tes = intervals_is_overlap(blocks, [tes_block])
        
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
            if anno_start_condon:
                start_condon_region = self.start_condon()
                if start_condon_region is not None:
                    start_condon = intervals_is_overlap(blocks, start_condon_region)
                stop_condon_region = self.stop_condon()
                if stop_condon_region is not None:
                    stop_condon = intervals_is_overlap(blocks, stop_condon_region)                    
        else:
            exons = intervals_is_overlap(blocks, self.exons())
            utr5_exons = False
            utr3_exons = False
            cds_exons = False
            start_condon = False
            stop_condon = False

        if anno_tss:
            if tss:
                anno.append("TSS")
        if anno_tes:
            if tes:
                anno.append("TES")
        if upstream:
            anno.append("Upstream")
        if self.is_protein_coding():
            if utr5_exons:
                anno.append("5UTR")
            if utr3_exons:
                anno.append("3UTR")
            if cds_exons:
                anno.append("CDS")
            if anno_start_condon:
                if start_condon:
                    anno.append("StartCondon")
            if anno_stop_condon:
                if stop_condon:
                    anno.append("StopCondon")
        elif exons:
            anno.append("Exon")
        if downstream:
            anno.append("Downstream")
        if not anno:
            anno.append("Intron")
        elif not {"5UTR", "3UTR", "CDS", "Exon"} & set(anno):
            if intervals_is_overlap(blocks, [(pos[1], pos[2])]):
                anno.append("Intron")

        return set(anno)
