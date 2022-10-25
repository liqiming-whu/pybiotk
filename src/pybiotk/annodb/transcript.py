# -*- coding: utf-8 -*-
from __future__ import annotations

import warnings
from typing import List, Tuple, Literal, Iterable, Optional, TYPE_CHECKING

from pybiotk.annodb import GFeature, GenomicAnnotation
from stream import flatten, window, skip_while, to_list

if TYPE_CHECKING:
    from pybiotk.io import TransInfo, Bed12, GTF


class Transcript(GFeature):
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
        exons: Iterable[Tuple[int, int]] = ()
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
        self._exons = list(exons)
        self._introns = None
        self._cds_exons = None
        self._utr5_exons = None
        self._utr3_exons = None
        self.id = self.transcript_id

    def __repr__(self) -> str:
        return f"{self.transcript_id}:{self.chrom}:{self.start}-{self.end}({self.strand})"

    __str__ = __repr__

    @classmethod
    def init_by_gtf(cls, gtf: GTF):
        return cls(gtf.transcript_id(), gtf.transcript_name(), gtf.transcript_type(),
                   gtf.gene_id(), gtf.gene_name(), gtf.gene_type(), gtf.chrom,
                   gtf.start - 1, gtf.end, gtf.strand, 0, 0, [(gtf.start-1, gtf.end)])

    def update(self, gtf: GTF):
        self.end = gtf.end
        self._exons.append((gtf.start-1, gtf.end))

    @classmethod
    def init_by_bed(cls, transinfo: TransInfo, bed12: Bed12):
        return cls(*transinfo.__slots__, bed12.chrom, bed12.start, bed12.end,
                   bed12.strand, bed12.thickStart, bed12.thickEnd, bed12.exons())

    def is_protein_coding(self) -> bool:
        if self.transcript_type == 'protein_coding':
            return True
        else:
            return False

    def exons(self) -> List[Tuple[int, int]]:
        return self._exons

    def introns(self) -> List[Tuple[int, int]]:
        if self._introns is None:
            if len(self._exons) <= 1:
                self._introns = []
            else:
                self._introns = [self.start, *(self.exons() | flatten), self.end] | window(2, 2) | skip_while(lambda x: x[0] == x[1]) | to_list
        return self._introns

    def tss_region(self, region: Tuple[int, int] = (-1000, 1000)) -> Tuple[int, ...]:
        if self.strand == '+':
            region = tuple(i+self.start for i in region)
        else:
            region = tuple(reversed(tuple(self.end-i for i in region)))
        return region

    def downstream(self, down: int = 5000) -> Tuple[int, int]:
        if self.strand == '+':
            downstream = (self.end, self.end+down)
        else:
            downstream = (self.start-down, self.start)
        return downstream

    def _classify_exons(self):
        if self.cds_start is not None and self.cds_end is not None:
            utr5_exons = []
            utr3_exons = []
            cds_exons = []
            for exon in self._exons:
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
            loginfo = f"exons: {self.exons()}\ncds_start: {self.cds_start}, cds_end: {self.cds_end}"
            if not cds_exons:
                warnings.warn(f"utr5_exons of {self.transcript_id} does not exist, please check:\n{loginfo}")
            self._utr5_exons = utr5_exons
            self._cds_exons = cds_exons
            self._utr3_exons = utr3_exons
        else:
            self._utr5_exons = []
            self._cds_exons = []
            self._utr3_exons = []

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

    def annotation(self, blocks: List[Tuple[int, int]], tss_region: Tuple[int, int] = (-1000, 1000), downstream: int = 3000) -> GenomicAnnotation:
        return GenomicAnnotation(self.transcript_id, self.gene_name, self.start, self.end, self.transcript_type, self.anno(blocks, tss_region, downstream))
