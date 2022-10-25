# -*- coding: utf-8 -*-
import warnings
from typing import List, Tuple, Literal, Iterable, Optional

from pybiotk.annodb.anno import GFeature, GenomicAnnotation
from pybiotk.annodb.transcript import Transcript
from pybiotk.intervals import merge_intervals
from stream import flatten, window, skip_while, to_list


class Gene(GFeature):
    def __init__(
        self,
        gene_id: Optional[str] = None,
        gene_name: Optional[str] = None,
        gene_type: Optional[str] = None,
        chrom: Optional[str] = None,
        start: int = 0,
        end: int = 0,
        strand: Optional[Literal['+', '-']] = None,
        cds_start: Optional[int] = None,
        cds_end: Optional[int] = None,
        exons: List[Tuple[int, int]] = None,
        transcripts: List[Transcript] = None,
    ):
        self.gene_id = gene_id
        self.gene_name = gene_name
        self.gene_type = gene_type
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.cds_start = int(cds_start) if cds_start is not None else cds_start
        self.cds_end = int(cds_end) if cds_end is not None else cds_end
        self._exons = exons if exons is not None else []
        self.transcripts = transcripts if transcripts is not None else []
        self._introns = None
        self._cds_exons = None
        self._utr5_exons = None
        self._utr3_exons = None
        self.id = self.gene_id

    def __repr__(self) -> str:
        return f"{self.gene_id}:{self.chrom}:{self.start}-{self.end}({self.strand})"

    __str__ = __repr__

    def add_transcript(self, transcript: Transcript):
        if not hasattr(self, "transcripts"):
            self.transcripts = []
        else:
            self.transcripts = list(self.transcripts)
        self.transcripts.append(transcript)

    def post_init(self):
        if self.transcripts:
            min_start = float("inf")
            max_end = 0
            min_cds_start = float("inf")
            max_cds_end = 0
            for transcript in self.transcripts:
                if self.gene_id is None:
                    self.gene_id = transcript.gene_id
                    self.gene_name = transcript.gene_name
                    self.gene_type = transcript.gene_type
                    self.chrom = transcript.chrom
                    self.strand = transcript.strand
                self._exons.extend(transcript.exons())
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

            self._exons = list(set(self._exons))
            del self.transcripts
            self.transcripts = []
            self.start = min_start
            self.end = max_end
            if min_cds_start == float("inf"):
                min_cds_start = None
            if max_cds_end == 0:
                max_cds_end = None
            self.cds_start = min_cds_start
            self.cds_end = max_cds_end
        else:
            warnings.warn(
                f"gene({self.gene_id})'s transcripts list is empty, nothing will be done. ")

    @classmethod
    def init_by_transcripts(cls, transcripts: Iterable[Transcript]):
        min_start = float("inf")
        max_end = 0
        min_cds_start = float("inf")
        max_cds_end = 0
        gene_id = None
        gene_name = None
        gene_type = None
        chrom = None
        strand = None
        exons = []
        for transcript in transcripts:
            if gene_id is None:
                gene_id = transcript.gene_id
                gene_name = transcript.gene_name
                gene_type = transcript.gene_type
                chrom = transcript.chrom
                strand = transcript.strand
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
        exons = list(set(exons))
        start = min_start
        end = max_end
        if min_cds_start == float("inf"):
            min_cds_start = None
        if max_cds_end == 0:
            max_cds_end = None
        cds_start = min_cds_start
        cds_end = max_cds_end
        return cls(gene_id, gene_name, gene_type, chrom, start, end, strand, cds_start, cds_end, exons)

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
        if not self._exons:
            self.post_init()
        self._exons = merge_intervals(self._exons)
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

    def downstream(self, down: int = 5000) -> Tuple[int, int]:
        if self.strand == '+':
            downstream = (self.end, self.end+down)
        else:
            downstream = (self.start-down, self.start)
        return downstream

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
        return GenomicAnnotation(self.gene_id, self.gene_name, self.start, self.end, self.gene_type, self.anno(blocks, tss_region, downstream))
