# -*- coding: utf-8 -*-
import re
from dataclasses import dataclass, field
from io import TextIOWrapper
from typing import List, Sequence, Tuple, Literal, Iterable, Iterator, Optional, Union, TextIO

from pybiotk.annodb import Transcript
from pybiotk.utils import logging
from pybiotk.io.bed import Bed6, Bed12, Intron, GeneInfo, TransInfo
from stream.pipe import Pipe, sort, drop_while, filter, kgroupby, groupby, apply, window


@dataclass
class GTF:
    seqname: str = field(default=None)
    source: str = field(default=None, repr=False)
    feature: str = field(default=None)
    start: int = field(default=0)
    end: int = field(default=0)
    score: str = field(default=None, repr=False)
    strand: Literal['+', '-'] = field(default='+')
    frame: str = field(default=None, repr=False)
    attributes: str = field(default=None, repr=False)
    chrom: str = field(init=False, repr=False)

    def __post_init__(self):
        self.chrom = self.seqname
        self.start = int(self.start)
        self.end = int(self.end)

    def __str__(self):
        return "\t".join(str(s) for s in list(self.__dict__.values())[:9])

    @staticmethod
    def parse_attributes(term: str, attributes: str) -> Optional[str]:
        parse_attr = re.compile(f'{term} "([^"]+)"')
        patterns = parse_attr.findall(attributes)

        if patterns:
            return patterns[0]

    def gene_type(self):
        for term in ['gene_type', 'gene_biotype']:
            attr = self.parse_attributes(term, self.attributes)
            if attr is not None:
                return attr

    def gene_id(self):
        return self.parse_attributes("gene_id", self.attributes)

    def gene_name(self):
        for term in ["gene_name", "gene"]:
            attr = self.parse_attributes(term, self.attributes)
            if attr is not None:
                return attr
        else:
            return self.gene_id()

    def transcript_type(self):
        for term in ["transcript_type", "transcript_biotype"]:
            attr = self.parse_attributes(term, self.attributes)
            if attr is not None:
                return attr

    def transcript_id(self):
        attr = self.parse_attributes("transcript_id", self.attributes)
        if attr is None:
            attr = self.gene_id()
        return attr

    def transcript_name(self):
        attr = self.parse_attributes("transcript_name", self.attributes)
        if attr is None:
            attr = self.gene_name()
        return attr

    def get_attribute(self, attribute):
        if attribute == "gene_type":
            return self.gene_type()
        elif attribute == "gene_name":
            return self.gene_name()
        elif attribute == "transcript_id":
            return self.transcript_id()
        elif attribute == "transcript_type":
            return self.transcript_type()
        elif attribute == "transcript_name":
            return self.transcript_name()
        else:
            return self.parse_attributes(attribute, self.attributes)


@Pipe
def to_GTF(iterable: Iterable[str]) -> Iterator[GTF]:
    for line in iterable:
        yield GTF(*line.rstrip("\r\n").split("\t"))


@Pipe
def to_Bed6(iterable: Iterable[Tuple[str, Tuple[GTF, ...]]]) -> Iterator[Bed6]:
    for (k, gtfs) in iterable:
        for gtf in gtfs:
            yield Bed6.init_by_gtf(gtf, k)


@Pipe
def to_GeneInfo(iterable: Iterable[Tuple[GTF, ...]]) -> Iterator[GeneInfo]:
    for gtfs in iterable:
        for gtf in gtfs:
            yield GeneInfo.init_by_gtf(gtf)


@Pipe
def to_TransInfo(iterable: Iterable[Tuple[GTF, ...]]) -> Iterator[TransInfo]:
    for group in iterable:
        gene_gtf: Sequence[GTF] = []
        trans_gtf: Sequence[GTF] = []
        for gtf in group:
            if gtf.feature == "gene":
                gene_gtf.append(gtf)
            else:
                trans_gtf.append(gtf)
        for gtf in trans_gtf:
            trans_info = TransInfo.init_by_gtf(gtf)
            if gene_gtf :
                if trans_info.gene_type is None:
                    trans_info.gene_type = gene_gtf[0].gene_type()
            yield trans_info


@Pipe
def to_Bed12(iterable: Iterable[Tuple[GTF, ...]], name: str = "transcript_id") -> Iterator[Bed12]:
    for group in iterable:
        gtfs: Iterator[GTF] = group | sort(key=lambda x: x.start)
        cds_exons: Sequence[GTF] = []
        bed: Optional[Bed12] = None
        last_gtf = None
        for gtf in gtfs:
            last_gtf = gtf
            if gtf.feature in {'CDS', 'stop_codon'}:
                cds_exons.append(gtf)
            elif gtf.feature == "exon":
                if bed is None:
                    bed = Bed12.init_by_gtf(gtf, gtf.get_attribute(name))
                else:
                    bed.update(gtf)
        if bed is None:
            logging.warning(f"No exon found for transcript: {last_gtf.get_attribute(name)}")
            continue
        if cds_exons:
            bed.thickStart = cds_exons[0].start - 1
            bed.thickEnd = cds_exons[-1].end
        else:
            bed.thickStart = bed.thickEnd = bed.end
        yield bed


@Pipe
def to_Intron(iterable: Iterable[Tuple[GTF, ...]], name: str = "transcript_id") -> Iterator[Intron]:
    for group in iterable:
        gtfs = group | sort(key=lambda x: x.start)
        for exons in gtfs | window:
            if len(exons) < 2:
                continue
            yield Intron.init_by_gtf(exons[0], exons[1], exons[0].get_attribute(name))


@Pipe
def to_Transcript(iterable: Iterable[Tuple[GTF, ...]]) -> Iterator[Transcript]:
    for group in iterable:
        gene_gtf: Sequence[GTF] = []
        exons_gtf: Sequence[GTF] = []
        for gtf in group:
            if gtf.feature == "gene":
                gene_gtf.append(gtf)
            else:
                exons_gtf.append(gtf)

        for group in exons_gtf | groupby(key=lambda x: x.transcript_id()):
            gtfs: Iterator[GTF] = group | sort(key=lambda x: x.start)
            cds_exons: Sequence[GTF] = []
            transcript: Optional[Transcript] = None
            for gtf in gtfs:
                if gtf.feature in {'CDS', 'stop_codon'}:
                    cds_exons.append(gtf)
                elif gtf.feature == "exon":
                    if transcript is None:
                        transcript = Transcript.init_by_gtf(gtf)
                    else:
                        transcript.update(gtf)
            if transcript is None:
                logging.warning(f"No exon found for transcript: {gtf.get_attribute('transcript_id')}")
                continue
            if cds_exons:
                cds_start = cds_exons[0].start - 1
                cds_end = cds_exons[-1].end
                transcript.cds_start = cds_start
                transcript.cds_end = cds_end
            else:
                transcript.cds_start = transcript.cds_end = None
    
            if gene_gtf:
                if transcript.gene_type is None:
                    transcript.gene_type = gene_gtf[0].gene_type()
            if transcript.transcript_type is None:
                transcript.transcript_type = transcript.gene_type

            yield transcript


class GtfFile:
    def __init__(self, filepath_or_buffer: Union[str, TextIO], comment: str = "#"):
        if isinstance(filepath_or_buffer, TextIOWrapper):
            self.filename: str = filepath_or_buffer.name
            self.gtf: TextIO = filepath_or_buffer
        else:
            self.filename: str = filepath_or_buffer
            self.gtf: TextIO = open(filepath_or_buffer)
        self.comment: str = comment
        self.gtf_list: Optional[Sequence[GTF]] = None
        self.ptr: int = 0
        self.closed: bool = False

    def __iter__(self) -> Iterator[GTF]:
        if self.gtf_list is not None:
            return iter(self.gtf_list)
        else:
            self.gtf_list = []
            iterator = self.gtf | drop_while(lambda x: x.startswith(self.comment)) | to_GTF | apply(lambda x: self.gtf_list.append(x))
            return iterator

    def __next__(self) -> GTF:
        if self.gtf_list is None:
            self.to_list()
        if self.ptr < len(self.gtf_list):
            gtf = self.gtf_list[self.ptr]
            self.ptr += 1
        else:
            raise StopIteration

        return gtf

    def close(self):
        if not self.closed:
            self.gtf.close()
            self.closed = True
            self.gtf_list = None
            self.ptr = 0

    def to_list(self) -> Optional[List[GTF]]:
        if self.gtf_list is None:
            _ = list(self)
        return self.gtf_list

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, traceback):
        self.close()

    def iter(self,
             features: Optional[Sequence[str]] = None,
             gene_types: Optional[Sequence[str]] = None,
             transcript_types: Optional[Sequence[str]] = None,
             transcript_ids: Optional[Sequence[str]] = None,
             transcript_names: Optional[Sequence[str]] = None,
             gene_ids: Optional[Sequence[str]] = None,
             gene_names: Optional[Sequence[str]] = None,
             ) -> Iterator[GTF]:

        iterator = iter(self)
        if features is not None:
            features = set(features)
            iterator = iterator | filter(lambda x: x.feature in features)
        if gene_types is not None:
            gene_types = set(gene_types)
            iterator = iterator | filter(
                lambda x: x.gene_type() in gene_types)
        if transcript_types is not None:
            transcript_types = set(transcript_types)
            iterator = iterator | filter(
                lambda x: x.transcript_type() in transcript_types)
        if transcript_ids is not None:
            transcript_ids = set(transcript_ids)
            iterator = iterator | filter(lambda x: x.transcript_id() in transcript_ids)
        if transcript_names is not None:
            transcript_names = set(transcript_names)
            iterator = iterator | filter(lambda x: x.transcript_name() in transcript_names)
        if gene_ids is not None:
            gene_ids = set(gene_ids)
            iterator = iterator | filter(lambda x: x.gene_id() in gene_ids)
        if gene_names is not None:
            gene_names = set(gene_names)
            iterator = iterator | filter(lambda x: x.gene_name() in gene_names)

        return iterator

    def iter_exon(self,
                  gene_types: Optional[Sequence[str]] = None,
                  transcript_types: Optional[Sequence[str]] = None,
                  transcript_ids: Optional[Sequence[str]] = None,
                  transcript_names: Optional[Sequence[str]] = None,
                  gene_ids: Optional[Sequence[str]] = None,
                  gene_names: Optional[Sequence[str]] = None,
                  ) -> Iterator[GTF]:
        return self.iter(features=["exon"],
                         gene_types=gene_types,
                         transcript_types=transcript_types,
                         transcript_ids=transcript_ids,
                         transcript_names=transcript_names,
                         gene_ids=gene_ids,
                         gene_names=gene_names
                         )

    def iter_gene(self,
                  gene_types: Optional[Sequence[str]] = None,
                  transcript_types: Optional[Sequence[str]] = None,
                  transcript_ids: Optional[Sequence[str]] = None,
                  transcript_names: Optional[Sequence[str]] = None,
                  gene_ids: Optional[Sequence[str]] = None,
                  gene_names: Optional[Sequence[str]] = None,
                  ) -> Iterator[GTF]:
        return self.iter(features=["gene"],
                         gene_types=gene_types,
                         transcript_types=transcript_types,
                         transcript_ids=transcript_ids,
                         transcript_names=transcript_names,
                         gene_ids=gene_ids,
                         gene_names=gene_names
                         )

    def iter_transcript(self,
                        gene_types: Optional[Sequence[str]] = None,
                        transcript_types: Optional[Sequence[str]] = None,
                        transcript_ids: Optional[Sequence[str]] = None,
                        transcript_names: Optional[Sequence[str]] = None,
                        gene_ids: Optional[Sequence[str]] = None,
                        gene_names: Optional[Sequence[str]] = None,
                        ) -> Iterator[GTF]:
        return self.iter(features=["transcript"],
                         gene_types=gene_types,
                         transcript_types=transcript_types,
                         transcript_ids=transcript_ids,
                         transcript_names=transcript_names,
                         gene_ids=gene_ids,
                         gene_names=gene_names
                         )

    def to_bed6(self, feature: Literal["exon", "gene"] = "gene",
                name: Literal["gene_id", "gene_name", "transcript_id", "transcript_name"] = "gene_name",
                gene_types: Optional[Sequence[str]] = None,
                transcript_types: Optional[Sequence[str]] = None,
                transcript_ids: Optional[Sequence[str]] = None,
                transcript_names: Optional[Sequence[str]] = None,
                gene_ids: Optional[Sequence[str]] = None,
                gene_names: Optional[Sequence[str]] = None
                ) -> Iterator[Bed6]:

        bed6 = self.iter(features=[feature],
                         gene_types=gene_types,
                         transcript_types=transcript_types,
                         transcript_ids=transcript_ids,
                         transcript_names=transcript_names,
                         gene_ids=gene_ids,
                         gene_names=gene_names
                         ) | kgroupby(lambda x: x.get_attribute(name)) | to_Bed6
        return bed6

    def to_gene_info(self,
                     gene_types: Optional[Sequence[str]] = None,
                     transcript_types: Optional[Sequence[str]] = None,
                     transcript_ids: Optional[Sequence[str]] = None,
                     transcript_names: Optional[Sequence[str]] = None,
                     gene_ids: Optional[Sequence[str]] = None,
                     gene_names: Optional[Sequence[str]] = None,
                     ) -> Iterator[GeneInfo]:
        gene_info = self.iter_gene(
            gene_types=gene_types,
            transcript_types=transcript_types,
            transcript_ids=transcript_ids,
            transcript_names=transcript_names,
            gene_ids=gene_ids,
            gene_names=gene_names
            ) | groupby(lambda x: x.gene_id()) | to_GeneInfo
        return gene_info

    def to_trans_info(self,
                      gene_types: Optional[Sequence[str]] = None,
                      transcript_types: Optional[Sequence[str]] = None,
                      transcript_ids: Optional[Sequence[str]] = None,
                      transcript_names: Optional[Sequence[str]] = None,
                      gene_ids: Optional[Sequence[str]] = None,
                      gene_names: Optional[Sequence[str]] = None,
                      ) -> Iterator[TransInfo]:
        trans_info = self.iter(
            features=["gene", "transcript"],
            gene_types=gene_types,
            transcript_types=transcript_types,
            transcript_ids=transcript_ids,
            transcript_names=transcript_names,
            gene_ids=gene_ids,
            gene_names=gene_names
            ) | groupby(lambda x: x.gene_id()) | to_TransInfo
        return trans_info

    def to_bed12(self, name: Literal["gene_id", "gene_name", "transcript_id", "transcript_name"] = "transcript_id",
                 gene_types: Optional[Sequence[str]] = None,
                 transcript_types: Optional[Sequence[str]] = None,
                 transcript_ids: Optional[Sequence[str]] = None,
                 transcript_names: Optional[Sequence[str]] = None,
                 gene_ids: Optional[Sequence[str]] = None,
                 gene_names: Optional[Sequence[str]] = None,
                 ) -> Iterator[Bed12]:
        bed12 = self.iter(features=["CDS", 'stop_codon', 'exon'],
                          gene_types=gene_types,
                          transcript_types=transcript_types,
                          transcript_ids=transcript_ids,
                          transcript_names=transcript_names,
                          gene_ids=gene_ids,
                          gene_names=gene_names
                          ) | groupby(lambda x: x.transcript_id()) | to_Bed12(name=name)
        return bed12

    def to_intron(self, name: Literal["gene_id", "gene_name", "transcript_id", "transcript_name"] = "transcript_id",
                  gene_types: Optional[Sequence[str]] = None,
                  transcript_types: Optional[Sequence[str]] = None,
                  transcript_ids: Optional[Sequence[str]] = None,
                  transcript_names: Optional[Sequence[str]] = None,
                  gene_ids: Optional[Sequence[str]] = None,
                  gene_names: Optional[Sequence[str]] = None,
                  ) -> Iterator[Intron]:
        intron = self.iter_exon(
            gene_types=gene_types,
            transcript_types=transcript_types,
            transcript_ids=transcript_ids,
            transcript_names=transcript_names,
            gene_ids=gene_ids,
            gene_names=gene_names
            ) | groupby(lambda x: x.transcript_id()) | to_Intron(name=name)
        return intron

    def to_transcript(self,
                      gene_types: Optional[Sequence[str]] = None,
                      transcript_types: Optional[Sequence[str]] = None,
                      transcript_ids: Optional[Sequence[str]] = None,
                      transcript_names: Optional[Sequence[str]] = None,
                      gene_ids: Optional[Sequence[str]] = None,
                      gene_names: Optional[Sequence[str]] = None,
                      ) -> Iterator[Transcript]:
        transcript = self.iter(features=["gene", "CDS", 'stop_codon', 'exon'],
                               gene_types=gene_types,
                               transcript_types=transcript_types,
                               transcript_ids=transcript_ids,
                               transcript_names=transcript_names,
                               gene_ids=gene_ids,
                               gene_names=gene_names
                               ) | groupby(lambda x: x.gene_id()) | to_Transcript
        return transcript
