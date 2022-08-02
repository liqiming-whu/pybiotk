# -*- coding: utf-8 -*-
import itertools
import random
import re
import sys
from typing import Dict, List, Sequence, Iterable, Iterator, Tuple, Literal, Callable, Optional

import pysam

from pybiotk.intervals import GRange
from pybiotk.utils import reverse_seq


class FastaFile(pysam.FastaFile):
    def __init__(self, /, *args, **kwargs):
        super().__init__()
        self.reference_dict: Optional[Dict[str, str]] = None

    def __iter__(self) -> Iterator:
        for reference in self.references:
            yield reference, self[reference]

    def stdout(self, referenceList: Optional[Sequence[str]] = None, wrap: bool = True):
        if not referenceList:
            referenceList = self.references

        for reference in referenceList:
            if reference not in self.references:
                continue
            sys.stdout.write(f">{reference}\n")
            sequence = self[reference]
            if wrap:
                seqlen = len(sequence)
                index = range(seqlen)
                for i in itertools.islice(index, None, None, 60):
                    sys.stdout.write(f"{sequence[i: i+60]}\n")
            else:
                sys.stdout.write(f"{sequence}\n")

    def load_into_dict(self) -> Dict[str, str]:
        self.reference_dict = dict((reference, seq) for reference, seq in self)
        return self.reference_dict

    def fetchs(self, reference: str, start: int, end: int, strand: Literal["+", "-"] = "+"):
        sequence = self.fetch(reference, int(start), int(end))
        if strand == "-":
            sequence = reverse_seq(sequence)
        return sequence

    @staticmethod
    def fetch_use_dict(reference_dict: Dict[str, str], reference: str, start: int, end: int, strand: Literal["+", "-"] = "+"):
        sequence = reference_dict[reference][int(start):int(end)]
        if strand == "-":
            sequence = reverse_seq(sequence)
        return sequence

    def dict_fetch(self, reference: str, start: int, end: int, strand: Literal["+", "-"] = "+"):
        assert self.reference_dict is not None
        return self.fetch_use_dict(self.reference_dict, reference, start, end, strand)


class GenomeFile(FastaFile):
    def __init__(self, /, *args, **kwargs):
        super().__init__()
        self.custom_chroms: List[str] = self.chroms_norm()

    get_chrom_length: Callable[[str], int] = FastaFile.get_reference_length
    chroms: List[str] = FastaFile.references

    def fetch_blocks(self, chrom: str, blocks: Iterable[Tuple[int, int]], strand: Literal["+", "-"] = "+") -> str:
        sequence = ""
        chromLen = self.get_chrom_length(chrom)
        for block in blocks:
            start = int(block[0]) if int(block[0]) > 0 else 0
            end = int(block[1]) if int(block[1]) < chromLen else chromLen
            if end <= start:
                continue
            sequence += self.fetch(chrom, start, end)
        seq = sequence if strand == '+' else reverse_seq(sequence)

        return seq

    def chroms_norm(self) -> List[str]:
        chroms = []
        for chrom in self.chroms:
            if chrom.startswith("chr") or re.match(r"\d.*", chrom):
                chroms.append(chrom)
        return chroms

    def random(self, length: int = 50, random_on_chroms: Optional[Sequence[str]] = None) -> Tuple[GRange, str]:
        strand = ["+", "-"][random.randint(0, 1)]
        chroms = self.custom_chroms
        if random_on_chroms:
            chroms = list(set(chroms) & set(random_on_chroms))
        chrom = chroms[random.randint(0, len(chroms)-1)]
        chrom_len = self.get_chrom_length(chrom)
        start = random.randint(0, chrom_len-length-1)
        end = start + 50
        if self.reference_dict:
            seq = self.dict_fetch(chrom, start, end, strand)
        else:
            seq = self.fetch_blocks(chrom, [(start, end)], strand).upper()
        gRange = GRange(chrom, start, end, strand)
        return gRange, seq

    def stdout(self, chromList: Optional[Sequence[str]] = None, wrap: bool = True):
        if not chromList:
            chromList = self.custom_chroms
        return super().stdout(chromList, wrap)
