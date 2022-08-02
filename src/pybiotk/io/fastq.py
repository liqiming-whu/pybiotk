# -*- coding: utf-8 -*-
import gzip
from typing import Iterator, Tuple, Optional

import pysam


class FastqFile(pysam.FastxFile):
    def __init__(self, /,*args, **kwargs):
        super().__init__()
        self.ptr: Optional[int] = None

    def to_fasta(self) -> Iterator[str]:
        for entry in self:
            yield f">{entry.name}\n{entry.sequence}"

    def uniq(self) -> Iterator[pysam.libcfaidx.FastxRecord]:
        self.ptr = 0
        sequence_set = set()
        for entry in self:
            self.ptr += 1
            if entry.sequence not in sequence_set:
                sequence_set.add(entry.sequence)
                yield entry

    def iter_len(self) -> Iterator[int]:
        for entry in self:
            yield len(entry.sequence)


class FastqPair:
    def __init__(self, read1: str, read2: str):
        self.filename1 = read1
        self.filename2 = read2
        self.read1 = FastqFile(read1)
        self.read2 = FastqFile(read2)
        self.ptr: Optional[int] = None

    def __iter__(self) -> Iterator[Tuple[pysam.libcfaidx.FastxRecord, ...]]:
        for entry1, entry2 in zip(self.read1, self.read2):
            yield entry1, entry2

    def uniq(self) -> Iterator[Tuple[pysam.libcfaidx.FastxRecord, ...]]:
        self.ptr = 0
        unique = set()
        for entry1, entry2 in zip(self.read1, self.read2):
            if not entry1.name == entry2.name:
                raise RuntimeError(f"{entry1.name} != {entry2.name}")
            self.ptr += 1
            whole_reads = entry1.sequence + entry2.sequence
            if whole_reads not in unique:
                unique.add(whole_reads)
                yield entry1, entry2

    def close(self):
        self.read1.close()
        self.read2.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, traceback):
        self.close()


class OpenFqGzip:
    def __init__(self, filename, mode="wb", compresslevel=9, encoding=None, errors=None, newline=None):
        self.filename = filename
        self.mode = mode
        self.compresslevel = compresslevel
        self.encoding = encoding
        self.errors = errors
        self.newline = newline
        self.fq = gzip.open(filename, mode, compresslevel, encoding, errors, newline)
        self.name = self.filename

    def write_entry(self, name: str, sequence: str, comment: Optional[str] = None, quality: Optional[str] = None):
        if quality is None:
            quality = "F"*len(sequence)
        if comment is None:
            self.fq.write(f"@{name}\n{sequence}\n+\n{quality}\n".encode("utf-8"))
        else:
            self.fq.write(f"@{name} {comment}\n{sequence}\n+\n{quality}\n".encode("utf-8"))

    def write_fastx_record(self, fq: pysam.libcfaidx.FastxRecord):
        self.write_entry(fq.name, fq.comment, fq.sequence, fq.quality)

    def write(self, string: str):
        self.fq.write(string)

    def close(self):
        self.fq.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, traceback):
        self.fq.close()
