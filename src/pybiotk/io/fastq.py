# -*- coding: utf-8 -*-
import gzip
import sys
from collections import defaultdict
from itertools import zip_longest
from typing import Iterator, Literal, Optional, Tuple

import pysam


class FastxFile(pysam.FastxFile):
    def __init__(self, /, *args, **kwargs):
        super().__init__()
        self.ptr: Optional[int] = None

    def to_fasta(self) -> Iterator[str]:
        for entry in self:
            yield f">{entry.name}\n{entry.sequence}"

    def to_fastq(self, path: str = "-"):
        with OpenFqGzip(path) as fq:
            for entry in self:
                fq.write_fastx_record(entry)

    def uniq(self, by: Literal["id", "name", "seq"] = "seq") -> Iterator[pysam.libcfaidx.FastxRecord]:
        if by not in {"id", "name", "seq"}:
            raise ValueError(f"unsupported uniqueness key: {by}")
        self.ptr = 0
        unique = set()
        for entry in self:
            self.ptr += 1
            if by == "seq":
                key = entry.sequence
            elif by == "id":
                key = entry.name
            else:
                key = (entry.name, entry.comment)
            if key not in unique:
                unique.add(key)
                yield entry

    def rename(self) -> Iterator[pysam.libcfaidx.FastxRecord]:
        self.ptr = 0
        used_names = set()
        next_suffix = defaultdict(lambda: 2)
        for entry in self:
            self.ptr += 1
            base_name = entry.name
            output_name = base_name
            if output_name in used_names:
                suffix = next_suffix[base_name]
                output_name = f"{base_name}_{suffix}"
                while output_name in used_names:
                    suffix += 1
                    output_name = f"{base_name}_{suffix}"
                next_suffix[base_name] = suffix + 1
            entry.name = output_name
            used_names.add(output_name)
            yield entry

    def iter_len(self) -> Iterator[int]:
        for entry in self:
            yield len(entry.sequence)


class FastqFile(FastxFile):...


class FastqPair:
    def __init__(self, read1: str, read2: str):
        self.filename1 = read1
        self.filename2 = read2
        self.read1 = FastqFile(read1)
        self.read2 = FastqFile(read2)
        self.ptr: Optional[int] = None

    def __iter__(self) -> Iterator[Tuple[pysam.libcfaidx.FastxRecord, ...]]:
        sentinel = object()
        for entry1, entry2 in zip_longest(self.read1, self.read2, fillvalue=sentinel):
            if entry1 is sentinel or entry2 is sentinel:
                raise RuntimeError("paired FASTQ files contain different numbers of reads")
            read1_name = entry1.name[:-2] if entry1.name.endswith("/1") else entry1.name
            read2_name = entry2.name[:-2] if entry2.name.endswith("/2") else entry2.name
            if read1_name != read2_name:
                raise RuntimeError(f"unmatched read names: {entry1.name} != {entry2.name}")
            yield entry1, entry2

    def uniq(self, by: Literal["id", "name", "seq"] = "seq") -> Iterator[Tuple[pysam.libcfaidx.FastxRecord, ...]]:
        if by not in {"id", "name", "seq"}:
            raise ValueError(f"unsupported uniqueness key: {by}")
        self.ptr = 0
        unique = set()
        for entry1, entry2 in self:
            self.ptr += 1
            if by == "seq":
                key = (entry1.sequence, entry2.sequence)
            elif by == "id":
                key = (entry1.name, entry2.name)
            else:
                key = (entry1.name, entry1.comment, entry2.name, entry2.comment)
            if key not in unique:
                unique.add(key)
                yield entry1, entry2

    def close(self):
        self.read1.close()
        self.read2.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, traceback):
        self.close()


class OpenFqGzip:
    def __init__(self, filename, mode="wb", compresslevel=4, encoding=None, errors=None, newline=None):
        self.filename = filename
        self.mode = mode
        self.compresslevel = compresslevel
        self.encoding = encoding
        self.errors = errors
        self.newline = newline
        if filename == "-":
            self.fq = sys.stdout
        else:
            self.fq = gzip.open(filename, mode, compresslevel, encoding, errors, newline)
        self.name = self.filename

    def _write(self, record: str):
        if self.filename == "-":
            self.fq.write(record)
        else:
            self.fq.write(record.encode("utf-8"))

    def write_entry(self, name: str, sequence: str, comment: Optional[str] = None, quality: Optional[str] = None):
        if quality is None:
            quality = "F"*len(sequence)
        if comment is None:
            self._write(f"@{name}\n{sequence}\n+\n{quality}\n")
        else:
            self._write(f"@{name} {comment}\n{sequence}\n+\n{quality}\n")

    def write_fastx_record(self, fq: pysam.libcfaidx.FastxRecord):
        self.write_entry(fq.name, fq.sequence, fq.comment, fq.quality)

    def write(self, string: str):
        self.fq.write(string)

    def close(self):
        self.fq.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, traceback):
        self.fq.close()
