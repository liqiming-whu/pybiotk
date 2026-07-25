#!/usr/bin/env python3
"""Rename FASTA/FASTQ records in single-end or paired-end mode.

``index`` is the default and assigns compact globally unique names while
preserving record order. Paired-end index mode reads both mates together and
assigns the same name to each pair using constant renaming memory. ``preserve``
keeps original names and adds a numeric suffix to repeated names.

.. warning::

   When ``index_base`` is set to 16 or 36, the generated read names contain
   hex-like or base36-like segments (e.g. ``read_1e00d8``).  Tools that sort
   read names with *natural* ordering — most notably ``samtools sort -n`` —
   parse embedded digits numerically, producing a different order than ASCII
   lexicographic (Python-``str`` / ``samtools sort -N``).  Mismatched sort
   order can break downstream steps that rely on contiguous ``groupby``
   semantics (e.g. chimeric-read splitting).  Prefer the default decimal
   encoding (``index_base=10``); if compact names are required, pair them
   with ``samtools sort -N`` instead of ``-n``.
"""
import argparse
import sys
import time
from collections import defaultdict
from contextlib import nullcontext
from typing import Literal, Optional, Sequence, Tuple

from pybiotk.io import FastqPair, FastxFile, OpenFqGzip
from pybiotk.utils import configure_logging, get_logger
from pybiotk.utils import ignore


OutputFormat = Literal["fastq", "fasta"]
RenameMode = Literal["index", "preserve"]
_DIGITS = "0123456789abcdefghijklmnopqrstuvwxyz"


logger = get_logger(__name__)


def _format_index(index: int, base: int) -> str:
    if base not in {10, 16, 36}:
        raise ValueError("index_base must be 10, 16, or 36")
    digits = []
    while index:
        index, remainder = divmod(index, base)
        digits.append(_DIGITS[remainder])
    return "".join(reversed(digits)) or "0"


def _validate_options(outfmt: str, mode: str, compresslevel: int, prefix: Optional[str]) -> None:
    if outfmt not in {"fastq", "fasta"}:
        raise ValueError(f"unsupported output format: {outfmt}")
    if mode not in {"index", "preserve"}:
        raise ValueError(f"unsupported rename mode: {mode}")
    if not 0 <= compresslevel <= 9:
        raise ValueError("compresslevel must be between 0 and 9")
    if prefix is not None and (not prefix or any(character.isspace() for character in prefix)):
        raise ValueError("prefix must be non-empty and contain no whitespace")


def _output_stream(path: str, outfmt: OutputFormat, compresslevel: int):
    if outfmt == "fastq":
        return OpenFqGzip(path, compresslevel=compresslevel)
    if path == "-":
        return nullcontext(sys.stdout)
    return open(path, "w")


def _write_record(stream, record, outfmt: OutputFormat) -> None:
    if outfmt == "fastq":
        stream.write_fastx_record(record)
    else:
        stream.write(f">{record.name}\n{record.sequence}\n")


def _mate_suffix(name: str) -> Tuple[str, str]:
    if name.endswith(("/1", "/2")):
        return name[:-2], name[-2:]
    return name, ""


def fastx_rename(input_fq: str, output: str, outfmt: OutputFormat = "fastq",
                 mode: RenameMode = "index", prefix: Optional[str] = None, index_base: int = 10,
                 compresslevel: int = 4) -> None:
    _validate_options(outfmt, mode, compresslevel, prefix)
    _format_index(1, index_base)
    input_str = "stdin" if input_fq == "-" else input_fq
    logger.info(f"Processing {input_str}, mode {mode} ...")
    start = time.perf_counter()
    input_reads = 0
    with FastxFile(input_fq) as fqi, _output_stream(output, outfmt, compresslevel) as fqo:
        if mode == "preserve":
            records = fqi.rename()
        else:
            records = fqi
        for input_reads, fq in enumerate(records, start=1):
            if mode == "index":
                base = prefix if prefix else fq.name
                fq.name = f"{base}_{_format_index(input_reads, index_base)}"
            _write_record(fqo, fq, outfmt)
    logger.info(f"Processed {input_reads} reads in {time.perf_counter() - start:.2f} seconds.")


def fastx_rename_pair(read1_files: Sequence[str], read2_files: Sequence[str], output1: str,
                      output2: str, outfmt: OutputFormat = "fastq",
                      mode: RenameMode = "index", prefix: Optional[str] = None, index_base: int = 10,
                      compresslevel: int = 4) -> None:
    _validate_options(outfmt, mode, compresslevel, prefix)
    _format_index(1, index_base)
    if isinstance(read1_files, str):
        read1_files = [read1_files]
    if isinstance(read2_files, str):
        read2_files = [read2_files]
    if not read1_files or len(read1_files) != len(read2_files):
        raise ValueError("read1 and read2 must contain the same non-zero number of files")
    if output1 == "-" and output2 == "-":
        raise ValueError("paired-end outputs cannot both use stdout")

    start = time.perf_counter()
    pair_count = 0
    used_names = set()
    next_suffix = defaultdict(lambda: 2)
    with _output_stream(output1, outfmt, compresslevel) as out1, \
            _output_stream(output2, outfmt, compresslevel) as out2:
        for read1, read2 in zip(read1_files, read2_files):
            with FastqPair(read1, read2) as pairs:
                for fq1, fq2 in pairs:
                    pair_count += 1
                    if mode == "index":
                        if prefix:
                            base_name = prefix
                        else:
                            base_name = fq1.name[:-2] if fq1.name.endswith(("/1", "/2")) else fq1.name
                        name = f"{base_name}_{_format_index(pair_count, index_base)}"
                        fq1.name = name
                        fq2.name = name
                    else:
                        base1, suffix1 = _mate_suffix(fq1.name)
                        base2, suffix2 = _mate_suffix(fq2.name)
                        if base1 != base2:
                            raise RuntimeError(f"unmatched read names: {fq1.name} != {fq2.name}")
                        output_base = base1
                        if output_base in used_names:
                            suffix = next_suffix[base1]
                            output_base = f"{base1}_{suffix}"
                            while output_base in used_names:
                                suffix += 1
                                output_base = f"{base1}_{suffix}"
                            next_suffix[base1] = suffix + 1
                        fq1.name = output_base + suffix1
                        fq2.name = output_base + suffix2
                        used_names.add(output_base)
                    _write_record(out1, fq1, outfmt)
                    _write_record(out2, fq2, outfmt)
    logger.info(f"Processed {pair_count} read pairs in {time.perf_counter() - start:.2f} seconds.")


@ignore
def run():
    configure_logging(rich=True, force=True)
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", nargs="?", help="single-end FASTA/FASTQ input; stdin when piped.")
    parser.add_argument("-o", "--output", default="-", help="single-end output.")
    parser.add_argument("-1", "--read1", nargs="+", help="ordered R1 input files for paired-end mode.")
    parser.add_argument("-2", "--read2", nargs="+", help="ordered R2 input files for paired-end mode.")
    parser.add_argument("--output1", help="R1 output for paired-end mode.")
    parser.add_argument("--output2", help="R2 output for paired-end mode.")
    parser.add_argument("--outfmt", choices=("fastq", "fasta"), default="fastq",
                        help="output format.")
    parser.add_argument("--mode", choices=("index", "preserve"),
                        default="index", help="rename mode.")
    parser.add_argument("--prefix", default=None, help="name prefix used by index mode. "
                        "Defaults to the original FASTQ record name (mate suffixes /1 /2 are "
                        "stripped in paired-end mode).")
    parser.add_argument("--index-base", type=int, choices=(10, 16, 36), default=10,
                        help="numeric base used for compact index names.")
    parser.add_argument("--gzip-level", dest="compresslevel", type=int, choices=range(10), default=4,
                        help="gzip compression level for FASTQ output.")
    args = parser.parse_args()

    paired = args.read1 is not None or args.read2 is not None
    try:
        if paired:
            if args.input is not None or args.read1 is None or args.read2 is None:
                parser.error("paired-end mode requires --read1 and --read2 without positional input")
            if args.output1 is None or args.output2 is None:
                parser.error("paired-end mode requires --output1 and --output2")
            fastx_rename_pair(args.read1, args.read2, args.output1, args.output2, args.outfmt,
                              args.mode, args.prefix, args.index_base, args.compresslevel)
        else:
            if args.input is None and not sys.stdin.isatty():
                args.input = "-"
            if args.input is None:
                parser.error("single-end input is required")
            if args.output1 is not None or args.output2 is not None:
                parser.error("--output1 and --output2 are only valid in paired-end mode")
            fastx_rename(args.input, args.output, args.outfmt, args.mode, args.prefix,
                         args.index_base, args.compresslevel)
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    run()
