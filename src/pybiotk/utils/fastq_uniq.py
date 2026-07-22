#!/usr/bin/env python3
"""Remove duplicate FASTQ records after length filtering.

The default ``digest`` key mode stores one 128-bit BLAKE2b digest per eligible
read or read pair. For paired-end data, both mates are encoded into one digest
with explicit field boundaries. ``exact`` mode stores the original values and
is available when absolutely collision-free comparison is required.

The collision estimates below use the birthday bound and assume one uniformly
distributed key per unique read pair. "Ordinary hash" means a 64-bit hash such
as the value previously stored with ``numpy.int64(hash(...))``. Python's hash
is also randomized between processes, whereas BLAKE2b digests are stable.

=================  ============  =========================  ========================
Unique read pairs  Key count     128-bit digest collision   64-bit hash collision
=================  ============  =========================  ========================
500,000,000        500,000,000   3.6734e-22                 6.7534e-03 (0.6753%)
1,000,000,000      1,000,000,000 1.4694e-21                 2.6741e-02 (2.6741%)
=================  ============  =========================  ========================

Approximate memory usage for exact paired-end 2x150-nt sequence keys versus
128-bit digest keys on 64-bit CPython is shown below. The estimate assumes two
191-byte ASCII string objects plus a 64-byte tuple for ``exact``, a 49-byte
``bytes`` object for ``digest``, a 32-byte ``numpy.int64`` object for the
64-bit hash, and about 42 bytes of set-table overhead per entry. Values use
decimal GB and exclude allocator fragmentation and temporary memory during set
resizing.

================  ===================  =================  ===================
Key mode          Bytes/read pair      500M read pairs    1B read pairs
================  ===================  =================  ===================
exact             about 488            about 244 GB       about 488 GB
digest            about 91             about 45.5 GB      about 91 GB
64-bit hash       about 74             about 37 GB        about 74 GB
saving            about 397 (81.4%)    about 198.5 GB     about 397 GB
================  ===================  =================  ===================
"""
import argparse
import hashlib
import sys
import time
from typing import Hashable, Literal, Optional, Sequence, Tuple, Union

from pybiotk.io import FastqFile, FastqPair, OpenFqGzip
from pybiotk.utils import configure_logging, get_logger


UniqBy = Literal["seq", "id", "name"]
KeyMode = Literal["digest", "exact"]


logger = get_logger(__name__)


def _record_parts(fq, by: UniqBy) -> Tuple[Optional[str], ...]:
    if by == "seq":
        return (fq.sequence,)
    if by == "id":
        return (fq.name,)
    return fq.name, fq.comment


def _make_key(parts: Tuple[Optional[str], ...], key_mode: KeyMode) -> Hashable:
    if key_mode == "exact":
        return parts[0] if len(parts) == 1 else parts
    digest = hashlib.blake2b(digest_size=16)
    for part in parts:
        if part is None:
            digest.update(b"\0")
            continue
        value = part.encode()
        digest.update(b"\1")
        digest.update(len(value).to_bytes(8, "big"))
        digest.update(value)
    return digest.digest()


def _validate_options(input_files: Sequence[str], output_files: Sequence[str], min_len: int,
                      by: str, key_mode: str, compresslevel: int) -> None:
    if len(input_files) not in {1, 2}:
        raise ValueError("exactly one or two input FASTQ files are required")
    if len(output_files) != len(input_files):
        raise ValueError("the number of output files must match the number of input files")
    if min_len < 0:
        raise ValueError("min_len must be greater than or equal to zero")
    if by not in {"seq", "id", "name"}:
        raise ValueError(f"unsupported uniqueness key: {by}")
    if key_mode not in {"digest", "exact"}:
        raise ValueError(f"unsupported key mode: {key_mode}")
    if not 0 <= compresslevel <= 9:
        raise ValueError("compresslevel must be between 0 and 9")


def _log_summary(label: str, input_reads: int, too_short_reads: int, duplicate_reads: int,
                 output_reads: int, min_len: int) -> None:
    eligible_reads = input_reads - too_short_reads
    duplication_rate = duplicate_reads * 100 / eligible_reads if eligible_reads else 0.0
    logger.info(
        f"result summary:\n{label} too short (<{min_len}nt): {too_short_reads}\n"
        f"input {label}: {input_reads}\neligible {label}: {eligible_reads}\n"
        f"duplicate {label}: {duplicate_reads}\noutput {label}: {output_reads}\n"
        f"duplication rate: {duplication_rate:.2f}%"
    )


def single_end(input_fq: str, output: str, min_len: int = 15, by: UniqBy = "seq",
               key_mode: KeyMode = "digest", compresslevel: int = 4) -> None:
    logger.info(f"Single end mode, by {by}, key mode {key_mode} ...")
    input_reads = 0
    too_short_reads = 0
    duplicate_reads = 0
    output_reads = 0
    unique = set()
    with FastqFile(input_fq) as fqi, OpenFqGzip(output, compresslevel=compresslevel) as fqo:
        for fq in fqi:
            input_reads += 1
            if len(fq.sequence) < min_len:
                too_short_reads += 1
                continue
            key = _make_key(_record_parts(fq, by), key_mode)
            if key in unique:
                duplicate_reads += 1
                continue
            unique.add(key)
            fqo.write_fastx_record(fq)
            output_reads += 1
    _log_summary("reads", input_reads, too_short_reads, duplicate_reads, output_reads, min_len)


def pair_end(input_r1: str, input_r2: str, output_r1: str, output_r2: str,
             min_len: int = 15, by: UniqBy = "seq", key_mode: KeyMode = "digest",
             compresslevel: int = 4) -> None:
    logger.info(f"Pair end mode, by {by}, key mode {key_mode} ...")
    input_reads = 0
    too_short_reads = 0
    duplicate_reads = 0
    output_reads = 0
    unique = set()

    with FastqPair(input_r1, input_r2) as pairfq, \
            OpenFqGzip(output_r1, compresslevel=compresslevel) as r1, \
            OpenFqGzip(output_r2, compresslevel=compresslevel) as r2:
        for fq1, fq2 in pairfq:
            input_reads += 1
            if min(len(fq1.sequence), len(fq2.sequence)) < min_len:
                too_short_reads += 1
                continue
            parts = _record_parts(fq1, by) + _record_parts(fq2, by)
            key = _make_key(parts, key_mode)
            if key in unique:
                duplicate_reads += 1
                continue
            unique.add(key)
            r1.write_fastx_record(fq1)
            r2.write_fastx_record(fq2)
            output_reads += 1
    _log_summary("read pairs", input_reads, too_short_reads, duplicate_reads, output_reads, min_len)


def main(input_files: Union[Sequence[str], str] = "-", output_files: Union[Sequence[str], str] = "-",
         min_len: int = 15, by: UniqBy = "seq", key_mode: KeyMode = "digest",
         compresslevel: int = 4) -> None:
    if isinstance(input_files, str):
        input_files = [input_files]
    if isinstance(output_files, str):
        output_files = [output_files]

    _validate_options(input_files, output_files, min_len, by, key_mode, compresslevel)
    input_str = " ".join(input_files)
    output_str = " ".join(output_files)
    if input_str == "-":
        input_str = "stdin"
    logger.info(f"Processing {input_str} ...")

    start = time.perf_counter()
    if len(input_files) == 1:
        single_end(input_files[0], output_files[0], min_len, by, key_mode, compresslevel)
    else:
        pair_end(*input_files[:2], *output_files[:2], min_len=min_len, by=by,
                 key_mode=key_mode, compresslevel=compresslevel)
    end = time.perf_counter()
    if output_str == "-":
        logger.info(f"task finished in {end-start:.2f}s")
    else:
        logger.info(f"task finished in {end-start:.2f}s, output saved in {output_str}.")


def run():
    configure_logging(rich=True, force=True)
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(dest="input", type=str, nargs="*", default=(None if sys.stdin.isatty() else "-"), help="input fastq files. [stdin]")
    parser.add_argument("-o", "--output", dest="output", type=str, nargs="+", default="-", help="output fastq files. [stdout]")
    parser.add_argument("-m", "--min-len", dest="min_len", type=int, default=15, help="min length.")
    parser.add_argument("--key-mode", choices=("digest", "exact"), default="digest",
                        help="store 128-bit digests or exact values as uniqueness keys.")
    parser.add_argument("--gzip-level", dest="compresslevel", type=int, choices=range(10), default=4,
                        help="gzip compression level for output files.")
    by_group = parser.add_mutually_exclusive_group()
    by_group.add_argument("-i", "--by-id", dest="id", action="store_true", help="by id instead of seq.")
    by_group.add_argument("-n", "--by-name", dest="name", action="store_true", help="by full name instead of just id.")
    args = parser.parse_args()
    if not args.input:
        parser.print_help()
        sys.exit(1)
    by = "seq"
    if args.id:
        by = "id"
    if args.name:
        by = "name"
    try:
        main(args.input, args.output, args.min_len, by, args.key_mode, args.compresslevel)
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    run()
