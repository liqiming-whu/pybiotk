# pybiotk

`pybiotk` is a Python toolkit and command-line collection for common
bioinformatics data-processing tasks, including FASTA/FASTQ, BAM, BED, GTF,
BigWig, annotation, and interval operations.

## Installation

Install the latest release from PyPI:

```bash
pip install pybiotk
```

Install from Gitee:

```bash
git clone https://gitee.com/liqiming_whu/pybiotk.git
cd pybiotk
pip install .
```

Install from GitHub:

```bash
git clone https://github.com/liqiming-whu/pybiotk.git
cd pybiotk
pip install .
```

The current release is **1.3.5**. See [CHANGELOG.rst](CHANGELOG.rst) for the
release notes.

## FASTQ deduplication

`fastq_uniq` removes duplicate reads after length filtering. The default
`digest` mode stores stable 128-bit BLAKE2b digests to reduce memory usage.
Use `--key-mode exact` when exact, collision-free key comparison is required.

Single-end:

```bash
fastq_uniq reads.fq.gz -o reads.unique.fq.gz \
  --key-mode digest --gzip-level 4
```

Paired-end:

```bash
fastq_uniq R1.fq.gz R2.fq.gz \
  -o R1.unique.fq.gz R2.unique.fq.gz \
  --key-mode digest --gzip-level 4
```

Uniqueness can be evaluated by sequence (default), read ID (`--by-id`), or
full read name (`--by-name`). Paired inputs must have the same number of
records and matching mate names.

## FASTA/FASTQ renaming

`fastx_rename` uses compact base-36 indexes by default and preserves input
order. Base 10 and hexadecimal indexes are also available through
`--index-base`.

Single-end or independently streamed mates:

```bash
zcat sub1_R1.fq.gz sub2_R1.fq.gz | \
  fastx_rename -o R1.fq.gz --mode index --gzip-level 4

zcat sub1_R2.fq.gz sub2_R2.fq.gz | \
  fastx_rename -o R2.fq.gz --mode index --gzip-level 4
```

For guaranteed matching mate identifiers, use paired-end mode. Input files
are processed in the order supplied, with one continuous index across lanes:

```bash
fastx_rename \
  --read1 sub1_R1.fq.gz sub2_R1.fq.gz \
  --read2 sub1_R2.fq.gz sub2_R2.fq.gz \
  --output1 R1.fq.gz --output2 R2.fq.gz \
  --mode index --index-base 36 --gzip-level 4
```

Use `--mode preserve` to retain original names and suffix only duplicate
names. The `index` mode uses names such as `read_1`, without adding `|`.

## Logging API

Library modules can use the centralized standard logging configuration:

```python
from pybiotk.utils import get_logger

logger = get_logger(__name__)
logger.info("Processing started")
```

Command-line applications can enable Rich-formatted stderr output:

```python
from pybiotk.utils import configure_logging

configure_logging(rich=True, force=True)
```

The legacy import `from pybiotk.utils import logging` remains supported.

## Commands

Installed console commands include:

- Format conversion: `gtf2bed`, `bed2bedgraph`, `fq2fasta`, `fa2fastq`,
  `bam2fastx`, `bampe_order_by_name`
- FASTA/FASTQ processing: `fastq_uniq`, `fastx_rename`, `fastq_join`,
  `fasta_filter`, `reverse_fastx`, `seq_random`
- BAM and annotation: `bam_random`, `pyanno`, `infer_experiment`,
  `rna_fragment_size`, `reference_count`
- Table and genomic utilities: `read_tables`, `merge_row`, `rmats_filter`,
  `count_normalize`, `genomefetcher`, `bigwigfetcher`, `summary_log`
- Transcript and plotting utilities: `merge_transcript`, `metaplot`

Run any command with `--help` for its complete options, for example:

```bash
fastq_uniq --help
fastx_rename --help
pyanno --help
```

## Building from source

Build an sdist and platform wheel in an isolated environment:

```bash
pyproject-build
```

Release versions are derived from Git tags through `setuptools_scm`.
