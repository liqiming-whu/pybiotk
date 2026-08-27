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

The current release is **1.3.6**. See [CHANGELOG.rst](CHANGELOG.rst) for the
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

`fastx_rename` uses short decimal names such as `read_1` by default and
preserves input order. Base 16 and base 36 indexes are also available through
`--index-base`. Use `--use-original-name` to use each original FASTQ record
name as the prefix instead of `read`.

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

To retain each original record name as the index prefix:

```bash
fastx_rename \
  --read1 sub1_R1.fq.gz sub2_R1.fq.gz \
  --read2 sub1_R2.fq.gz sub2_R2.fq.gz \
  --output1 R1.fq.gz --output2 R2.fq.gz \
  --mode index --use-original-name
```

Use `--mode preserve` to retain original names and suffix only duplicate
names. The `index` mode uses names such as `read_1`, without adding `|`.

## Genomic annotation with `pyanno`

`pyanno` annotates BAM alignments or BED intervals against a GTF gene model.
The input type is inferred from the `.bam` or `.bed*` filename suffix. Output
is a tab-separated table containing the input coordinates, annotation class,
gene coordinates, gene name, gene or transcript ID, and gene type.

Annotate a BAM at transcript level:

```bash
pyanno -i aligned.bam -g genes.gtf -o aligned.annotation.tsv
```

Annotate paired-end fragments from a name-sorted BAM with the common
first-strand RNA-seq rule:

```bash
pyanno -i aligned.name_sorted.bam -g genes.gtf \
  -o fragments.annotation.tsv \
  --level gene --pair --ordered_by_name --strand \
  --rule '1+-,1-+,2++,2--'
```

Annotate stranded BED intervals and customize the regulatory regions:

```bash
pyanno -i regions.bed6 -g genes.gtf -o regions.annotation.tsv \
  --strand --tss_region -3000 0 --downstream 3000 \
  --tss --tes --start_condon --stop_condon
```

Important options:

- `--level transcript|gene` selects the annotation level.
- `--strand` requires compatible feature and read/interval strands.
- `--rule` describes how mapped BAM reads encode the originating RNA strand.
- `--pair` annotates paired-end fragments rather than individual alignments.
- `--ordered_by_name` enables streaming pair iteration for a name-sorted BAM.
- `--tss_region START END` and `--downstream LENGTH` control regulatory
  annotation windows.

## GTF conversion with `gtf2bed`

`gtf2bed` converts a GTF file into BED12 transcripts, BED6 features, introns,
or gene/transcript information tables. It can also filter records by gene
type, transcript type, IDs, or names.

Create transcript BED12 using transcript IDs as names:

```bash
gtf2bed genes.gtf --outfmt bed12 --name transcript_id \
  -o transcripts.bed12
```

Create a gene-level BED6 file using gene names:

```bash
gtf2bed genes.gtf --outfmt bed6 --feature gene --name gene_name \
  -o genes.bed6
```

Extract introns from protein-coding transcripts:

```bash
gtf2bed genes.gtf --outfmt intron --name transcript_id \
  --transcript_types protein_coding -o introns.bed6
```

The supported output formats are `bed12`, `bed6`, `intron`, `gene_info`, and
`trans_info`. GTF input and converted output can also be streamed:

```bash
zcat genes.gtf.gz | gtf2bed --outfmt bed12 > transcripts.bed12
```

## BAM to FASTA/FASTQ with `bam2fastx`

`bam2fastx` detects single-end or paired-end BAM input and reconstructs FASTA
or compressed FASTQ records:

```bash
bam2fastx aligned.bam -o recovered --outfmt fastq
```

For paired-end input, it writes paired and unpaired outputs separately, such
as `recovered.R1.fq.gz`, `recovered.R2.fq.gz`, and the corresponding
`unpaired` files. Use name-sorted input for streaming paired reads:

```bash
bam2fastx aligned.name_sorted.bam -o recovered --bamtype PE \
  --outfmt fastq --ordered_by_name
```

## Sequence extraction with `genomefetcher`

Extract one stranded genomic interval (quote the location to protect shell
parentheses):

```bash
genomefetcher -f genome.fa -l 'chr1:100000-101000(+)' -o region.fa
```

Extract spliced transcript exons from a GTF:

```bash
genomefetcher -f genome.fa -g genes.gtf --regions exons \
  --transcript_ids ENST00000000000 -o transcript.fa
```

Other region modes include `all`, `introns`, `first_exon`, `last_exon`,
`first_intron`, `last_intron`, `5utr`, `3utr`, and `cds`. Add `--separate` to
write multi-block regions separately, or `--no-sequence` to output selected
coordinates without sequences.

## RNA-seq library inference with `infer_experiment`

Estimate whether a mapped RNA-seq library is single-end or paired-end and
which strand rule best explains its alignments:

```bash
infer_experiment aligned.bam -g genes.gtf --mapq 30
```

The reported strand rules can be passed directly to tools such as `pyanno`.
`infer_experiment` can also filter alignments by a selected rule:

```bash
infer_experiment aligned.bam -g genes.gtf \
  --filter '1+-,1-+,2++,2--' --outbam stranded.bam
```

## Additional useful CLI tools

Convert one or more BED files into combined bedGraph coverage:

```bash
bed2bedgraph input.bed > coverage.bedgraph
bed2bedgraph lane1.bed.gz lane2.bed.gz --header > coverage.bedgraph
```

The built-in command accepts unsorted BED input and uses an endpoint sweep
rather than expanding intervals base by base. For very large BED files that
are grouped in chromosome order, a compiled bedtools alternative is:

```bash
bedtools genomecov -i input.bed -bg -g chrom.sizes
```

Generate `chrom.sizes` from the same reference FASTA used for alignment. The
script reuses an existing `.fai` index or creates one with `samtools faidx`:

```bash
scripts/get_chrom_length.sh reference.fa chrom.sizes
```

Count mapped reads, or paired-end fragments, by reference:

```bash
reference_count aligned.bam -o reference_counts.tsv
reference_count aligned.name_sorted.bam -o fragment_counts.tsv \
  --pair --ordered_by_name
```

Calculate mapped fragment lengths and save a length-frequency table:

```bash
rna_fragment_size aligned.name_sorted.bam \
  --ordered_by_name -s fragment_sizes.tsv > fragments.tsv
```

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

List every command installed with pybiotk and its short description:

```bash
pybiotk
```

The list is read from the installed package metadata, so it automatically
includes newly registered tools. Installed console commands include:

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

## Building a source distribution

Build the release source archive in an isolated environment:

```bash
pyproject-build --sdist
```

Release versions are derived from Git tags through `setuptools_scm`.
