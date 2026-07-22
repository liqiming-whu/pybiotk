#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
    echo "Usage: merge_subseq.sh <threads> <input.fa> <outdir>" >&2
    exit 1
fi

threads=$1
fasta=$2
outdir=$3

if [[ ! $threads =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: threads must be a positive integer" >&2
    exit 1
fi
if [[ ! -f $fasta ]]; then
    echo "Error: input FASTA not found: $fasta" >&2
    exit 1
fi

mkdir -p "$outdir/bowtie_index"

bowtie-build --threads "$threads" "$fasta" "$outdir/bowtie_index/ref"
bowtie -p "$threads" -x "$outdir/bowtie_index/ref" "$fasta" -f -a -v 0 --norc -S \
    2> "$outdir/mapped.log" |
    samtools view - |
    awk '$1 != $3 { print $1 "\t" $3 }' > "$outdir/overlap.reads.txt"

subseq_analysis -f "$fasta" -o "$outdir" -r "$outdir/overlap.reads.txt" \
    > "$outdir/collapse_seq.log"
