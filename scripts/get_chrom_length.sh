#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: get_chrom_length.sh <reference.fa> <output.chrom.sizes>" >&2
    exit 1
fi

reference=$1
output=$2

if [[ ! -f $reference ]]; then
    echo "Error: reference FASTA not found: $reference" >&2
    exit 1
fi

if [[ ! -f ${reference}.fai ]]; then
    samtools faidx "$reference"
fi

cut -f 1,2 "${reference}.fai" > "$output"
