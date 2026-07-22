#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: fasta_u2t.sh <input.fa[.gz]|->" >&2
    exit 1
fi

input=$1

if [[ $input == - ]]; then
    cat
elif [[ $input == *.gz ]]; then
    gzip -cd -- "$input"
else
    cat -- "$input"
fi | awk '/^>/ { print; next } { gsub(/U/, "T"); gsub(/u/, "t"); print }'
