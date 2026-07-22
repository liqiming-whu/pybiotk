#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: fasta_len.sh <input.fa[.gz]|->" >&2
    exit 1
fi

input=$1

read_fasta() {
    if [[ $input == - ]]; then
        cat
    elif [[ $input == *.gz ]]; then
        gzip -cd -- "$input"
    else
        cat -- "$input"
    fi
}

read_fasta |
    awk '
        /^>/ {
            if (seen) print sequence_length
            sequence_length = 0
            seen = 1
            next
        }
        seen {
            gsub(/[[:space:]]/, "")
            sequence_length += length($0)
        }
        END {
            if (seen) print sequence_length
        }
    ' |
    sort -n |
    uniq -c |
    awk 'BEGIN { OFS="\t"; print "Length", "Count" } { print $2, $1 }'
