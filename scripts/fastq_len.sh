#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
    echo "Usage: fastq_len.sh <input.fq[.gz]|->" >&2
    exit 1
fi

input=$1

if [[ $input == - ]]; then
    cat
elif [[ $input == *.gz ]]; then
    gzip -cd -- "$input"
else
    cat -- "$input"
fi |
    awk '
        NR % 4 == 2 { print length($0) }
        END {
            if (NR % 4 != 0) {
                print "Error: incomplete four-line FASTQ record" > "/dev/stderr"
                exit 1
            }
        }
    ' |
    sort -n |
    uniq -c |
    awk 'BEGIN { OFS="\t"; print "Length", "Count" } { print $2, $1 }'
