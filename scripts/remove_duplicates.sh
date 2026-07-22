#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 5 ]]; then
    echo "Usage: remove_duplicates.sh <threads> <read1.fq.gz> <read2.fq.gz> <out1.fq.gz> <out2.fq.gz>" >&2
    exit 1
fi

threads=$1
read1=$2
read2=$3
out1=$4
out2=$5

if [[ ! $threads =~ ^[1-9][0-9]*$ ]]; then
    echo "Error: threads must be a positive integer" >&2
    exit 1
fi

tempdir=$(mktemp -d "${out1}.temp.XXXXXXXX")
trap 'rm -rf "$tempdir"' EXIT

decompressed1="$tempdir/R1.decompressed.fastq"
decompressed2="$tempdir/R2.decompressed.fastq"
uncompressed1="$tempdir/R1.deduplicated.fastq"
uncompressed2="$tempdir/R2.deduplicated.fastq"

echo "Decompressing $read1 and $read2 ..." >&2
pigz -p "$threads" -d -c -- "$read1" > "$decompressed1"
pigz -p "$threads" -d -c -- "$read2" > "$decompressed2"

lines1=$(wc -l < "$decompressed1")
lines2=$(wc -l < "$decompressed2")
if (( lines1 % 4 != 0 || lines2 % 4 != 0 )); then
    echo "Error: input files must contain complete four-line FASTQ records" >&2
    exit 1
fi
if (( lines1 != lines2 )); then
    echo "Error: read1 and read2 contain different numbers of records" >&2
    exit 1
fi

input_num=$((lines1 / 4))
echo "Input read pairs: $input_num" >&2

paste "$decompressed1" "$decompressed2" |
    awk 'BEGIN { OFS="\t" } { printf "%s%s", $0, (NR % 4 == 0 ? ORS : OFS) }' |
    LC_ALL=C sort -T "$tempdir" -t $'\t' -k3,4 -u |
    awk -F '\t' -v out1="$uncompressed1" -v out2="$uncompressed2" '
        {
            print $1 "\n" $3 "\n" $5 "\n" $7 > out1
            print $2 "\n" $4 "\n" $6 "\n" $8 > out2
        }
    '

if [[ -f $uncompressed1 ]]; then
    output_lines=$(wc -l < "$uncompressed1")
else
    : > "$uncompressed1"
    : > "$uncompressed2"
    output_lines=0
fi
output_num=$((output_lines / 4))

echo "Output read pairs: $output_num" >&2
awk -v input="$input_num" -v output="$output_num" '
    BEGIN {
        rate = input == 0 ? 0 : (input - output) * 100 / input
        printf "Duplication rate: %.2f%%\n", rate > "/dev/stderr"
    }
'

pigz -4 -p "$threads" -c -- "$uncompressed1" > "$out1"
pigz -4 -p "$threads" -c -- "$uncompressed2" > "$out2"
