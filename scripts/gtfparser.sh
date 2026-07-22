#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 2 ]]; then
    echo "Usage: gtfparser.sh <annotation.gtf> <output-prefix>" >&2
    exit 1
fi

gtf=$1
prefix=$2

if [[ ! -f $gtf ]]; then
    echo "Error: GTF file not found: $gtf" >&2
    exit 1
fi

tempdir=$(mktemp -d "${TMPDIR:-/tmp}/gtfparser.XXXXXXXX")
trap 'rm -rf "$tempdir"' EXIT

gtfToGenePred "$gtf" "$tempdir/annotation.genePred"
genePredToBed "$tempdir/annotation.genePred" "${prefix}.bed12"

awk -F "\t" -v OFS="\t" \
    -v gene_output="${prefix}.gene.bed" \
    -v transcript_output="${prefix}.transcript.bed" \
    -v exon_output="${prefix}.exon.bed" '
    function attribute(name, fields, count, i, value) {
        count = split($9, fields, ";")
        for (i = 1; i <= count; i++) {
            value = fields[i]
            sub(/^[[:space:]]+/, "", value)
            if (value ~ ("^" name "[[:space:]]+")) {
                sub("^" name "[[:space:]]+", "", value)
                gsub(/^"|"$/, "", value)
                return value
            }
        }
        return "."
    }
    BEGIN {
        printf "%s", "" > gene_output
        close(gene_output)
        printf "%s", "" > transcript_output
        close(transcript_output)
        printf "%s", "" > exon_output
        close(exon_output)
    }
    $3 == "gene" {
        print $1, $4-1, $5, attribute("gene_id"), $6, $7 > gene_output
    }
    $3 == "transcript" {
        print $1, $4-1, $5, attribute("transcript_id"), $6, $7 > transcript_output
    }
    $3 == "exon" {
        print $1, $4-1, $5, attribute("transcript_id"), $6, $7 > exon_output
    }
' "$gtf"

awk -v OFS="\t" '{ print $4, $2, $3, $1, $5, $6 }' "${prefix}.transcript.bed" > "$tempdir/transcript.reverse.bed"
awk -v OFS="\t" '{ print $4, $2, $3, $1, $5, $6 }' "${prefix}.exon.bed" > "$tempdir/exon.reverse.bed"
bedtools subtract -a "$tempdir/transcript.reverse.bed" -b "$tempdir/exon.reverse.bed" -s |
    awk -v OFS="\t" '{ print $4, $2, $3, $1, $5, $6 }' > "${prefix}.intron.bed"
