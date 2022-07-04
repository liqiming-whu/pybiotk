#!/bin/sh
if [ $# != 3 ];then
    echo "Usage: merge_subseq.sh threads input.fa outdir"
    exit 1;
fi

threads=$1
fasta=$2
outdir=$3

mkdir -p $outdir/bowtie_index

bowtie-build --threads $threads $fasta $outdir/bowtie_index/ref

bowtie -p $threads -x $outdir/bowtie_index/ref $fasta -f -a -v 0 --norc -S $outdir/mapped.sam 2> $outdir/mapped.log

samtools view $outdir/mapped.sam | awk '$1 != $3' | cut -f 1,3 > $outdir/overlap.reads.txt
rm $outdir/mapped.sam

subseq_analysis -f $fasta -o $outdir -r $outdir/overlap.reads.txt > $outdir/collapse_seq.log
