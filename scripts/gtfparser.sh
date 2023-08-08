#！/bin/sh
if [ $# != 2 ];then
    echo "Usage: gtfparser.sh gencode.v32.annotation.gtf hg38"
    exit 1;
fi

GTF=$1
PREFIX=$2

## gtf2bed12
gtfToGenePred $GTF $PREFIX.genePhred && genePredToBed $PREFIX.genePhred $PREFIX.bed12 && rm $PREFIX.genePhred


## gtf2intron
awk -v OFS="\t" '{if($3=="gene"){match($0,/gene_id "(\S*)"/,a);print $1,$4-1,$5,a[1],$6,$7}}' $GTF > $PREFIX.gene.bed
awk -v OFS="\t" '{if($3=="transcript"){match($0,/transcript_id "(\S*)"/,a);print $1,$4-1,$5,a[1],$6,$7}}' $GTF > $PREFIX.transcript.bed
awk -v OFS="\t" '{if($3=="exon"){match($0,/transcript_id "(\S*)"/,a);print $1,$4-1,$5,a[1],$6,$7}}' $GTF > $PREFIX.exon.bed
awk -v OFS="\t" '{print $4,$2,$3,$1,$5,$6}' $PREFIX.transcript.bed > transcript.reverse.bed && \
awk -v OFS="\t" '{print $4,$2,$3,$1,$5,$6}' $PREFIX.exon.bed > exon.reverse.bed && \
bedtools subtract -a transcript.reverse.bed -b exon.reverse.bed -s | awk -v OFS="\t" '{print $4,$2,$3,$1,$5,$6}' > $PREFIX.intron.bed && \
rm transcript.reverse.bed exon.reverse.bed
