#！/bin/sh
if [ $# != 2 ];then
    echo "Usage: get_chrom_length.sh GRCh38.p13.genome.fa hg38"
    exit 1;
fi

samtools faidx $1
cut -f 1,2 $1.fai > $2.chrom.sizes
