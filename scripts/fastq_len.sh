#！/bin/sh
if [ $# != 1 ];then
    echo "Usage: fastq_len.sh input.fa"
    exit 1;
fi

zcat $1 | awk '{if(NR%4==2)print length}' | sort | uniq -c | awk 'BEGIN{OFS="\t";print "Length\tCount"}{print $2,$1}'
