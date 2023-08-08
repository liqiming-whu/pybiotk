#！/bin/sh
if [ $# != 1 ];then
    echo "Usage: fasta_len.sh input.fa"
    exit 1;
fi

awk '{if(!/>/)print length}' $1 | sort | uniq -c | awk 'BEGIN{OFS="\t";print "Length\tCount"}{print $2,$1}'
