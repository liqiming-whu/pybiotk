#！/bin/sh
if [ $# != 1 ];then
    echo "Usage: fasta_u2t.sh input.fa"
    exit 1;
fi

awk '{if(/^>/)print $1; else {gsub("U", "T");print $0}}' $1