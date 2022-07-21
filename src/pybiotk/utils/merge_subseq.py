#!/usr/bin/env python3
import argparse
import sys

import pysam


def merge_sub(seq1, seq2):
    if seq1.sequence.find(seq2.sequence) != -1:
        seq1.name = seq1.name + ";" + seq2.name
        return seq1
    elif seq2.sequence.find(seq1.sequence) != -1:
        seq2.name = seq1.name + ";" + seq2.name
        return seq2
    else:
        return None


def overlaps(a):
    i = 0
    while i < len(a):
        j = 0
        while j < len(a) and i < len(a):
            if i == j:
                j += 1
                continue
            if a[i].sequence.find(a[j].sequence) != -1:
                a[i].name = a[i].name + ";" + a[j].name
                a.pop(j)
                j -= 1
            elif a[j].sequence.find(a[i].sequence) != -1:
                a[j].name = a[j].name + ";" + a[j].name
                a[i] = a[j]
                a.pop(j)
                j -= 1
            j += 1
        i += 1
    return a


def run_method1(inputfa):
    merge_list = []
    with pysam.FastxFile(inputfa) as fh:
        i = 0
        for entry in fh:
            i += 1
            sys.stderr.write(f"\rprocess {i} reads")
            if not merge_list:
                merge_list.append(entry)
                continue
            collapse_seq = entry
            filter_seq_list = []
            for i, seq in enumerate(merge_list):
                merge_seq = merge_sub(seq, collapse_seq)
                if merge_seq:
                    collapse_seq = merge_seq
                else:
                    filter_seq_list.append(merge_list[i])

            merge_list = filter_seq_list
            merge_list.append(collapse_seq)
        sys.stderr.write("\n")
        for i in merge_list:
            print(f">{i.name}\n{i.sequence}")


def run_method2(inputfa):
    with pysam.FastxFile(inputfa) as fh:
        results = overlaps(list(fh))

    for i in results:
        print(f">{i.name}\n{i.sequence}")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(dest="fasta", type=str, help="fasta")
    parser.add_argument("-m", dest="method", type=int, default=1, help="method1 or method2")

    args = parser.parse_args()
    if not args.method == 1:
        run_method2(args.fasta)
    else:
        run_method1(args.fasta)


if __name__ == "__main__":
    run()
