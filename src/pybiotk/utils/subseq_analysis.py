#!/usr/bin/env python3
import argparse
import os
from collections import defaultdict

import pysam


def main(originfa, path, outdir, more_x="C"):
    if not os.path.exists(path):
        os.makedirs(path)
    reads_dict = defaultdict(list)
    a = []
    b = []
    with open(path) as reads:
        for line in reads:
            fileds = line.split()
            reads_dict[fileds[1]].append(fileds[0])
            a.append(fileds[0])
            b.append(fileds[1])
    a = set(a)
    b = set(b)

    newset = b - a
    print(f"contained only by other sequences: {len(a - b)}")
    print(f"is contained by some sequences, and contains some sequences: {len(a & b)}")
    print(f"only contains other sequences: {len(newset)}")

    collaspe_dict = {}

    for read in reads_dict:
        if read not in newset:
            continue
        subs = []
        pirs = reads_dict[read].copy()
        for p in pirs:
            subs.append(p)
            if p in reads_dict:
                subs.extend(reads_dict[p])
        subs = set(subs)
        collaspe_dict[read] = list(subs)

    total = []
    with open(os.path.join(outdir, "collapse_info.tsv"), "w") as cof:
        for key, value in collaspe_dict.items():
            total.append(key)
            total.extend(value)
            cof.write(f"{key}\t{len(value)}\t{','.join(value)}\n")

    col_set = set(total)

    print(f"total overlaped reads: {len(col_set)}")

    fa = pysam.FastxFile(originfa)
    with open(os.path.join(outdir, "collapse.fa"), "w") as foa, open(os.path.join(outdir, "uniq.fa"), "w") as uoa:
        for entry in fa:
            if entry.name not in col_set:
                foa.write(f">{entry.name}\n{entry.sequence}\n")
                uoa.write(f">{entry.name}\n{entry.sequence}\n")
            elif entry.name in newset:
                foa.write(f">{entry.name}\n{entry.sequence}\n")

    rna_fetch = pysam.FastaFile(originfa)
    with open(path) as reads, open(os.path.join(outdir, "isoform_info.tsv"), "w") as isof:
        for line in reads:
            fileds = line.split()
            pir1 = fileds[0]
            pir2 = fileds[1]
            pir1s = rna_fetch.fetch(pir1)
            pir2s = rna_fetch.fetch(pir2)
            short = len(pir2s) - len(pir1s)
            pos = pir2s.find(pir1s)
            assert pos >= 0

            if pos == 0:
                info = "3isoform"
            else:
                if pos + len(pir1s) == len(pir2s):
                    info = "5isoform"
                else:
                    info = "53isoform"
            isof.write(f"{pir1}\t{pir2}\t{info}\t{short}\t{pir1s}\t{pir2s}\n")

    with open(path) as reads, open(os.path.join(outdir, f"more_{more_x}_info.tsv"), "w") as c:
        for line in reads:
            info = ""
            fileds = line.split()
            pir1 = fileds[0]
            pir2 = fileds[1]
            pir1s = rna_fetch.fetch(pir1)
            pir2s = rna_fetch.fetch(pir2)
            if len(pir2s) - len(pir1s) > 2:
                continue
            if more_x + pir1s == pir2s:
                info = "L" + more_x
            if pir1s + more_x == pir2s:
                info = "R" + more_x
            if more_x + pir1s + more_x == pir2s:
                info = "LR" + more_x
            if info:
                c.write(f"{pir1}\t{pir2}\t{info}\t{pir1s}\t{pir2s}\n")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", dest="originfa", type=str, required=True, help="origin fasta.")
    parser.add_argument("-o", dest="outdir", type=str, required=True, help="output dir.")
    parser.add_argument("-r", dest="reads_overlap", type=str, required=True, help="overlaped reads file, create by bowtie and shell scripts.")
    parser.add_argument("--more_x", dest="more_x", type=str, default="C", choices=("A", "T", "G", "C"), help="more_x analysis.")
    args = parser.parse_args()
    main(args.originfa, args.reads_overlap, args.outdir, args.reads_overlap)


if __name__ == "__main__":
    run()
