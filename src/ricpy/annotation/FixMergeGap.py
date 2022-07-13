#!/usr/bin/env python3
import os
import sys
import time
import argparse
import pandas as pd
from pybiotk.utils import logging, blocks_len, write_table
from pybiotk.io import GenomeFile, Openbed
from pybiotk.intervals import GRangeTree, GRange, merge_intervals


def main(known_JS: str, annofile:str, outannofile: str, outrna:str, genome_fa: str):
    start = time.perf_counter()
    genome = GenomeFile(genome_fa)
    annodf = pd.read_table(annofile, header=0)
    outannolist = []
    rna = open(outrna, "w")
    logging.info("start to load knownJS tree...")
    js_tree = GRangeTree()
    with Openbed(known_JS) as bedfile:
        for bed in bedfile:
            rec = GRange(chrom=bed.chrom, start=bed.start, end=bed.end)
            js_tree.add_range(rec)
    logging.info("load completed.")

    for anno in annodf.iterrows():
        read_name = anno["seqname"]
        strand = anno["strand"]
        chrom = anno["chrom"]
        block = eval(anno["blocks"])
        block_len = len(block)
        chimeric_JS_list = None
        if block_len > 1:
            JS = []
            for indx in range(1, block_len):
                JS.append((block[indx-1][1], block[indx][0]))
            new_JS_li = []
            for js in JS:
                js_rec = GRange(chrom, js[0], js[1])
                overlap_js = js_tree.find(
                    js_rec.chrom, js_rec.start, js_rec.end)
                if not any([x == js_rec for x in overlap_js]):
                    new_JS_li.append((js_rec.start, js_rec.end))
            if new_JS_li:
                chimeric_JS_list = new_JS_li
        if chimeric_JS_list is None:
            outannolist.append(anno)
            seq = genome.fetch(chrom, block, strand)
        else:
            if blocks_len(chimeric_JS_list) > 150:
                sys.stderr.write(f"{read_name} merged reads gap length > 150.\n")
                continue
            new_blocks = merge_intervals(block + chimeric_JS_list)
            seq = genome.fetch(chrom, new_blocks, strand)
            new_anno = anno.copy()
            new_anno["blocks"] = str(new_blocks)
            outannolist.append(new_anno)
        rna.write(f">{read_name}\n{seq}\n")
    rna.close()
    outdf = pd.DataFrame(outannolist)
    write_table(outdf, outannofile)
    end = time.perf_counter()
    logging.info(f"task finished in {end-start:.2f}s.")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-g", dest="genome",
                        type=str, required=True, help="genome fasta")
    parser.add_argument("-a", dest="anno",
                        type=str, required=True, help="rna anno file")
    parser.add_argument("-j", dest="js",
                        type=str, required=True, help="js file")
    parser.add_argument("-f", dest="fasta", default=os.devnull,
                        type=str, help="mapped reads")
    parser.add_argument("-o", dest="rna_anno", default=os.devnull,
                        type=str, help="rna anno")
    args = parser.parse_args()

    main(args.js, args.anno, args.rna_anno, args.fasta, args.genome)


if __name__ == "__main__":
    run()
