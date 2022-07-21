#!/usr/bin/env python3
import argparse
import os
import sys
from collections import defaultdict
from itertools import groupby

import pandas as pd

from pybiotk.utils import split_discontinuous, ignore


def main(bed_list, header):
    if header:
        name = os.path.splitext(os.path.basename(bed_list[0]))[0]
        print(f'track type=bedGraph name="{name}" visibility=full color=0,51,51 AutoScale=on alwaysZero=on maxHeightPixels=50:50:50')
    df_list = []
    for bedfile in bed_list:
        df_list.append(pd.read_table(bedfile, header=None))
    df = pd.concat(df_list)
    bed = df.sort_values(by=0).reset_index(drop=True)
    for chrom, group in groupby(bed.iterrows(), lambda x: x[1][0]):
        values = defaultdict(int)
        for _, row in group:
            start = row[1]
            end = row[2]
            for pos in range(start, end):
                values[int(pos)] += 1
        values_items = sorted(values.items(), key=lambda x: x[0])
        for value, g in groupby(values_items, lambda x: x[1]):
            pos = [i[0] for i in g]
            for pos_list in split_discontinuous(pos):
                sys.stdout.write(f"{chrom}\t{pos_list[0]}\t{pos_list[-1]}\t{value}\n")


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str, nargs="+",
                        help="Input bedfiles")
    parser.add_argument("--header", dest="header", action="store_true",
                        help="Add a header to the output file")
    args = parser.parse_args()
    main(args.input, args.header)


if __name__ == "__main__":
    run()
