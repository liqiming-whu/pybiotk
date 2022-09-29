#!/usr/bin/env python3
"""
A simple tool for joining and filtering tables
"""
import argparse
import csv
import os
import sys

import pandas as pd
from pybiotk.utils import ignore, read_table, write_table


def main(table_list, outfile, namefile, noheader, column, delimiter=None, exclude=False, contains=False):
    df_list = []
    for table in table_list:
        header = None if noheader else 0
        df = read_table(table, header=header, dtype=str, comment="#")
        if namefile is not None:
            names = set(j for i in namefile for j in i.split())
            if names:
                if exclude:
                    if contains:
                        df = df.loc[~df.iloc[:, column].str.contains("|".join(names))]
                    else:
                        if delimiter is not None:
                            df = df.loc[~df.iloc[:, column].str.split(delimiter).apply(lambda x: any(set(x) & names))]
                        else:
                            df = df.loc[~df.iloc[:, column].isin(names)]
                else:
                    if contains:
                        df = df.loc[df.iloc[:, column].str.contains("|".join(names))]
                    else:
                        if delimiter is not None:
                            df = df.loc[df.iloc[:, column].str.split(delimiter).apply(lambda x: any(set(x) & names))]
                        else:
                            df = df.loc[df.iloc[:, column].isin(names)]
        df_list.append(df)
    if len(df_list) == 1:
        out_df = df_list[0]
    else:
        out_df = pd.concat(df_list)
    if os.path.splitext(outfile)[1] == ".xlsx":
        write_table(out_df, outfile, header=not noheader)
    else:
        write_table(out_df, outfile, header=not noheader, quoting=csv.QUOTE_NONE)


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str, nargs="+",
                        help="input tables.")
    parser.add_argument('-o', dest='output', type=str,
                        default="-", help="output file name [stdout]")
    parser.add_argument('-n', dest="namefile", type=argparse.FileType('r'), default=(None if sys.stdin.isatty() else sys.stdin),
                        help="whose name is listed in FILE|stdin")
    parser.add_argument('-H', "--noheader", dest="noheader", action="store_true", help="if noheader")
    parser.add_argument('-c', dest="column", type=int, default=0, help="name column")
    parser.add_argument('-e', dest="exclude", action="store_true", help="name is not listed in FILE|stdin")
    parser.add_argument('-d', dest="delimiter", type=str, default=None, help="name is separated by delimiter in input tables")
    parser.add_argument('--contains', dest="contains", action="store_true", help="contains one of a substrings listed in FILE|stdin.")
    args = parser.parse_args()
    main(args.input, args.output, args.namefile, args.noheader, args.column, args.delimiter, args.exclude, args.contains)


if __name__ == "__main__":
    run()
