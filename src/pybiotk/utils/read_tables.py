#!/usr/bin/env python3
"""
A simple tool for joining and filtering tables
"""
import argparse
import csv
import os
import re
import sys

import pandas as pd
from pybiotk.utils import configure_logging
from pybiotk.utils import ignore, read_table, write_table


def _load_names(namefile):
    if namefile == "-":
        namefile = sys.stdin
    if isinstance(namefile, (str, os.PathLike)):
        with open(namefile) as stream:
            return {name for line in stream for name in line.split()}
    return {name for line in namefile for name in line.split()}


def _matching_mask(series, names, delimiter=None, pattern=None):
    if not names:
        return pd.Series(False, index=series.index, dtype=bool)
    if pattern is not None:
        return series.str.contains(pattern, na=False)
    if delimiter is not None:
        return (
            series.fillna("")
            .str.split(delimiter, regex=False)
            .apply(lambda fields: not names.isdisjoint(fields))
        )
    return series.isin(names)


def main(table_list, outfile, namefile=None, noheader=False, column=0, delimiter=None, exclude=False, contains=False):
    names = _load_names(namefile) if namefile is not None else None
    pattern = None
    if contains and names:
        pattern = re.compile("|".join(re.escape(name) for name in sorted(names)))

    df_list = []
    for table in table_list:
        header = None if noheader else 0
        df = read_table(table, header=header, dtype=str, comment="#")
        if names is not None:
            if not -df.shape[1] <= column < df.shape[1]:
                raise ValueError(
                    f"column index {column} is out of range for table {table!r} "
                    f"with {df.shape[1]} columns"
                )
            mask = _matching_mask(df.iloc[:, column], names, delimiter, pattern)
            df = df.loc[~mask if exclude else mask]
        df_list.append(df)
    if not df_list:
        raise ValueError("table_list must contain at least one input table")
    out_df = df_list[0] if len(df_list) == 1 else pd.concat(df_list, ignore_index=True)
    if os.path.splitext(outfile)[1] == ".xlsx":
        write_table(out_df, outfile, header=not noheader)
    else:
        write_table(out_df, outfile, header=not noheader, quoting=csv.QUOTE_NONE)


@ignore
def run():
    configure_logging(rich=True, force=True)
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str, nargs="+",
                        help="input tables.")
    parser.add_argument('-o', dest='output', type=str,
                        default="-", help="output file name [stdout]")
    parser.add_argument('-n', dest="namefile", type=str, default=(None if sys.stdin.isatty() else "-"),
                        help="whose name is listed in FILE|stdin")
    parser.add_argument('-H', "--noheader", dest="noheader", action="store_true", help="if noheader")
    parser.add_argument('-c', dest="column", type=int, default=0, help="name column")
    parser.add_argument('-e', dest="exclude", action="store_true", help="name is not listed in FILE|stdin")
    parser.add_argument('-d', dest="delimiter", type=str, default=None, help="name is separated by delimiter in input tables")
    parser.add_argument('--contains', dest="contains", action="store_true", help="contains one of a substrings listed in FILE|stdin.")
    args = parser.parse_args()
    if args.namefile == "-":
        args.namefile = sys.stdin
    main(args.input, args.output, args.namefile, args.noheader, args.column, args.delimiter, args.exclude, args.contains)


if __name__ == "__main__":
    run()
