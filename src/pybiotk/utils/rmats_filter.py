#!/usr/bin/env python3
import argparse
import sys

import pandas as pd


def main(inf, outf, up_out, down_out, topnumber, p_value=None, readnumber=None, fdr_value=None, inclevel=None):
    df = pd.read_table(inf)
    if p_value:
        df = df[df['PValue'] < p_value]
    if fdr_value:
        df = df[df['FDR'] < fdr_value]
    if inclevel:
        df = df[abs(df['IncLevelDifference']) > inclevel]
    if readnumber:
        df['AvgSAMPLE1'] = (df['IJC_SAMPLE_1'] + df['SJC_SAMPLE_1']) / 2
        df['AvgSAMPLE2'] = (df['IJC_SAMPLE_2'] + df['SJC_SAMPLE_2']) / 2
        df = df[(df['AvgSAMPLE1'] > readnumber) & (df['AvgSAMPLE2'] > readnumber)]
        df = df.drop('AvgSAMPLE2', axis=1)
        df = df.drop('AvgSAMPLE1', axis=1)

    df['sort_key'] = df.IncLevelDifference.abs()

    df = df.sort_values(by=['sort_key', 'FDR'], ascending=[False, True])
    df = df.drop('sort_key', axis=1)

    if topnumber:
        df = df[:topnumber]
    if outf:
        df.to_csv(outf, sep="\t", index=False)

    if up_out:
        up_df = df[df['IncLevelDifference'] > 0]
        up_df.to_csv(up_out, sep="\t", index=False)

    if down_out:
        down_df = df[df['IncLevelDifference'] < 0]
        down_df.to_csv(down_out, sep="\t", index=False)

    sys.stderr.write(f"Outfile {outf}: {df.shape[0]} genes.\n")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input", dest="input_file", type=str,
                        required=True, help="Input file gnerated by rmats")
    parser.add_argument("-o", "--output", dest="output_file",
                        type=str, default=None, help="Output file")
    parser.add_argument("-o1", "--upout", dest="up_out",
                        type=str, default=None, help="Up output file")
    parser.add_argument("-o2", "--downout", dest="down_out",
                        type=str, default=None, help="Down output file") 
    parser.add_argument("-p", "--pvalue", dest="pvalue",
                        type=float, default=None, help="Pvalue cutoff")
    parser.add_argument("-f", "--fdr", dest="fdr",
                        type=float, default=None, help="FDR value cutoff")
    parser.add_argument("-l", "--inclevel", dest="inclevel",
                        type=float, default=0.01, help="IncLevel Difference cutoff")
    parser.add_argument("-r", "--readnumber", dest="readnumber",
                        type=int, default=10, help="""
                        Average junction (inclusion junction and skipping junction) reads number in either group cutoff""")
    parser.add_argument("-n", "--number", dest="topnumber",
                        type=int, default=None, help="How many records to save")

    args = parser.parse_args()

    main(args.input_file, args.output_file, args.up_out, args.down_out,
         args.topnumber, args.pvalue, args.readnumber, args.fdr, args.inclevel)


if __name__ == "__main__":
    run()
