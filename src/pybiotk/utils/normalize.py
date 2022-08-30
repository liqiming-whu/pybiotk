#!/usr/bin/env python3
import argparse
import os
import sys
from typing import Sequence, Literal

import pandas as pd

from pybiotk.utils import write_table, ignore


def check_for_nonnumeric(pd_series=None):
    if pd.to_numeric(pd_series, errors='coerce').isna().sum() == 0:
        return 0
    else:
        return 1


class Norm:
    def __init__(self, df: pd.DataFrame):
        self.df = df.dropna()

    def cpm(self):
        df = self.df
        # check for non-numeric values
        for i in df.columns:
            assert check_for_nonnumeric(df[i]) == 0, \
                'dataframe contains non-numeric values in {} column'.format(i)
        self.lib_size = df.sum()
        self.cpm_norm = (df * 1e6) / df.sum()

        return self.cpm_norm

    def rpkm(self, gl="Length"):
        df = self.df
        assert gl is not None, "Provide column name for gene length in bp"
        # check for non-numeric values
        for i in df.columns:
            assert check_for_nonnumeric(df[i]) == 0, \
                'dataframe contains non-numeric values in {} column'.format(i)
        self.rpkm_norm = (df.div(df[gl], axis=0) * 1e9) / df.sum()
        self.rpkm_norm = self.rpkm_norm.drop([gl], axis=1)

        return self.rpkm_norm

    def tpm(self, gl="Length"):
        df = self.df
        assert gl is not None, "Provide column name for gene length in bp"
        # check for non-numeric values
        for i in df.columns:
            assert check_for_nonnumeric(df[i]) == 0, \
                'dataframe contains non-numeric values in {} column'.format(i)
        # gene length must be in bp
        self.a = df.div(df[gl], axis=0) * 1e3
        self.tpm_norm = (self.a * 1e6) / self.a.sum()
        self.tpm_norm = self.tpm_norm.drop([gl], axis=1)

        return self.tpm_norm


def read_data(path: str, cols: Sequence[int] = (0, 5, 6)) -> pd.DataFrame:
    df = pd.read_table(path, comment="#")
    df = df.iloc[:, cols]

    if len(cols) == 3:
        samples = [os.path.basename(i).split(".")[0] for i in df.columns[2:]]
        df.columns = ["Gene", "Length"] + samples
    else:
        samples = [os.path.basename(i).split(".")[0] for i in df.columns[1:]]
        df.columns = ["Gene"] + samples
    df = df.set_index("Gene")

    return df


def main(path: str, cols: Sequence[int], method: Literal["CPM", "RPKM", "TPM"], outpath):
    df = read_data(path, cols)

    norm = Norm(df)

    if method == "CPM":
        if "Length" in df.columns:
            norm.df = norm.df.drop(['Length'], axis=1)
        data = norm.cpm()
    elif method == "RPKM":
        data = norm.rpkm()
    elif method == "TPM":
        data = norm.tpm()
    else:
        raise RuntimeError(f"{method} normalize not support!")

    data = data[data.apply(sum, axis=1) != 0]
    data = data.round(2)
    data = data.sort_values(by=data.columns[0], ascending=False)
    write_table(data, outpath, index=True)


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str, help="featureCount out file.")
    parser.add_argument("-c", dest="columns", type=int, nargs="+", default=[0, 6, 5], help="geneid, count, length(required by RPKM and TPM) column, 0 based.")
    parser.add_argument("-m", dest="method", type=str, default="CPM", choices=("CPM", "RPKM", "TPM"), help="normalize method, [CPM|RPKM|TPM].")
    parser.add_argument("-o", dest="output", type=argparse.FileType('w'), default=sys.stdout, help="output file.")

    args = parser.parse_args()
    main(args.input, args.columns, args.method, args.output)


if __name__ == "__main__":
    run()
