#!/usr/bin/env python3
import argparse
import os
import sys
from typing import Optional, Sequence, TextIO

from pybiotk.utils import read_table, ignore


def main(input: str = "-", outfa: Optional[TextIO] = None, outgtf: Optional[str] = None, seq_ids: Sequence[str] = None):
    ercc = read_table(input)
    outgtf = os.devnull if outgtf is None else outgtf
    seq_ids = set(seq_ids) if seq_ids is not None else False
    with outfa or sys.stdout as fa, open(outgtf, "w") as gtf:
        for _, row in ercc.iterrows():
            id = row["ERCC_ID"]
            seq = row["Sequence"]
            if seq_ids:
                if id.lstrip("ERCC-") not in seq_ids:
                    continue
            outfa.write(f">{id}\n{seq}\n")
            gtf.write(f"""{id}\tERCC\tgene\t1\t{len(seq)}\t.\t+\t.\tgene_id "{id}"; gene_version "1"; gene_name "{id}"; gene_source "ERCC"; gene_biotype "Spike-In";\n""")
            gtf.write(f"""{id}\tERCC\ttranscript\t1\t{len(seq)}\t.\t+\t.\tgene_id "{id}"; gene_version "1"; transcript_id "{id}"; transcript_version "1"; gene_name "{id}"; gene_source "ERCC"; gene_biotype "Spike-In"; transcript_name "{id}"; transcript_source "ERCC"; transcript_biotype "Spike-In";\n""")
            gtf.write(f"""{id}\tERCC\texon\t1\t{len(seq)}\t.\t+\t.\tgene_id "{id}"; gene_version "1"; transcript_id "{id}"; transcript_version "1"; exon_number "1"; gene_name "{id}"; gene_source "ERCC"; gene_biotype "Spike-In"; transcript_name "{id}"; transcript_source "ERCC"; transcript_biotype "Spike-In"; exon_id "{id}"; exon_version "1";\n""")


@ignore
def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(dest="input", type=str, nargs="?", default=(None if sys.stdin.isatty() else "-"), help="""ERCC table(https://www.thermofisher.com): "ERCC Controls Annotation: ERCC RNA Spike-In Control Mixes (English)". [stdin]""")
    parser.add_argument('-o', '--outfa', dest='outfa', type=argparse.FileType('w'), default=sys.stdout, help="output fasta.")
    parser.add_argument('-g', '--outgtf', dest="outgtf", type=str, default=None, help="output gtf.")
    parser.add_argument("-s", "--seq-ids", dest="seq_ids", nargs="+", default=None, help="ERCC ids.")
    args = parser.parse_args()
    if args.input is None:
        parser.print_help()
        sys.exit(1)
    main(**vars(args))
