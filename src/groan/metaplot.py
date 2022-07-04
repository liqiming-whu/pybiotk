#!/usr/bin/env python3
import os
import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as st
from typing import List, Literal, Dict, Optional
from stream import count
from pybiotk.io import Openbed
from pybiotk.annodb import MergedTranscript
from pybiotk.utils import logging
from groan.task import scale_regions_task, reference_point_task
from groan.merge_transcript import load_gene


def bed2merge_transcript(bedpath: str) -> Dict[str, List[MergedTranscript]]:
    a = []
    b = []
    with Openbed(bedpath) as bedfile:
        for bed in bedfile:
            transcript = MergedTranscript(gene_name=bed.name,
                                          chrom=bed.chrom,
                                          start=bed.start,
                                          end=bed.end,
                                          strand=bed.strand,
                                          before=50000,
                                          after=50000)
            if bed.strand == '-':
                b.append(transcript)
            else:
                a.append(transcript)
    return {'+': a, '-': b}


def point_metaplot(gene: Dict[str, List[MergedTranscript]],
                   numpy_file: str,
                   out_plot: str,
                   bw_fwd: List[str],
                   group: Optional[List[str]] = None,
                   loci: Literal['TES', 'TSS'] = 'TES',
                   method: Literal['reads', 'coverage'] = 'reads',
                   upStream: int = 2000, downStream: int = 5000,
                   bins: int = 10, title=''):

    with open(numpy_file, 'rb') as f:
        a = np.load(f)
    _, ax = plt.subplots()
    x = np.arange(upStream+downStream) if method == 'reads' else np.arange(upStream+downStream-bins+1)
    names = [os.path.basename(x).split(".")[0] for x in bw_fwd]
    if group is not None:
        assert len(group) == len(bw_fwd)
        group_d = {}
        for n, g in enumerate(group):
            if g not in group_d:
                group_d[g] = [a[n]]
            else:
                group_d[g].append(a[n])
        for g in group_d:
            data = group_d[g]
            if len(data) == 1:
                ax.plot(x, data[0], label=g)
            else:
                data = np.array(group_d[g])
                y = data.mean(axis=0)
                ax.plot(x, y, label=g)
                # low_CI_bound, high_CI_bound = st.t.interval(0.95, df=len(data)-1, loc=y, scale=st.sem(data))
                # ax.fill_between(x, low_CI_bound, high_CI_bound, color="grey", alpha=0.3)
                sem = st.sem(data)
                ax.fill_between(x, y-sem, y+sem, color="grey", alpha=0.3)
    else:
        for n, j in enumerate(a):
            ax.plot(x, j, label=names[n])
    n = (gene['+'] | count) + (gene['-'] | count)
    if upStream >= 1000:
        ax.set_xticks([0, upStream, upStream+downStream])
        ax.set_xticklabels([f'-{upStream/1000}', loci, f'+{downStream/1000}'])
    else:
        ax.set_xticks([upStream, upStream+downStream])
        ax.set_xticklabels([loci, f'+{downStream/1000}'])
    ax.set_xlabel('Genomic region [Kb]')
    ax.set_ylabel(f'Average density\n[{method}]')
    ax.set_title(f'n={n}   '+title)
    plt.legend(loc=1, frameon=False)
    plt.tight_layout()
    plt.savefig(out_plot)


def body_metaplot(gene: Dict[str, List[MergedTranscript]],
                  numpy_file: str, out_file: str, bw_fwd: List[str],
                  group: Optional[List[str]] = None,
                  upStream: int = 2000, downStream: int = 2000,
                  length: int = 5000, bins: int = 10):
    with open(numpy_file, 'rb') as f:
        a = np.load(f)
    _, ax = plt.subplots()
    if group is not None:
        assert len(group) == len(bw_fwd)
        group_d = {}
        for n, g in enumerate(group):
            if g not in group_d:
                group_d[g] = [a[n]]
            else:
                group_d[g].append(a[n])
        x = np.arange((upStream+downStream+length)/bins)
        for g in group_d:
            y = np.array(group_d[g])
            ax.plot(x, y.mean(axis=0), label=g)
            ax.fill_between(x, y.min(axis=0), y.max(axis=0), color="grey", alpha=0.3)
    else:
        names = [os.path.basename(x).split(".")[0] for x in bw_fwd]
        for n, j in enumerate(a):
            ax.plot(np.arange((upStream+downStream+length)/bins), j, label=names[n])
    n = (gene['+'] | count) + (gene['-'] | count)
    ax.set_xticks([0, upStream/bins, (length+upStream)/bins, (upStream+downStream+length)/bins])
    ax.set_xticklabels([f'-{upStream/1000}', 'TSS', 'TES', f'+{downStream/1000}'])
    ax.set_xlabel('Genomic region [Kb]')
    ax.set_ylabel('Average density\n[reads]')
    ax.set_title(f'n={n}   ')
    # handles, labels = ax.get_legend_handles_labels()
    # labels = legend
    # ax.legend(handles[::len(handles)//2], labels)
    plt.legend(loc=1, frameon=False)
    plt.tight_layout()
    plt.savefig(out_file)


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers(dest='subparser_name')

    parser_s = subparsers.add_parser("scale-regions")
    parser_s.add_argument("--load", dest="load_pickle", type=str, default=None, help="load pickle file.")
    parser_s.add_argument("--reference", dest="reference", type=str, default=None, help="reference bedfile")
    parser_s.add_argument("-f", "--fwd", dest="fwd", nargs="+", required=True, type=str,
                          help="input forward bigwigs or collapsed bigwigs.")
    parser_s.add_argument("-r", "--rev", dest="rev", nargs="+", default=None, type=str,
                          help="input reverse bigwigs or None.")
    parser_s.add_argument("-bs", dest="bins", default=10, type=int, help="bins")
    parser_s.add_argument("-g", "--group", dest="group", nargs="+", type=str, default=None, help="groups")
    parser_s.add_argument("-o", dest="output", type=str, required=True, help="outfig")
    parser_s.add_argument("-b", "--upStream", dest="upStream", type=int, default=2000, help="upStream cutoff.")
    parser_s.add_argument("-a", "--downStream", dest="downStream", type=int, default=2000, help="downStream cutoff.")
    parser_s.add_argument("-m", "--length", dest="length", type=int, default=5000, help="gene body scale length.")

    parser_p = subparsers.add_parser("reference-point")
    parser_p.add_argument("--load", dest="load_pickle", type=str, default=None, help="load pickle file.")
    parser_p.add_argument("--reference", dest="reference", type=str, default=None, help="reference bedfile")
    parser_p.add_argument("-f", "--fwd", dest="fwd", nargs="+", required=True, type=str,
                          help="input forward bigwigs or collapsed bigwigs.")
    parser_p.add_argument("-r", "--rev", dest="rev", nargs="+", default=None, type=str,
                          help="input reverse bigwigs or None.")
    parser_p.add_argument("-bs", dest="bins", default=10, type=int, help="bins")
    parser_p.add_argument("-g", dest="group", nargs="+", type=str, default=None, help="groups")
    parser_p.add_argument("-o", dest="output", type=str, required=True, help="outfig")
    parser_p.add_argument("-b", "--upStream", dest="upStream", type=int, default=2000, help="upStream cutoff.")
    parser_p.add_argument("-a", "--downStream", dest="downStream", type=int, default=10000, help="downStream cutoff.")
    parser_p.add_argument("-m", "--method", dest="method", type=str, default="reads", choices=["reads", "coverage"],
                          help="calculte method.")
    parser_p.add_argument("-l", "--loci", dest="loci", type=str, default="TES", choices=["TES", "TSS"],
                          help="reference point.")

    return parser


def run():
    parser = parse_args()
    args = parser.parse_args()
    start = time.perf_counter()
    if args.subparser_name == "scale-regions":
        if args.reference is None and args.load_pickle is None:
            args = parser.parse_args(['-h'])
        if args.reference is not None:
            gene = bed2merge_transcript(args.reference)
        else:
            gene = load_gene(args.load_pickle)
        outfig = args.output
        outnp = os.path.splitext(outfig)[0] + ".np"
        logging.info("choose scale_regions mode.")
        logging.info("start to calculte np matrix...")
        scale_regions_task(outnp, gene, args.fwd, args.rev, args.upStream, args.downStream, args.length, args.bins)
        logging.info("start to plot gene body...")
        body_metaplot(gene, outnp, outfig, args.fwd, args.group, args.upStream, args.downStream, args.length, args.bins)
        logging.info(f"figure saved in {outfig}.")
    elif args.subparser_name == "reference-point":
        if args.reference is None and args.load_pickle is None:
            args = parser.parse_args(['-h'])
        if args.reference is not None:
            gene = bed2merge_transcript(args.reference)
        else:
            gene = load_gene(args.load_pickle)
        outfig = args.output
        outnp = os.path.splitext(outfig)[0] + ".np"
        logging.info("choose reference-point mode.")
        logging.info("start to calculte np matrix...")
        reference_point_task(outnp, gene, args.fwd, args.rev, args.loci, args.method, args.upStream, args.downStream, args.bins)
        logging.info(f"start to plot {args.loci} {args.method} ...")
        point_metaplot(gene, outnp, outfig, args.fwd, args.group, args.loci, args.method, args.upStream, args.downStream, args.bins)
    else:
        args = parser.parse_args(["-h"])
    end = time.perf_counter()
    logging.info(f"task finished in {end-start:.2f}s.")


if __name__ == '__main__':
    run()
