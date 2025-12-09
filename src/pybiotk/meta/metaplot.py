#!/usr/bin/env python3
import os
import sys
import time
import argparse
import itertools
import numpy as np
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['svg.fonttype'] = 'none'
import matplotlib.pyplot as plt
from matplotlib import colors as pltcolors
from matplotlib.font_manager import FontProperties
import scipy.stats as st
from scipy.signal import savgol_filter
from typing import List, Literal, Dict, Optional
from stream import count
from pybiotk.io import Openbed
from pybiotk.annodb import MergedTranscript
from pybiotk.utils import logging
from pybiotk.meta.task import scale_regions_task, reference_point_task
from pybiotk.meta.merge_transcript import load_gene


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


def smooth_window(size: int, bins: int = 10) -> int:
    window = size // bins
    if window % 2 == 0:
        window += 1
    if window < 5:
        window = 5
    return window


def cm2inch(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)


def point_metaplot(gene: Dict[str, List[MergedTranscript]],
                   numpy_file: str,
                   out_plot: str,
                   bw_fwd: List[str],
                   group: Optional[List[str]] = None,
                   labels: Optional[List[str]] = None,
                   loci: Literal['TES', 'TSS'] = 'TES',
                   method: Literal['reads', 'coverage'] = 'reads',
                   upStream: int = 2000, downStream: int = 5000,
                   bins: int = 10,
                   ylab: str = "Average density",
                   plotsem: bool = False,
                   smooth: bool = False,
                   smooth_k: int = 5
                   ):
    plot_height = 7
    plot_width = 11
    color_list = ["#3b5076", "#4cb6c4", "#1f9c81",  "#eb9277", "#db5342"]
    plt.rcParams['font.size'] = 8.0
    font_p = FontProperties()
    font_p.set_size('small')
    with open(numpy_file, 'rb') as f:
        a = np.load(f)
    fig = plt.figure(figsize=cm2inch(plot_width, plot_height))
    ax = fig.add_subplot()
    x = np.arange((upStream+downStream)//bins) if method == 'reads' else np.arange(upStream+downStream-bins+1)
    if labels is not None:
        assert len(labels) == len(bw_fwd)
        names = labels
    else:
        names = [os.path.basename(x).split(".")[0] for x in bw_fwd]
    if group is not None:
        assert len(group) == len(bw_fwd)
        group_d = {}
        for n, g in enumerate(group):
            if g not in group_d:
                group_d[g] = [a[n]]
            else:
                group_d[g].append(a[n])
        numlines = len(group_d)
        if numlines == 2:
            color_list = ["#3b5076", "#db5342"]
        if numlines == 3:
            color_list = ["#3b5076", "#eb9277", "#db5342"]
        if len(color_list) < numlines:
            logging.warning("\nThe given list of colors is too small, "
                     "at least {} colors are needed\n".format(numlines))
            color_list = None
        if color_list is None:
            cmap_plot = plt.get_cmap('jet')
            color_list = cmap_plot(np.arange(numlines, dtype=float) / float(numlines))
        for color in color_list:
            if not pltcolors.is_color_like(color):
                sys.exit("\nThe color name {} is not valid. Check "
                         "the name or try with a html hex string "
                         "for example #eeff22".format(color))
        for i, g in enumerate(group_d):
            data = group_d[g]
            color = color_list[i]
            if isinstance(color, np.ndarray):
                color = pltcolors.to_hex(color, keep_alpha=True)
            if len(data) == 1:
                y = data[0]
                ysize = len(y)
                if smooth and ysize > 10:
                    y = savgol_filter(y, smooth_window(ysize), smooth_k)
                ax.plot(x, y, label=g, color=color, alpha=1.0, linewidth=1)
            else:
                data = np.array(group_d[g])
                y = data.mean(axis=0)
                ysize = len(y)
                if smooth and ysize > 10:
                    y = savgol_filter(y, smooth_window(ysize), smooth_k)
                ax.plot(x, y, label=g, color=color, alpha=1.0, linewidth=1)
                if plotsem:
                    # low_CI_bound, high_CI_bound = st.t.interval(0.95, df=len(data)-1, loc=y, scale=st.sem(data))
                    # ax.fill_between(x, low_CI_bound, high_CI_bound, color="grey", alpha=0.3)
                    sem = st.sem(data)
                    ax.fill_between(x, y-sem, y+sem, color="grey", alpha=0.3)
    else:
        numlines = len(a)
        if numlines == 2:
            color_list = ["#3b5076", "#db5342"]
        if numlines == 3:
            color_list = ["#3b5076", "#eb9277", "#db5342"]
        if len(color_list) < numlines:
            logging.warning("\nThe given list of colors is too small, "
                            "at least {} colors are needed\n".format(numlines))
            color_list = None
        if color_list is None:
            cmap_plot = plt.get_cmap('jet')
            color_list = cmap_plot(np.arange(numlines, dtype=float) / float(numlines))
        for color in color_list:
            if not pltcolors.is_color_like(color):
                sys.exit("\nThe color name {} is not valid. Check "
                         "the name or try with a html hex string "
                         "for example #eeff22".format(color))
        for n, j in enumerate(a):
            color = color_list[n]
            if isinstance(color, np.ndarray):
                color = pltcolors.to_hex(color, keep_alpha=True)
            jsize = len(j)
            if smooth and jsize > 10:
                j = savgol_filter(j, smooth_window(jsize), smooth_k)
            ax.plot(x, j, label=names[n], color=color, alpha=1.0, linewidth=1)
    ax.set_xlim(0, max(x))
    n = (gene['+'] | count) + (gene['-'] | count)
    ticks = []
    labels = []
    k = upStream
    while True:
        k += 5000
        if k >= (upStream + downStream):
            break
        ticks.append(k)
        labels.append(f"+{int((k-upStream)/1000)}")
    ticks.append(upStream+downStream)
    ticks = [i/bins for i in ticks]
    labels.append(f"+{int(downStream/1000)}(kb)")
    if upStream >= 1000:
        ax.set_xticks(list(itertools.chain([0, upStream/bins], ticks)))
        ax.set_xticklabels(list(itertools.chain([f'-{int(upStream/1000)}', loci], labels)))
    else:
        ax.set_xticks(list(itertools.chain([upStream/bins], ticks)))
        ax.set_xticklabels(list(itertools.chain([loci], labels)))
    ax.set_xlabel(f'Genomic region downstream {loci}')
    ax.set_ylabel(f'{ylab}')
    ax.set_title(f'n={n}   ')
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    plt.legend(loc=1, frameon=False, markerscale=0.5)
    plt.tight_layout()
    plt.savefig(out_plot, dpi=300)
    plt.close()


def body_metaplot(gene: Dict[str, List[MergedTranscript]],
                  numpy_file: str,
                  out_plot: str,
                  bw_fwd: List[str],
                  group: Optional[List[str]] = None,
                  labels: Optional[List[str]] = None,
                  upStream: int = 2000,
                  downStream: int = 2000,
                  length: int = 5000,
                  bins: int = 10,
                  ylab: str = 'Average density',
                  plotsem: bool = False,
                  smooth: bool = False,
                  smooth_k: int = 5):
    plot_height = 7
    plot_width = 11
    color_list = ["#3b5076", "#4cb6c4", "#1f9c81",  "#eb9277", "#db5342"]
    plt.rcParams['font.size'] = 8.0
    font_p = FontProperties()
    font_p.set_size('small')
    with open(numpy_file, 'rb') as f:
        a = np.load(f)
    fig = plt.figure(figsize=cm2inch(plot_width, plot_height))
    ax = fig.add_subplot()
    x = np.arange((upStream+downStream+length)//bins)
    if labels is not None:
        assert len(labels) == len(bw_fwd)
        names = labels
    else:
        names = [os.path.basename(x).split(".")[0] for x in bw_fwd]
    if group is not None:
        assert len(group) == len(bw_fwd)
        group_d = {}
        for n, g in enumerate(group):
            if g not in group_d:
                group_d[g] = [a[n]]
            else:
                group_d[g].append(a[n])
        numlines = len(group_d)
        if numlines == 2:
            color_list = ["#3b5076", "#db5342"]
        if numlines == 3:
            color_list = ["#3b5076", "#eb9277", "#db5342"]
        if color_list is None:
            cmap_plot = plt.get_cmap('jet')
            color_list = cmap_plot(np.arange(numlines, dtype=float) / float(numlines))
        if len(color_list) < numlines:
            sys.exit("\nThe given list of colors is too small, "
                     "at least {} colors are needed\n".format(numlines))
        for color in color_list:
            if not pltcolors.is_color_like(color):
                sys.exit("\nThe color name {} is not valid. Check "
                         "the name or try with a html hex string "
                         "for example #eeff22".format(color))
        for i, g in enumerate(group_d):
            data = group_d[g]
            color = color_list[i]
            if isinstance(color, np.ndarray):
                color = pltcolors.to_hex(color, keep_alpha=True)
            if len(data) == 1:
                y = data[0]
                ysize = len(y)
                if smooth and ysize > 10:
                    y = savgol_filter(y, smooth_window(ysize), smooth_k)
                ax.plot(x, y, label=g, color=color, alpha=1.0, linewidth=1)
            else:
                data = np.array(group_d[g])
                y = data.mean(axis=0)
                ysize = len(y)
                if smooth and ysize > 10:
                    y = savgol_filter(y, smooth_window(ysize), smooth_k)
                ax.plot(x, y, label=g, color=color, alpha=1.0, linewidth=1)
                if plotsem:
                    # low_CI_bound, high_CI_bound = st.t.interval(0.95, df=len(data)-1, loc=y, scale=st.sem(data))
                    # ax.fill_between(x, low_CI_bound, high_CI_bound, color="grey", alpha=0.3)
                    sem = st.sem(data)
                    ax.fill_between(x, y-sem, y+sem, color="grey", alpha=0.3)
    else:
        numlines = len(a)
        if numlines == 2:
            color_list = ["#3b5076", "#db5342"]
        if numlines == 3:
            color_list = ["#3b5076", "#eb9277", "#db5342"]
        if len(color_list) < numlines:
            logging.warning("\nThe given list of colors is too small, "
                     "at least {} colors are needed\n".format(numlines))
            color_list = None
        if color_list is None:
            cmap_plot = plt.get_cmap('jet')
            color_list = cmap_plot(np.arange(numlines, dtype=float) / float(numlines))
        for color in color_list:
            if not pltcolors.is_color_like(color):
                sys.exit("\nThe color name {} is not valid. Check "
                         "the name or try with a html hex string "
                         "for example #eeff22".format(color))
        for n, j in enumerate(a):
            color = color_list[n]
            if isinstance(color, np.ndarray):
                color = pltcolors.to_hex(color, keep_alpha=True)
            jsize = len(j)
            if smooth and jsize > 10:
                j = savgol_filter(j, smooth_window(jsize), smooth_k)
            ax.plot(x, j, label=names[n], color=color, alpha=1.0, linewidth=1)
    ax.set_xlim(0, max(x))
    n = (gene['+'] | count) + (gene['-'] | count)
    skip = length/(3*bins)
    ax.set_xticks([0, upStream/bins, upStream/bins+skip, upStream/bins+2*skip, (length+upStream)/bins, (upStream+downStream+length)/bins])
    ax.set_xticklabels([f'-{int(upStream/1000)}', 'TSS', '33%', '66%', 'TES', f'+{int(downStream/1000)}(kb)'])
    ax.set_xlabel('Genomic region gene body')
    ax.set_ylabel(f'{ylab}')
    ax.set_title(f'n={n}   ')
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # handles, labels = ax.get_legend_handles_labels()
    # labels = legend
    # ax.legend(handles[::len(handles)//2], labels)
    plt.legend(loc=1, frameon=False, markerscale=0.5)
    plt.tight_layout()
    plt.savefig(out_plot, dpi=300)
    plt.close()


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
    parser_s.add_argument("--labels", dest="labels", nargs="+", type=str, default=None, help="labels")
    parser_s.add_argument("-o", dest="output", type=str, required=True, help="outfig")
    parser_s.add_argument("-b", "--upStream", dest="upStream", type=int, default=2000, help="upStream cutoff.")
    parser_s.add_argument("-a", "--downStream", dest="downStream", type=int, default=2000, help="downStream cutoff.")
    parser_s.add_argument("-m", "--length", dest="length", type=int, default=5000, help="gene body scale length.")
    parser_s.add_argument("-y", "--ylab", dest="ylab", type=str, default="Average density", help="y label.")
    parser_s.add_argument("--sem", dest="sem", action="store_true", help="plot sem.")
    parser_s.add_argument("--smooth", dest="smooth", action="store_true", help="plot smooth.")
    parser_s.add_argument("--smooth-k", dest="smooth_k", type=int, default=5, help="smooth use k-th order polynomial fit.")

    parser_p = subparsers.add_parser("reference-point")
    parser_p.add_argument("--load", dest="load_pickle", type=str, default=None, help="load pickle file.")
    parser_p.add_argument("--reference", dest="reference", type=str, default=None, help="reference bedfile")
    parser_p.add_argument("-f", "--fwd", dest="fwd", nargs="+", required=True, type=str,
                          help="input forward bigwigs or collapsed bigwigs.")
    parser_p.add_argument("-r", "--rev", dest="rev", nargs="+", default=None, type=str,
                          help="input reverse bigwigs or None.")
    parser_p.add_argument("-bs", dest="bins", default=10, type=int, help="bins")
    parser_p.add_argument("-g", dest="group", nargs="+", type=str, default=None, help="groups")
    parser_p.add_argument("--labels", dest="labels", nargs="+", type=str, default=None, help="labels")
    parser_p.add_argument("-o", dest="output", type=str, required=True, help="outfig")
    parser_p.add_argument("-b", "--upStream", dest="upStream", type=int, default=2000, help="upStream cutoff.")
    parser_p.add_argument("-a", "--downStream", dest="downStream", type=int, default=10000, help="downStream cutoff.")
    parser_p.add_argument("-m", "--method", dest="method", type=str, default="reads", choices=["reads", "coverage"],
                          help="calculte method.")
    parser_p.add_argument("-l", "--loci", dest="loci", type=str, default="TES", choices=["TES", "TSS"],
                          help="reference point.")
    parser_p.add_argument("-y", "--ylab", dest="ylab", type=str, default="Average density", help="y label.")
    parser_p.add_argument("--sem", dest="sem", action="store_true", help="plot sem.")
    parser_p.add_argument("--smooth", dest="smooth", action="store_true", help="plot smooth.")
    parser_p.add_argument("--smooth-k", dest="smooth_k", type=int, default=5, help="smooth use k-th order polynomial fit.")

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
        body_metaplot(gene, outnp, outfig, args.fwd, args.group, args.labels, args.upStream, args.downStream, args.length, args.bins, args.ylab, args.sem, args.smooth, args.smooth_k)
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
        point_metaplot(gene, outnp, outfig, args.fwd, args.group, args.labels, args.loci, args.method, args.upStream, args.downStream, args.bins, args.ylab, args.sem, args.smooth, args.smooth_k)
    else:
        args = parser.parse_args(["-h"])
    end = time.perf_counter()
    logging.info(f"task finished in {end-start:.2f}s.")


if __name__ == '__main__':
    run()
