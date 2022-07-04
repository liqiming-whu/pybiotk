#!/usr/bin/env python3
"""
GroAn: a tool for GroSeq Analysis.
"""
import os
import sys
import argparse
import time
from pprint import pprint
from typing import Dict, List, Literal, Optional, Union
from stream.pipe import mkdir, cd
from pybiotk.utils import logging
from groan.merge_transcript import merge_transcripts, MergedTranscript, load_gene
from groan.filter_gene import filter_gene
from groan.expression_qc import expression_qc
from groan.signal_continuity_qc import signal_continuity_qc
from groan.task import scale_regions_task, reference_point_task
from groan.metaplot import body_metaplot, point_metaplot
from groan.read_through_stats import readthrough_stats, readthrough_interpolation_coverage, plot_readthrough_stats


def parse_args():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-f", dest="fwd", nargs="+", required=True, help="forward bw files.")
    parser.add_argument("-r", dest="rev", nargs="+", help="reverse bw files.")
    parser.add_argument("-g", dest="gtf", type=str, required=True, help="gtf file.")
    parser.add_argument("-s", "--species", dest="species", type=str, default=None, choices=["hg38", "mm10"],
                        help="use default chrom_lenth file: hg38 or mm10. override --chrom_lenth option.")
    parser.add_argument("--chrom_length", dest="chrom_length", type=str, default=None, help="chrom_length file.")
    parser.add_argument("-o", dest="outdir", type=str, default="results", help="outdir.")
    parser.add_argument("-m", dest="method", type=str, choices=["coverage", "reads"], default="reads", help="statics method.")
    parser.add_argument("-v", dest="value", type=float, default=0.25, help="minimum coverage or reads value.")
    parser.add_argument("-l", "--loci", dest="loci", type=str, default="tes", choices=["body", "tss", "tes"], help="calc in which loci")
    parser.add_argument("-d", "--distance", dest="distance", type=int, default=3000, help="tss/tes +- distance.")
    parser.add_argument("-bs", dest="bins", default=10, type=int, help="bins")
    parser.add_argument("--filter_gene_by", dest="filter_method", type=str, default="mean", choices=["value", "top", "mean", "quantile"], help="filter method")
    parser.add_argument("--filter_signal_continuity", dest="filter_signal", type=str, default=None, choices=["value", "top", "mean", "quantile"], help="filter signal method")
    parser.add_argument("--filter_signal_value", dest="filter_signal_value", type=float, default=0.75, help="filter signal value")
    parser.add_argument("--group", dest="group", nargs="+", type=str, default=None, help="samples group.")
    parser.add_argument("--gene_length", dest="gene_length", type=int, default=5000, help="gene_length cutoff.")
    parser.add_argument("--downStream", dest="downStream", type=int, default=10000, help="downStream cutoff.")
    parser.add_argument("--downStart", dest="downStart", type=int, default=0, help="downStream distance from TES.")
    args = parser.parse_args()
    return args


def run_merge_transcripts(gtf:str, outdir:str, species: Optional[str] = None, chrom_length: Optional[str] = None):
    mkdir(outdir)
    cd(outdir)
    logging.info("Convert gtf to merged transcripts, make sure the gtf is sorted.")
    start = time.perf_counter()
    merged_transcripts = merge_transcripts(gtf, "merged_transcripts.bed", "merge_tanscripts.pickle",
                      escape_gene_name_startswith=["MIR"], chrom_length_file=chrom_length, species=species)
    end = time.perf_counter()
    logging.info(f"Convert completed, use {end-start} s.")
    return merged_transcripts


def run_filter_gene(merged_transcripts: Dict[str, List[MergedTranscript]],
                    outdir: str, gene_length=5000, downStream=10000):
    mkdir(outdir)
    cd(outdir)
    logging.info("Start to filter gene ...")
    start = time.perf_counter()
    filtered_gene = filter_gene(merged_transcripts,
                                "filtered_gene.bed",
                                "filtered_gene.pickle",
                                gene_length=gene_length,
                                downStream=downStream,
                                gene_types=["protein_coding", "lncRNA"])
    end = time.perf_counter()
    logging.info(f"Filter gene completed, use {end-start} s.")
    return filtered_gene


def run_expression_qc(filtered_gene: Dict[str, List[MergedTranscript]],
                      bw_fwd: List[str], bw_rev: List[str], outdir: str,
                      loci: Literal["body", "tss", "tes"] = "body",
                      method: Literal["coverage", "reads"] = "coverage",
                      distance: int = 3000, value: Union[int, float] = 0.25,
                      filter_method: Literal["value", "top", "mean", "quantile"] = "mean"):
    mkdir(outdir)
    cd(outdir)
    logging.info("Start expression QC ...")
    start = time.perf_counter()
    filtered_high_expressed_gene = expression_qc(
        filtered_gene, "expression.xls",
        bw_fwd, bw_rev, "filtered_high_expressed_genes.pickle",
        loci=loci,
        method=method,
        distance=distance,
        value=value,
        filter_method=filter_method)
    end = time.perf_counter()
    logging.info(f"expression QC completed, use {end-start} s.")
    return filtered_high_expressed_gene


def run_signal_continuity_qc(filtered_gene: Dict[str, List[MergedTranscript]],
                             bw_fwd: List[str], bw_rev: List[str], outdir: str,
                             downStream: int = 10000, downStart: int=0, value: Union[int, float] = 0.75,
                             filter_method: Literal["value", "top", "mean", "quantile"] = "quantile"):
    mkdir(outdir)
    cd(outdir)
    logging.info("Start signal continuity QC ...")
    start = time.perf_counter()
    filtered_high_signal_continuity_gene = signal_continuity_qc(
        filtered_gene, "signal_continuity.xls",
        bw_fwd, bw_rev, "filtered_high_signal_continuity_genes.pickle",
        downStream=downStream,
        shift=downStart,
        value=value,
        filter_method=filter_method)
    end = time.perf_counter()
    logging.info(f"signal continuity QC completed, use {end-start} s.")
    return filtered_high_signal_continuity_gene


def run_metaplot(filtered_high_expressed_gene: Dict[str, List[MergedTranscript]],
                 outdir: str, bw_fwd: List[str], bw_rev: List[str],
                 downStream: int = 10000, bins=10, group: Optional[List[str]] = None):
    mkdir(outdir)
    cd(outdir)
    logging.info("Start metaplot task ...")
    start = time.perf_counter()
    logging.info("start to plot gene_body -2k, 5k, +2k ...")
    if os.path.exists("gene_body.np"):
        logging.info("using exists file: gene_body.np")
    else:
        scale_regions_task("gene_body.np", filtered_high_expressed_gene, bw_fwd, bw_rev, 2000, 2000, 5000, bins)
    body_metaplot(filtered_high_expressed_gene, "gene_body.np", "gene_body.pdf", bw_fwd, group, 2000, 2000, 5000, bins)
    logging.info(f"start to plot TES coverage -2k, TES, +{downStream/1000:.0f}k ...")
    if os.path.exists("tes.coverage.np"):
        logging.info("using exists file: tes.coverage.np")
    else:
        reference_point_task("tes.coverage.np", filtered_high_expressed_gene, bw_fwd, bw_rev, "TES", "coverage", 2000, downStream, bins)
    point_metaplot(filtered_high_expressed_gene, "tes.coverage.np", "tes.coverage.pdf", bw_fwd, group, "TES", "coverage", 2000, downStream, bins)
    logging.info(f"start to plot TES reads -2k, TES, +{downStream/1000:.0f}k ...")
    if os.path.exists("tes.reads.np"):
        logging.info("using exists file: tes.reads.np")
    else:
        reference_point_task("tes.reads.np", filtered_high_expressed_gene, bw_fwd, bw_rev, "TES", "reads", 2000, downStream, bins)
    point_metaplot(filtered_high_expressed_gene, "tes.reads.np", "tes.reads.pdf", bw_fwd, group, "TES", "reads", 2000, downStream, bins)
    logging.info(f"start to plot TES reads TES, +{downStream/1000:.0f}k ...")
    if os.path.exists("tes_start.reads.np"):
        logging.info("using exists file: tes_start.reads.np")
    else:
        reference_point_task("tes_start.reads.np", filtered_high_expressed_gene, bw_fwd, bw_rev, "TES", "reads", 0, downStream, bins)
    point_metaplot(filtered_high_expressed_gene, "tes_start.reads.np", "tes_start.reads.pdf", bw_fwd, group, "TES", "reads", 0, downStream, bins)
    end = time.perf_counter()
    logging.info(f"metaplot task compelted, figure saved in gene_body.pdf, tes.coverage.pdf, tes.reads.pdf, tes_start.reads.pdf, use {end-start} s.")


def run_readthrough(filtered_high_expressed_gene: Dict[str, List[MergedTranscript]],
                    outdir: str, bw_fwd: List[str], bw_rev: List[str],
                    downStream: int = 10000, downStart: int = 0):
    mkdir(outdir)
    cd(outdir)
    logging.info("Start readthrough task ...")
    start = time.perf_counter()
    if os.path.exists("readthrough_coverage.xls"):
        logging.info("using exists file: readthrough_coverage.xls")
    else:
        readthrough_stats(filtered_high_expressed_gene, "readthrough_coverage.xls", bw_fwd, bw_rev,
                          downStream=downStream, shift=downStart, method="coverage")
    plot_readthrough_stats("readthrough_coverage.xls", "readthrough_coverage.pdf")
    if os.path.exists("readthrough_reads.xls"):
        logging.info("using exists file: readthrough_reads.xls")
    else:
        readthrough_stats(filtered_high_expressed_gene, "readthrough_reads.xls", bw_fwd, bw_rev,
                          downStream=downStream, shift=downStart, method="reads")
    plot_readthrough_stats("readthrough_reads.xls", "readthrough_reads.pdf")
    end = time.perf_counter()
    logging.info(f"readthrough task compelted, use {end-start} s.")

    logging.info("Start readthrough interpolation coverage task ...")
    start = time.perf_counter()
    if os.path.exists("readthrough_interpolation_coverage.xls"):
        logging.info("using exists file: readthrough_interpolation_coverage.xls")
    else:
        readthrough_interpolation_coverage(filtered_high_expressed_gene, "readthrough_interpolation_coverage.xls", bw_fwd, bw_rev,
                                           downStream=downStream, shift=downStart)
    plot_readthrough_stats("readthrough_interpolation_coverage.xls", "readthrough_interpolation_coverage.pdf")
    end = time.perf_counter()
    logging.info(f"readthrough interpolation coverage task compelted, use {end-start} s.")


def main():
    import __main__
    __main__.MergedTranscript = MergedTranscript
    args = parse_args()
    logging.info("Start GroAn ...")
    logging.info("check parameters: ")
    pprint(args.__dict__, stream=sys.stderr)
    gtf = os.path.abspath(args.gtf)
    fwd = [os.path.abspath(x) for x in args.fwd]
    method = args.method
    species = args.species
    chrom_length = args.chrom_length
    value = args.value
    loci = args.loci
    distance = args.distance
    bins = args.bins
    filter_method = args.filter_method
    filter_signal = args.filter_signal
    filter_signal_value = args.filter_signal_value
    group = args.group
    gene_length = args.gene_length
    downStream = args.downStream
    downStart = args.downStart
    if args.rev is not None:
        rev = [os.path.abspath(x) for x in args.rev]
        assert len(fwd) == len(rev)
    else:
        rev = args.rev

    outdir = os.path.abspath(args.outdir)
    mkdir(outdir)
    cd(outdir)
    merged_transcripts_dir = "merged_transcripts"
    filtered_gene_dir = "filtered_gene"
    qc_dir = "qc"
    metaplot_dir = "metaplot"
    readthrough_dir = "read_through"
    if os.path.exists(os.path.join(merged_transcripts_dir, "merge_tanscripts.pickle")):
        merged_transcripts = load_gene(os.path.join(merged_transcripts_dir, "merge_tanscripts.pickle"))
        logging.info("skip merge transcript, using exist file: merge_tanscripts.pickle.")
    else:
        merged_transcripts = run_merge_transcripts(gtf, merged_transcripts_dir, species, chrom_length)
    cd(outdir)

    if os.path.exists(os.path.join(filtered_gene_dir, "filtered_gene.pickle")):
        filtered_gene = load_gene(os.path.join(filtered_gene_dir, "filtered_gene.pickle"))
        logging.info("skip filter gene, using exist file: filter_gene.pickle")
    else:
        filtered_gene = run_filter_gene(merged_transcripts, filtered_gene_dir, gene_length, downStream)
    cd(outdir)

    filtered_high_expressed_gene = run_expression_qc(filtered_gene, fwd, rev, qc_dir, loci=loci, method=method, distance=distance, value=value, filter_method=filter_method)
    cd(outdir)
    
    if filter_signal is not None:
        filtered_high_expressed_gene = run_signal_continuity_qc(filtered_high_expressed_gene, fwd, rev, qc_dir, downStream, downStart, filter_signal_value, filter_signal)
        cd(outdir)
    run_metaplot(filtered_high_expressed_gene, metaplot_dir, fwd, rev, downStream=downStream, group=group, bins=bins)
    cd(outdir)
    run_readthrough(filtered_high_expressed_gene, readthrough_dir, fwd, rev, downStream=downStream, downStart=downStart)


if __name__ == '__main__':
    main()
