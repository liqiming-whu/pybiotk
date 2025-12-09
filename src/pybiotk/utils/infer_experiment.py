#!/usr/bin/env python3
'''=================================================================================================
Infer RNA-seq experiment design from SAM/BAM file. This module will determine if the RNA-seq
experiment is:
1) pair-end or single-end
2) if experiment is strand-specific, how reads were stranded.
 * For pair-end RNA-seq, there are two different ways to strand reads:
  i) 1++,1--,2+-,2-+
     read1 mapped to '+' strand indicates parental gene on '+' strand
     read1 mapped to '-' strand indicates parental gene on '-' strand
     read2 mapped to '+' strand indicates parental gene on '-' strand
     read2 mapped to '-' strand indicates parental gene on '+' strand
  ii) 1+-,1-+,2++,2--
     read1 mapped to '+' strand indicates parental gene on '-' strand
     read1 mapped to '-' strand indicates parental gene on '+' strand
     read2 mapped to '+' strand indicates parental gene on '+' strand
     read2 mapped to '-' strand indicates parental gene on '-' strand
 * For single-end RNA-seq, there are two different ways to strand reads:
  i) ++,--
     read mapped to '+' strand indicates parental gene on '+' strand
     read mapped to '-' strand indicates parental gene on '-' strand
  ii) +-,-+
     read mapped to '+' strand indicates parental gene on '-' strand
     read mapped to '-' strand indicates parental gene on '+' strand		

 NOTE:
   You don't need to know the RNA sequencing protocol before mapping your reads to the reference
   genome. Mapping your RNA-seq reads as if they were non-strand specific, this script can
   "guess" how RNA-seq reads were stranded.
================================================================================================='''

import os
import sys
import argparse
import collections

from pybiotk.io import GtfFile, Bam
from pybiotk.bx.intersection import IntervalTree
from pybiotk.utils import logging


def configure_experiment(samfile, gtf_file, sample_size, q_cut = 30, filter_strandness = None, outbamfile = os.devnull):
    '''Given a BAM/SAM file, this function will try to guess the RNA-seq experiment:
    1) single-end or pair-end
    2) strand_specific or not
    3) if it is strand-specific, what's the strand_ness of the protocol
    '''

    sam = Bam(samfile)
    outsam = Bam(outbamfile, 'wb', template=sam)
    filter_strandness = set(filter_strandness.split(',')) if filter_strandness else None
    count = 0
    p_strandness = collections.defaultdict(int)
    s_strandness = collections.defaultdict(int)
    # load reference gene model
    gene_ranges = {}
    logging.info("Reading reference " + gtf_file + ' ...')
    try:
        for bed in GtfFile(gtf_file).to_bed12():
            chrom = bed.chrom
            txStart = bed.thickStart
            txEnd = bed.thickEnd
            strand = bed.strand
            if chrom not in gene_ranges:
                gene_ranges[chrom] = IntervalTree()
            gene_ranges[chrom].add(txStart, txEnd, strand)
        assert gene_ranges
    except AssertionError:
        for gtf in GtfFile(gtf_file).iter_gene():
            chrom = gtf.chrom
            start = gtf.start
            end = gtf.end
            strand = gtf.strand
            if chrom not in gene_ranges:
                gene_ranges[chrom] = IntervalTree()
            gene_ranges[chrom].add(start, end, strand)
    assert gene_ranges
    logging.info("Done!")

    # read SAM/BAM file
    logging.info("Loading SAM/BAM file ... ")
    for aligned_read in sam:
        if aligned_read.is_qcfail:              # skip low quanlity
            continue
        if aligned_read.is_duplicate:           # skip duplicate read
            continue
        if aligned_read.is_secondary:           # skip non primary hit
            continue
        if aligned_read.is_unmapped:            # skip unmap read
            continue
        if aligned_read.mapping_quality < q_cut:
            continue

        chrom = aligned_read.reference_name
        if aligned_read.is_paired:
            if aligned_read.is_read1:
                read_id = '1'
            if aligned_read.is_read2:
                read_id = '2'
            if aligned_read.is_reverse:
                map_strand = '-'
            else:
                map_strand = '+'
            readStart = aligned_read.pos
            readEnd = readStart + aligned_read.qlen
            if chrom in gene_ranges:
                tmp = set(gene_ranges[chrom].find(readStart, readEnd))
                if len(tmp) == 0:
                    continue
                strand_from_gene = ':'.join(tmp)
                strandness = read_id + map_strand + strand_from_gene
                p_strandness[strandness] += 1        
                count += 1
        else:
            if aligned_read.is_reverse:
                map_strand = '-'
            else:
                map_strand = '+'
            readStart = aligned_read.reference_start
            readEnd = readStart + aligned_read.query_alignment_length
            if chrom in gene_ranges:
                tmp = set(gene_ranges[chrom].find(readStart, readEnd))
                if len(tmp) == 0:
                    if filter_strandness is not None:
                        outsam.write(aligned_read)
                    continue
                strand_from_gene = ':'.join(tmp)
                strandness = map_strand + strand_from_gene
                s_strandness[strandness] += 1
                count += 1
        if filter_strandness is None:
            if count >= sample_size:
                break
        else:
            if strandness in filter_strandness:
                outsam.write(aligned_read)

    outsam.close()
    logging.info("Finished")
    if filter_strandness is None:
        logging.info(f"Total {count} usable reads were sampled")
    else:
        logging.info("filter mode, all reads processed")

    protocol = "unknown"
    strandness = None
    spec1 = 0.0
    spec2 = 0.0
    other = 0.0
    if len(p_strandness) > 0 and len(s_strandness) == 0:
        protocol = "PairEnd"
        spec1 = (p_strandness['1++'] + p_strandness['1--'] + p_strandness['2+-'] + p_strandness['2-+'])/float(sum(p_strandness.values()))
        spec2 = (p_strandness['1+-'] + p_strandness['1-+'] + p_strandness['2++'] + p_strandness['2--'])/float(sum(p_strandness.values()))
        other = 1-spec1-spec2
    elif len(s_strandness) > 0 and len(p_strandness) == 0:
        protocol = "SingleEnd"
        spec1 = (s_strandness['++'] + s_strandness['--'])/float(sum(s_strandness.values()))
        spec2 = (s_strandness['+-'] + s_strandness['-+'])/float(sum(s_strandness.values()))
        other = 1-spec1-spec2
    else:
        protocol = "Mixture"
        spec1 = "NA"
        spec2 = "NA"
        other = "NA"

    if other < 0:
        other = 0.0
    if protocol == "PairEnd":
        print("\n\nThis is PairEnd Data")
        print(f"Fraction of reads failed to determine: {other:.4f}")
        print(f"Fraction of reads explained by \"1++,1--,2+-,2-+\": {spec1:.4f}")
        print(f"Fraction of reads explained by \"1+-,1-+,2++,2--\": {spec2:.4f}")
    elif protocol == "SingleEnd":
        print("\n\nThis is SingleEnd Data")
        print(f"Fraction of reads failed to determine: {other:.4f}")
        print(f"Fraction of reads explained by \"++,--\": {spec1:.4f}")
        print(f"Fraction of reads explained by \"+-,-+\": {spec2:.4f}")
    else:
        print("Unknown Data type")


def run():
    parser = argparse.ArgumentParser(
       description=__doc__,
       formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument(dest="input_file", type=str, nargs='?',
                        default=(None if sys.stdin.isatty() else "-"), help="Input alignment file in SAM or BAM format.[stdin]")
    parser.add_argument("-g", "--gtf", dest="gtf_file", type=str, required=True, help="Reference gene model in gtf fomat.")
    parser.add_argument("-s", "--sample-size", dest="sample_size",default=200000, help="Number of reads sampled from SAM/BAM file.")
    parser.add_argument("-q", "--mapq", dest="map_qual", type=int, default=30, help="Minimum mapping quality (phred scaled) for an alignment to be considered as \"uniquely mapped\".")
    parser.add_argument("-f", "--filter", dest="filter_strandness", type=str, default=None, choices=('1++,1--,2+-,2-+', '1+-,1-+,2++,2--', '++,--', '+-,-+'),
                        help="Filter reads with specified strandness, e.g. 1++,1--,2+-,2-+,1+-,1-+,2++,2--, when specified, --sample_size option will be ignored.")
    parser.add_argument("-o", "--outbam", dest="output_bam", type=str, default=os.devnull, help="output file.")
    args = parser.parse_args()

    if args.input_file is None:
        parser.print_help()
        sys.exit(1)
    configure_experiment(args.input_file, args.gtf_file, args.sample_size, args.map_qual, args.filter_strandness, args.output_bam)


if __name__ == "__main__":
    run()
