#!/usr/bin/env python3
import os
import argparse
import time
import pysam
from pybiotk.io import Bam, Openbed, OpenFqGzip, check_bam_type, BamType
from pybiotk.intervals import GRangeTree, GRange, merge_intervals
from pybiotk.utils import cigartuples2blocks, logging


class NoChimAlign:
    def __init__(self, read: pysam.AlignedSegment):
        self.read = read
        self.name = read.query_name.replace(":", "_")
        self.sequence = read.get_forward_sequence()
        try:
            self.quality = pysam.qualities_to_qualitystring(
                read.get_forward_qualities())
        except TypeError:
            self.quality = 'F' * len(self.sequence)
        if not self.quality:
            self.quality = 'F' * len(self.sequence)
        self.cigar = read.cigartuples
        self.chrom = read.reference_name
        self.chromStart = read.reference_start
        self.blocks = None
        self.chimeric_JS_list = None

    def is_NoChimeric(self, knownJS: GRangeTree):
        if self.read.is_qcfail:
            return False
        if self.read.is_unmapped:
            return False

        soft_clip_len = 0
        for cigar_tuple in self.cigar:
            if cigar_tuple[0] == 4:
                soft_clip_len += cigar_tuple[1]

        if soft_clip_len >= 25:
            return False
        block = cigartuples2blocks(self.chromStart, self.cigar)
        self.blocks = block
        block_len = len(block)
        if block_len > 1:
            JS = []
            for indx in range(1, block_len):
                JS.append((block[indx-1][1], block[indx][0]))
            new_JS_li = []
            for js in JS:
                js_rec = GRange(self.chrom, js[0], js[1])
                overlap_js = knownJS.find(
                    js_rec.chrom, js_rec.start, js_rec.end)
                if not any([x == js_rec for x in overlap_js]):
                    new_JS_li.append((js_rec.start, js_rec.end))
            if not new_JS_li:
                return True
            else:
                self.chimeric_JS_list = new_JS_li
                return False
        else:
            return True

    def fix_chimeric_js(self):
        if self.blocks is not None and self.chimeric_JS_list is not None:
            return merge_intervals(self.blocks+self.chimeric_JS_list)


def main(bamfile: str, known_JS: str, fastq: str, junction_info: str = os.devnull):
    start = time.perf_counter()
    assert check_bam_type(bamfile) is BamType.SE
    bam = Bam(bamfile)
    logging.info("start to load knownJS tree...")
    js_tree = GRangeTree()
    with Openbed(known_JS) as bedfile:
        for bed in bedfile:
            rec = GRange(chrom=bed.chrom, start=bed.start, end=bed.end)
            js_tree.add_range(rec)
    logging.info("load completed.")
    logging.info(f"start to filter bamfile: {bamfile}")
    reads_count = no_js = 0
    with OpenFqGzip(fastq) as fq, open(junction_info, "w") as info:
        info.write("read_name\tchrom\tstrand\tblocks\tchimeric_js\n")
        for read in bam.iter():
            reads_count += 1
            align = NoChimAlign(read)
            if not align.is_NoChimeric(js_tree):
                fq.write_entry(name=align.name, sequence=align.sequence, quality=align.quality)
                info.write(f"{align.name}\t{align.chrom}\t{align.strand}\t{align.blocks}\t{align.chimeric_JS_list}\n")
            else:
                no_js += 1
    logging.info(f"total: {reads_count} reads. no junction reads mapped: {no_js}. rate: {no_js/reads_count:.2f}.")
    end = time.perf_counter()
    logging.info(f"task finished in {end-start:.2f}s.")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str, help="input bam.")
    parser.add_argument("-j", dest="knownjs", type=str, required=True, help="input knownjs table.")
    parser.add_argument("-o", dest="output", type=str, default=os.devnull, help="output fastq.")
    parser.add_argument("--js_info", dest="js_info", type=str, default=os.devnull, help="junction info.")
    args = parser.parse_args()
    main(args.input, args.knownjs, args.output, args.js_info)


if __name__ == "__main__":
    run()
