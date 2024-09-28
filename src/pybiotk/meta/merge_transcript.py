import argparse
import pickle
import time
from typing import Optional, Dict, List, Sequence
from pybiotk.io import GtfFile
from pybiotk.annodb import MergedTranscript, merge_transcripts_groupby_strand
from pybiotk.utils import load_chrom_length_dict, default_hg38_chrom_length_dict, default_mm10_chrom_length_dict, logging


def merge_transcripts(gtfpath: str,
                      merged_bed: Optional[str] = None,
                      dump: Optional[str] = None,
                      escape_gene_types: Sequence[str] = [],
                      escape_gene_name_startswith: Sequence[str] = [],
                      chrom_length_file: Optional[str] = None,
                      remove_strand_overlap: bool = False,
                      species: Optional[str] = None,
                      ) -> Dict[str, List[MergedTranscript]]:
    if species is not None:
        if species == "hg38":
            chrom_length_dict = default_hg38_chrom_length_dict()
        elif species == "mm10":
            chrom_length_dict = default_mm10_chrom_length_dict()
    elif chrom_length_file is not None:
        chrom_length_dict = load_chrom_length_dict(chrom_length_file)
    else:
        chrom_length_dict = None
    start = time.perf_counter()
    logging.info("start to merge transcripts...")
    with GtfFile(gtfpath) as gtf:
        merged_transcripts = merge_transcripts_groupby_strand(
            gtf,
            escape_gene_types=escape_gene_types,
            escape_gene_name_startswith=tuple(escape_gene_name_startswith),
            chrom_length_dict=chrom_length_dict,
            savebed=merged_bed,
            remove_strand_overlap=remove_strand_overlap
            )
    end = time.perf_counter()
    logging.info(f"saved {len(merged_transcripts['+']) + len(merged_transcripts['-'])} merged transcripts.")
    logging.info(f"task finished in {end-start:.2f}s.")
    if dump is not None:
        with open(dump, "wb") as f:
            pickle.dump(merged_transcripts, f)
    return merged_transcripts


def load_gene(picklefile):
    with open(picklefile, "rb") as f:
        merged_transcripts = pickle.load(f)
    return merged_transcripts


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-g", dest="gtf", type=str, required=True, help="sorted gtf file.")
    parser.add_argument("-b", dest="bed", type=str, default=None, help="output merged bed file.")
    parser.add_argument("-p", dest="pickle", type=str, default=None, help="output pickle dump file.")
    parser.add_argument("--escape_gene_types", dest="escape_gene_types", nargs="+", type=str, 
                        default=[], help="escape gene_types.")
    parser.add_argument("--escape_gene_name_startswith", dest="escape_gene_name_startswith", nargs="+",
                        type=str, default=[], help="escape gene_name startswith.")
    parser.add_argument("--remove_strand_overlap", dest="remove_strand_overlap", action="store_true",
                        help="remove strand overlap.")
    parser.add_argument("-s", "--species", dest="species", type=str, default=None, choices=["hg38", "mm10"],
                        help="use default chrom_lenth file: hg38 or mm10. override --chrom_lenth option.")
    parser.add_argument("--chrom_length", dest="chrom_length", type=str, default=None,
                        help="chrom_length file.")
    args = parser.parse_args()
    merge_transcripts(args.gtf, args.bed, args.pickle,
                      args.escape_gene_types, args.escape_gene_name_startswith,
                      args.chrom_length, args.remove_strand_overlap, args.species)


if __name__ == "__main__":
    run()
