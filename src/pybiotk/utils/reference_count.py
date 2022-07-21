import argparse
from collections import defaultdict

from pybiotk.io import Bam, BamPE, BamType, check_bam_type


def main(filename: str, output: str, count_fragments: bool = False, ordered_by_name: bool = False):
    reference_count_dict = defaultdict(int)
    bamtype = check_bam_type(filename)
    if bamtype is BamType.PE and count_fragments:
        with BamPE(filename, ordered_by_name=ordered_by_name) as bam:
            for read1, read2 in bam.iter_pair(properly_paired=False):
                if read1 is not None and read2 is not None:
                    if read1.reference_name == read2.reference_name:
                        reference_count_dict[read1.reference_name] += 1
                    elif not read1.reference_name == read2.reference_name:
                        reference_count_dict[read1.reference_name] += 1
                        reference_count_dict[read2.reference_name] += 1
                if read1 is not None and read2 is None:
                    reference_count_dict[read1.reference_name] += 1
                if read2 is not None and read1 is None:
                    reference_count_dict[read2.reference_name] += 1
    else:
        with Bam(filename) as bam:
            for read in bam.iter_mapped():
                reference_count_dict[read.reference_name] += 1

    with open(output, 'w') as outfile:
        outfile.write(f"Geneid\t{filename}\n")
        for reference in reference_count_dict:
            outfile.write(f"{reference}\t{reference_count_dict[reference]}\n")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("input", type=str, help="input bam.")
    parser.add_argument("-o", dest="output", type=str, required=True,
                        help="output file.")
    parser.add_argument("-p", "--pair", dest="pair", action="store_true",
                        help="count fragments instead of reads.")
    parser.add_argument("--ordered_by_name", dest="ordered_by_name", action="store_true",
                        help="if input bam is ordered by name.")

    args = parser.parse_args()
    main(args.input, args.output, args.pair, args.ordered_by_name)
