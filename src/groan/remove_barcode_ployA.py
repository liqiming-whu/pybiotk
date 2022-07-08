import os
import re
import argparse
import pysam
from collections import defaultdict
from pybiotk.io import OpenFqGzip
from pybiotk.utils import reverse_seq

def parse_args():
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    base_group = parser.add_argument_group("Base")
    base_group.add_argument("--R1", type=str, dest="R1", metavar="R1.fastq.gz", required=True, help="input R1")
    base_group.add_argument("--R2", type=str, dest="R2", metavar="R2.fastq.gz", required=True, help="input R2")
    base_group.add_argument("--out-R1", type=str, dest="out_R1", metavar="out_R1.fastq.gz", required=True, help="output R1")
    base_group.add_argument("--out-R2", type=str, dest="out_R2", metavar="out_R2.fastq.gz", required=True, help="output R2")
    base_group.add_argument("--out-stat", type=str, dest="out_stat", metavar="stat.tsv", default=os.devnull, help="stats file")
    base_group.add_argument("--min-A-num", type=int, dest="min_A_num", metavar="min_A_num", required=False, default=10, help="min_A_num")
    base_group.add_argument("--barcode-len", type=int, dest="barcode_len", metavar="barcode_len", required=False, default=8, help="barcode length")
    base_group.add_argument("--min-seq-len", type=int, dest="min_seq_len", metavar="min_seq_len", required=False, default=15, help="min seq length")
    return parser.parse_args()

def main():
    args = parse_args()
    f_R1 = args.R1
    f_R2 = args.R2
    f_out_R1 = OpenFqGzip(args.out_R1)
    f_out_R2 = OpenFqGzip(args.out_R2, "w")
    f_out_stat = args.out_stat
    min_A_num = args.min_A_num
    barcode_len = args.barcode_len
    min_seq_len = args.min_seq_len

    half_A_num = int(min_A_num/2)
    half_T_pattern = f"T{{{half_A_num},}}"
    A_pattern = f"A{{{min_A_num},}}"

    R1_fq = pysam.FastxFile(f_R1)
    R2_fq = pysam.FastxFile(f_R2)
    polyA_cnt = defaultdict(int)  # max polyA len, polyA num pass, polyT num pass, Barcode pass, Seq len pass

    for R1_rec, R2_rec in zip(R1_fq, R2_fq):
        assert R1_rec.name == R2_rec.name
        name = R1_rec.name
        R1_seq = R1_rec.sequence.upper()
        R2_seq = R2_rec.sequence.upper()
        R1_quality = R1_rec.quality
        R2_quality = R2_rec.quality
        polyT_indx = sorted(re.finditer(half_T_pattern, R2_seq), key=lambda x: x.start() - x.end())
        if not polyT_indx:
            stat = ("NA", "NA", "NA", "NA", "NA")
            polyA_cnt[stat] += 1
            continue
        polyT_start = polyT_indx[0].start()
        polyT_end = polyT_indx[0].end()
        polyT_len = polyT_end - polyT_start
        polyT_num = len(list(filter(lambda x: (x.end() - x.start()) >= min_A_num, polyT_indx)))
        polyA_indx = list(re.finditer(A_pattern, R1_seq))
        polyA_num = len(polyA_indx)
        if (polyT_num!=1) or (polyA_num>1):
            stat = (polyT_len, polyA_num<=1, polyT_num==1, "NA", "NA")
            polyA_cnt[stat] += 1
            continue
        if (R2_seq[barcode_len:polyT_end].count("T") + 1) < (polyT_end-barcode_len):
            stat = (polyT_len, True, False, "NA", "NA")
            polyA_cnt[stat] += 1
            continue

        T_barcode = R2_seq[:barcode_len]
        barcode_seq = reverse_seq(T_barcode)
        if (len(barcode_seq) < barcode_len) or barcode_seq.count("A") or barcode_seq.count("N"):
            stat = (polyT_len, True, True, False, "NA")
            polyA_cnt[stat] += 1
            continue

        ## trim R2
        R2_new_seq = R2_seq[polyT_end:]
        R2_new_quality = R2_quality[polyT_end:]
        assert len(R2_new_seq) == len(R2_new_quality), name
        new_name = f"{name}_{barcode_seq}"

        ## trim R1
        if polyA_indx:
            R1_new_seq = R1_seq[:polyA_indx[0].start()]
            R1_new_quality = R1_quality[:polyA_indx[0].start()]
            assert len(R1_new_seq) == len(R1_new_quality), name
        else:
            ## no polyA
            R1_new_seq = R1_seq.rstrip("A")
            R1_new_quality = R1_quality[:len(R1_new_seq)]
            assert len(R1_new_seq) == len(R1_new_quality), name

        if (len(R1_new_seq) < min_seq_len) or (len(R2_new_seq) < min_seq_len):
            stat = (polyT_len, True, True, True, False)
            polyA_cnt[stat] += 1
            continue

        stat = (polyT_len, True, True, True, True)
        polyA_cnt[stat] += 1

        f_out_R1.write_entry(name=new_name, sequence=R1_new_seq, quality=R1_new_quality)
        f_out_R2.write_entry(name=new_name, sequence=R2_new_seq, quality=R2_new_quality)

    with open(f_out_stat, "w") as f:
        header = "polyALen\tR1PaNumPass\tR2PaNumPass\tBarcodePass\tInsertLenPass\tReadNum\n"
        f.write(header)
        for key, val in sorted(polyA_cnt.items(), key=lambda x: -x[1]):
            data = list(key) + [val]
            f.write("\t".join(list(map(str, data)))+"\n")

    f_out_R1.close()
    f_out_R2.close()


def run():
    main()


if __name__ == "__main__":
    run()
