"""
Joins two paired-end reads on the overlapping ends. need fastq-join 1.3.1.
"""
#!/usr/bin/env python3
import os
import argparse
import subprocess
from typing import Literal
from pybiotk.io import FastqPair, OpenFqGzip
from pybiotk.utils import logging, reverse_seq
from pybiotk.utils.reverse_fastx import reverse_fastx


def fastq_join(fq1: str, fq2: str, outprefix: str, threads: int = 1, 
               save_as: Literal["read1", "read2"] = "read1", collapse: bool = False):
    fq1_uncom = fq1 + ".uncom.fq"
    fq2_uncom = fq2 + ".uncom.fq"
    logging.info("start decompressing fastq ...")
    try:
        subprocess.check_call(f"pigz -p {threads} -d -c {fq1} > {fq1_uncom}", shell=True)
        subprocess.check_call(f"pigz -p {threads} -d -c {fq2} > {fq2_uncom}", shell=True)
    except subprocess.CalledProcessError:
        logging.warning(f"pigz is not installed, use gzip instead.")
        subprocess.check_call(f"gzip -d -c {fq1} > {fq1_uncom}", shell=True)
        subprocess.check_call(f"gzip -d -c {fq2} > {fq2_uncom}", shell=True)
    
    un1_uncom = outprefix + ".unmerge_R1.fq"
    un2_uncom = outprefix + ".unmerge_R2.fq"
    join_uncom = outprefix + ".merge.fq"
    
    logging.info("start join fastq ...")
    try:
        subprocess.check_call(f"fastq-join {fq1_uncom} {fq2_uncom} -o {un1_uncom} -o {un2_uncom} -o {join_uncom}")
    except subprocess.CalledProcessError:
        logging.error(f"fastq-join is not installed.")
        raise
    if os.path.exists(join_uncom):
        logging.info("join completed, remove decompressed fastq")
        os.remove(fq1_uncom)
        os.remove(fq2_uncom)
    
    logging.info("start compressing fastq ...")
    try:
        subprocess.check_call(f"pigz -p {threads} {un1_uncom}", shell=True)
        subprocess.check_call(f"pigz -p {threads} {un2_uncom}", shell=True)
        subprocess.check_call(f"pigz -p {threads} {join_uncom}", shell=True)
    except subprocess.CalledProcessError:
        logging.warning(f"pigz is not installed, use gzip instead.")
        subprocess.check_call(f"gzip {un1_uncom}", shell=True)
        subprocess.check_call(f"gzip {un2_uncom}", shell=True)
        subprocess.check_call(f"gzip {join_uncom}", shell=True)
    
    un1 = un1_uncom + ".gz"
    un2 = un2_uncom + ".gz"
    join = join_uncom + ".gz"
    
    if save_as == "read2":
        logging.info(f"start to reverse sequence ...")
        rev_join = join_uncom + ".rev.gz"
        reverse_fastx(join, rev_join)
        os.remove(join)
        os.rename(rev_join, join)
    
    if collapse:
        un_collapse = outprefix + ".unmerge.fq.gz"
        with OpenFqGzip(un_collapse) as fq, FastqPair(un1, un2) as fp:
            for entry1, entry2 in fp:
                name1 = entry1.name + "_R1"
                sequence1 = entry1.sequence
                comment1 = entry1.comment
                quality1 = entry1.quality
                name2 = entry2.name + "_R2"
                sequence2 = entry2.sequence
                comment2 = entry2.comment
                quality2 = entry2.quality
                
                if save_as == "read2":
                    sequence1 = reverse_seq(sequence1)
                    quality1 = "".join(reversed(quality1))
                else:
                    sequence2 = reverse_seq(sequence2)
                    quality2 = "".join(reversed(quality2))
                fq.write_entry(name1, sequence1, comment1, quality1)
                fq.write_entry(name2, sequence2, comment2, quality2)
        merge_collapse = outprefix + ".collapse.fq.gz"
        subprocess.check_call(f"cat {join} {un_collapse} > {merge_collapse}", shell=True)
        if os.path.exists(merge_collapse):
            logging.info(f"collapse completed, remove temp fastq files")
            os.remove(un1)
            os.remove(un2)
            os.remove(join)
            os.remove(un_collapse)


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-1", dest="fq1", type=str, required=True, help="R1 fastq.")
    parser.add_argument("-2", dest="fq2", type=str, required=True, help="R2 fastq.")
    parser.add_argument("-o", dest="outprefix", type=str, default="fastq", help="output file prefix.")
    parser.add_argument("-p", dest="threads", type=int, default=1, help="use pgzip")
    parser.add_argument("--save_as", dest="save_as", default="read1", choices=["read1", "read2"], help="save as read1 or read2.")
    parser.add_argument("--collapse", dest="collapse", action="store_true", help="collapse merged and unmerged fastq.")
    
    args = parser.parse_args()
    fastq_join(args.fq1, args.fq2, args.outputprefix, args.threads, args.save_as, args.collapse)


if __name__ == "__main__":
    run()
