#!/usr/bin/env python3
import os
import time
import argparse
import pysam
import pandas as pd
from typing import NamedTuple, Literal, Tuple
from stream import cmdin, mkdir
from pybiotk.utils import write_table, logging


class RNAhybridDotBracket(NamedTuple):
    query_structure: str
    target_structure: str
    mfe: str
    pvalue: str
    postion: Literal['left', 'right'] = 'left'
    
    def __str__(self):
        structure = f"{self.query_structure}&{self.target_structure}" if self.postion == "left" else f"{self.target_structure}&{self.query_structure}"
        return structure


def parse_pair(q_paired: str, q_unpaired: str, t_paired: str, t_unpaired: str) -> Tuple[str, ...]:
    q_p = ''
    q_u = ''
    t_p = ''
    t_u = ''
    i = 0
    while i < len(q_paired):
        if q_paired[i] != " " or q_unpaired[i] != " ":
            q_p += q_paired[i]
            q_u += q_unpaired[i]
        i += 1

    i = 0
    while i < len(t_paired):
        if t_paired[i] != " " or t_unpaired[i] != " ":
            t_p += t_paired[i]
            t_u += t_unpaired[i]
        i += 1

    return "".join(reversed(q_p)), "".join(reversed(q_u)), t_p, t_u


def rnahybrid2dotbracket(res: str, q_len: int, t_len: int, postion: Literal['left', 'right']) -> RNAhybridDotBracket:
    q = ''
    t = ''
    fields = res.split(":")
    assert int(fields[1]) == t_len
    assert int(fields[3]) == q_len
    pos = int(fields[6])
    mfe = fields[4]
    p_value = fields[5]
    q_paired, q_unpaired, t_paired, t_unpaired = parse_pair(fields[9], fields[10], fields[8], fields[7])
    assert len(q_paired) == len(q_unpaired), res
    assert len(t_paired) == len(t_unpaired), res
    if not q_paired:
        return q_len*".", t_len*".", mfe, p_value
    assert len(q_paired) == q_len, res
    for i in range(q_len):
        if q_paired[i] == " ":
            q += "."
        elif postion == 'left':
            q += "("
        elif postion == 'right':
            q += ")"
        else:
            q += "("
    for i in range(t_len):
        if i < pos - 1 or i >= pos -1 + len(t_paired):
            t += '.'
            continue
        if t_paired[i-pos+1] == " ":
            t += '.'
        elif postion == 'left':
            t += ')'
        elif postion == 'right':
            t += '('
        else:
            t += ')'

    return RNAhybridDotBracket(q, t, mfe, p_value, postion)


def main(query_fasta: str, target_fasta: str, annofile: str, chimeric_out1:str, chimeric_out2: str, tempdir: str):
    logging.info("start to predict rna structure...")
    start = time.perf_counter()
    mkdir(tempdir)
    query = pysam.FastaFile(query_fasta)
    target = pysam.FastxFile(target_fasta)
    anno = pd.read_table(annofile, header=0, index_col=0)
    chimeric = open(chimeric_out1, "w")
    
    df_list = []
    for entry in target:
        name = entry.name
        target_sequence = entry.sequence
        strand = anno.loc[name]["strand"]
        gene_name = anno.loc[name]["geneName"]
        namefileds = name.split("|")
        read_name = namefileds[0]
        ref_name = namefileds[1]
        ref_loc = namefileds[2]
        if ref_loc == 'mid':
            query_sequence = query.fetch(name.rstrip("|L").rstrip("|R"))
            LpCp = namefileds[3]
            RpCp = namefileds[4]
            if len(namefileds) >= 5:
                target_loc = namefileds[5]
            else:
                target_loc = None
        else:
            query_sequence = query.fetch(name)
            pCp = namefileds[3]
        
        q_fa = os.path.join(tempdir, "query.fa")
        t_fa = os.path.join(tempdir, "target.fa")
        with open(q_fa, "w") as q, open(t_fa, "w") as t:
            q.write(f">{ref_name}\n{query_sequence}\n")
            t.write(f">{read_name}\n{target_sequence}\n")
        cmd = f'RNAhybrid -b 1 -c -s 3utr_human -n 36 -t {t_fa} -q {q_fa}'
        res = list(cmdin(cmd))
        os.remove(q_fa)
        os.remove(t_fa)
        loc = "right" if target_loc == 'L' else "left"
        if ref_loc == 'mid':
            dotbracket = rnahybrid2dotbracket(res[0], len(query_sequence), len(target_sequence), loc)
        else:
            dotbracket = rnahybrid2dotbracket(res[0], len(query_sequence), len(target_sequence), ref_loc)
        seqname = name + "|" + strand
        query_structure = dotbracket.query_structure
        target_structure = dotbracket.target_structure
        mfe = dotbracket.mfe
        pvalue = dotbracket.pvalue
        if ref_loc == 'mid':
            if not target_loc:
                if RpCp == 'withoutRpCp':
                    seq = query_sequence + target_sequence
                    description = len(query_sequence)*'P' + len(target_sequence)*'M'
                    description_simple = f"{len(query_sequence)}P{len(target_sequence)}M"
                    structure = query_structure + target_structure
                elif RpCp == 'withRpCp':
                    seq = query_sequence + 'C' + target_sequence
                    description = len(query_sequence)*'P' + '.' + len(target_sequence)*'M'
                    description_simple = f"{len(query_sequence)}P.{len(target_sequence)}M"
                    structure = query_structure + '.' + target_structure
                if LpCp == 'withLpCp':
                    seq = 'C' + seq
                    description = '.' + description
                    description_simple = '.' + description_simple
                    structure = '.' + structure
            elif target_loc == 'R':
                if RpCp == 'withoutRpCp':
                    seq = query_sequence + target_sequence
                    description = len(query_sequence)*'P' + len(target_sequence)*'M'
                    description_simple = f"{len(query_sequence)}P{len(target_sequence)}M"
                    structure = query_structure + target_structure
                elif RpCp == 'withRpCp':
                    seq = query_sequence + 'C' + target_sequence
                    description = len(query_sequence)*'P' + '.' + len(target_sequence)*'M'
                    description_simple = f"{len(query_sequence)}P.{len(target_sequence)}M"
                    structure = query_structure + '.' + target_structure
                if LpCp == 'withLpCp':
                    seq = 'C' + seq
                    description = '.' + description
                    description_simple = '.' + description_simple
                    structure = '.' + structure
            elif target_loc == 'L':
                if LpCp == 'withoutLpCp':
                    seq = target_sequence + query_sequence
                    description = len(target_sequence)*'M' + len(query_sequence)*'P'
                    description_simple = f"{len(target_sequence)}M{len(query_sequence)}P"
                    structure = target_structure + query_structure
                elif LpCp == 'withLpCp':
                    seq = target_sequence + 'C' + query_sequence
                    description = len(target_sequence)*'M' + '.' + len(query_sequence)*'P'
                    description_simple = f"{len(target_sequence)}M.{len(query_sequence)}P"
                    structure = target_structure + '.' + query_structure
                if RpCp == 'withRpCp':
                    seq = seq + 'C'
                    description = description + '.'
                    description_simple = description_simple + '.'
                    structure = structure + '.'
            data = {
                "read": read_name,
                "ref_name": ref_name,
                "target_gene": gene_name,
                "strand": strand,
                "ref_location": ref_loc,
                "LpCp": LpCp,
                "RpCp": RpCp,
                "ref_sequence": query_sequence,
                "ref_len": len(query_sequence),
                "target_sequence": target_sequence,
                "target_len": len(target_sequence),
                "ref_structure": query_structure,
                "target_structure": target_structure,
                "mfe": mfe,
                "pvalue": pvalue,
                "sequence": seq,
                "sequence_length": len(seq),
                "description": description_simple,
                "structure": structure
                }
        else:
            if ref_loc == 'left':
                if pCp == 'withoutpCp':
                    seq = query_sequence + target_sequence
                    description = len(query_sequence)*'P' + len(target_sequence)*'M'
                    description_simple = f"{len(query_sequence)}P{len(target_sequence)}M"
                    structure = query_structure + target_structure
                elif pCp == 'withpCp':
                    seq = query_sequence + 'C' + target_sequence
                    description = len(query_sequence)*'P' + '.' + len(target_sequence)*'M'
                    description_simple = f"{len(query_sequence)}P.{len(target_sequence)}M"
                    structure = query_structure + '.' + target_structure
            elif ref_loc == 'right':
                if pCp == 'withoutpCp':
                    seq = target_sequence + query_sequence
                    description = len(target_sequence)*'M' + len(query_sequence)*'P'
                    description_simple = f"{len(target_sequence)}M{len(query_sequence)}P"
                    structure = target_structure + query_structure
                elif pCp == 'withpCp':
                    seq = target_sequence + 'C' + query_sequence
                    description = len(target_sequence)*'M' + '.' + len(query_sequence)*'P'
                    description_simple = f"{len(target_sequence)}M.{len(query_sequence)}P"
                    structure = target_structure + '.' + query_structure
            data = {
                "read": read_name,
                "ref_name": ref_name,
                "target_gene": gene_name,
                "strand": strand,
                "ref_location": ref_loc,
                "pCp": pCp,
                "ref_sequence": query_sequence,
                "ref_len": len(query_sequence),
                "target_sequence": target_sequence,
                "target_len": len(target_sequence),
                "ref_structure": query_structure,
                "target_structure": target_structure,
                "mfe": mfe,
                "pvalue": pvalue,
                "sequence": seq,
                "sequence_length": len(seq),
                "description": description_simple,
                "structure": structure
                }
        df_list.append(data)
        chimeric.write(f"{seqname}\n{seq}\n{description}\n{structure}\n")
    chimeric.close()
    df = pd.DataFrame(df_list)
    write_table(df, chimeric_out2)
    end = time.perf_counter()
    logging.info(f"task finished in {end-start:.2f}s.")


def run():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-q", dest="query", type=str, required=True, help="query fasta.")
    parser.add_argument("-t", dest="target", type=str, default=True, help="target fasta.")
    parser.add_argument("-a", dest="anno", type=str, default=True, help="annofile.")
    parser.add_argument("--outfile", dest="outfile", type=str, default=os.devnull, help="out 4 lines file.")
    parser.add_argument("-o", "--outtable", dest="outtable", type=str, default=os.devnull, help="output table")
    args = parser.parse_args()
    
    main(args.query, args.target, args.anno, args.outfile, args.outtable)
    

if __name__ == "__main__":
    run()
