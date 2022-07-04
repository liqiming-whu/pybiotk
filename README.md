# pybiotk: A python toolkit for bioinformatics analysis

## Install

    python setup.py install

## Modules

```
console_scripts =
    gtf2bed = pybiotk.convert.gtf2bed:run
    bed2bedgraph = pybiotk.convert.bed2bedgraph:run
    fq2fasta = pybiotk.convert.fq2fasta:run
    bam2fastx = pybiotk.convert.bam2fastx:run
    bampe_order_by_name = pybiotk.convert.bampe_order_by_name:run
    gtf_filter = pybiotk.utils.gtf_filter:run
    fasta_filter = pybiotk.utils.fasta_filter:run
    fastq_uniq = pybiotk.utils.fastq_uniq:run
    seq_random = pybiotk.utils.seq_random:run
    merge_row = pybiotk.utils.merge_row:run
    read_tables = pybiotk.utils.read_tables:run
    rmats_filter = pybiotk.utils.rmats_filter:run
    count_normalize = pybiotk.utils.normalize:run
    reference_count = pybiotk.utils.reference_count:run
    pyanno = pybiotk.utils.pyanno:run
    rna_fragment_size = pybiotk.utils.fragment_size:run
    merge_subseq = pybiotk.utils.merge_subseq:run
    subseq_analysis = pybiotk.utils.subseq_analysis:run
    merge_transcript = groan.merge_transcript:run
    filter_gene = groan.filter_gene:run
    expression_qc = groan.expression_qc:run
    signal_continuity_qc = groan.signal_continuity_qc:run
    readthrough_stats = groan.read_through_stats:run
    readthrough_peak = groan.read_through_peak:run
    GroAn = groan.GroAn:main
    metaplot = groan.metaplot:run
    summary_log = pybiotk.utils.summary_log:run
    peak_qc = groan.peak_qc:run
    remove_no_chimeric = ricpy.chimeric.NoChimAlign:run
```

## Usage

```
usage: GroAn [-h] -f FWD [FWD ...] [-r REV [REV ...]] -g GTF [-o OUTDIR] [-m {coverage,reads}] [-v VALUE] [-l {body,tss,tes}] [-d DISTANCE]
             [--filter_gene_by {value,top,mean,quantile}] [--filter_signal_continuity {value,top,mean,quantile}] [--filter_signal_value FILTER_SIGNAL_VALUE]
             [--group GROUP [GROUP ...]] [--gene_length GENE_LENGTH] [--downStream DOWNSTREAM] [--downStart DOWNSTART]

GroAn: a tool for GroSeq Analysis.

optional arguments:
  -h, --help            show this help message and exit
  -f FWD [FWD ...]      forward bw files. (default: None)
  -r REV [REV ...]      reverse bw files. (default: None)
  -g GTF                gtf file. (default: None)
  -o OUTDIR             outdir. (default: results)
  -m {coverage,reads}   statics method. (default: reads)
  -v VALUE              minimum coverage or reads value. (default: 0.25)
  -l {body,tss,tes}, --loci {body,tss,tes}
                        calc in which loci (default: tes)
  -d DISTANCE, --distance DISTANCE
                        tss/tes +- distance. (default: 3000)
  --filter_gene_by {value,top,mean,quantile}
                        filter method (default: mean)
  --filter_signal_continuity {value,top,mean,quantile}
                        filter signal method (default: None)
  --filter_signal_value FILTER_SIGNAL_VALUE
                        filter signal value (default: 0.75)
  --group GROUP [GROUP ...]
                        samples group. (default: None)
  --gene_length GENE_LENGTH
                        gene_length cutoff. (default: 5000)
  --downStream DOWNSTREAM
                        downStream cutoff. (default: 10000)
  --downStart DOWNSTART
                        downStream distance from TES. (default: 0)
```

```
usage: pyanno [-h] -i INPUT -o OUTPUT -g GTF [-l {transcript,gene}] [--tss_region TSS_REGION [TSS_REGION ...]] [--downstream DOWNSTREAM] [-s]
              [--rule {1+-,1-+,2++,2--,1++,1--,2+-,2-+,+-,-+,++,--}] [-p] [--ordered_by_name]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input file, bam or bed. The file type will bed inferred from the filename suffix ['*.bam', '*.bed']. (default: None)
  -o OUTPUT, --output OUTPUT
                        output file name. (default: None)
  -g GTF, --gtf GTF     gtf file download from Genecode, or a sorted gtf file. (default: None)
  -l {transcript,gene}, --level {transcript,gene}
                        annotation level, transcript or gene. (default: transcript)
  --tss_region TSS_REGION [TSS_REGION ...]
                        choose region from tss. (default: [-1000, 1000])
  --downstream DOWNSTREAM
                        downstream length from tes. (default: 3000)
  -s, --strand          require same strandedness. (default: False)
  --rule {1+-,1-+,2++,2--,1++,1--,2+-,2-+,+-,-+,++,--}
                        how read(s) were stranded during sequencing. only for bam. (default: 1+-,1-+,2++,2--)
  -p, --pair            annotate fragments instead of reads. (default: False)
  --ordered_by_name     if input bam is ordered by name, only for pair-end bam. (default: False)
```
