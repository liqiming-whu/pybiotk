# pybiotk: A python toolkit for bioinformatics analysis

## Install

### use PyPi(official)

```
pip install pybiotk
```

**An older version may be installed**

### or

```
git clone https://gitee.com/liqiming_whu/pybiotk.git
cd pybiotk
pip install .
```

## Modules

```
console_scripts =
    gtf2bed = pybiotk.convert.gtf2bed:run
    bed2bedgraph = pybiotk.convert.bed2bedgraph:run
    fq2fasta = pybiotk.convert.fq2fasta:run
    fa2fastq = pybiotk.convert.fa2fastq:run
    bam2fastx = pybiotk.convert.bam2fastx:run
    bampe_order_by_name = pybiotk.convert.bampe_order_by_name:run
    bam_random = pybiotk.utils.bam_random:run
    gtf_filter = pybiotk.utils.gtf_filter:run
    fasta_filter = pybiotk.utils.fasta_filter:run
    fastq_uniq = pybiotk.utils.fastq_uniq:run
    fastq_join = pybiotk.utils.fastq_join:run
    fastx_rename = pybiotk.utils.fastx_rename:run
    genomefetcher = pybiotk.utils.genomefetcher:run
    bigwigfetcher = pybiotk.utils.bigwigfetcher:run
    reverse_fastx = pybiotk.utils.reverse_fastx:run
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
    summary_log = pybiotk.utils.summary_log:run
    ercc_parser = pybiotk.utils.ercc_parser:run
    infer_experiment = pybiotk.utils.infer_experiment:run
```

## Usage

```
usage: pyanno [-h] -i INPUT [-o OUTPUT] -g GTF [-l {transcript,gene}] [--tss_region TSS_REGION [TSS_REGION ...]] [--tss_region_name TSS_REGION_NAME] [--downstream DOWNSTREAM] [--downstream_name DOWNSTREAM_NAME] [--tss] [--tes] [--start_condon] [--stop_condon] [-s]
              [--rule {1+-,1-+,2++,2--,1++,1--,2+-,2-+,+-,-+,++,--}] [-p] [--ordered_by_name]

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input file, bam or bed. The file type will be inferred from the filename suffix ['*.bam', '*.bed']. (default: None)
  -o OUTPUT, --output OUTPUT
                        output file name. [stdout] (default: -)
  -g GTF, --gtf GTF     gtf file download from Genecode, or a sorted gtf file. (default: None)
  -l {transcript,gene}, --level {transcript,gene}
                        annotation level, transcript or gene. (default: transcript)
  --tss_region TSS_REGION [TSS_REGION ...]
                        choose region from tss. (default: [-3000, 0])
  --tss_region_name TSS_REGION_NAME
                        tss region name. (default: Upstream)
  --downstream DOWNSTREAM
                        downstream length from tes. (default: 3000)
  --downstream_name DOWNSTREAM_NAME
                        downstream name. (default: Downstream)
  --tss                 annotate tss. (default: False)
  --tes                 annotate tes. (default: False)
  --start_condon        annotate start condon. (default: False)
  --stop_condon         annotate stop condon. (default: False)
  -s, --strand          require same strandedness. (default: False)
  --rule {1+-,1-+,2++,2--,1++,1--,2+-,2-+,+-,-+,++,--}
                        how read(s) were stranded during sequencing. only for bam. (default: 1+-,1-+,2++,2--)
  -p, --pair            annotate fragments instead of reads. (default: False)
  --ordered_by_name     if input bam is ordered by name, only for pair-end bam. (default: False)

```


