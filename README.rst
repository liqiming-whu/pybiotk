.. These are examples of badges you might want to add to your README:
   please update the URLs accordingly

    .. image:: https://api.cirrus-ci.com/github/<USER>/pybiotk.svg?branch=main
        :alt: Built Status
        :target: https://cirrus-ci.com/github/<USER>/pybiotk
    .. image:: https://readthedocs.org/projects/pybiotk/badge/?version=latest
        :alt: ReadTheDocs
        :target: https://pybiotk.readthedocs.io/en/stable/
    .. image:: https://img.shields.io/coveralls/github/<USER>/pybiotk/main.svg
        :alt: Coveralls
        :target: https://coveralls.io/r/<USER>/pybiotk
    .. image:: https://img.shields.io/pypi/v/pybiotk.svg
        :alt: PyPI-Server
        :target: https://pypi.org/project/pybiotk/
    .. image:: https://img.shields.io/conda/vn/conda-forge/pybiotk.svg
        :alt: Conda-Forge
        :target: https://anaconda.org/conda-forge/pybiotk
    .. image:: https://pepy.tech/badge/pybiotk/month
        :alt: Monthly Downloads
        :target: https://pepy.tech/project/pybiotk
    .. image:: https://img.shields.io/twitter/url/http/shields.io.svg?style=social&label=Twitter
        :alt: Twitter
        :target: https://twitter.com/pybiotk

.. image:: https://img.shields.io/badge/-PyScaffold-005CA0?logo=pyscaffold
    :alt: Project generated with PyScaffold
    :target: https://pyscaffold.org/

|

=======
pybiotk
=======


pybiotk: A python toolkit for bioinformatics analysis


Install
====

python setup.py install

Modules
====

gtf2bed = pybiotk.convert.gtf2bed:run

bed2bedgraph = pybiotk.convert.bed2bedgraph:run

fq2fasta = pybiotk.convert.fq2fasta:run

bam2fastx = pybiotk.convert.bam2fastx:run

bampe_order_by_name = pybiotk.convert.bampe_order_by_name:run

bam_random = pybiotk.utils.bam_random:run

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


Usage
====


usage: pyanno [-h] -i INPUT -o OUTPUT -g GTF [-l {transcript,gene}] [--tss_region TSS_REGION [TSS_REGION ...]] [--downstream DOWNSTREAM] [-s]
              [--rule {1+-,1-+,2++,2--,1++,1--,2+-,2-+,+-,-+,++,--}] [-p] [--ordered_by_name]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        input file, bam or bed. The file type will be inferred from the filename suffix ['*.bam', '*.bed']. (default: None)
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


.. _pyscaffold-notes:

Note
====

This project has been set up using PyScaffold 4.2.3. For details and usage
information on PyScaffold see https://pyscaffold.org/.
