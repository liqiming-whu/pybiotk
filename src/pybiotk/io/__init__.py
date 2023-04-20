from .bed import *
from .gtf import *
try:
    from .bw import *
except ModuleNotFoundError:
    import warnings
    warnings.warn("pyBigWig may not be installed correctly, the bw module cannot be used.")
    pass
try:
    from .bam import *
    from .fasta import *
    from .fastq import *
except ModuleNotFoundError:
    import warnings
    warnings.warn("pysam may not be installed correctly, the bam,fastq,fasta module cannot be used.")
    pass
