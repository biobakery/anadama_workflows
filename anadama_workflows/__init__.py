from .wgs import metaphlan2, humann
from .alignment import bowtie2_align
from .general import (
    extract, 
    fastq_split,
    sequence_convert
)
from .sixteen import (
    write_map, 
    demultiplex, 
    pick_otus_closed_ref, 
    merge_otu_tables
)
from .association import (
    maaslin,
    biom2tsv
)

all = (
    metaphlan2, humann,
    bowtie2_align, 
    extract, fastq_split, sequence_convert,
    write_map, demultiplex, pick_otus_closed_ref, merge_otu_tables,
    maaslin, biom2tsv
)


