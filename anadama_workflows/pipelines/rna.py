from os.path import basename, abspath

from anadama import util
from anadama.pipelines import Pipeline

from .. import subread
from .. import samtools

from . import SampleMetadataMixin
from . import (
    maybe_convert_to_fastq,
    infer_pairs,
    _to_merged
)

class RNAPipeline(Pipeline, SampleMetadataMixin):
    """
    Pipeline for analyzing RNA-seq. Produces read count tables for all
    samples.

    Steps:

    * For each sequence set:

      - Convert sequences into paired and single fastq files
      - Align sequences to a genome (default GRCh38/hg38)
      - Combine with annotations and calculate read-counts


    Workflows used:

      - :py:func:`anadama_workflows.general.sequence_convert`
      - :py:func:`anadama_workflows.general.pe_split`
      - :py:func:`anadama_workflows.subread.align`
      - :py:func:`anadama_workflows.subread.featureCounts`

    """

    name = "RNA"

    products = {
        "sample_metadata"     : list(),
        "raw_seq_files"       : list(),
        "paired_fastq_files"  : list(),
        "unpaired_fastq_files": list(),
        "align_sams"          : list(),
        "count_tables"        : list()
    }

    default_options = {
        'infer_pairs':         {
            'infer': True
        },
        'sequence_convert': { },
        'to_paired_fastq' : { },
        'subread_align'   : { },
        'featureCounts'   : { }
    }
    
    def __init__(self, sample_metadata, 
                 raw_seq_files=list(),
                 paired_fastq_files=list(), 
                 unpaired_fastq_files=list(), 
                 align_sams=list(),
                 products_dir=str(), 
                 workflow_options=dict(), 
                 *args, **kwargs):
        """Initialize the pipeline. """

        super(RNAPipeline, self).__init__(*args, **kwargs)

        if not products_dir:
            products_dir = settings.workflows.product_directory
        self.products_dir = abspath(products_dir)

        self.options = self.default_options.copy()
        self.options.update(workflow_options)

        self.add_products(
            sample_metadata      = sample_metadata,
            raw_seq_files        = raw_seq_files,
            paired_fastq_files   = paired_fastq_files,
            unpaired_fastq_files = unpaired_fastq_files,
            align_sams           = align_sams,
            count_tables         = list()
        )
        
        def _default_metadata():
            cls = namedtuple("Sample", ['SampleID'])
            return [ cls(basename(f)) for f in raw_seq_files ]
        self._unpack_metadata(default = _default_metadata)


    def _configure(self):
        if self.options['infer_pairs'].get('infer'):
            paired, notpaired = infer_pairs(self.raw_seq_files)
            self.raw_seq_files = paired + notpaired

        maybe_tasks = list()
        for maybe_pair in self.raw_seq_files:
            is_pair = type(maybe_pair) in (tuple, list)
            if is_pair:
                pair, tasks = maybe_convert_to_fastq(maybe_pair,
                                                     self.products_dir)
                self.paired_fastq_files.append(pair)
                maybe_tasks.extend(tasks)
            elif util.guess_seq_filetype(maybe_pair) == 'bam':
                prefix = util.new_file( util.rmext(basename(maybe_pair)),
                                        basedir=self.products_dir )
                t = samtools.to_paired_fastq(maybe_pair, prefix)
                paired, single = t['targets'][:2], t['targets'][2]
                self.paired_fastq_files.append(paired)
                self.unpaired_fastq_files.append(single)
                maybe_tasks.append(t)
            else:
                single, tasks = maybe_convert_to_fastq([maybe_pair],
                                                       self.products_dir)
                self.unpaired_fastq_files.append(single[0])
                maybe_tasks.extend(tasks)

        for task in maybe_tasks:
            yield task
                
        for pair in self.paired_fastq_files:
            align_sam = util.new_file( 
                _to_merged(basename(pair[0]), tag="align"),
                basedir=self.products_dir 
            )
            align_sam += ".sam"
            self.align_sams.append(align_sam)
            yield subread.align(
                pair, align_sam,
                self.options.get('subread_align', dict())
            )

        for single in self.unpaired_fastq_files:
            align_sam = util.new_file( util.addtag(basename(single), "align"),
                                       basedir=self.products_dir )
            align_sam += ".sam"
            self.align_sams.append(align_sam)
            yield subread.align(
                single, align_sam,
                self.options.get('subread_align', dict())
            )
        
        for align_sam in self.align_sams:
            count_table = util.new_file( 
                util.addtag(basename(align_sam), "count"),
                basedir=self.products_dir
            )
            self.count_tables.append(count_table)
            yield subread.featureCounts(
                [align_sam], count_table,
                self.options.get('featureCounts', dict())
            )



        
