import os
from os.path import basename
from collections import namedtuple

from anadama import util
from anadama.pipelines import Pipeline

from .. import settings
from .. import general, wgs, alignment

from . import SampleFilterMixin, SampleMetadataMixin
from . import maybe_stitch, infer_pairs


class WGSPipeline(Pipeline, SampleFilterMixin, SampleMetadataMixin):

    """Pipeline for analyzing whole metagenome shotgun sequence data.
    Produces taxonomic profiles with metaphlan2 and gene, pathway
    lists with HUMAnN

    Steps:

    * For each sequence set:

      - Convert, join paired end reads, and concatenate into a single fastq file
      - Filter sequences for length > 60 bases
      - Perform taxonomic profiling with metaphlan2
      - Align the fastq file agains the KEGG proks reduced dataset with
          bowtie2
      - Infer pathway, gene lists with HUMAnN


    Workflows used:

    * :py:func:`anadama_workflows.general.sequence_convert`
    * :py:func:`anadama_workflows.general.fastq_join`
    * :py:func:`anadama_workflows.wgs.metaphlan2`
    * :py:func:`anadama_workflows.wgs.humann2`

    """

    name = "WGS"
    products = {
        "sample_metadata"          : list(),
        "raw_seq_files"            : list(),
        "intermediate_fastq_files" : list(),
        "alignment_result_files"   : list(),
        "metaphlan_results"        : list(),
        "otu_tables"               : list(),
    }

    default_options = {
        'infer_pairs':         {
            'infer': True
        },
        'sequence_convert': { },
        'decontaminate':    { },
        'metaphlan2':       { },
        'bowtie2_align':    { },
        'humann':           { }
    }

    def __init__(self, sample_metadata,
                 raw_seq_files=list(), 
                 intermediate_fastq_files=list(),
                 decontaminated_fastq_files=list(),
                 products_dir=str(),
                 workflow_options=dict(),
                 *args, **kwargs):
        """Initialize the pipeline.
        
        :keyword sample_metadata: List of namedtuples or String; Sample-level 
                                  metadata. This can either be in a file of 
                                  tab-separated-values formatted akin to 
                                  http://qiime.org/tutorials/tutorial.html#mapping-file-tab-delimited-txt, 
                                  or in a list of python's namedtuples 
                                  (https://docs.python.org/2.7/library/collections.html#collections.namedtuple). 
                                  In either case, the sample name should be the
                                  first attribute of each sample.
        :keyword raw_seq_files: List of strings; File paths to raw WGS 
                                sequence reads. raw_seq_files can be a list
                                of pairs, if you want to stitch paired-end
                                reads.
        :keyword intermediate_fastq_files: List of strings; List of files to 
                                           be fed into metaphlan for taxonomic
                                           profiling and bowtie2 for alignment
        :keyword products_dir: String; Directory path for where outputs will 
                               be saved.
        :keyword workflow_options: Dictionary; **opts to be fed into the 
                                  respective workflow functions.

        """
        super(WGSPipeline, self).__init__(*args, **kwargs)

        if not products_dir:
            products_dir = settings.workflows.product_directory
        self.products_dir = os.path.realpath(products_dir)

        self.options = self.default_options.copy()
        self.options.update(workflow_options)

        self.add_products(
            sample_metadata            = sample_metadata,
            raw_seq_files              = raw_seq_files,
            intermediate_fastq_files   = intermediate_fastq_files,
            decontaminated_fastq_files = decontaminated_fastq_files,
            metaphlan_results          = list(),
            otu_tables                 = list()
        )

        def _default_metadata():
            cls = namedtuple("Sample", ['SampleID'])
            return [ cls(basename(f)) for f in raw_seq_files ]
        self._unpack_metadata(default = _default_metadata)


    def _configure(self):
        if self.options['infer_pairs'].get('infer'):
            paired, notpaired = infer_pairs(self.raw_seq_files)
            self.raw_seq_files = paired + notpaired

        self.raw_seq_files, _, maybe_tasks = maybe_stitch(self.raw_seq_files,
                                                          self.products_dir)
        for t in maybe_tasks:
            yield t

        for file_ in self.raw_seq_files:
            fastq_file = util.new_file( basename(file_)+"_filtered.fastq",
                                        basedir=self.products_dir )
            yield general.sequence_convert(
                [file_], fastq_file, 
                **self.options.get('sequence_convert', dict())
            )
            self.intermediate_fastq_files.append(fastq_file)

            name_base = os.path.join(self.products_dir,
                                     basename(file_))
            name_base = util.rmext(name_base)
            task_dict = wgs.knead_data([fastq_file], name_base).next()
            decontaminated_fastq = task_dict['targets'][0]
            self.decontaminated_fastq_files.append(decontaminated_fastq)
            yield task_dict

            metaphlan_file = util.new_file( 
                basename(file_)+".metaphlan2.pcl",
                basedir=self.products_dir )
            otu_table = metaphlan_file.replace('.pcl', '.biom')
            yield wgs.metaphlan2(
                [decontaminated_fastq], output_file=metaphlan_file,
                biom=otu_table,
                # first index is for first item in list of samples
                # second index is to get the sample id from the sample
                sample_id=self._filter_samples_for_file(self.sample_metadata,
                                                        fastq_file)[0][0],
                **self.options.get('metaphlan2', dict())
            )
            self.metaphlan_results.append(metaphlan_file)
            self.otu_tables.append(otu_table)

            # Finally, HUMAnN all alignment files
            humann_output_dir = fastq_file+"_humann"
            yield wgs.humann2(
                decontaminated_fastq, humann_output_dir, 
                **self.options.get('humann', dict())
            )
            
