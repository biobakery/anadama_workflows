import os
from os.path import basename

from anadama import util
from anadama.pipelines import Pipeline

from .. import settings
from .. import general, wgs, alignment

from . import SampleFilterMixin


class WGSPipeline(Pipeline, SampleFilterMixin):

    """Pipeline for analyzing whole metagenome shotgun sequence data.
    Produces taxonomic profiles with metaphlan2 and gene, pathway
    lists with HUMAnN

    Steps:

    * For each sequence set:

      - Convert the sequences and concatenate into a single fastq file
      - Filter sequences for length > 60 bases
      - Perform taxonomic profiling with metaphlan2
      - Align the fastq file agains the KEGG proks reduced dataset with
          bowtie2
      - Infer pathway, gene lists with HUMAnN


    Workflows used:

    * anadama_workflows.general.sequence_convert
    * anadama_workflows.wgs.metaphlan2
    * anadama_workflows.alignment.bowtie2
    * anadama_workflows.wgs.humann

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

    def __init__(self, sample_metadata,
                 raw_seq_files=list(), 
                 intermediate_fastq_files=list(),
                 alignment_result_files=list(),
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
                                of lists, if you want to aggregate results by
                                something other than by file
        :keyword intermediate_fastq_files: List of strings; List of files to 
                                           be fed into metaphlan for taxonomic
                                           profiling and bowtie2 for alignment
        :keyword alignment_result_files: List of strings; List of files to be
                                         fed into HUMAnN for gene, pathway 
                                         inference.
        :keyword products_dir: String; Directory path for where outputs will 
                               be saved.
        :keyword workflow_options: Dictionary; **opts to be fed into the 
                                  respective workflow functions.

        """
        super(WGSPipeline, self).__init__(*args, **kwargs)

        if not products_dir:
            products_dir = settings.workflows.product_directory
        self.products_dir = os.path.realpath(products_dir)

        self.options = {
            'sequence_convert': {
                'lenfilters_list': [ '>=60' ]
            }, 
            'metaphlan2':       { },
            'bowtie2_align':    { },
            'humann':           { }
        }
        self.options.update(workflow_options)

        self.add_products(
            sample_metadata          = sample_metadata,
            raw_seq_files            = raw_seq_files,
            intermediate_fastq_files = intermediate_fastq_files,
            alignment_result_files   = alignment_result_files,
            metaphlan_results        = list(),
            otu_tables               = list()
        )


    def _configure(self):
        # Convert all raw files into fastq files; run them through
        # metaphlan2 to get community profiles
        for files in self.raw_seq_files:
            if type(files) is str:
                files = [files]
            fastq_file = util.new_file( basename(files[0])+"_merged.fastq",
                                        basedir=self.products_dir )
            yield general.sequence_convert(
                files, fastq_file, 
                **self.options.get('sequence_convert', dict())
            )
            self.intermediate_fastq_files.append(fastq_file)

            metaphlan_file = util.new_file( 
                basename(files[0])+"_merged.metaphlan2.pcl",
                basedir=self.products_dir )
            otu_table = metaphlan_file.replace('.pcl', '.biom')
            yield wgs.metaphlan2(
                [fastq_file], output_file=metaphlan_file,
                biom=otu_table,
                # first index is for first item in list of samples
                # second index is to get the sample id from the sample
                sample_id=self._filter_samples_for_file(self.sample_metadata,
                                                        fastq_file)[0][0],
                **self.options.get('metaphlan2', dict())
            )
            self.metaphlan_results.append(metaphlan_file)
            self.otu_tables.append(otu_table)

        # Next align the fastq files against the kegg proks reduced db
        for fastq_file in self.intermediate_fastq_files:
            alignment_file = fastq_file+".sam"
            self.alignment_result_files.append(alignment_file)
            yield alignment.bowtie2_align(
                [fastq_file], alignment_file,
                **self.options.get('bowtie2_align', dict())
            )

        # Finally, HUMAnN all alignment files
        for alignment_file in self.alignment_result_files:
            # create a directory for the humann environment
            # util.new_file does that for us
            new_basedir = alignment_file+"_humann"
            scons_fname = util.new_file('SConstruct', basedir=new_basedir)
            yield wgs.humann(
                [alignment_file], workdir=new_basedir,
                **self.options.get('humann', dict())
            )

