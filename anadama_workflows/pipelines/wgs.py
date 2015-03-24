import os
from os.path import basename
from collections import namedtuple

from anadama import util
from anadama.pipelines import Pipeline

from .. import settings
from .. import general, wgs, alignment

from . import SampleFilterMixin, SampleMetadataMixin
from . import (
    infer_pairs,
    split_pairs,
    _to_merged,
)


class WGSPipeline(Pipeline, SampleFilterMixin, SampleMetadataMixin):

    """Pipeline for analyzing whole metagenome shotgun sequence data.
    Produces taxonomic profiles with metaphlan2 and gene, pathway
    lists with HUMAnN

    Steps:

    * For each sequence set:

      - Convert sequences and concatenate into a single fastq file
      - Filter sequences for human contaminants with knead_data
      - Perform taxonomic profiling with metaphlan2
      - Infer pathway and gene lists with HUMAnN v2


    Workflows used:

    * :py:func:`anadama_workflows.general.sequence_convert`
    * :py:func:`anadama_workflows.wgs.knead_data`
    * :py:func:`anadama_workflows.wgs.metaphlan2`
    * :py:func:`anadama_workflows.wgs.humann2`

    """

    name = "WGS"
    products = {
        "sample_metadata"            : list(),
        "raw_seq_files"              : list(),
        "intermediate_fastq_files"   : list(),
        "decontaminated_fastq_files" : list(),
        "metaphlan_results"          : list(),
        "otu_tables"                 : list(),
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
        self.products_dir = os.path.abspath(products_dir)

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

        self.sequence_attrs = ("raw_seq_files",
                               "intermediate_fastq_files",
                               "decontaminated_fastq_files")

        def _default_metadata():
            cls = namedtuple("Sample", ['SampleID'])
            return [ cls(basename(util.rmext(f, all=True)))
                     for f in raw_seq_files ]
        self._unpack_metadata(default = _default_metadata)


    def _configure(self):
        for attr in self.sequence_attrs:
            seq_set = getattr(self, attr)

            if self.options['infer_pairs'].get('infer'):
                paired, notpaired = infer_pairs(seq_set)
                seq_set = paired + notpaired

            seq_set, maybe_tasks = maybe_concatenate(seq_set, self.products_dir)
            setattr(self, attr, seq_set)
            for t in maybe_tasks:
                yield t

        for file_ in self.raw_seq_files:
            if util.guess_seq_filetype(file_) != "fastq":
                fastq_file = util.new_file( basename(file_)+"_filtered.fastq",
                                            basedir=self.products_dir )
                yield general.sequence_convert(
                    [file_], fastq_file, 
                    **self.options.get('sequence_convert', dict())
                )
            else:
                fastq_file = file_
            self.intermediate_fastq_files.append(fastq_file)
                

        for fastq_file in self.intermediate_fastq_files:
            name_base = os.path.join(self.products_dir,
                                     basename(fastq_file))
            name_base = util.rmext(name_base, all=True)
            task_dict = wgs.knead_data([fastq_file], name_base).next()
            decontaminated_fastq = first_half(task_dict['targets'])
            self.decontaminated_fastq_files.extend(decontaminated_fastq)
            yield task_dict

        for d_fastq in self.decontaminated_fastq_files:
            metaphlan_file = util.new_file(
                basename(d_fastq)+".metaphlan2.pcl",
                basedir=self.products_dir )
            otu_table = metaphlan_file.replace('.pcl', '.biom')
            import pdb; pdb.set_trace()
            yield wgs.metaphlan2(
                [d_fastq], output_file=metaphlan_file,
                biom=otu_table,
                # first index is for first item in list of samples
                # second index is to get the sample id from the sample
                sample_id=self._filter_samples_for_file(self.sample_metadata,
                                                        d_fastq)[0][0],
                **self.options.get('metaphlan2', dict())
            )
            self.metaphlan_results.append(metaphlan_file)
            self.otu_tables.append(otu_table)

            # Finally, HUMAnN all alignment files
            humann_output_dir = d_fastq+"_humann"
            yield wgs.humann2( d_fastq, humann_output_dir, 
                               **self.options.get('humann', dict()) )
            


def maybe_concatenate(maybe_pairs, products_dir):
    pairs, singles = split_pairs(maybe_pairs)
    tasks = list()

    if not pairs:
        return singles, tasks

    for pair in pairs:
        catted_fname = util.new_file( 
            _to_merged(pair[0], tag="cat", strip_ext=False),
            basedir=products_dir )

        simply_cat = all( util.guess_seq_filetype(s) in ('fastq' , 'fasta')
                          for s in pair )
        if simply_cat:
            tasks.append(general.cat(pair, catted_fname))
        else:
            tasks.append(general.sequence_convert(pair, catted_fname))

        singles.append(catted_fname)

    return singles, tasks


def first_half(list_):
    n = len(list_)
    return list_[:n/2]
