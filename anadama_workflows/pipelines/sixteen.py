import os
from os.path import join
from operator import itemgetter
from collections import Counter
from itertools import groupby

from anadama import util
from anadama.pipelines import Pipeline

from .. import settings
from .. import general, sixteen, biom

from . import SampleFilterMixin


class SixteenSPipeline(Pipeline, SampleFilterMixin):

    """Pipeline for analyzing 16S data.

    Steps:

    * Decompress any compressed sequences
    * Aggregate samples by SampleID.
    * For each sample:
    
      - Write map.txt entries for that SampleID to its own map.txt file
      - Convert all sequence files for that SampleID into fasta and qual
      - Demultiplex and quality filter
      - Perform closed reference OTU picking against greengenes
      - Infer genes, pathways with picrust

    Workflows used:

    * anadama_workflows.general.extract
    * anadama_workflows.sixteen.write_map
    * anadama_workflows.general.fastq_split
    * anadama_workflows.sixteen.demultiplex
    * anadama_workflows.sixteen.pick_otus_closed_ref
    * anadama_workflows.sixteen.picrust
    """

    name = "16S"
    products = {
        "sample_metadata"     : list(),
        "raw_seq_files"       : list(),
        "demuxed_fasta_files" : list(),
        "otu_tables"          : list()
    }

    def __init__(self,
                 sample_metadata,
                 raw_seq_files=list(),
                 demuxed_fasta_files=list(), # assumed to be QC'd
                 otu_tables=list(),
                 products_dir=str(),
                 workflow_options=dict(),
                 *args, **kwargs):

        """Initialize the pipeline.

        :param sample_metadata: String or list of namedtuples; sample metadata 
                                as deserialized by 
                                ``anadama.util.deserialize_map_file``. If a 
                                string is given, the string is treated as a path
                                to the qiime-formatted map.txt for all samples. 
                                For more information about sample-level 
                                metadata, refer to the qiime documentation: 
                                http://qiime.org/tutorials/tutorial.html#mapping-file-tab-delimited-txt
        :keyword raw_seq_files: List of strings; File paths to raw 16S sequence
                                files. Supported sequence file formats are 
                                whatever biopython recognizes.
        :keyword demuxed_fasta_files: List of strings; File paths to 
                                      demultiplexed fasta files. Assumed to be
                                      quality-checked and ready to be used in
                                      OTU picking.
        :keyword otu_tables: List of strings; File paths to produced OTU tables
                             ready to be fed into picrust for gene, pathway 
                             inference.
        :keyword products_dir: String; Directory path for where outputs will 
                               be saved.
        :keyword workflow_options: Dictionary; **opts to be fed into the 
                                  respective workflow functions.
        """

        super(SixteenSPipeline, self).__init__(*args, **kwargs)

        self.add_products(sample_metadata     = sample_metadata,
                          raw_seq_files       = raw_seq_files,
                          demuxed_fasta_files = demuxed_fasta_files,
                          otu_tables          = otu_tables)

        if type(sample_metadata) is str:
            self.sample_metadata = util.deserialize_map_file(samples)

        if not products_dir:
            products_dir = settings.workflows.product_directory
        self.products_dir = os.path.realpath(products_dir)

        self.options = {
            'write_map':            { },
            'fastq_split':          { },
            'demultiplex':          {
                'qiime_opts': { 
                    'M': '2'    
                }
            },
            'pick_otus_closed_ref': { },
            'picrust':              { },
        }
        self.options.update(workflow_options)
        

    def _configure(self):
        # ensure all files are decompressed
        compressed_files = util.filter_compressed(self.raw_seq_files)
        if compressed_files:
            self.raw_seq_files = [
                os.path.splitext(f)[0] if util.is_compressed(f) else f
                for f in self.raw_seq_files
            ]
            yield general.extract(compressed_files)

        # split to fasta, qual, map.txt triplets and demultiplex
        if self.raw_seq_files:
            firstitem = itemgetter(0)
            self.sample_metadata = sorted(self.sample_metadata, key=firstitem)
            for sample_id, sample_group in groupby(self.sample_metadata, 
                                                   firstitem):
                sample_dir = join(self.products_dir, sample_id)
                sample_group = list(sample_group)
                map_fname = util.new_file("map.txt", basedir=sample_dir)
                yield sixteen.write_map(
                    sample_group, sample_dir, 
                    **self.options.get('write_map', dict())
                )

                files_list = self._filter_files_for_sample(
                    self.raw_seq_files, sample_group)
                fasta_fname = util.new_file(sample_id+".fa", 
                                            basedir=sample_dir)
                qual_fname = util.new_file(sample_id+".qual", 
                                           basedir=sample_dir)
                yield general.fastq_split(
                    files_list, fasta_fname, qual_fname, 
                    **self.options.get('fastq_split', dict())
                )

                qiime_opts = self.options['demultiplex'].pop('qiime_opts', {})
                if 'barcode-type' not in qiime_opts:
                    qiime_opts['barcode-type'] = self._determine_barcode_type(
                        sample_group)
                demuxed_fname = util.new_file("seqs.fna", basedir=sample_dir)
                yield sixteen.demultiplex(
                    map_fname, fasta_fname, qual_fname, demuxed_fname,
                    qiime_opts=qiime_opts,
                    **self.options.get('demultiplex', dict())
                )
                self.demuxed_fasta_files.append(demuxed_fname)

        # do closed reference otu picking
        for fasta_fname in self.demuxed_fasta_files:
            otu_dir = join(os.path.dirname(fasta_fname), "otus")
            yield sixteen.pick_otus_closed_ref(
                input_fname=fasta_fname, output_dir=otu_dir,
                **self.options.get('pick_otus_closed_ref', dict())
            )
            self.otu_tables.append(join(otu_dir, "otu_table.biom"))

        # convert biom file to tsv
        for otu_table in self.otu_tables:
            tsv_filename = otu_table+".tsv"
            yield biom.biom_to_tsv(
                otu_table, 
                tsv_filename)

        # infer genes and pathways with picrust
        for otu_table in self.otu_tables:
            yield sixteen.picrust(
                otu_table, 
                **self.options.get('picrust', dict())
            )

    @staticmethod
    def _determine_barcode_type(sample_group):
        lengths = ( len(s.BarcodeSequence) for s in sample_group )
        length_histogram = Counter(lengths)

        if len(length_histogram) > 1:
            return "variable_length"
        else:
            bcode_len = length_histogram.most_common(1)[0][0]
            return str(bcode_len)


