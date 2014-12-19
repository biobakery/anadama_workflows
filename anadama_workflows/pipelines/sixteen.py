import os
import re
from os.path import join, basename
from operator import itemgetter
from collections import Counter
from itertools import groupby

from anadama import util
from anadama.pipelines import Pipeline

from .. import settings
from .. import general, sixteen, biom

from . import SampleFilterMixin, SampleMetadataMixin, maybe_stitch, maybe_decompress


class SixteenSPipeline(Pipeline, SampleFilterMixin, SampleMetadataMixin):

    """Pipeline for analyzing 16S data.

    Steps:

      * Decompress any compressed sequences
      * Aggregate samples by SampleID.
      * For each sample:
    
        * Write map.txt entries for that SampleID to its own map.txt file
        * Convert all sequence files for that SampleID into fasta and qual
        * Demultiplex and quality filter
        * Perform closed reference OTU picking against greengenes
        * Infer genes, pathways with picrust

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
        "barcode_seq_files"   : list(),
        "demuxed_fasta_files" : list(),
        "otu_tables"          : list()
    }

    def __init__(self,
                 sample_metadata,
                 raw_seq_files=list(),
                 barcode_seq_files=list(),
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
                                files. ``raw_seq_files`` can be a list
                                of pairs, if you want to stitch paired-end
                                reads. Supported sequence file formats are 
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
                          barcode_seq_files   = barcode_seq_files,
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
            'demultiplex_illumina': { },
            'pick_otus_closed_ref': { },
            'picrust':              { },
        }
        self.options.update(workflow_options)
        

    def _configure(self):
        # ensure all files are decompressed
        self.raw_seq_files, compressed_files = maybe_decompress(
                self.raw_seq_files)
        self.barcode_seq_files, compressed_files2 = maybe_decompress(
                self.barcode_seq_files)
        for compressed_file in compressed_files + compressed_files2:
            yield general.extract(compressed_file)

        # possibly stitch paired reads, demultiplex, and quality filter
        if self.raw_seq_files:
            self.raw_seq_files, maybe_tasks = maybe_stitch(self.raw_seq_files,
                                                           self.products_dir)
            for t in maybe_tasks:
                yield t

            if self.barcode_seq_files:
                # must be an illumina run, then
                demuxed, tasks = self.split_illumina_style(
                    self.raw_seq_files, self.barcode_seq_files)
            else:
                # no barcode files, assume 454 route
                demuxed, tasks = self.split_454_style(self.raw_seq_files)
            self.demuxed_fasta_files.extend(demuxed)
            for t in tasks:
                yield t

        # do closed reference otu picking
        for fasta_fname in self.demuxed_fasta_files:
            dirname = util.rmext(os.path.basename(fasta_fname))
            otu_dir = join(os.path.dirname(fasta_fname), dirname+"_otus")
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


    def split_illumina_style(self, seqfiles_to_split, barcode_seqfiles):
        demuxed, tasks = list(), list()
        bcode_pairs = zip(sorted(seqfiles_to_split), sorted(barcode_seqfiles))
        map_fname = self._get_or_create_sample_metadata()
        for seqfile, bcode_file in bcode_pairs:
            sample_dir = join(self.products_dir, basename(seqfile)+"_split")
            outfile = util.new_file("seqs.fna", basedir=sample_dir) 
            demuxed.append(outfile)
            tasks.append( sixteen.demultiplex_illumina(
                [seqfile], [bcode_file], map_fname, outfile,
                qiime_opts=self.options.get("demultiplex_illumina",dict())
            ) )

        return demuxed, tasks
            

    def split_454_style(self, seqfiles_to_split):
        tasks, demuxed = list(), list()
        firstitem = itemgetter(0)
        self.sample_metadata = sorted(self.sample_metadata, key=firstitem)
        for sample_id, sample_group in groupby(self.sample_metadata, 
                                               firstitem):
            sample_dir = join(self.products_dir, sample_id)
            sample_group = list(sample_group)
            map_fname = util.new_file("map.txt", basedir=sample_dir)
            tasks.append( sixteen.write_map(
                sample_group, sample_dir, 
                **self.options.get('write_map', dict())
            ) )

            files_list = self._filter_files_for_sample(
                seqfiles_to_split, sample_group)
            fasta_fname = util.new_file(sample_id+".fa", 
                                        basedir=sample_dir)
            qual_fname = util.new_file(sample_id+".qual", 
                                       basedir=sample_dir)
            tasks.append( general.fastq_split(
                files_list, fasta_fname, qual_fname, 
                **self.options.get('fastq_split', dict())
            ) )

            qiime_opts = self.options['demultiplex'].pop('qiime_opts', {})
            if 'barcode-type' not in qiime_opts:
                qiime_opts['barcode-type'] = self._determine_barcode_type(
                    sample_group)
            demuxed_fname = util.new_file("seqs.fna", basedir=sample_dir)
            tasks.append( sixteen.demultiplex(
                map_fname, fasta_fname, qual_fname, demuxed_fname,
                qiime_opts=qiime_opts,
                **self.options.get('demultiplex', dict())
            ))
            demuxed.append(demuxed_fname)

        return demuxed, tasks


    @staticmethod
    def _determine_barcode_type(sample_group):
        lengths = ( len(s.BarcodeSequence) for s in sample_group )
        length_histogram = Counter(lengths)

        if len(length_histogram) > 1:
            return "variable_length"
        else:
            bcode_len = length_histogram.most_common(1)[0][0]
            return str(bcode_len)


