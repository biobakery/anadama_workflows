import os
import re
from os.path import join, basename
from operator import itemgetter
from collections import Counter, namedtuple
from itertools import groupby, izip_longest

from anadama import util
from anadama.pipelines import Pipeline

from .. import settings
from .. import general, sixteen, biom, usearch

from . import (
    SampleFilterMixin, 
    SampleMetadataMixin, 
    infer_pairs,
    split_pairs,
    maybe_convert_to_fastq,
    _to_merged,
    maybe_decompress
)

firstitem = itemgetter(0)

class DemultiplexMixin(object):
            
    def split_illumina_style(self, seqfiles_to_split, barcode_seqfiles):
        demuxed, tasks = list(), list()
        bcode_pairs = zip(seqfiles_to_split, barcode_seqfiles)

        options = self.options.get("demultiplex_illumina", dict())
        if 'barcode_type' not in options:
            options['barcode_type'] = _determine_barcode_type(
                self.sample_metadata)

        do_groupby = options.pop("group_by_sampleid", False)
        for seqfile, bcode_file in bcode_pairs:
            sample_dir = join(self.products_dir, basename(seqfile)+"_split")

            map_fname = util.new_file("map.txt", basedir=sample_dir)
            sample_group = self._filter_samples_for_file(
                self.sample_metadata, seqfile,
                key=lambda val: val.Run_accession)
            tasks.append( sixteen.write_map(
                sample_group, sample_dir, 
                **self.options.get('write_map', dict())
            ) )

            outfile = util.new_file(
                util.rmext(basename(seqfile))+"_demuxed.fna",
                basedir=sample_dir
            ) 
            tasks.append( sixteen.demultiplex_illumina(
                [seqfile], [bcode_file], map_fname, outfile,
                qiime_opts=options
            ) )
            demuxed.append(outfile)

        if do_groupby:
            output_dir = join(self.products_dir, "demuxed_by-sampleid")
            sample_ids = [ s[0] for s in sample_group ]
            task_dict = general.group_by_sampleid(
                demuxed, output_dir, sample_ids
                )
            demuxed = task_dict['targets']
            tasks.append(task_dict)

        return demuxed, tasks
            

    def split_454_style(self, seqfiles_to_split):
        tasks, demuxed = list(), list()
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
                qiime_opts['barcode-type'] = _determine_barcode_type(
                    sample_group)
            demuxed_fname = util.addtag(fasta_fname, "demuxed")
            tasks.append( sixteen.demultiplex(
                map_fname, fasta_fname, qual_fname, demuxed_fname,
                qiime_opts=qiime_opts,
                **self.options.get('demultiplex', dict())
            ))
            demuxed.append(demuxed_fname)

        return demuxed, tasks



class SixteenSPipeline(Pipeline, DemultiplexMixin, SampleFilterMixin,
                       SampleMetadataMixin):

    """Pipeline for analyzing 16S data.

    Steps:

      * Decompress any compressed sequences
      * Paired end reads are stitched.
      * Aggregate samples by SampleID.
      * For each sample:
    
        * Write map.txt entries for that SampleID to its own map.txt file
        * Convert all sequence files for that SampleID into fasta and qual
        * Demultiplex and quality filter
        * Perform closed reference OTU picking against greengenes
        * Infer genes, pathways with picrust

    Workflows used:

      * :py:func:`anadama_workflows.general.extract`
      * :py:func:`anadama_workflows.general.sequence_convert`
      * :py:func:`anadama_workflows.general.fastq_join`
      * :py:func:`anadama_workflows.sixteen.write_map`
      * :py:func:`anadama_workflows.general.fastq_split`
      * :py:func:`anadama_workflows.usearch.filter`
      * :py:func:`anadama_workflows.sixteen.demultiplex`
      * :py:func:`anadama_workflows.sixteen.pick_otus_closed_ref`
      * :py:func:`anadama_workflows.sixteen.picrust`
    """

    name = "16S"
    products = {
        "sample_metadata"         : list(),
        "raw_seq_files"           : list(),
        "barcode_seq_files"       : list(),
        "raw_demuxed_fastq_files" : list(),
        "demuxed_fasta_files"     : list(),
        "otu_tables"              : list()
    }

    default_options = {
        'infer_pairs':         {
            'infer': True
        },
        'write_map':            { },
        'fastq_split':          { },
        'demultiplex':          {
            'qiime_opts': { 
                'M': '2'    
            }
        },
        'fastq_filter':        {
            'fastq_minlen':    200,
            'fastq_truncqual': 25,
            'do_mangle':       True,
            'mangle_to':       None,
        },
        'demultiplex_illumina': { },
        'pick_otus_closed_ref': { },
        'picrust':              { },
    }

    workflows = {
        'infer_pairs':          None,
        'write_map':            None,
        'fastq_split':          general.fastq_split,
        'fastq_filter':         usearch.filter,
        'demultiplex':          sixteen.demultiplex,
        'demultiplex_illumina': sixteen.demultiplex_illumina,
        'pick_otus_closed_ref': sixteen.pick_otus_closed_ref,
        'picrust':              sixteen.picrust
    }

    def __init__(self,
                 sample_metadata=list(),
                 raw_seq_files=list(),
                 barcode_seq_files=list(),
                 raw_demuxed_fastq_files=list(), # not QC'd
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

        :keyword barcode_seq_files: List of strings; File paths to
        barcode fastq files. Typically, an Illumina sequencer will
        give sequences and barcodes in separate files. This is where
        those Illumina barcode files go.

        :keyword raw_demuxed_fastq_files: List of strings; File paths
        to already demuxed but unfiltered sequence files. Adjust the
        fastq_filter options to change how these fastq files are
        filtered.

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

        self.add_products(
            sample_metadata         = sample_metadata,
            raw_seq_files           = raw_seq_files,
            barcode_seq_files       = barcode_seq_files,
            raw_demuxed_fastq_files = raw_demuxed_fastq_files,
            demuxed_fasta_files     = demuxed_fasta_files,
            otu_tables              = otu_tables
        )


        maybe_seqs = (self.raw_seq_files, self.demuxed_fasta_files)
        def _default_metadata():
            cls = namedtuple("Sample", ['SampleID'])
            for seq_list in filter(bool, maybe_seqs):
                # return the first one that contains filenames
                return [ cls(basename(util.rmext(f, all=True)))
                         for f in seq_list ]

        self._unpack_metadata(default = _default_metadata)

        if not products_dir:
            products_dir = settings.workflows.product_directory
        self.products_dir = os.path.abspath(products_dir)

        self.options = self.default_options.copy()
        self.options.update(workflow_options)
        

    def _handle_raw_seqs(self):
        attrs = ("raw_seq_files", "barcode_seq_files",
                 "raw_demuxed_fastq_files")
        for attr in attrs:
            seqs, maybe_tasks = maybe_decompress(getattr(self, attr),
                                                 self.products_dir)
            setattr(self, attr, seqs)
            yield maybe_tasks

        paired_demuxed, single_demuxed = list(), list()
        if self.options['infer_pairs'].get('infer'):
            paired, notpaired = infer_pairs(self.raw_seq_files)
            self.raw_seq_files = paired + notpaired
            paired_demuxed, single_demuxed = infer_pairs(
                self.raw_demuxed_fastq_files)

        packed = maybe_stitch(
            self.raw_seq_files,
            self.products_dir,
            barcode_files=self.barcode_seq_files
        )
        self.raw_seq_files, self.barcode_seq_files, maybe_tasks = packed
        yield maybe_tasks

        if paired_demuxed:
            singles, _, maybe_tasks = maybe_stitch(paired_demuxed,
                                                   self.products_dir)
            yield maybe_tasks
            self.raw_demuxed_fastq_files = singles
        self.raw_demuxed_fastq_files += single_demuxed

        for t in self._process_raw_demuxed_fastq_files():
            yield t


    def _process_raw_demuxed_fastq_files(self):
        for fname in self.raw_demuxed_fastq_files:
            filtered_fname = util.addtag(fname, "filtered")
            opts = self.options.get('fastq_filter', {})
            opts['mangle_to'] = self._filter_samples_for_file(
                self.sample_metadata, fname)[0][0]
            yield usearch.filter(fname, filtered_fname, **opts)
            self.demuxed_fasta_files.append(filtered_fname)


    def _demultiplex(self):
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



    def _configure(self):
        if self.raw_seq_files or self.raw_demuxed_fastq_files:
            for task in self._handle_raw_seqs():
                yield task
        if self.raw_seq_files:
            for task in self._demultiplex():
                yield task

        # ensure all files are decompressed
        # possibly stitch paired reads, demultiplex, and quality filter
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

            

def _determine_barcode_type(sample_group, use_most_common=False):
    lengths = ( len(s.BarcodeSequence) for s in sample_group )
    length_histogram = Counter(lengths)

    if len(length_histogram) == 1 or use_most_common:
        bcode_len = length_histogram.most_common(1)[0][0]
        bcode_len = str(bcode_len)
    else:
        return "variable_length"

    if bcode_len == "12":
        bcode_len = "golay_12"

    return bcode_len



def maybe_stitch(maybe_pairs, products_dir, 
                 barcode_files=list(), drop_unpaired=False):
    pairs, singles = split_pairs(maybe_pairs)
    tasks = list()
    barcodes = list()

    if not pairs:
        return singles, barcode_files, tasks

    pairs = sorted(pairs, key=firstitem)
    barcode_files = sorted(barcode_files)
    for pair, maybe_barcode in izip_longest(pairs, barcode_files):
        (forward, reverse), maybe_tasks = maybe_convert_to_fastq(
            pair, products_dir)
        tasks.extend(maybe_tasks)
        output = util.new_file( 
            _to_merged(forward),
            basedir=products_dir 
        )
        singles.append(output)
        if maybe_barcode and drop_unpaired:
            tasks.append(
                general.fastq_join(forward, reverse, output, 
                                   options={'drop_unpaired': drop_unpaired})
            )
            filtered_barcode = util.new_file(
                util.addtag(maybe_barcode, "filtered"),
                basedir=products_dir
            )
            pairtask = general.sequence_pair(
                maybe_barcode, output,
                outfname1=filtered_barcode,
                options={"inner_join": "right"}
            )
            barcodes.append(filtered_barcode)
            tasks.append(pairtask)
        else:
            tasks.append(
                general.fastq_join(forward, reverse, output,
                                   maybe_barcode,
                                   {'drop_unpaired': drop_unpaired})
            )
            barcodes.append(maybe_barcode)

    return singles, barcodes, tasks


