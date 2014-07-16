import os
from os.path import join, basename
from operator import itemgetter
from collections import Counter
from itertools import chain, groupby

from anadama import util
from anadama.pipelines import Pipeline

from . import settings
from . import general, sixteen, wgs, alignment

class SixteenSPipeline(Pipeline):

    def __init__(self,
                 sample_metadata, # list or possibly string (handled later)
                 raw_seq_files=list(),
                 demuxed_fasta_files=list(), # assumed to be QC'd
                 otu_tables=list(),
                 products_dir=str(),
                 pathway_options=dict(),
                 *args, **kwargs):
        self.sample_metadata = sample_metadata
        self.raw_seq_files = raw_seq_files
        self.demuxed_fasta_files = demuxed_fasta_files
        self.otu_tables = otu_tables

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
            'merge_otu_tables':     {
                'name': 'all_otu_tables_merged.biom'
            }
        }
        self.options.update(pathway_options)
        
        super(SixteenSPipeline, self).__init__(*args, **kwargs)


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
        firstitem = itemgetter(0)
        self.sample_metadata = sorted(self.sample_metadata, key=firstitem)
        for sample_id, sample_group in groupby(self.sample_metadata, firstitem):
            sample_dir = join(self.products_dir, sample_id)
            sample_group = list(sample_group)
            map_fname = util.new_file("map.txt", basedir=sample_dir)
            yield sixteen.write_map(sample_group, sample_dir, 
                                    **self.options.get('write_map', dict()))

            files_list = self._filter_files_for_sample(
                self.raw_seq_files, sample_group)
            fasta_fname = util.new_file(sample_id+".fa", basedir=sample_dir)
            qual_fname = util.new_file(sample_id+".qual", basedir=sample_dir)
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

        # infer genes and pathways with picrust
        for otu_table in self.otu_tables:
            yield sixteen.picrust(
                otu_table, 
                **self.options.get('picrust', dict())
            )

        # now merge all otus together
        if self.otu_tables:
            yield sixteen.merge_otu_tables(
                self.otu_tables,
                output_dir=self.products_dir,
                **self.options.get('merge_otu_tables', dict())
            )


    @staticmethod
    def _filter_files_for_sample(files_list, sample_group):
        try:
            return [ 
                f for f in files_list
                if any(s.Run_accession in f for s in sample_group) 
            ]
        except AttributeError:
            return files_list

    
    @staticmethod
    def _determine_barcode_type(sample_group):
        lengths = ( len(s.BarcodeSequence) for s in sample_group )
        length_histogram = Counter(lengths)

        if len(length_histogram) > 1:
            return "variable_length"
        else:
            bcode_len = length_histogram.most_common(1)[0][0]
            return str(bcode_len)



class WGSPipeline(Pipeline):

    def __init__(self, raw_seq_files=list(), 
                 intermediate_fastq_files=list(),
                 alignment_result_files=list(),
                 products_dir=str(),
                 pathway_options=dict(),
                 *args, **kwargs):
        """raw_seq_files can be a list of lists, if you want to aggregate
        results by something other than by file

        """
        self.raw_seq_files = raw_seq_files
        self.intermediate_fastq_files = intermediate_fastq_files
        self.alignment_result_files = alignment_result_files
        if not products_dir:
            products_dir = settings.workflows.product_directory
        self.products_dir = os.path.realpath(products_dir)

        self.options = {
            'sequence_convert': { }, # no defaults for sequence_convert
            'metaphlan2':       { },
            'bowtie2_align':    { },
            'humann':           { }
        }
        self.options.update(pathway_options)


        super(WGSPipeline, self).__init__(*args, **kwargs)
        

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
            yield wgs.metaphlan2(
                [fastq_file],
                **self.options.get('metaphlan2', dict())
            )
            
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
            sample_dir = util.new_file('SConstruct', basedir=new_basedir)
            yield wgs.humann(
                [alignment_file], workdir=sample_dir,
                **self.options.get('humann', dict())
            )


