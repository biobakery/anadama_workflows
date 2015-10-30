import os
from os.path import join

from anadama import util

from .sixteen import SixteenSPipeline

from .. import biom
from .. import general
from .. import sixteen
from .. import usearch
from ..usearch import pick_otus_closed_ref, truncate



class Usearch16SPipeline(SixteenSPipeline):
    """USEARCH-based 16S pipeline"""

    default_options = {
        'infer_pairs':         {
            'infer': True
        },
        'write_map':            { },
        'fastq_split':          { },
        'fastq_filter':         {
            'fastq_minlen':     200,
            'fastq_truncqual':  25,
            'do_mangle':        True,
            'mangle_to':        None,
        },
        'demultiplex':          {
            'qiime_opts': { 
                'M': '2'    
            }
        },
        'demultiplex_illumina': {
            "group_by_sampleid": True
        },
        'pick_otus_closed_ref': {
            "usearch_closed_opts": {
                "strand": "both",
                "id": "0.90"
            }
        },
        'picrust':              { },
    }

    workflows = {
        'infer_pairs': None,
        'write_map':   None,
        'truncate':  truncate,
        'fastq_split': general.fastq_split,
        'fastq_filter': usearch.filter,
        'demultiplex': sixteen.demultiplex,
        'demultiplex_illumina': sixteen.demultiplex_illumina,
        'pick_otus_closed_ref': pick_otus_closed_ref,
        'picrust':              sixteen.picrust
    }

    def _handle_raw_seqs_and_demultiplex(self):
        """ Process raw sequence files and demultiplex """
        if self.raw_seq_files or self.raw_demuxed_fastq_files:
            for task in self._handle_raw_seqs():
                yield task
        if self.raw_seq_files:
            for task in self._demultiplex():
                yield task


class Usearch32_16SPipeline(Usearch16SPipeline):
    """USEARCH32bit-based 16S pipeline"""

    def _configure(self):
        self._handle_raw_seqs_and_demultiplex()

        for fasta_fname in self.demuxed_fasta_files:
            otu_table = util.rmext(fasta_fname)+"_tax.biom"
            otu_table = join(self.products_dir, os.path.basename(otu_table))
            yield pick_otus_closed_ref(
                fasta_fname, otu_table,
                **self.options.get('pick_otus_closed_ref', dict())
            )
            self.otu_tables.append(otu_table)

        # infer genes and pathways with picrust
        for otu_table in self.otu_tables:
            yield sixteen.picrust(
                otu_table, 
                **self.options.get('picrust', dict())
            )

class Usearch64_16SPipeline(Usearch16SPipeline):
    """USEARCH64bit-based 16S pipeline"""

    def _configure(self):
        self._handle_raw_seqs_and_demultiplex()

        # merge all of the demultiplexed files into a single file
        merged_fasta = util.new_file("all_samples.fa", basedir=self.products_dir)
        yield general.cat(self.demuxed_fasta_files, merged_fasta)

        otu_table_biom = util.new_file("all_samples_otu_tax.biom", basedir=self.products_dir)
        otu_table_tsv = util.new_file("all_samples_otu_tax.tsv", basedir=self.products_dir)
        # run closed reference picking
        yield pick_otus_closed_ref(
            merged_fasta, otu_table_biom,
            out_tsv=otu_table_tsv,
            **self.options.get('pick_otus_closed_ref', dict())
        )

        # infer genes and pathways with picrust
        yield sixteen.picrust(
            otu_table_biom,
            **self.options.get('picrust', dict())
        )

