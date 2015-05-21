import os
from os.path import join

from anadama import util

from .sixteen import SixteenSPipeline

from .. import biom
from .. import sixteen
from ..usearch import pick_otus_closed_ref, truncate


class Usearch16SPipeline(SixteenSPipeline):
    """USEARCH-based 16S pipeline"""

    default_options = {
        'infer_pairs':         {
            'infer': True
        },
        'write_map':            { },
        'truncate':             {
            "trunclen": "215"
        },
        'fastq_split':          { },
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
                "id": "0.97"
            }
        },
        'picrust':              { },
    }


    def _configure(self):
        if self.raw_seq_files:
            for task in self._handle_raw_seqs():
                yield task
            for task in self._demultiplex():
                yield task

        for fasta_fname in self.demuxed_fasta_files:
            # truncate to a uniform length first
            base = util.new_file(os.path.basename(fasta_fname),
                                 basedir=self.products_dir)
            truncated = util.addtag(base, "truncated")
            yield truncate( fasta_fname, fasta_out=truncated,
                            **self.options.get("truncate", {}) )
            fasta_fname = truncated

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

