import os
from os.path import join

from anadama import util

from .sixteen import SixteenSPipeline

from .. import biom
from .. import sixteen
from ..usearch import pick_otus_closed_ref


class Usearch16SPipeline(SixteenSPipeline):
    """USEARCH-based 16S pipeline"""

    def _configure(self):
        if self.raw_seq_files:
            tasks = self._handle_raw_seqs()
            for task in tasks:
                yield task

        for fasta_fname in self.demuxed_fasta_files:
            otu_table = util.rmext(fasta_fname)+".biom"
            otu_table = join(self.products_dir, os.path.basename(otu_table))
            yield pick_otus_closed_ref(
                fasta_fname, otu_table,
                **self.options.get('pick_otus_closed_ref', dict())
            )
            self.otu_tables.append(otu_table)

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

