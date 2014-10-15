import os
from os.path import join, basename
from operator import itemgetter
from collections import Counter
from itertools import chain, groupby

from anadama import util
from anadama.pipelines import Pipeline

from . import settings
from . import general, sixteen, wgs
from . import alignment, visualization, association, biom


def _filter(func, iterable):
    try:
        filtered = filter(func, iterable)
    except AttributeError:
        return iterable
    return filtered or iterable

class SampleFilterMixin(object):

    @staticmethod
    def _filter_files_for_sample(files_list, sample_group, 
                                 key=lambda val: getattr(val, "Run_accession")):
        return _filter(lambda f: any(key(s) in f for s in sample_group),
                       files_list )

    @staticmethod
    def _filter_samples_for_file(sample_group, file_, 
                                 key=lambda val: val[0]):
        return _filter(lambda sample: key(sample) in file_,
                       sample_group)



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

      * Finally, merge all OTU tables

    Workflows used:
      * anadama_workflows.general.extract
      * anadama_workflows.sixteen.write_map
      * anadama_workflows.general.fastq_split
      * anadama_workflows.sixteen.demultiplex
      * anadama_workflows.sixteen.pick_otus_closed_ref
      * anadama_workflows.sixteen.picrust
      * anadama_workflows.sixteen.merge_otu_tables
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



class WGSPipeline(Pipeline, SampleFilterMixin):

    """Pipeline for analyzing whole metagenome shotgun sequence data.
    Produces taxonomic profiles with metaphlan2 and gene, pathway
    lists with HUMAnN

    Steps:
      * For each sequence set:

        - Convert the sequences and concatenate into a single fastq file
        - Filter sequences for length > 60 bases
        - Perform taxonomic profiling with metaphlan2
        - Align the fastq file agains the KEGG proks reduced dataset with bowtie2
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


class VisualizationPipeline(Pipeline):
    
    name = "Visualization"

    products = {
        'sample_metadata'  : str(),
        'otu_tables'       : list(),
        'merged_otu_tables': list(),
        'pcl_files'        : list()
    }

    def __init__(self, sample_metadata,
                 otu_tables=list(),
                 pcl_files=list(),
                 merged_otu_tables=list(),
                 workflow_options=dict(),
                 products_dir=str(),
                 *args, **kwargs):
        super(VisualizationPipeline, self).__init__(*args, **kwargs)

        self.options = {
            'stacked_bar_chart': { }
        }
        self.options.update(workflow_options)
        
        self.add_products(sample_metadata   = sample_metadata,
                          otu_tables        = otu_tables,
                          merged_otu_tables = merged_otu_tables,
                          pcl_files         = pcl_files)

        if not products_dir:
            products_dir = settings.workflows.product_directory
        self.products_dir = os.path.realpath(products_dir)


    @property
    def _inferred_sample_metadata_fname(self):
        if self.otu_tables:
            base_input = self.otu_tables[0]
        elif self.pcl_files:
            base_input = self.pcl_files[0]
        else:
            raise ValueError("Unable to infer map.txt file location"
                             " because pipeline inputs are empty")

        dir_ = os.path.dirname(base_input)
        return os.path.join(dir_, "map.txt")


    def _get_or_create_sample_metadata(self):
        if type(self.sample_metadata) is not str:
            sample_metadata_fname = self._inferred_sample_metadata_fname
            try:
                util.serialize_map_file(self.sample_metadata, 
                                        sample_metadata_fname)
            except IndexError as e:
                raise ValueError("The provided sample metadata is not in list"
                                 " format, nor is it a string. Sample_metadata"
                                 " should either be a string for a map.txt"
                                 " filename or a list of namedtuples"
                                 " representing the sample metadata")
            return sample_metadata_fname
        else:
            if not os.path.exists(self.sample_metadata):
                raise ValueError("The provided sample metadata file "
                                 "does not exist: "+self.sample_metadata)
            else:
                return self.sample_metadata
                

    def _configure(self):

        if self.otu_tables:
            merged_name = util.addtag(self.otu_tables[0], "merged")
            merged_file = util.new_file(merged_name, basedir=self.products_dir)
            yield sixteen.merge_otu_tables(
                self.otu_tables,
                name=merged_file
            )
            meta_biom_name = util.addtag(merged_file, "meta")
            yield biom.add_metadata(
                merged_file, meta_biom_name, 
                self._get_or_create_sample_metadata()
            )
            self.merged_otu_tables.append(meta_biom_name)

        for otu_table in self.merged_otu_tables:
            barchart_path = util.new_file(
                otu_table+"_barcharts", basedir=self.products_dir)
            yield visualization.stacked_bar_chart(otu_table, barchart_path)

            tsv_filename = otu_table+".tsv"
            yield association.biom_to_tsv(otu_table, tsv_filename)
            nice_tsv_filename = util.addtag(tsv_filename, 'maaslin')
            yield association.qiime_to_maaslin(tsv_filename, nice_tsv_filename)
            pcl_filename = otu_table+".pcl"
            yield association.merge_otu_metadata(
                nice_tsv_filename, 
                self._get_or_create_sample_metadata(),
                pcl_filename
            )
            self.pcl_files.append(pcl_filename)

        for pcl_file in self.pcl_files:
            yield visualization.breadcrumbs_pcoa_plot(
                pcl_file, pcl_file+"_pcoa_plot.png",
                CoordinatesMatrix = pcl_file+"_pcoa_coords.txt"
            )

                
        
            
