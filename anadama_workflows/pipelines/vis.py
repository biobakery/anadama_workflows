import os

from anadama import util
from anadama.pipelines import Pipeline

from .. import settings
from .. import sixteen, biom, association, visualization

from . import SampleMetadataMixin

class VisualizationPipeline(Pipeline, SampleMetadataMixin):
    """Pipeline for visualizing taxonomic profiles.  This pipeline can be
    run as is or it can be appended as an augmentation to existing
    pipelines. See the `anadama.pipelines` module documentation for
    help on using the `append` method.

    Steps: 

      * Merge all biom-formatted OTU-tables into a single biom table
      * Add all sample-level metadata to the merged OTU table
      * Convert all merged OTU-tables into pcl files
      * Create stacked bar charts and taxa abundance summaries
      * Create PCoA plots and charts

    Workflows used:

      * :py:func:`anadama_workflows.sixteen.merge_otu_tables`
      * :py:func:`anadama_workflows.biom.add_metadata`
      * :py:func:`anadama_workflows.visualization.stacked_bar_chart`
      * :py:func:`anadama_workflows.biom.to_tsv`
      * :py:func:`anadama_workflows.association.qiime_to_maaslin`
      * :py:func:`anadama_workflows.association.merge_otu_metadata`
      * :py:func:`anadama_workflows.visualization.breadcrumbs_pcoa_plot`

    """


    name = "Visualization"

    products = {
        'sample_metadata'  : str(),
        'otu_tables'       : list(),
        'merged_otu_tables': list(),
        'pcl_files'        : list()
    }

    default_options = {
        'stacked_bar_chart':     { },
        'breadcrumbs_pcoa_plot': {
            "meta"       : True,
            "id"         : True,
            "noShape"    : True,
        },
    }

    workflows = {
        'stacked_bar_chart':     visualization.stacked_bar_chart,
        'breadcrumbs_pcoa_plot': visualization.breadcrumbs_pcoa_plot
    }

    def __init__(self, sample_metadata,
                 otu_tables=list(),
                 pcl_files=list(),
                 merged_otu_tables=list(),
                 workflow_options=dict(),
                 products_dir=str(),
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
        :keyword otu_tables: List of Strings; file paths to biom-formatted 
                             OTU-tables to be visualized.
        :keyword pcl_files: List of Strings; file paths to pcl-formatted 
                            OTU-tables to be visualized.
        :keyword merged_otu_tables: List of Strings; file paths to any already 
                                    merged, already metadata-enriched, 
                                    biom-formatted OTU tables
        :keyword workflow_options: Dictionary; **opts to be fed into the 
                                  respective workflow functions.
        :keyword products_dir: String; Directory path for where outputs will 
                               be saved.

        """

        super(VisualizationPipeline, self).__init__(*args, **kwargs)

        self.options = self.default_options.copy()
        self.options.update(workflow_options)
        
        self.add_products(sample_metadata   = sample_metadata,
                          otu_tables        = otu_tables,
                          merged_otu_tables = merged_otu_tables,
                          pcl_files         = pcl_files)

        if not products_dir:
            products_dir = settings.workflows.product_directory
        self.products_dir = os.path.abspath(products_dir)
                

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
            yield visualization.stacked_bar_chart(
                otu_table, barchart_path,
                **self.options.get('stacked_bar_chart', {}))

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
                CoordinatesMatrix = pcl_file+"_pcoa_coords.txt",
                **self.options.get('breadcrumbs_pcoa_plot', {})
            )
