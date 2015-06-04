import os
import re

from anadama.action import CmdAction
from anadama.decorators import requires
from anadama.util import addtag, dict_to_cmd_opts

@requires(binaries=["summarize_taxa_through_plots.py"],
          version_methods=["print_qiime_config.py "
                           "| awk '/QIIME library version/{print $NF;}'"])
def stacked_bar_chart(biom_fname, output_dir, qiime_opts=dict()):
    """Workflow to produce stacked bar charts of biom-formatted taxonomic
    profiles using QIIME's `summarize_taxa_through_plots.py`.

    :param biom_fname: String; the file name of a single biom-formatted otu 
                       table or taxonomic profile to be visualized.
    :param output_dir: String; the full path to a directory wherein the 
                       summary plots and charts will be placed
    :keyword qiime_opts: Dictionary; A dictionary of command line options to be
                         passed to the wrapped 
                         summarize_taxa_through_plots.py script. No - or -- 
                         flags are necessary; the correct - or --t flags are 
                         inferred based on the length of the option. For 
                         boolean options, use the key/value pattern of 
                         { "my-option": "" }.

    External dependencies
      - Qiime 1.8.0: https://github.com/qiime/qiime-deploy

    """

    cmd = ("summarize_taxa_through_plots.py "
           "-i {} -o {} ".format(biom_fname, output_dir))
    
    default_opts = {
        "force": ""
    }
    default_opts.update(qiime_opts)
    
    opts = dict_to_cmd_opts(default_opts)
    cmd += opts

    target = os.path.join(
        output_dir,
        addtag( os.path.basename(biom_fname), "L2")
    )

    yield {
        "name"     : "stacked_bar_chart: "+output_dir,
        "actions"  : [cmd],
        "file_dep" : [biom_fname],
        "targets"  : [target]
    }


@requires(binaries=["scriptPcoa.py"],
          version_methods=["md5sum $(which scriptPcoa.py) "
                           "| awk '{print $1;}'"])
def breadcrumbs_pcoa_plot(pcl_fname, output_plot_fname, **opts):
    """Use breadcrumbs `scriptPcoa.py` script to produce principal
    coordinate plots of pcl files.

    :param pcl_fname: String; file name of the pcl-formatted taxonomic profile
                      to visualize via `scriptPcoa.py`.
    :param output_plot_fname: String; file name of the resulting image file.
    :keyword **opts: Any additional keyword arguments are passed to 
                     `scriptPcoa.py`  as command line flags. By default, 
                     it passes `meta=None`, `id=None` and `noShape=None`, 
                     which are converted into `--meta`, `--id`, and 
                     `--noShape`, respectively.

    External dependencies
      - Breadcrumbs: https://bitbucket.org/biobakery/breadcrumbs

    """

    pcoa_cmd = ("scriptPcoa.py ")

    default_opts = {
        "meta"       : True,
        "id"         : True,
        "noShape"    : True,
        "outputFile" : output_plot_fname
    }
    default_opts.update(opts)

    def sample_id(fname):
        id_ = str()
        with open(fname) as f:
            for line in f:
                if line.startswith("#"):
                    id_ = line.split('\t')[0]
                    continue
                else:
                    return id_ or line.split('\t')[0]

    def last_meta_name(fname):
        prev_line = str()
        with open(fname) as f:
            for line in f:
                if re.search(r'[Bb]acteria|[Aa]rchaea.*\s+\d', line):
                    return prev_line.split('\t')[0]
                prev_line = line

            return prev_line.split('\t')[0]

    def run(pcoa_cmd=pcoa_cmd):
        if default_opts['meta'] is True or not default_opts['meta']:
            default_opts['meta'] = last_meta_name(pcl_fname)
        if default_opts['id'] is True or not default_opts['id']:
            default_opts['id'] = sample_id(pcl_fname)
        pcoa_cmd += dict_to_cmd_opts(default_opts)
        pcoa_cmd += " "+pcl_fname+" "
        return CmdAction(pcoa_cmd, verbose=True).execute()

    targets = [output_plot_fname]
    if 'CoordinatesMatrix' in default_opts:
        targets.append(default_opts['CoordinatesMatrix'])

    yield {
        "name": "breadcrumbs_pcoa_plot: "+output_plot_fname,
        "actions": [run],
        "file_dep": [pcl_fname],
        "targets": targets
    }

