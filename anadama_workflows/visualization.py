import os
import re

from anadama.action import CmdAction
from anadama.decorators import requires
from anadama.util import addtag, dict_to_cmd_opts

@requires(binaries=["qiime_cmd"])
def stacked_bar_chart(biom_fname, output_dir, qiime_opts=dict()):
    cmd = ("qiime_cmd summarize_taxa_through_plots.py "
           "-i {} -o {} ".format(biom_fname, output_dir))
    
    default_opts = {
        "force": ""
    }
    default_opts.update(qiime_opts)
    
    opts = dict_to_cmd_opts(default_opts)
    cmd += opts

    yield {
        "name"     : "stacked_bar_chart: "+output_dir,
        "actions"  : [cmd],
        "file_dep" : [biom_fname],
        "targets"  : [os.path.join(output_dir, addtag(biom_fname, "L1"))]
    }


@requires(binaries=["scriptPcoa.py"])
def breadcrumbs_pcoa_plot(pcl_fname, output_plot_fname, **opts):
    pcoa_cmd = ("scriptPcoa.py ")

    default_opts = {
        "meta"       : None,
        "id"         : None,
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
        default_opts['meta'] = default_opts['meta'] or last_meta_name(pcl_fname)
        default_opts['id'] = default_opts['id'] or sample_id(pcl_fname)
        pcoa_cmd += dict_to_cmd_opts(default_opts)
        pcoa_cmd += " "+pcl_fname+" "
        return CmdAction(pcoa_cmd, verbose=True).execute()


    yield {
        "name": "breadcrumbs_pcoa_plot: "+output_plot_fname,
        "actions": [run],
        "file_dep": [pcl_fname],
        "targets": [output_plot_fname]
    }

