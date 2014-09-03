import os

from anadama.decorators import requires
from anadama.util import addtag, dict_to_cmd_opts

@requires(binaries=["qiime_cmd"])
def stacked_bar_chart(biom_fname, output_dir, qiime_opts=dict()):
    cmd = ("qiime_cmd summarize_taxa_through_plots.py "
           "-i {} -o {} ".format(biom_fname, output_dir))
    
    opts = dict_to_cmd_opts(qiime_opts)
    cmd += opts

    yield {
        "name": "stacked_bar_chart: "+output_dir,
        "actions": [cmd],
        "file_dep": [biom_fname],
        "targets": [os.path.join(output_dir, addtag(biom_fname, "L1"))]
    }
