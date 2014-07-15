from anadama.util import dict_to_cmd_opts, addext, new_file

from . import (
    settings
)

def bowtie2_align(infiles_list, output_file, **opts):
    all_opts = { # defaults in here
        "reference_db": settings.workflows.alignment.kegg_bowtie2_db,
        "threads": 2,
    }
    all_opts.update(opts)

    cmd = ("bowtie2 "
           + " -x "+all_opts.pop('reference_db')
           + " -p "+str(all_opts.pop('threads'))
           + " -U "+",".join(infiles_list)
           + " --no-head"
           + " --very-sensitive"
           + " "+dict_to_cmd_opts(all_opts)
           + " > "+output_file)
    
    return {
        "name": "bowtie2_align:"+output_file,
        "actions": [cmd],
        "file_dep": infiles_list,
        "targets": [output_file]
    }
