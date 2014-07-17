from anadama.util import dict_to_cmd_opts, addext, new_file

from . import (
    settings
)

def bowtie2_align(infiles_list, output_file, **opts):
    """Workflow to use bowtie2 to map a list of input sequence files
    against a bowtie2 database. Additional keyword options are used
    directly as bowtie2 command-line flags.

    :param infiles_list: List of strings; File paths to input search
                         queries as sequences in fastq format
    :param output_file: String; File path to the search results, in 
                        sam format.
    :keyword reference_db: String; File path to the bowtie2 reference 
                           db basename. Fed immediately into bowtie2's
                           -x option.
    :keyword threads: String or int; Number of threads to use when 
                      performing the mapping. Uses bowtie2's -p option.


    External dependencies:
      - Bowtie2 2.2.1: http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

    Resource utilization:
      - Ram: 2.0-3.0G
      - CPU: 1 core; > 1 core depending on 'threads' option
    """

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
