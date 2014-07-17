import os.path

from anadama.util import addtag, addext, guess_seq_filetype
from anadama.util import biopython_to_metaphlan, dict_to_cmd_opts, new_file

from . import (
    starters,
    settings
)


def humann(infiles_list, workdir):
    """Workflow to find pathway and gene lists from homology search and/or
    genome mapping results with HUMAnN.

    :param infiles_list: List of strings; Paths to files to be fed
                         into HUMAnN.
    :param workdir: String; Directory path to where a HUMAnN environment 
                    will be created. Input files are softlinked into this 
                    directory's 'input' subdirectory. All products can be
                    found in the ``workdir``'s 'output' subdirectory.

    External dependencies
      - HUMAnN 0.99b: https://bitbucket.org/biobakery/humann

    Resource utilization:
      - Ram: 18-32G

    """
    _join = lambda f, label: os.path.join(workdir, "output", f+label)
    targets = []
    for infile in infiles_list:
        targets.extend([
            _join(infile, "_04a-hit-keg-mpm-cop-nul-nve-nve-xpe.txt"),
            _join(infile, "_04a-hit-keg-mpt-cop-nul-nve-nve-xpe.txt"),
            _join(infile, "_04b-hit-keg-mpm-cop-nul-nve-nve.txt"),
            _join(infile, "_04b-hit-keg-mpt-cop-nul-nve-nve.txt"),
        ])

    humann_input_dir = os.path.join(workdir, "input")

    return {
        "name": "humann:"+targets[0],
        "file_dep": infiles_list,
        "targets": targets,
        "actions": [
            "cd "+workdir+"; humann_init.py",
            "rm "+os.path.join(humann_input_dir, "*"),
            "ln -s %s %s"%(" ".join(infiles_list), humann_input_dir),
            "cd "+workdir+"; scons"
        ]
    }
        

def metaphlan2(files_list, **opts):
    """Workflow to perform taxonomic profiling from whole metagenome
    shotgun sequences. Additional keyword options are used directly as
    bowtie2 command-line flags.

    :param files_list: List of strings; File paths to input sequences,
                       in fastq format.
    
    External dependencies
      - Metaphlan2 @tip: https://bitbucket.org/biobakery/metaphlan2

    Resource utilization:
      - Ram: 1.5-3.0G

    """

    infiles_list  = files_list
    outfile       = new_file(addext(files_list[0], "metaphlan2"))
    bowtie2out    = new_file(addext(files_list[0], "bowtie2out.txt"))

    all_opts = { 'bt2_ps'   : 'very-sensitive',
                 'bowtie2db': settings.workflows.metaphlan2.bowtie2db,
                 'mpa_pkl'  : settings.workflows.metaphlan2.mpa_pkl,
    }
    all_opts.update(opts)
    all_opts = dict_to_cmd_opts(all_opts)

    seqtype = guess_seq_filetype(infiles_list[0])
    if seqtype not in ('fasta', 'fastq'):
        raise ValueError("Need sequences in fasta or fastq format")

    cmd = starters.cat(infiles_list, guess_from=infiles_list[0])
    cmd += (" | metaphlan2.py"
            + " --bowtie2out="+bowtie2out
            + " --output_file="+outfile
            + " --input_type="+biopython_to_metaphlan[seqtype]
            + " "+all_opts )

    return dict(name     = "metaphlan2:"+outfile,
                actions  = [cmd],
                file_dep = infiles_list,
                targets  = [outfile, bowtie2out])
