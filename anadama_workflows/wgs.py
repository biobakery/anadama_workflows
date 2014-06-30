import os.path

from anadama.util import addtag, addext, guess_seq_filetype
from anadama.util import biopython_to_metaphlan, dict_to_cmd_opts, new_file

from . import (
    starters,
    settings
)


def humann(infiles_list, workdir):
    
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
    """Workflow to transform WGS sequences to OTU tables, now with doit
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
