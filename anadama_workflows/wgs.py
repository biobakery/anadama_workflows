import os.path

from anadama.decorators import requires
from anadama.util import addtag, addext, guess_seq_filetype
from anadama.util import biopython_to_metaphlan, dict_to_cmd_opts, new_file

from . import (
    starters,
    settings
)

@requires(binaries=['humann_init.py'],
          version_methods=["apt-cache show humann "
                           "| awk '/Version: /{ print $NF; }'"])
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
            "cd "+workdir+"; test -f SConstruct || humann_init.py",
            ("cd "+humann_input_dir+"; "
             "ls "+humann_input_dir+" | grep -v 'dat$' | xargs rm"),
            "ln -s %s %s"%(" ".join(infiles_list), humann_input_dir),
            "cd "+workdir+"; scons"
        ]
    }
        
@requires(binaries=['metaphlan2.py'], 
          version_methods=['metaphlan2.py --version'])
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
    def_outfile       = new_file(addext(files_list[0], "metaphlan2"))
    def_bowtie2out    = new_file(addext(files_list[0], "bowtie2out.txt"))

    seqtype = guess_seq_filetype(infiles_list[0])
    if seqtype not in ('fasta', 'fastq'):
        raise ValueError("Need sequences in fasta or fastq format")

    all_opts = { 'bt2_ps'      : 'very-sensitive',
                 'bowtie2db'   : settings.workflows.metaphlan2.bowtie2db,
                 'mpa_pkl'     : settings.workflows.metaphlan2.mpa_pkl,
                 "bowtie2out"  : def_bowtie2out,
                 "output_file" : def_outfile,
                 "input_type"  : biopython_to_metaphlan[seqtype]  }

    all_opts.update(opts)
    
    cmd = starters.cat(infiles_list, guess_from=infiles_list[0])
    cmd += (" | metaphlan2.py"
            + " "+dict_to_cmd_opts(all_opts) )

    targets = [all_opts['output_file'], all_opts['bowtie2out']]
    if 'biom' in opts:
        targets.append(opts['biom'])

    return dict(name     = "metaphlan2:"+all_opts['output_file'],
                actions  = [cmd],
                file_dep = infiles_list,
                targets  = targets )
