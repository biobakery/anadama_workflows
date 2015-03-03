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
    """Workflow to find pathway and gene lists from homology search
    and/or genome mapping results with HUMAnN.

    :param infiles_list: List of strings; Paths to files to be fed
                         into HUMAnN.  :param workdir: String;
                         Directory path to where a HUMAnN environment
                         will be created. Input files are softlinked
                         into this directory's 'input'
                         subdirectory. All products can be found in
                         the ``workdir``'s 'output' subdirectory.

    External dependencies - HUMAnN 0.99b:
    https://bitbucket.org/biobakery/humann

    Resource utilization: - Ram: 18-32G

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
        

@requires(binaries=['humann2.py'],
          version_methods=['pip freeze | grep humann2lib'])
def humann2(seqfile_in, output_dir, **opts):
    """Workflow to find pathway and gene lists grouped by organism from
    raw whole genome shotgun reads.

    Additional keywords are interpreted as command line options to be
    passed to the wrapped knead_data.py script.  No - or -- flags are
    necessary; the correct - or --t flags are inferred based on the
    length of the option.  For boolean options, use the key/value
    pattern of { "my-option": "" }.

    :param seqfile_in: String; Paths to file to be fed into HUMAnN2.

    :param output_dir: String; Directory path to where a HUMAnN2
      deposits its results.

    External dependencies:

      - `HUMAnN2 <https://bitbucket.org/biobakery/humann2>`_


    Resource utilization: 

      - Ram: 4-6G
      - Time: 50 hrs (fifty)

    """

    default_opts = {
        "input"  : seqfile_in,
        "output" : os.path.abspath(output_dir),
        
        "uniref"             : settings.workflows.humann2.uniref_path,
        "chocophlan"         : settings.workflows.humann2.chocophlan_path,
        "pathways_databases" : settings.workflows.humann2.pathways_databases,
        "o_log"              : os.path.join(output_dir, "humann2_log.txt"),
        "output_format"      : "tsv"
    }
    default_opts.update(opts)
    opts_str = dict_to_cmd_opts(default_opts, sep=" ")

    cmd = "humann2.py " + opts_str

    suffix = default_opts['output_format']
    def _join(s):
        file, ext = os.path.splitext(seqfile_in)
        file = file + '.' + suffix
        return new_file(addtag(file, s), basedir=default_opts['output'])

    targets = map(_join, ("genefamilies", 
                          "pathcoverage", 
                          "pathabundance"))

    return {
        "name"     : "humann2:"+output_dir,
        "file_dep" : [seqfile_in],
        "targets"  : targets,
        "actions"  : [cmd]
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

    seqtype = guess_seq_filetype(files_list[0])
    if seqtype not in ('fasta', 'fastq'):
        raise ValueError("Need sequences in fasta or fastq format")

    def_base = opts.get("output_file") or files_list[0]
    def_outfile    = new_file(addext(def_base, "metaphlan2"))
    def_bowtie2out = new_file(addext(def_base, "bowtie2out.txt"))

    all_opts = { 'bt2_ps'      : 'very-sensitive',
                 'bowtie2db'   : settings.workflows.metaphlan2.bowtie2db,
                 'mpa_pkl'     : settings.workflows.metaphlan2.mpa_pkl,
                 "bowtie2out"  : def_bowtie2out,
                 "output_file" : def_outfile,
                 "input_type"  : biopython_to_metaphlan[seqtype]  }

    all_opts.update(opts)
    
    cmd = starters.cat(files_list, guess_from=files_list[0])
    cmd += (" | metaphlan2.py"
            + " "+dict_to_cmd_opts(all_opts) )

    targets = [all_opts['output_file'], all_opts['bowtie2out']]
    if 'biom' in opts:
        targets.append(opts['biom'])

    return dict(name     = "metaphlan2:"+all_opts['output_file'],
                actions  = [cmd],
                file_dep = files_list,
                targets  = targets )


@requires(binaries=['knead_data.py', 'bowtie2'],           
          version_methods=["pip freeze | grep knead_datalib",
                           "bowtie2 --version  |head"])
def knead_data(infiles, output_basestr, **opts):
    """Workflow to sanitize host data and otherwise quality filter
    metagenomic reads. Input sequences are mapped against a host
    database using bowtie2; any sequences that map back to the host
    database are discarded.

    Additional keywords are interpreted as command line options to be
    passed to the wrapped knead_data.py script.  No - or -- flags are
    necessary; the correct - or --t flags are inferred based on the
    length of the option.  For boolean options, use the key/value
    pattern of { "my-option": "" }.

    :param infiles: Iterable of strings; File path to the input
      sequences. Should be either a one-length or two-length
      iterable. Two-length iterables are treated as paired-end data.

    :param output_basestr: String; Path to the directory and base
      filename where the output cleaned sequences will be saved.

    
    External dependencies:
      - `knead_data.py <https://bitbucket.org/biobakery/kneaddata>`_
      - `bowtie2 <http://bowtie-bio.sourceforge.net/index.shtml>`_

    Resource utilization:
      - RAM: 3 G

    """
    
    default_opts = {
        "output-prefix": output_basestr,
        "reference-db": settings.workflows.knead.reference_db,
    }
    default_opts.update(opts)
    
    db_base = os.path.basename(settings.workflows.knead.reference_db)
    def _targ(num_tag=None):
        if num_tag:
            base = "_".join([output_basestr, db_base, "clean", num_tag])
        else:
            base = "_".join([output_basestr, db_base, "clean"])
        return base+".fastq"

    infiles_list = list(infiles)
    if len(infiles_list) > 1:
        one, two = infiles_list
        default_opts['1'] = one
        default_opts['2'] = two
        targets = [ _targ("1"), _targ("2") ]
    else:
        default_opts['1'] = infiles_list[0]
        targets = [ _targ() ]

    cmd = "knead_data.py " + dict_to_cmd_opts(default_opts)

    return {
        "name": "knead_data:"+output_basestr,
        "targets": targets,
        "file_dep": infiles_list,
        "actions": [cmd]
    }
