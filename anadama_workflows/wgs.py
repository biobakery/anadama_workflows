import os
from math import log
from glob import glob

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
        
_humann2_default_dbs = None
def _get_humann2_dbs(opts_d):
    global _humann2_default_dbs
    if ("chocophlan" not in opts_d or "uniref" not in opts_d) \
       and not _humann2_default_dbs:
        from re import findall
        from subprocess import check_output
        _humann2_default_dbs = findall(r'DEFAULT: (.+uniref|.+chocophlan)]',
                                       check_output(["humann2", "-h"]))

    if "chocophlan" not in opts_d and "uniref" not in opts_d:
        return _humann2_default_dbs
    elif "chocophlan" not in opts_d and 'uniref' in opts_d:
        return _humann2_default_dbs[0], opts_d['uniref']
    elif "chocophlan" in opts_d and 'uniref' not in opts_d:
        return opts_d['chocophlan'], _humann2_default_dbs[1]
    else: # both uniref and chocophlan are in opts_d
        return opts_d['chocophlan'], opts_d['uniref']


@requires(binaries=['humann2'],
          version_methods=['humann2 --version'])
def humann2(seqfile_in, output_dir, scratch=None, **opts):
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

      - `HUMAnN2 v0.1.9 <https://bitbucket.org/biobakery/humann2>`_


    Resource utilization: 

      - Ram: 4-6G
      - Time: 1 hr

    """

    default_opts = {
        "input"  : seqfile_in,
        "output" : os.path.abspath(output_dir),
        
        "o-log"              : os.path.join(output_dir, "humann2_log.txt"),
        "memory-use"         : "minimum",
        "log-level"          : "INFO",
        "remove-temp-output" : True,
        "output-format"      : "tsv"
    }
    default_opts.update(opts)

    suffix = default_opts['output-format']
    def _join(s):
        file, ext = os.path.splitext(seqfile_in)
        file = file + '.' + suffix
        return new_file(addtag(file, s), basedir=default_opts['output'])

    targets = map(_join, ("genefamilies", 
                          "pathcoverage", 
                          "pathabundance"))

    if scratch:
        old_out = default_opts['output']
        default_opts.pop('output', None)
        dbs = _get_humann2_dbs(default_opts)
        default_opts.pop('chocophlan', None), default_opts.pop('uniref', None)
        cmd = "humann2 " + dict_to_cmd_opts(default_opts, longsep=" ")
        actions = [
            """ tdir=$(mktemp -d -p {sdir});
                cd ${{tdir}}; 
                mkdir -pv ${{tdir}}/dbs;
                cp -rv {dbs} ${{tdir}}/dbs/;
                {humann2} --output ${{tdir}} \
                          --chocophlan ${{tdir}}/dbs/chocophlan \
                          --uniref ${{tdir}}/dbs/uniref;
                mv -iv ${{tdir}}/*.* {final_out};
                rm -rvf ${{tdir}};
            """.format(sdir=scratch, dbs=" ".join(dbs),
                       humann2=cmd, final_out=old_out)
        ]
    else:
        actions = ["humann2 " + dict_to_cmd_opts(default_opts, longsep=" ")]


    def _perfhint(task):
        threads = default_opts.get('threads', 1)
        insize = os.stat(seqfile_in).st_size
        est_reads = insize/100/10e5
        return "{n} estimated mem={mem} time={time}, threads={threads}".format(
            n=task.name,
            mem=5 + (3.5*log(est_reads)),
            time=3.5 + ((2*est_reads)/threads),
            threads=threads)

        
    return {
        "name"     : "humann2:"+output_dir,
        "file_dep" : [seqfile_in],
        "targets"  : targets,
        "actions"  : actions,
        "title"    : _perfhint
    }


@requires(binaries=['metaphlan2.py'], 
          version_methods=['metaphlan2.py --version'])
def metaphlan2(files_list, scratch=None, **opts):
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
    def_base = opts.get("output_file") or files_list[0]
    all_opts = { 'bt2_ps'      : 'very-sensitive',
                 'bowtie2db'   : settings.workflows.metaphlan2.bowtie2db,
                 'mpa_pkl'     : settings.workflows.metaphlan2.mpa_pkl,
                 "bowtie2out"  : new_file(addext(def_base, "bowtie2out.txt")),
                 "output_file" : new_file(addext(def_base, "metaphlan2")) }
    all_opts.update(opts)
    
    if 'input_type' not in all_opts:
        guessed = guess_seq_filetype(files_list[0])
        if guessed not in ('fasta', 'fastq'):
            raise ValueError("Need sequences in fasta or fastq format, "
                             "or provide keyword 'input_type'")
        all_opts['input_type'] = biopython_to_metaphlan[guessed]


    targets = [all_opts['output_file'], all_opts['bowtie2out']]
    if 'biom' in opts:
        targets.append(opts['biom'])

    cmd = starters.cat(files_list, guess_from=files_list[0])
    if scratch:
        db, pkl = all_opts['bowtie2db'], all_opts['mpa_pkl']
        all_opts.pop('bowtie2db', None), all_opts.pop('mpa_pkl', None)
        dbbase, pklbase = map(os.path.basename, (db, pkl))
        cmd += (" | metaphlan2.py"
                + " "+dict_to_cmd_opts(all_opts) )
        actions = [
            """ tdir=$(mktemp -d -p {sdir});
                cd ${{tdir}};
                mkdir -pv ${{tdir}}/dbs;
                cp {db}* {pkl} ${{tdir}}/dbs;
                {cmd} --mpa_pkl ${{tdir}}/dbs/{pklbase} \
                      --bowtie2db ${{tdir}}/dbs/{dbbase};
                rm -rvf ${{tdir}};
            """.format(sdir=scratch, pkl=pkl, db=db, cmd=cmd,
                       pklbase=pklbase, dbbase=dbbase)
        ]
    else:
        cmd += (" | metaphlan2.py"
                + " "+dict_to_cmd_opts(all_opts) )
        actions = [cmd]
    
    def _perfhint(task):
        threads = all_opts.get('nproc', 1)
        insize = sum(os.stat(f).st_size for f in files_list)
        return "{n} estimated mem={mem} time={time}, threads={threads}".format(
            n=task.name,
            mem=1.5*1024,
            time=15+(insize/1.2e9/(threads)),
            threads=threads)


    return dict(name     = "metaphlan2:"+all_opts['output_file'],
                actions  = actions,
                file_dep = files_list,
                targets  = targets,
                title    = _perfhint,)



@requires(binaries=['knead_data.py', 'bowtie2'],           
          version_methods=["pip freeze | grep knead_datalib",
                           "bowtie2 --version  |head"])
def knead_data(infiles, output_basestr, scratch=None, **opts):
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
      - RAM: 4 G

    """
    
    path, base = os.path.split(output_basestr)
    default_opts = {
        "output-prefix": base,
        "output-dir": output_basestr+"_knead",
        "reference-db": settings.workflows.knead.reference_db,
        "strategy": "memory",
        "logging": "WARNING"
    }
    default_opts.update(opts)
    
    db_bases = map(os.path.basename, default_opts['reference-db'])
    def _targets(nums=[None]):
        outdir = default_opts['output-dir']
        prefix = default_opts['output-prefix']
        yield os.path.join(outdir, prefix+".fastq")
        for num in nums:
            for db_base in db_bases:
                to_join = [prefix, db_base, num, "contam.fastq"]
                n = "_".join(filter(bool, to_join))
                yield os.path.join(outdir, n)

    if type(infiles) in (unicode, str):
        infiles_list = [infiles]
    else:
        infiles_list = list(infiles)

    if len(infiles_list) > 1:
        one, two = infiles_list
        default_opts['1'] = one
        default_opts['2'] = two
        targets = list(_targets(nums=[1,2]))
    else:
        default_opts['1'] = infiles_list[0]
        targets = list(_targets())

    if scratch:
        db_patterns = " ".join(s+"*" for s in default_opts['reference-db'])
        db_printf_cmds = " ".join([
            '--reference-db "${{tdir}}/dbs/{}"'.format(db)
            for db in db_bases
        ])
        default_opts.pop("reference-db", None)
        knead = "knead_data.py " + dict_to_cmd_opts(default_opts)
        cmd = """ tdir=$(mktemp -d -p {sdir});
                  cd ${{tdir}};
                  mkdir -pv ${{tdir}}/dbs;
                  cp {db_patterns} ${{tdir}}/dbs/;
                  {knead} {db_printf_cmds};
                  rm -rvf ${{tdir}};
        """.format(sdir=scratch, db_patterns=db_patterns,
                  knead=knead, db_printf_cmds=db_printf_cmds)
    else:
        cmd = "knead_data.py " + dict_to_cmd_opts(default_opts)

    def _perfhint(task):
        threads = default_opts.get('threads', 1)
        insize = sum(os.stat(f).st_size for f in infiles)
        dbsize = sum(os.stat(f).st_size
                     for pat in default_opts['reference-db']
                     for f in glob(pat+"*"))
        return "{n} estimated mem={mem} time={time}, threads={threads}".format(
            n=task.name,
            mem=dbsize/1024/1024,
            time=60+(insize/9e8/(threads)),
            threads=threads)

        
    return {
        "name": "knead_data:"+output_basestr,
        "targets": targets,
        "file_dep": infiles_list,
        "actions": [cmd],
        "title": _perfhint,
    }
