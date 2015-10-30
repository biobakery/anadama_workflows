import os
from os.path import join
from itertools import imap
from collections import Counter
from operator import itemgetter

from anadama import strategies
from anadama.action import CmdAction
from anadama.decorators import requires
from anadama.strategies import if_exists_run
from anadama.util import addtag, rmext, dict_to_cmd_opts

from . import settings
from .sixteen import assign_taxonomy

snd = itemgetter(1)
statsum = lambda fs: sum(os.stat(f).st_size for f in fs)

def usearch_dict_flags(opts_dict):
    opts = [ " -%s %s " %(key, val) 
             for key, val in opts_dict.iteritems() ]
    return " ".join(opts)

def usearch_rusage(input_seqs, time_multiplier=1, threads=1):
    def _titlefunc(task):
        msg = task.name+(" Estimated mem={mem:.2f} "
                         "time={time:.0f} threads={threads:.0f}")
        return msg.format(
            mem=100 + (statsum(task.file_dep)/1024/1024.),
            time=10 + (statsum(input_seqs)*5e-7*time_multiplier),
            threads=threads
        )
    return _titlefunc


from itertools import imap
from collections import Counter
from operator import itemgetter

snd = itemgetter(1)


class util:
    @staticmethod
    def fasta_sequences(seqs_f):
        id = None
        lines = iter(seqs_f)
        line = next(lines).strip()
        while True:
            if line.startswith(">"):
                id = line.split(None, 1)[0]
                id = id.replace(">", "")
                line = next(lines).strip()
            else:
                seq = str()
                while not line.startswith(">"):
                    seq += line
                    line = next(lines).strip()
                yield (id, seq)
                id, seq = str(), str()


    @staticmethod
    def hist(fname):
        with open(fname) as f:
            seqs = util.fasta_sequences(f)
            cnts = Counter( imap(len, imap(snd, seqs)) )
        return sorted(cnts.iteritems(), reverse=True)


    @staticmethod
    def cumsum(it):
        sm = 0
        for item in it:
            sm += item
            yield sm


    @staticmethod
    def cutoff(fname, cutoff_ratio=.95):
        h = util.hist(fname)
        read_length_counts = map(snd, h)
        read_depth = sum(read_length_counts)
        for i, sm in enumerate(util.cumsum(read_length_counts)):
            if float(sm)/read_depth > cutoff_ratio:
                return h[i][0]



@requires(binaries=['usearch7', 'sequence_pair'],
          version_methods=["usearch7 -version"
                           "pip freeze | grep anadama_workflows"])
def stitch(input_fastq_pair, output_fastq, verbose=True, 
                   remove_tempfiles=True, **opts):
    """Workflow to stitch together paried reads using USEARCH version
    7. The USEARCH binary should be named usearch7 in order for this
    workflow to operate.  USEARCH version 7 expects the two sequence
    read files to be completely aligned and will fail if this
    condition is not met. The `usearch_stitch` workflow tries to
    surmount this 'feature' by using the following rubric:

    1. Run usearch 7
    2. If 1. fails, re-sort the two sequence read files and drop unpaired
       reads using the `sequence_pair` script. Run usearch 7 again.
    3. If 2. fails, just copy the first sequence file to the `output_fastq` 
       file.

    :param input_fastq_pair: 2-Tuple of strings; file names of Input sequence
                             read files. Forward reads should be the first 
                             item, while reverse reads should be the second 
                             item in the tuple
    :param output_fastq: String; name of resulting stitched fastq file
    :keyword verbose: Boolean; if true, print commands at runtime as they 
                      are executed. Defaults to true.
    :keyword remove_tempfiles: Boolean; if true, removes any files produced 
                               by `sequence_pair`
    :keyword **opts: Any additional keyword arguments are passed to usearch7 
                     as command line flags. By default, it passes 
                     `fastq_truncqual=10` and `fastq_maxdiffs=3` as 
                     '-fastq_truncqual 10' and '-fastq_maxdiffs 3', 
                     respectively.

    External dependencies
      - USEARCH v7: http://www.drive5.com/usearch
      - sequence_pair: Sequence pairing script that should come preinstalled 
        with the `anadama_workflows` python module.

    """

    stitch_cmd = ("usearch7"+ 
                  " -fastq_mergepairs {r1}"
                  " -reverse {r2}"
                  " {opts_str}"
                  " -fastqout {output_fastq}")

    default_options = {
        "fastq_truncqual": "10",
        "fastq_maxdiffs": "3"
    }
    default_options.update(opts)
    
    opts_str = usearch_dict_flags(default_options)

    pair_cmd = ("sequence_pair -f fastq -t fastq"
                " -1 {r1out} -2 {r2out}"
                " {r1} {r2}")

    def run():
        r1, r2 = input_fastq_pair
        paired = lambda s: s.replace('.fastq', '_paired.fastq')
        ret = strategies.backup(
            (CmdAction(stitch_cmd.format(r1=r1, r2=r2, opts_str=opts_str,
                                         output_fastq=output_fastq),
                       verbose=verbose),
             strategies.Group(
                 CmdAction(pair_cmd.format(r1=r1, r2=r2,
                                           r1out=paired(r1),
                                           r2out=paired(r2)), 
                           verbose=verbose),
                 CmdAction(stitch_cmd.format(r1=paired(r1),
                                             r2=paired(r2),
                                             opts_str=opts_str,
                                             output_fastq=output_fastq),
                           verbose=verbose)),
             CmdAction('cp {r1} {out}'.format(r1=r1, out=output_fastq),
                       verbose=verbose)),
        )
            
        if remove_tempfiles:
            for f in (paired(r1), paired(r2)):
                if os.path.exists(f):
                    os.remove(f)
        return ret


    return { "name"    : "usearch_stitch: "+output_fastq,
             "actions" : [run],
             "file_dep": input_fastq_pair,
             "targets" : [output_fastq],
             "title"   : usearch_rusage(input_fastq_pair)}


@requires(binaries=['usearch7'],
          version_methods=["usearch7 -version"])
def filter(input_fastq, output_fasta, verbose=True, do_mangle=False,
           mangle_to=None, **opts):
    """Filter a fastq file, outputting sequences as fasta, using USEARCH
    version 7. The USEARCH binary should be named usearch7 in order
    for this workflow to operate.

    :param input_fastq: String; file name of a single fastq file to be filtered.
    :param output_fasta: String; name of resulting filtered fasta file

    :keyword verbose: Boolean; if true, print commands at runtime as
    they are executed. Defaults to true.

    :keyword do_mangle: Boolean; if true, mangle the header for each
    sequence to only include a string (defaults to the name of the
    output file minus the file extension) and the number of the
    sequence. Sequence headers will look something like myfastaseqs_1.

    :keyword mangle_to: String; The base string to use when mangling
    sequence headers. Defaults to the name of the output file minus
    the file extension.

    :keyword **opts: Any additional keyword arguments are passed to
    usearch7 as command line flags. By default, it passes
    `fastq_minlen=200` and `fastq_truncqual=25` as '-fastq_minlen 200'
    and '-fastq_truncqual 25', respectively.

    External dependencies
      - USEARCH v7 or greater http://www.drive5.com/usearch

    """
    
    def _maybe_mangle():
        if not do_mangle:
            return
        if not os.path.exists(output_fasta) or os.stat(output_fasta).st_size <1:
            return
        m = rmext(output_fasta, all=True) if mangle_to is False else mangle_to
        cmd = "sequence_convert -m {m} -f fasta -t fasta {o} > {o}.tmp".format(
            m=m, o=output_fasta)
        CmdAction(cmd, verbose=verbose).execute()
        CmdAction("mv {o}.tmp {o}".format(o=output_fasta),
                  verbose=verbose).execute()
        

    cmd = ("usearch7"+
           " -fastq_filter "+input_fastq+
           " -fastaout "+output_fasta)

    default_options = {
        "fastq_minlen": "200",
        "fastq_truncqual": "25"
    }
    default_options.update(opts)

    cmd += usearch_dict_flags(default_options)

    def run():
        if os.stat(input_fastq).st_size > 1:
            ret = CmdAction(cmd, verbose=verbose).execute()
            if ret is None or not issubclass(type(ret), Exception):
                _maybe_mangle()
        else:
            open(output_fasta, 'w').close()
        

    return { "name"     : "usearch_filter: "+output_fasta,
             "actions"  : [run],
             "file_dep" : [input_fastq],
             "targets"  : [output_fasta],
             "title"   : usearch_rusage([input_fastq])}


@requires(binaries=["usearch8"],
          version_methods="usearch8 -version")
def truncate(fastx_in, fasta_out=None, fastq_out=None, **opts):
    default_opts = dict([
        ("trunclen", "215"),
    ]+list(opts.items()))

    if fasta_out:
        out = fasta_out
        default_opts['fastaout'] = fasta_out
    elif fastq_out:
        out = fastq_out
        default_opts['fastqout'] = fastq_out
    else:
        raise ValueError("must provide either fasta_out or fastq_out arguments")

    cmd = ("usearch8 -fastx_truncate "+fastx_in+
           " "+usearch_dict_flags(default_opts))

    return { "name": "usearch_truncate: "+out,
             "actions": [cmd],
             "targets": [out],
             "file_dep": [fastx_in],
             "title"   : usearch_rusage([fastx_in]) }


@requires(binaries=["usearch8", "uclust_otutable"],
          version_methods=["usearch8 -version",
                           "pip freeze | fgrep anadama_workflows"])
def pick_denovo_otus(fasta_in, otutab_out, keep_tempfiles=False,
                     strand="plus", log_file=None, resume=False,
                     quiet=False, chimera_standard=None,
                     truncate_opts={}, derep_opts={}, sort_opts={},
                     cluster_opts={}, chimera_opts={}, map_opts={}):

    opts = dict(input=fasta_in, output=otutab_out, print_cmd=True)
    if 'db' in chimera_opts:
        s = chimera_opts['db']
    elif bool(chimera_standard) is True:
        s = chimera_standard
    else:
        s = settings.workflows.usearch.chimera_gold_standard
    opts['chimera_standard'] = s
    opts['tmp_dir'] = otutab_out+"_usearch"
    kvopts = [("truncate_opts", truncate_opts), ("derep_opts", derep_opts),
              ("sort_opts", sort_opts),         ("cluster_opts", cluster_opts),
              ("chimera_opts", chimera_opts), ("map_opts", map_opts)]
    for name, value in kvopts:
        if value:
            s = " ".join('='.join(pair) for pair in value.iteritems())
            opts[name] = "'" + s + "'"

    if not log_file:
        log_file = opts['tmp_dir']+".log"
    opts['log_file'] = log_file
    
    cmd = "usearch_denovo_otus "+dict_to_cmd_opts(opts)
    targets = [otutab_out, join(opts['tmp_dir'], "nonchimeric.fa"), log_file]
    def _run():
        ret = CmdAction(cmd).execute()
        if ret is None or not issubclass(type(ret), Exception):
            if not keep_tempfiles:
                for f in os.listdir(opts['tmp_dir']):
                    if f != "nonchimeric.fa" and \
                       os.path.isfile(join(opts['tmp_dir'], f)):
                        os.remove(join(opts['tmp_dir'], f))
        else:
            for t in targets:
                if not os.path.exists(t):
                    open(t, 'w').close()
        return ret
                        
            

    return { "name"     : "usearch_pick_denovo_otus: "+otutab_out,
             "actions"  : [run],
             "file_dep" : [fasta_in],
             "targets"  : targets,
             "title"   : usearch_rusage([fasta_in]) }


@requires(binaries=["usearch8", "biom"])
def pick_otus_closed_ref(in_fasta, out_biom, out_tsv=None,
                         non_chimeric_otu_seqs=None,
                         denovo_otu_txt=None,
                         sample_metadata_fname=None,
                         taxonomy_fname=None, ref_fasta=None,
                         keep_tempfiles=False, strand='plus',
                         chimera_standard=None, log_file=None,
                         resume=False, quiet=False, tmp_folder=None,
                         usearch_closed_opts={}, denovo_opts={}):
    
    if not tmp_folder:
        tmp_folder = in_fasta+"_usearch"
    if not log_file:
        log_file = in_fasta + "_usearch.log"
    if not taxonomy_fname:
        taxonomy_fname = settings.workflows.sixteen.otu_taxonomy
    if not ref_fasta:
        ref_fasta = settings.workflows.usearch.otu_db
    if not out_tsv:
        out_tsv = out_biom+".tsv"

    if 'db' in denovo_opts.get('chimera_opts', {}):
        chimera = denovo_opts['chimera_opts']['db']
    elif 'chimera_standard' in denovo_opts:
        chimera = denovo_opts['chimera_standard']
    elif bool(chimera_standard) is True:
        chimera = chimera_standard
    else:
        chimera = settings.workflows.usearch.chimera_gold_standard

    opts = dict(input=in_fasta, output=out_tsv,
                taxonomy=taxonomy_fname, reference=ref_fasta,
                strand=strand, chimera_standard=chimera,
                quiet=quiet, print_cmd=True, log_file=log_file,
                denovo_otu_table=denovo_otu_txt, resume=resume,
                keep_tempfiles=True, tmp_dir=tmp_folder,
                otu_sequences=non_chimeric_otu_seqs)

    kvopts = list(denovo_opts.items())+[("closed_opts",usearch_closed_opts)]
    for name, value in kvopts:
        if value:
            s = " ".join('='.join(pair) for pair in value.iteritems())
            opts[name] = "'" + s + "'"

    usearch_cmd = "uclust_closed_otus "+dict_to_cmd_opts(opts)

    biom_cmd = ("biom convert -i "+out_tsv+" -o "+out_biom+
                " --table-type='OTU Table' --process-obs-metadata=taxonomy"+
                " --output-metadata-id=taxonomy")

    if sample_metadata_fname:
        biom_cmd += " --sample-metadata-fp=" + sample_metadata_fname

    def _run():
        ret = CmdAction(usearch_cmd).execute()
        if ret is None or not issubclass(type(ret), Exception):
            ret = CmdAction(biom_cmd).execute()
            if not keep_tempfiles:
                for f in os.listdir(opts['tmp_dir']):
                    if f != "nonchimeric.fa" and \
                       os.path.isfile(join(opts['tmp_dir'], f)):
                        os.remove(join(opts['tmp_dir'], f))
        else:
            for t in targets:
                if not os.path.exists(t):
                    open(t, 'w').close()
        return ret

    file_dep = [in_fasta, taxonomy_fname, ref_fasta, chimera]
    targets = [out_biom, out_tsv,
               join(opts['tmp_dir'], "nonchimeric.fa"), log_file]
    yield { "name": "usearch_pick_otus_closed_ref: "+out_biom,
            "targets": targets,
            "actions": [_run],
            "file_dep": file_dep,
            "title"   : usearch_rusage(file_dep) }
    

