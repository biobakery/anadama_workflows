import os
from os.path import join
from itertools import imap
from collections import Counter
from operator import itemgetter

from anadama import strategies
from anadama.action import CmdAction
from anadama.decorators import requires
from anadama.strategies import if_exists_run
from anadama.util import addtag, rmext

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
def pick_denovo_otus(fasta_in, otutab_out, remove_tempfiles=True,
                     strand="plus", truncate_opts={}, derep_opts={},
                     sort_opts={}, cluster_opts={}, chimera_opts={},
                     map_opts={}):

    tmpfolder = otutab_out+"_usearch"
    mkdir_cmd = "test -d "+tmpfolder+" || mkdir "+tmpfolder

    trunc_out = join(tmpfolder, "truncated.fa")
    def _truncate_cmd():
        trunclen = util.cutoff(fasta_in)
        truncate_opts['trunclen'] = str(trunclen)
        truncate_opts['fastaout'] = trunc_out
        return CmdAction("usearch8 -fastx_truncate "+fasta_in+
                         " "+usearch_dict_flags(truncate_opts),
                         verbose=True
        ).execute()


    default_derep_opts = dict([
        ("sizeout", ""),
    ]+list(derep_opts.items()))
    
    derep_out = join(tmpfolder, "derep.fa")
    derep_cmd = ("usearch8 -derep_fulllength "+trunc_out+
                 " -fastaout "+derep_out+
                 " "+usearch_dict_flags(default_derep_opts))

    default_sort_opts = dict([
        ("minsize","2"),
    ]+list(sort_opts.items()))

    sort_out = join(tmpfolder, "sorted.fa")
    sort_cmd = ("usearch8 -sortbysize "+derep_out+" -fastaout "+sort_out+
                " "+usearch_dict_flags(default_sort_opts))

    default_cluster_opts = dict([
        ("relabel", "OTU_"),
        ("sizein", ""),
        ("sizeout", ""),
    ]+list(cluster_opts.items()))
                
    cluster_otus_out = join(tmpfolder, "otus.fa")
    cluster_otus_log = join(tmpfolder, "cluster_results.txt")
    cluster_cmd = ("usearch8 -cluster_otus "+sort_out+
                   " -otus "+cluster_otus_out+
                   " -uparseout "+cluster_otus_log+
                   " "+usearch_dict_flags(default_cluster_opts))

    default_chimera_opts = dict([
        ("db", settings.workflows.usearch.chimera_gold_standard),
        ("strand", strand),
    ]+list(chimera_opts.items()))

    chimera_out = join(tmpfolder, "nonchimeric.fa")
    chimera_cmd = ("usearch8 -uchime_ref "+cluster_otus_out+
                   " -nonchimeras "+chimera_out+
                   " "+usearch_dict_flags(default_chimera_opts))

    default_map_opts = dict([
        ("id", "0.97"),
        ("strand", strand),
    ]+list(map_opts.items()))

    map_out = join(tmpfolder, "mapping_results.uc")
    map_cmd = ("usearch8 -usearch_global "+fasta_in+
               " -db "+chimera_out+
               " -uc "+map_out+
               " "+usearch_dict_flags(default_map_opts))

    otu_cmd = ("uclust_otutable %s > %s"%(map_out, otutab_out))

    cleanup_cmd = "rm " + " ".join((derep_out, sort_out,
                                    cluster_otus_out, map_out))

    actions = [ mkdir_cmd, _truncate_cmd, derep_cmd, sort_cmd,
                cluster_cmd, chimera_cmd, map_cmd, otu_cmd ]
    if remove_tempfiles:
        actions.append(cleanup_cmd)

    return { "name"     : "usearch_pick_denovo_otus: "+otutab_out,
             "actions"  : actions,
             "file_dep" : [fasta_in],
             "targets"  : [otutab_out, chimera_out],
             "title"   : usearch_rusage([fasta_in]) }


@requires(binaries=["usearch8", "biom"])
def pick_otus_closed_ref(in_fasta, out_biom,
                         out_tsv=None,
                         non_chimeric_otu_seqs=None,
                         denovo_otu_txt=None,
                         sample_metadata_fname=None,
                         taxonomy_fname=None,
                         ref_fasta=None,
                         remove_tempfiles=True,
                         usearch_closed_opts={},
                         denovo_opts={}):


    if not taxonomy_fname:
        taxonomy_fname = settings.workflows.sixteen.otu_taxonomy
    if not ref_fasta:
        ref_fasta = settings.workflows.usearch.otu_db
    if not out_tsv:
        out_tsv = out_biom+".tsv"

    if denovo_otu_txt and non_chimeric_otu_seqs:
        denovo_otutab = denovo_otu_txt
        nonchimera = non_chimeric_otu_seqs
    else:
        denovo_otutab = in_fasta+".otus.txt"
        usearchfolder = denovo_otutab+"_usearch"
        nonchimera = join(usearchfolder, "nonchimeric.fa")
        denovo_opts.pop("remove_tempfiles", None)
        task_dict = next(iter(
            pick_denovo_otus(in_fasta, denovo_otutab,
                             remove_tempfiles=remove_tempfiles,
                             **denovo_opts)
        ))
        yield task_dict

    default_u_opts = dict([
        ("strand", "both"),
        ("id", "0.97"),
    ]+list(usearch_closed_opts.items()))

    tmpfolder = in_fasta+"_closedref"
    mkdir_cmd = "test -d "+tmpfolder+" || mkdir "+tmpfolder

    closed_out = join(tmpfolder, "closed.uc")
    usearch_cmd = ("usearch8 -usearch_global "+nonchimera+
                   " -db "+ref_fasta+
                   " -uc "+closed_out+
                   " -top_hit_only "+usearch_dict_flags(default_u_opts))

    def _fmt_otutab():
        import re
        from operator import add
        from collections import defaultdict

        def fields(fname, get_idxs=None):
            get = lambda item: item
            if get_idxs:
                get = lambda item: [item[idx] for idx in get_idxs]
            with open(fname) as f:
                for line in f:
                    line = line.strip().split('\t')
                    if line:
                        yield get(line)

        idx = dict(fields(taxonomy_fname))
        otu_rows = fields(denovo_otutab)
        otu_header = next(otu_rows)
        otu_idx = dict([ (row[0], row) for row in otu_rows ])

        default = lambda: [0 for _ in otu_header[1:]]
        output_dict = defaultdict(default)
        for query, target in fields(closed_out, get_idxs=(-2, -1)):
            otu_id = re.search("OTU_(\d+)", query).group(1)
            if otu_id not in otu_idx:
                continue
            abd = map(int, otu_idx[otu_id][1:])
            taxy = idx.get(target, "Unclassified")
            if taxy == "Unclassified":
                otu_id = "0"
            else:
                otu_id = target
            current = output_dict[(otu_id, taxy)]
            output_dict[(otu_id, taxy)] = map(add, current, abd)
                    
        with open(out_tsv, 'w') as out_f:
            print >> out_f, "\t".join(list(otu_header)+["taxonomy"])
            for (otu_id, taxy), abd in output_dict.iteritems():
                abd = map(str, abd)
                print >> out_f, "\t".join([otu_id]+abd+[taxy])

    biom_cmd = ("biom convert -i "+out_tsv+" -o "+out_biom+
                " --table-type='OTU Table' --process-obs-metadata=taxonomy"+
                " --output-metadata-id=taxonomy")

    if sample_metadata_fname:
        biom_cmd += " --sample-metadata-fp=" + sample_metadata_fname

    actions = [mkdir_cmd, usearch_cmd, _fmt_otutab, biom_cmd]
    if remove_tempfiles:
        cmd = "rm -rf "+tmpfolder
        actions.append(cmd)

    file_dep = [in_fasta, denovo_otutab, nonchimera, ref_fasta]
    yield { "name": "usearch_pick_otus_closed_ref: "+out_biom,
            "targets": [out_biom],
            "actions": actions,
            "file_dep": [in_fasta, denovo_otutab, nonchimera],
            "title"   : usearch_rusage(file_dep)}
    

