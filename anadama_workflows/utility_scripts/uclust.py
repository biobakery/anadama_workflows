import re
import os
import sys
import subprocess
import tempfile
import optparse
from operator import add
from itertools import ifilterfalse
from itertools import count
from collections import defaultdict
from collections import namedtuple
from os.path import join

from ..usearch import util
from ..usearch import usearch_dict_flags as dict_flags


class parse_otutable(object):
    """Goes from a usearch mapping results file to an OTU table, 
    -> samples \/ OTUs 

    Assumes that search query sequences are from qiime-formatted
    demultiplexed sequence files. That is, the sequence label is composed
    of the sampleid along with other items, separated by underscores '_'.
    The first item before the underscore is assumed to be the sample id.

    """

    @staticmethod
    def parsetarget(target_str):
        return re.search(r'OTU_(\d+)', target_str).group(1)

    @staticmethod
    def fields(uc_fname):
        with open(uc_fname) as f:
            for line in f:
                if line.startswith("H"):
                    fields = line.split('\t')
                    query, target = fields[8], fields[9]
                    sample_id = query.split("_", 1)[0]
                    otu_id = parse_otutable.parsetarget(target)
                    yield sample_id, otu_id

    @staticmethod
    def get(d, items, default=0):
        return [d.get(item, default) for item in items]

    @staticmethod
    def output(d, sample_ids):
        sample_ids = list(sample_ids)
        print "\t".join(["OTUId"]+list(sample_ids))
        for otu_id, hits_dict in d.iteritems():
            hits = parse_otutable.get(hits_dict, sample_ids, default=0)
            print "\t".join([otu_id]+map(str, hits))


    @staticmethod
    def main(uc_fname=None):
        if not uc_fname:
            uc_fname = sys.argv[1]
        table_dict = defaultdict(lambda:defaultdict(int))
        sample_ids = set()
        for sample_id, otu_id in parse_otutable.fields(uc_fname):
            sample_ids.add(sample_id)
            table_dict[otu_id][sample_id] += 1

        parse_otutable.output(table_dict, sample_ids)

def parse_otu_table():
    ret = parse_otutable.main(*sys.argv[1:])
    sys.exit(ret)

class ShellException(IOError):
    pass

def sh(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True,
       **kwargs):
    proc = subprocess.Popen(cmd, stdout=stdout, stderr=stderr,
                            shell=shell,**kwargs)
    ret = proc.communicate()
    if proc.returncode:
        raise ShellException("Command `{}' failed. \nOut: {}\nErr: {}".format(
            cmd, ret[0], ret[1]))
    return ret


Step = namedtuple("Step", "idx cmd targets needs note")
class ExecutionPlan(object):
    def __init__(self, resume=False, quiet=False,
                 cmd_stdout=sys.stderr, cmd_stderr=sys.stderr,
                 report_cmd=False,
                 report_f = lambda s: sys.stderr.write(s)):
        self.resume = resume
        self.stdout = cmd_stdout
        self.stderr = cmd_stderr
        self.quiet = quiet
        self._report = report_f
        self._report_cmd = report_cmd
        self.steps = list()
        self._cntr = count(0)

    def step(self, cmd, targets, needs, note=""):
        self.steps.append(Step(next(self._cntr), cmd, targets, needs, note))

    def _pct(self, step):
        return (float(step.idx)/len(self.steps))*100

    def _msg(self, step, msg):
        msg =  "[step #{i:02d} - {pct:.2f}% ] {msg} -- {note}.".format(
            i=step.idx+1, pct=self._pct(step), msg=msg, note=step.note
        )
        if self._report_cmd:
            msg += " "+step.cmd if type(step.cmd) is str else str(step.cmd)
        self._report(msg+"\n")
        return msg

    def _skip(self, step):
        if not self.quiet:
            self._msg(step, "Skipped")

    def _quit(self, step, override=False):
        if not self.quiet and not override:
            return self._msg(step, "Error! Quitting.")

    def _runfunc(self, step):
        try:
            step.cmd()
        except Exception as e:
            self._quit(step._replace(note=e.message), override=True)
            return False
        else:
            return True
        
    def _runsh(self, step):
        try:
            if self.quiet:
                with open(os.devnull, 'w') as null_f:
                    sh(step.cmd, stdout=null_f, stderr=null_f)
            else:
                sh(step.cmd, stdout=self.stdout, stderr=self.stderr)
        except ShellException as e:
            self._quit(step._replace(note=e.message), override=True)
            return False
        else:
            return True

    def run(self, step):
        if not self.quiet:
            self._msg(step, "Running")
        if type(step.cmd) is str:
            return self._runsh(step)
        else:
            return self._runfunc(step)

    def go(self):
        for step in self.steps:
            missing = list(ifilterfalse(os.path.exists, step.needs))
            if missing:
                msg = self._quit(step)
                raise IOError(msg)
            elif all(os.path.exists(t) for t in step.targets) and self.resume:
                self._skip(step)
                continue
            self.run(step)
            


def pick_denovo_otus(execution_plan, fasta_in, otutab_out,
                     chimera_gold_standard_fname, tmp_folder=None,
                     remove_tempfiles=True, strand="plus",
                     truncate_opts={}, derep_opts={}, sort_opts={},
                     cluster_opts={}, chimera_opts={}, map_opts={}):

    if not tmp_folder:
        tmp_folder = tempfile.mkdtemp(prefix="./"+otutab_out)
    else:
        if not os.path.exists(tmp_folder):
            os.mkdir(tmp_folder)

    plan = execution_plan

    # truncate 
    trunc_out = join(tmp_folder, "truncated.fa")
    if not truncate_opts.get("trunclen"):
        truncate_opts['trunclen'] = util.cutoff(fasta_in)
    truncate_opts['fastaout'] = trunc_out
    cmd = "usearch8 -fastx_truncate "+fasta_in+" "+dict_flags(truncate_opts)
    plan.step(cmd, [trunc_out], [fasta_in],
              note="Truncating to uniform length")

    # dereplicate
    default_derep_opts = dict([
        ("sizeout", ""),
    ]+list(derep_opts.items()))
    
    derep_out = join(tmp_folder, "derep.fa")
    derep_cmd = ("usearch8 -derep_fulllength "+trunc_out+
                 " -fastaout "+derep_out+
                 " "+dict_flags(default_derep_opts))
    plan.step(derep_cmd, [derep_out], [trunc_out],
              note="Dereplicating reads")

    # sort
    default_sort_opts = dict([
        ("minsize","2"),
    ]+list(sort_opts.items()))

    sort_out = join(tmp_folder, "sorted.fa")
    sort_cmd = ("usearch8 -sortbysize "+derep_out+" -fastaout "+sort_out+
                " "+dict_flags(default_sort_opts))
    plan.step(sort_cmd, [sort_out], [derep_out], note="Sorting reads")

    # cluster
    default_cluster_opts = dict([
        ("relabel", "OTU_"),
        ("sizein", ""),
        ("sizeout", ""),
    ]+list(cluster_opts.items()))
                
    cluster_otus_out = join(tmp_folder, "otus.fa")
    cluster_otus_log = join(tmp_folder, "cluster_results.txt")
    cluster_cmd = ("usearch8 -cluster_otus "+sort_out+
                   " -otus "+cluster_otus_out+
                   " -uparseout "+cluster_otus_log+
                   " "+dict_flags(default_cluster_opts))
    plan.step(cluster_cmd, [cluster_otus_out, cluster_otus_log],
              [sort_out], note="Cluster OTUs")

    # chimera filter
    default_chimera_opts = dict([
        ("db", chimera_gold_standard_fname),
        ("strand", strand),
    ]+list(chimera_opts.items()))

    chimera_out = join(tmp_folder, "nonchimeric.fa")
    chimera_cmd = ("usearch8 -uchime_ref "+cluster_otus_out+
                   " -nonchimeras "+chimera_out+
                   " "+dict_flags(default_chimera_opts))
    plan.step(chimera_cmd, [chimera_out],
              [default_chimera_opts['db'], cluster_otus_out],
              note="Remove chimeric sequences")

    # count otus
    default_map_opts = dict([
        ("id", "0.97"),
        ("strand", strand),
    ]+list(map_opts.items()))

    map_out = join(tmp_folder, "mapping_results.uc")
    map_cmd = ("usearch8 -usearch_global "+fasta_in+
               " -db "+chimera_out+
               " -uc "+map_out+
               " "+dict_flags(default_map_opts))
    plan.step(map_cmd, [map_out], [chimera_out], note="Count OTUs")

    # format the table
    otu_cmd = ("uclust_otutable %s > %s"%(map_out, otutab_out))
    plan.step(otu_cmd, [otutab_out], [map_out], note="Tabulate OTU counts")

    if remove_tempfiles:
        to_rm = [derep_out, sort_out, cluster_otus_out, map_out]
        cleanup_cmd = "rm " + " ".join(to_rm)
        plan.step(cleanup_cmd, [], to_rm, "Remove temporary files")
        
    return plan



def format_otu_table(taxonomy_fname, denovo_otutab, closed_out, out_tsv):
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
    

def pick_otus_closed_ref(execution_plan, in_fasta, out_tsv,
                         taxonomy_fname, ref_fasta, chimera_std,
                         non_chimeric_otu_seqs=None,
                         denovo_otu_txt=None,
                         remove_tempfiles=True, tmp_folder=None,
                         usearch_closed_opts={}, denovo_opts={}):

    if denovo_otu_txt and non_chimeric_otu_seqs:
        denovo_otutab = denovo_otu_txt
        nonchimera = non_chimeric_otu_seqs
        plan = execution_plan
        if not tmp_folder:
            tmpfolder = tempfile.mkdtemp(prefix="./"+tmp_folder)
        else:
            tmpfolder = tmp_folder
            if not os.path.exists(tmp_folder):
                os.mkdir(tmp_folder)
    else:
        denovo_otutab = in_fasta+".otus.txt"
        nonchimera = join(tmp_folder, "nonchimeric.fa")
        plan = pick_denovo_otus(execution_plan, in_fasta,
                                denovo_otutab, chimera_std,
                                tmp_folder=tmp_folder,
                                remove_tempfiles=False, **denovo_opts)
        tmpfolder = tmp_folder

    default_u_opts = dict([
        ("strand", "both"),
        ("id", "0.97"),
    ]+list(usearch_closed_opts.items()))
    closed_out = join(tmpfolder, "closed.uc")
    usearch_cmd = ("usearch8 -usearch_global "+nonchimera+
                   " -db "+ref_fasta+
                   " -uc "+closed_out+
                   " -top_hit_only "+dict_flags(default_u_opts))
    plan.step(usearch_cmd, [closed_out], [nonchimera, ref_fasta],
              note="Find OTU names")

    format_f = lambda: format_otu_table(taxonomy_fname, denovo_otutab,
                                        closed_out, out_tsv)
    plan.step(format_f, [out_tsv],
              [taxonomy_fname, denovo_otutab, closed_out],
              note="Tabulate named OTUs")

    if remove_tempfiles:
        plan.step("rm -rf "+tmpfolder, [], [],
                  note="Removing temporary files")

    return plan


global_opts = [
    optparse.make_option("-i", '--input', action="store", type="string",
                         help="16S sequences input (fasta format)"),
    optparse.make_option("-o", '--output', action="store", type="string",
                         default=None, help="OTU Table output (tsv format)."),
    optparse.make_option("-d", "--tmp_dir", default=None, dest="tmp_folder",
                         help="Possibly useful output files are stored here"),
    optparse.make_option("--keep_tempfiles", default=True,
                         dest="remove_tempfiles",
                         action="store_false", help="Removes the tmp_dir"),
    optparse.make_option("--resume", default=False, action="store_true",
                         help="Resume from intermediate files"),
    optparse.make_option('-q', "--quiet", default=False, action="store_true",
                         help="Shush!"),
    optparse.make_option('--print_cmd', action="store_true", default=False,
                         help="Print what commands I run"),
    optparse.make_option('-l', '--log_file', help="Write logs to this file"),
    optparse.make_option("-c", '--chimera_standard', action="store",
                         type="string", 
                         help=("Fasta DNA sequence file that contains "
                               "gold standard chimera-free sequences")),
    optparse.make_option('--strand', default="plus", type="choice",
                         choices=["plus", "both"],
                         help=("What strand should I search for"
                               " chimeric sequences? 'plus' or 'both'")),
]

kv_opts = dict(
    truncate_opts = optparse.make_option(
        "--truncate_opts",
        help=("key=value options passed to underlying usearch binary")),
    derep_opts = optparse.make_option(
        "--derep_opts",
        help=("key=value options passed to underlying usearch binary")),
    sort_opts = optparse.make_option(
        "--sort_opts",
        help=("key=value options passed to underlying usearch binary")),
    cluster_opts = optparse.make_option(
        "--cluster_opts",
        help=("key=value options passed to underlying usearch binary")),
    chimera_opts = optparse.make_option(
        "--chimera_opts",
        help=("key=value options passed to underlying usearch binary")),
    map_opts = optparse.make_option(
        "--map_opts",
        help=("key=value options passed to underlying usearch binary")),
)

def parse_kv_opts(opts):
    for key in kv_opts.iterkeys():
        s = getattr(opts, key, "") 
        if not s:
            yield key, {}
        else:
            yield key, dict([pair.split("=", 1) for pair in s.split()])
        

def denovo_cli():
    HELP = "Denovo OTUS from 16S sequences in a fasta file"
    required_opts = ("input", "output", "chimera_standard")
    options = global_opts + kv_opts.values()

    parser = optparse.OptionParser(option_list=options, usage=HELP)
    opts, _ = parser.parse_args()

    reqd = ( (v, getattr(opts, v, False)) for v in required_opts )
    missing_reqd = [ name for name, val in reqd if not val]
    if missing_reqd:
        msg = "Missing required options: "+", ".join(missing_reqd)
        print >> sys.stderr, msg
        parser.print_usage()
        sys.exit(1)
    
    kv_parsed = dict(parse_kv_opts(opts))

    if opts.log_file:
        log_f = open(opts.log_file, 'w')
        plan = ExecutionPlan(resume=opts.resume, quiet=False,
                             cmd_stdout=log_f, cmd_stderr=log_f,
                             report_cmd=opts.print_cmd, report_f=lambda s:
                             log_f.write(s+"\n") )
    else:
        plan = ExecutionPlan(resume=opts.resume, quiet=opts.quiet,
                             report_cmd=opts.print_cmd)

    plan = pick_denovo_otus(plan, opts.input, opts.output,
                            opts.chimera_standard,
                            tmp_folder=opts.tmp_folder,
                            remove_tempfiles=opts.remove_tempfiles,
                            strand=opts.strand, **kv_parsed)
    try:
        plan.go()
    except (ShellException, IOError) as e:
        print >> sys.stderr, str(e)
        sys.exit(1)


def closed_cli():
    HELP = "Denovo OTUS from 16S sequences in a fasta file"
    required_opts = ("input", "output", "taxonomy",
                     "reference", "chimera_standard")
    options = global_opts + [
        optparse.make_option("-t", '--taxonomy',
                             help=("Tab separated file connecting reference "
                                   "sequence IDs (left column) and names "
                                   "(right column).")),
        optparse.make_option("-r", '--reference',
                             help=("OTUs are compared to this fasta or "
                                   "usearch-indexed sequence set for "
                                   "identity.")),
        optparse.make_option("--denovo_otu_table",
                             help=("assign names to taxonomic features "
                                   "in this precomputed otu table")),
        optparse.make_option("--otu_sequences",
                             help=("Use these sequences as queries when "
                                   "searching for names of taxonomic "
                                   "features")),
        optparse.make_option(
            "--closed_opts", default="",
            help=("key=value options passed to underlying usearch binary")),
        
    ] + kv_opts.values()

    parser = optparse.OptionParser(option_list=options, usage=HELP)
    opts, _ = parser.parse_args()

    reqd = ( (v, getattr(opts, v, False)) for v in required_opts )
    missing_reqd = [ name for name, val in reqd if not val]
    if missing_reqd:
        msg = "Missing required options: "+", ".join(missing_reqd)
        print >> sys.stderr, msg
        parser.print_usage()
        sys.exit(1)
    
    denovo_opts = dict(parse_kv_opts(opts))
    denovo_opts['strand'] = opts.strand
    closed_opts = dict([pair.split("=", 1)
                        for pair in opts.closed_opts.split()])

    if opts.log_file:
        log_f = open(opts.log_file, 'w')
        plan = ExecutionPlan(resume=opts.resume, quiet=False,
                             cmd_stdout=log_f, cmd_stderr=log_f,
                             report_cmd=opts.print_cmd, report_f=lambda s:
                             log_f.write(s+"\n") )
    else:
        plan = ExecutionPlan(resume=opts.resume, quiet=opts.quiet,
                             report_cmd=opts.print_cmd)

    plan = pick_otus_closed_ref(plan, opts.input, opts.output,
                                opts.taxonomy, opts.reference,
                                opts.chimera_standard,
                                tmp_folder=opts.tmp_folder,
                                remove_tempfiles=opts.remove_tempfiles,
                                denovo_opts=denovo_opts,
                                usearch_closed_opts=closed_opts)
    try:
        plan.go()
    except (ShellException, IOError) as e:
        print >> sys.stderr, str(e)
        sys.exit(1)

