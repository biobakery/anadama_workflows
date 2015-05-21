import os
from os.path import join

from anadama import strategies
from anadama.action import CmdAction
from anadama.decorators import requires
from anadama.strategies import if_exists_run
from anadama.util import addtag

from . import settings
from .sixteen import assign_taxonomy

def usearch_dict_flags(opts_dict):
    opts = [ " -%s %s " %(key, val) 
             for key, val in opts_dict.iteritems() ]
    return " ".join(opts)


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
             "targets" : [output_fastq] }


@requires(binaries=['usearch7'],
          version_methods=["usearch7 -version"])
def filter(input_fastq, output_fasta, verbose=True, **opts):
    """Filter a fastq file, outputting sequences as fasta, using USEARCH
    version 7. The USEARCH binary should be named usearch7 in order
    for this workflow to operate.

    :param input_fastq: String; file name of a single fastq file to be filtered.
    :param output_fasta: String; name of resulting filtered fasta file
    :keyword verbose: Boolean; if true, print commands at runtime as they 
                      are executed. Defaults to true.
    :keyword **opts: Any additional keyword arguments are passed to usearch7 
                     as command line flags. By default, it passes 
                     `fastq_minlen=200` and `fastq_truncqual=25` as 
                     '-fastq_minlen 200' and '-fastq_truncqual 25', 
                     respectively.

    External dependencies
      - USEARCH v7 http://www.drive5.com/usearch

    """
    

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
            return CmdAction(cmd, verbose=verbose).execute()
        else:
            open(output_fasta, 'w').close()
        

    return { "name"     : "usearch_filter: "+output_fasta,
             "actions"  : [run],
             "file_dep" : [input_fastq],
             "targets"  : [output_fasta] }


@requires(binaries=["usearch8", "uclust_otutable"],
          version_methods=["usearch8 -version",
                           "pip freeze | fgrep anadama_workflows"])
def pick_denovo_otus(fasta_in, otutab_out, remove_tempfiles=True,
                     strand="plus", derep_opts={}, sort_opts={},
                     cluster_opts={}, chimera_opts={}, map_opts={}):

    tmpfolder = otutab_out+"_usearch"
    mkdir_cmd = "test -d "+tmpfolder+" || mkdir "+tmpfolder

    default_derep_opts = dict([
        ("sizeout", ""),
    ]+list(derep_opts.items()))
    
    derep_out = join(tmpfolder, "derep.fa")
    derep_cmd = ("usearch8 -derep_fulllength "+fasta_in+
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
                                    cluster_otus_out, chimera_out, map_out))

    actions = [ mkdir_cmd, derep_cmd, sort_cmd, cluster_cmd,
                chimera_cmd, map_cmd, otu_cmd ]
    if remove_tempfiles:
        actions.append(cleanup_cmd)

    return { "name"     : "usearch_pick_denovo_otus: "+otutab_out,
             "actions"  : actions,
             "file_dep" : [fasta_in],
             "targets"  : [otutab_out] }


@requires(binaries=["usearch8", "assign_taxonomy.py", "biom"])
def pick_otus_closed_ref(in_fasta, out_biom,
                         sample_metadata_fname=None,
                         remove_tempfiles=True, qiime_opts={}):

    denovo_otutab = in_fasta+".otus.txt"
    task_dict = next(iter(
        pick_denovo_otus(in_fasta, denovo_otutab,
                         remove_tempfiles=False)
    ))
    yield task_dict

    tmpfolder = in_fasta+"_otus"
    usearchfolder = denovo_otutab+"_usearch"
    nonchimera = join(usearchfolder, "nonchimeric.fa")
    task_dict = next(iter(
        assign_taxonomy(nonchimera, tmpfolder, qiime_opts=qiime_opts)
    ))
    yield task_dict

    tax_fname = task_dict["targets"][0]
    formatted_tax_fname = addtag(tax_fname, "fmt")
    def _format_tax():
        import re
        with open(tax_fname) as f_in, open(formatted_tax_fname, 'w') as f_out:
            print >> f_out, "\t".join(("#OTUID", "taxonomy"))
            for line in f_in:
                fields = line.split('\t')
                otu_id = re.search(r'OTU_(\d+)', fields[0]).group(1)
                print >> f_out, "\t".join((otu_id, fields[1]))
                
    tmpbiom = join(tmpfolder, "tmp.biom")
    biom_cmd_a = ("biom convert -i "+denovo_otutab+
                  " -o "+tmpbiom+" --table-type='OTU Table'")
    biom_cmd_b = ("biom add-metadata -i "+tmpbiom+
                  " --observation-metadata="+formatted_tax_fname+
                  " -o "+out_biom)

    if sample_metadata_fname:
        biom_cmd_b += " --sample-metadata-fp=" + sample_metadata_fname

    actions = [ _format_tax, biom_cmd_a, biom_cmd_b ]
    if remove_tempfiles:
        cmd = "rm -rf "+tmpfolder+" "+usearchfolder+" "+denovo_otutab
        actions.append(cmd)

    yield { "name": "usearch_pick_otus_closed_ref: "+out_biom,
            "targets": [out_biom],
            "actions": actions,
            "file_dep": [in_fasta] }
    
