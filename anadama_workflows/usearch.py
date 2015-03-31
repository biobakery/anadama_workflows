import os

from doit.action import CmdAction

from anadama import strategies
from anadama.decorators import requires


@requires(binaries=['usearch7', 'sequence_pair'],
          version_methods=["usearch7 -version"
                           "pip freeze | grep anadama_workflows"])
def usearch_stitch(input_fastq_pair, output_fastq, verbose=True, 
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
    
    opts = [ " -%s %s " %(key, val) 
             for key, val in default_options.iteritems() ]
    opts_str = " ".join(opts)

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


    yield { "name"    : "usearch_stitch: "+output_fastq,
            "actions" : [run],
            "file_dep": input_fastq_pair,
            "targets" : [output_fastq] }


@requires(binaries=['usearch7'],
          version_methods=["usearch7 -version"])
def usearch_filter(input_fastq, output_fasta, verbose=True, **opts):
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
    
    opts = [ " -%s %s " %(key, val) 
             for key, val in default_options.iteritems() ]
    cmd += " ".join(opts)

    def run():
        if os.stat(input_fastq).st_size > 1:
            return CmdAction(cmd, verbose=verbose).execute()
        else:
            open(output_fasta, 'w').close()
        

    yield { "name"     : "usearch_filter: "+output_fasta,
            "actions"  : [run],
            "file_dep" : [input_fastq],
            "targets"  : [output_fasta] }

