import os

from doit.exceptions import TaskError, TaskFailed

from anadama.action import CmdAction
from anadama.decorators import requires


@requires(binaries=['usearch7', 'sequence_pair'])
def usearch_stitch(input_fastq_pair, output_fastq, verbose=True, 
                   remove_tempfiles=True, **opts):
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
        ret = CmdAction(stitch_cmd.format(r1=r1, r2=r2,
                                   opts_str=opts_str,
                                   output_fastq=output_fastq),
                        verbose=verbose).execute()
        if type(ret) in (TaskError, TaskFailed):
            paired = lambda s: s.replace('.fastq', '_paired.fastq')
            CmdAction(pair_cmd.format(r1=r1, r2=r2,
                                      r1out=paired(r1),
                                      r2out=paired(r2)), 
                      verbose=verbose).execute()
            ret = CmdAction(stitch_cmd.format(r1=paired(r1),
                                              r2=paired(r2),
                                              opts_str=opts_str,
                                              output_fastq=output_fastq),
                            verbose=verbose).execute()
            if remove_tempfiles:
                os.remove(paired(r1))
                os.remove(paired(r2))
        if type(ret) in (TaskError, TaskFailed):
            ret = CmdAction('cp {r1} {out}'.format(r1=r1, out=output_fastq),
                            verbose=verbose).execute()
        else:
            return ret


    yield { "name"    : "usearch_stitch: "+output_fastq,
            "actions" : [run],
            "file_dep": input_fastq_pair,
            "targets" : [output_fastq] }


@requires(binaries=['usearch7'])
def usearch_filter(input_fastq, output_fasta, verbose=True, **opts):
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

