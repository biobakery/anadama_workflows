import os

from anadama.util import dict_to_cmd_opts
from anadama.decorators import requires
from anadama.action import CmdAction

from anadama_workflows import settings

@requires(binaries=['subread-align'],
          version_methods=["subread-align -v 2>&1 | awk '/Sub/{print $2;}'"])
def align(maybe_paired_fastq, output_sam, options=dict()):
    opts = {
        "unique": "",
        "hamming": "",
        "index": settings.workflows.subread.index
    }
    opts.update(options)
    opts['output'] = output_sam
    if type(maybe_paired_fastq) in (tuple, list):
        opts['read'], opts['read2'] = maybe_paired_fastq
        deps = maybe_paired_fastq
    else:
        opts['read'] = maybe_paired_fastq
        deps = [maybe_paired_fastq]

    cmd = "subread-align "+dict_to_cmd_opts(opts)

    def run():
        if any(os.stat(f).st_size < 1 for f in deps):
            open(output_sam, 'w').close()
        else:
            return CmdAction(cmd, verbose=True).execute()

    return { "name": "subread_align: "+output_sam,
             "actions": [run],
             "file_dep": deps,
             "targets": [output_sam] }


@requires(binaries=['featureCounts'],
          version_methods=["featureCounts -v 2>&1 | awk '/fea/{print $2;}'"])
def featureCounts(input_sams, output_table, options=dict()):
    opts = {
        "a": settings.workflows.subread.annotations,
    }
    opts.update(options)
    opts['o'] = output_table

    cmd = ("featureCounts"
           +" "+dict_to_cmd_opts(opts)
           +" ")
    
    def run():
        files = [f for f in input_sams
                 if os.path.exists(f) and
                 os.stat(f).st_size > 0 ]
        if files:
            return CmdAction(cmd+" ".join(files), verbose=True).execute()
        else:
            open(output_table, 'w').close()


    return { "name": "featureCounts: "+output_table,
             "file_dep": input_sams,
             "targets": [output_table],
             "actions": [run] }

