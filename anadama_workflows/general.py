"""General purpose workflows"""

import os
import mimetypes

from anadama.util import addext, guess_seq_filetype, new_file

from . import ( 
    starters
)

def extract(files_list):
    actions = list()
    targets = list()
    for fname in files_list:
        _, min_file_type = mimetypes.guess_type(fname)
        if min_file_type == 'gzip':
            actions.append( "gunzip "+fname )
            targets.append( os.path.splittext(fname)[0] )
        if min_file_type == 'bzip2':
            actions.append( "bunzip2 "+fname )
            targets.append( os.path.splittext(fname)[0] )
        else:
            pass

    return {
        "name": "decompress:"+infiles_list[0],
        "targets": targets,
        "actions": actions,
        "file_dep": infiles_list
    }


def fastq_split(files_list, fasta_fname, qual_fname,
                reverse_complement=False, trim=4, from_format=None):
    if not from_format:
        seqtype = guess_seq_filetype(files_list[0])
    else:
        seqype = from_format
    cmd = ("fastq_split"+
           " --fasta_out="+fasta_fname+
           " --qual_out="+qual_fname+
           " --format="+seqtype+
           " --trim="+str(trim))

    if reverse_complement:
        cmd += " -r"

    cmd += " "+" ".join(files_list)

    return {
        "name": "fastq_split:"+files_list[0],
        "actions": [cmd],
        "file_dep": files_list,
        "targets": [fasta_fname, qual_fname]
    }


def sequence_convert(files_list, output_file=None, format_to="fastq"):
    if not output_file:
        output_file = files_list[0] + "_merged."+format_to

    if not from_format:
        from_format = guess_seq_filetype(files_list[0])

    cmd = ("sequence_convert"
           + " --format="+from_format
           + " --to="+format_to )

    if reverse_complement:
        cmd += " --reverse_complement"

    cmd += ( " "+" ".join(files_list)
             + " > "+output_file)

    return {
        "name": "sequence_convert_to_%s: %s..."%(format_to, files_list[0]),
        "actions": [cmd],
        "file_dep": files_list,
        "targets": [output_file]
    }


###
# Example workflow function for returning multiple tasks
# 
# def myworkflow(somefiles):
#     stuff = _magic()

#     for item in stuff:
#         yield {
#             "name": item.name,
#             ...
#         }
