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


def fastq_split(files_list, name, dir, reverse_complement=False, trim=4):
    fasta_fname = new_file(addext(name, "fa"), basedir=dir)
    qual_fname = new_file(addext(name, "qual"), basedir=dir)

    seqtype = guess_seq_filetype(files_list[0])
    cmd = ("fastq_split"+
           " --fasta_out="+fasta_fname+
           " --qual_out="+qual_fname+
           " --format="+seqtype+
           " --trim="+str(trim))

    if reverse_complement:
        cmd += " -r"

    cmd += " "+" ".join(files_list)

    return {
        "name": "fastq_split:"+name,
        "actions": [cmd],
        "file_dep": files_list,
        "targets": [fasta_fname, qual_fname]
    }


def sequence_convert(files_list, output_file=None, format_to="fastq"):
    if not output_file:
        output_file = files_list[0] + "_merged."+format_to

    cmd = ("sequence_convert"
           + " --format="+guess_seq_filetype(files_list[0])
           + " --to="+format_to
           + " "+" ".join(files_list)
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
