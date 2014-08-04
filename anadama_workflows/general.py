"""General purpose workflows"""

import os
import mimetypes

from anadama.util import addext, guess_seq_filetype, new_file
from anadama.decorators import requires

from . import ( 
    starters
)

@requires(binaries=['gzip', 'bzip2', 'gunzip', 'bunzip2'])
def extract(files_list):
    """Workflow for converting a list of input files from their zipped to
    their unzipped equivalent.

    :param files_list: List; The input files to decompress

    External dependencies:
      - gunzip: should come with gzip
      - bunzip2: Should come with the bzip2 package

    """
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


@requires(binaries=['fastq_split'])
def fastq_split(files_list, fasta_fname, qual_fname,
                reverse_complement=False, trim=4, from_format=None):
    """ Workflow for concatenating and converting a list of sequence files
    into a fasta file and a qual file. 

    :param files_list: List; List of input files
    :param fasta_fname: String; File name for output fasta file
    :param qual_fname: String; File name for output qual file
    :keyword reverse_complement: Boolean; Set to True if the resulting 
                                 sequence files should be the reverse 
                                 complement of the input sequences
    :keyword trim: Integer; trim these number of sequence items from the 
                   start of the sequence 
    :keyword from_format: String; biopython-recognized string to convert
                                  the sequence from. If not specified, we 
                                  guess with ``guess_seq_filetype``

    External dependencies:
      - fastq_split: python script that should come pre-installed with
        the anadama_workflows module

    """

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
        "name": "fastq_split:"+fasta_fname,
        "actions": [cmd],
        "file_dep": files_list,
        "targets": [fasta_fname, qual_fname]
    }

@requires(binaries=['sequence_convert'])
def sequence_convert(files_list, output_file=None, reverse_complement=False,
                     from_format=None, format_to="fastq", lenfilters_list=list()):
    """ Workflow for converting between sequence file formats.

    :param files_list: List; List of input files
    :param output_file: String; File name for output file
    :keyword reverse_complement: Boolean; Set to True if the resulting 
                                 sequence file should be the reverse 
                                 complement of the input sequences
    :keyword from_format: String; biopython-recognized string to convert
                                  the sequence from. If not specified, we 
                                  guess with ``guess_seq_filetype``
    :keyword format_to: String; output file format as recognized by 
                        biopython
    :keyword lenfilters_list: List of strings; conditions for filtering 
                              sequences by length.  To keep all sequences 
                              longer than 60 chars, for example, use >60.

    External dependencies:
    - sequence_convert: python script that should come pre-installed with
      the anadama_workflows module
    
    """

    if not output_file:
        output_file = files_list[0] + "_merged."+format_to

    if not from_format:
        from_format = guess_seq_filetype(files_list[0])

    cmd = ("sequence_convert"
           + " --format="+from_format
           + " --to="+format_to
           + " ".join([" -n '%s'"%(s) for s in lenfilters_list]) )

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
