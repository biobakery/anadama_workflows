"""General purpose workflows"""

import os
import mimetypes

from anadama.util import addext, guess_seq_filetype, new_file
from anadama.util import dict_to_cmd_opts
from anadama.decorators import requires

from . import ( 
    starters
)

@requires(binaries=['gzip', 'bzip2'],
          version_methods=["gzip --version | head -1",
                           "bzip2 --version < /dev/null 2>&1 | head -1"])
def extract(fname_from, fname_to=None):
    """Workflow for converting an input file from their zipped to
    their unzipped equivalent.

    :param files_list: String; The input files to decompress

    :keyword fname_to: String; optional name of the resulting decompressed 
                       file. Defaults to the original file name, but with
                       the outermost file extension removed.

    External dependencies:
      - gunzip: should come with gzip
      - bunzip2: Should come with the bzip2 package

    """

    target = fname_to if fname_to else os.path.splitext(fname_from)[0]
    task = {
        "name": "decompress:"+fname_from,
        "targets": [target],
        "file_dep": [fname_from]
    }

    _, min_file_type = mimetypes.guess_type(fname_from)
    if min_file_type == 'gzip':
        task['actions'] = [ "gzip -d <"+fname_from+" > "+target ]
        return task
    elif min_file_type == 'bzip2':
        task['actions'] = [ "gzip -d <"+fname_from+" > "+target ]
        return task
    else:
        return None


@requires(binaries=['fastq_split'],
          version_methods=["pip freeze | grep anadama_workflows"])
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

    seqtype = from_format if from_format else guess_seq_filetype(files_list[0])

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


@requires(binaries=['sequence_convert'],
          version_methods=["pip freeze | grep anadama_workflows"])
def sequence_convert(files_list, output_file=None,
                     reverse_complement=False, from_format=None,
                     format_to="fastq", lenfilters_list=list()):
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

    if type(lenfilters_list) is not list:
        lenfilters_list = list(lenfilters_list)

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


@requires(binaries=['fastq-join'],
          version_methods=["md5sum `which fastq-join` | awk '{print $1;}'"])
def fastq_join(forward_fname, reverse_fname, output_file,
               reorder_to=None, options=dict()):
    """Workflow function for joining (aka stitching) paired-end fastq
    files with ea-utils' ``fastq-join``. If the ``drop_unpaired``
    option is set to True, unpaired forward reads are concatenated to
    the joined fastq file.

    :param forward_fname: String; file name for the forward reads.

    :param reverse_fname: String; file name for the reverse reads.

    :param output_file: String; file name for the the finished, joined reads.
   
    :param reorder_to: String; file name to reorder sequences against

    :param options: Dictionary; interpreted as command line options to
    be passed to the wrapped knead_data.py script.  No - or -- flags
    are necessary; the correct - or --t flags are inferred based on
    the length of the option.  For boolean options, use the key/value
    pattern of { "my-option": "" }.

    External Dependencies
      - fastq-join: Part of ea-utils 1.1.2-806 'https://drive.google.com/folderview?id=0B7KhouP0YeRAOTFWWGVFYkFSQjg&usp=sharing'

    """


    drop_unpaired = options.pop('drop_unpaired', False)

    default_opts = {
        "o": output_file
    }
    default_opts.update(options)
    opts = dict_to_cmd_opts(default_opts)

    cmd = "fastq-join "+opts+ " "+forward_fname+" "+reverse_fname

    if '%' in output_file:
        renamed_output = output_file.replace("%", "join")
    else:
        renamed_output = output_file+"join"

    actions = [cmd]
    if not drop_unpaired and reorder_to:
        unpaired_forward = renamed_output.replace("join", "un1")
        actions.append( "sequence_re-pair -f fastq -t fastq -b %s %s %s > %s"%(
            reorder_to, renamed_output, unpaired_forward, output_file) )
    else:
        actions.append("mv %s %s"%(renamed_output, output_file))


    return {
        "name": "fastq_join: %s..."%(output_file),
        "actions": actions,
        "file_dep": [forward_fname, reverse_fname],
        "targets": [output_file]
    }


def sequence_pair(seqfname1, seqfname2, 
                  outfname1="/dev/null", outfname2="/dev/null", 
                  from_format=None, format_to="fastq", 
                  options=dict()):

    extra_options = dict_to_cmd_opts(options)
    targets = [ target for target in [outfname1, outfname2] 
                if target != "/dev/null" ]

    if not from_format:
        from_format = guess_seq_filetype(targets[0])

    pair_cmd = ("sequence_pair"
                " -f {from_format} -t {format_to}"
                " -1 {r1out} -2 {r2out} "
    ).format(
        from_format=from_format, format_to=format_to,
        r1out=outfname1, r2out=outfname2
    )
    pair_cmd += extra_options + " {} {}".format(seqfname1, seqfname2)
    
    return {
        "name": "sequence_pair: %s %s"%(outfname1, outfname2),
        "actions" : [pair_cmd],
        "file_dep": [seqfname1, seqfname2],
        "targets": targets
    }



def cat(input_files, output_file):
    """Just concatenates files with the ``cat`` command"""

    return {
        "name": "cat: "+output_file,
        "actions": ["cat %s > %s" %(" ".join(input_files), output_file)],
        "file_dep": input_files,
        "targets": [output_file]
    }
