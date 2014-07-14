"""16S workflows"""

import os
import operator
import itertools
from os.path import join

# should probably import from anadama directly?
from anadama.action import CmdAction
from doit.exceptions import TaskError, TaskFailed

from anadama.util import (
    addext, 
    addtag,
    new_file,
    dict_to_cmd_opts, 
    guess_seq_filetype
)

from . import ( 
    settings
)

def write_map(sample_group, sample_dir):
    map_fname = new_file("map.txt", basedir=sample_dir)

    def _write(targets):
        with open(map_fname, 'w') as map_file:
            # print the headers first
            print >> map_file, "#"+"\t".join(sample_group[0]._fields)

            for _, samples_bycode in itertools.groupby(
                    sample_group, operator.attrgetter("BarcodeSequence")):
                # get the first (hopefully only) sample from the samples
                # grouped by ID then barcode. Ignore any other samples
                # under the same ID for the same barcode
                sample = samples_bycode.next()
                bcode = sample.BarcodeSequence
                # uniq-ify to make qiime happy
                sample = list(sample)
                sample[0] += "_" + bcode
                print >> map_file, "\t".join(sample)

    return {
        "name": "write_map:"+map_fname,
        "actions": [_write],
        "targets": [map_fname]
    }


def demultiplex(map_fname, fasta_fname, qual_fname, output_fname,
                qiime_opts={}):
    output_dir=os.path.dirname(output_fname)
    opts = dict_to_cmd_opts(qiime_opts)
    
    cmd = ("qiime_cmd split_libraries.py"+
           " --map="+map_fname+
           " --fasta="+fasta_fname+
           " --qual="+qual_fname+
           " --dir-prefix="+output_dir+
           " "+opts)

    return {
        "name": "demultiplex:"+fasta_fname,
        "actions": [cmd],
        "file_dep": [map_fname, fasta_fname, qual_fname],
        "targets": [output_fname]
    }


def pick_otus_closed_ref(input_fname, output_dir, verbose=False, qiime_opts={}):
    output_fname = new_file("otu_table.biom", basedir=output_dir)
    revcomp_fname = new_file(
        "revcomp.fna", basedir=os.path.dirname(input_fname))

    default_opts = {
        "taxonomy_fp": settings.workflows.sixteen.otu_taxonomy,
        "reference_fp": settings.workflows.sixteen.otu_refseq
    }
    default_opts.update(qiime_opts)
    opts = dict_to_cmd_opts(default_opts)

    cmd = ("qiime_cmd pick_closed_reference_otus.py"+
           " --input_fp={}"+
           " --output_dir="+output_dir+
           " -f"+
           " "+opts)

    revcomp_cmd = ("sequence_convert"+
                   " --format=fasta"+
                   " --to=fasta "+
                   " -r"+
                   " "+input_fname+
                   " > "+revcomp_fname)

    def run(targets):
        ret = CmdAction(cmd.format(input_fname), 
                        verbose=verbose).execute()
        if type(ret) in (TaskError, TaskFailed):
            CmdAction(revcomp_cmd).execute()
            return CmdAction(cmd.format(revcomp_fname), 
                             verbose=verbose).execute()
        else:
            return ret
        

    return {
        "name": "pick_otus_closed_ref:"+input_fname,
        "actions": [run],
        "targets": [output_fname],
        "file_dep": [input_fname]
    }


def merge_otu_tables(files_list, name, output_dir):
    output_file = new_file(name, basedir=output_dir)
    cmd = "qiime_cmd merge_otu_tables.py -i {filenames} -o {output}"
    cmd = cmd.format( filenames = ",".join(files_list), 
                      output    = output_file  )
    return {
        "name": "merge_otu_tables: "+os.path.basename(name),
        "actions": [cmd],
        "targets": [output_file],
        "file_dep": files_list
    }


def picrust(file, **opts):
    """Workflow to predict metagenome functional content from 16S OTU tables.

	Optional options:
	tab_in		use 1 if the input is a tabulated file (default:0)
	tab_out		use 1 if the output file is to be tabulated (default:0)
	gg_version	provide a string containing the greengenes version to be used (default:most recent version)
	t		option to use a different type of prediction (default:KO)
	with_confidence	use 1 to output confidence intervals (default:0)
	custom		specify a file containing a custom trait to predict metagenomes
    """
    norm_out = new_file(addtag(file, "normalized_otus"))
    predict_out = new_file(addtag(file, "picrust"))

    all_opts = { 'tab_in' : 0,#flag if input is a tabulated file
                 'tab_out' : 0, #flag if output file is to be tabulated
                 'gg_version' : '', #greengenes version
                 't' : '', #empty for KO, else can be COG or RFAM
                 'with_confidence' : 0, #flag if want to output confidence intervals
                 'custom' : '', #if not empty, use a custom trait to predict metagenomes, specified by a file
    }
    all_opts.update(opts)

    cmd1 = ("picrust_cmd normalize_by_copy_number.py "
            + "-i " + file
            + " -o " + norm_out)
    if all_opts['gg_version']:
        cmd1 += " -g " + all_opts['gg_version']
    if all_opts['tab_in']:
        cmd1 += " -f"

    cmd2 = ("picrust_cmd predict_metagenomes.py "
            + "-i " + norm_out
            + " -o " + predict_out)
    if all_opts['gg_version']:
        cmd2 += " -g " + all_opts['gg_version']
    if all_opts['tab_out']:
        cmd2 += " -f"
    if all_opts['t']:
        cmd2 += " -t " + all_opts['t']
    if all_opts['with_confidence']:
        cmd2 += " --with_confidence"
    if all_opts['custom']:
        cmd2 += " -c " + all_opts['custom']

    return dict(name = "picrust:"+predict_out,
                actions = [cmd1, cmd2],
                file_dep = [file],
                targets = [predict_out, norm_out])

