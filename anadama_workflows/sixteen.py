"""16S workflows"""

import os
import operator
import itertools
from os.path import join


from anadama.action import CmdAction, PythonAction
from anadama.decorators import requires
from anadama import strategies

from anadama.util import (
    addext, 
    rmext,
    addtag,
    new_file,
    dict_to_cmd_opts, 
    guess_seq_filetype
)

from . import ( 
    settings
)

def write_map(sample_group, sample_dir):
    """Workflow to write a new map.txt file from a list of samples.  The
    resultant map.txt file is always named 'map.txt' and is placed in
    the ``sample_dir`` directory
    
    :param sample_group: List of namedtuples; A list of samples as
                         deserialized by anadama.util.deserialize_map_file
    :param sample_dir: String; Directory path indicating where to write 
                       the map.txt file

    """


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


@requires(binaries=['qiime_cmd'],
          version_methods=["qiime_cmd print_qiime_config.py "
                           "| awk '/QIIME library version/{print $NF;}'"])
def demultiplex(map_fname, fasta_fname, qual_fname, output_fname,
                qiime_opts={}):
    """Workflow to demultiplex a barcoded set of 16S sequences from a
    single run. This workflow wraps the qiime split_libraries.py
    script. For information on what the split_libraries.py script
    does, check out the qiime documentation:
    - http://qiime.org/tutorials/tutorial.html#assign-samples-to-multiplex-reads
    - http://qiime.org/scripts/split_libraries.html

    :param map_fname: String; File path location of the map.txt metdata file
    :param fasta_fname: String; File path to the input, multiplex, fasta files
    :param qual_fname: String; File path to the qual file corresponding 
                       to ``fasta_fname``.
    :param output_fname: String; File path to where the demultiplexed reads 
                         will be saved in fasta format.
    :keyword qiime_opts: Dictionary; A dictionary of command line options to
                         be passed to the wrapped split_libraries.py script. 
                         No - or -- flags are necessary; the correct - or --t
                         flags are inferred based on the length of the option. 
                         For boolean options, use the key/value pattern 
                         of { "my-option": "" }.

    External dependencies:
      - Qiime 1.8.0: https://github.com/qiime/qiime-deploy

    """
    
    
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


@requires(binaries=['qiime_cmd'],
          version_methods=["qiime_cmd print_qiime_config.py "
                           "| awk '/QIIME library version/{print $NF;}'"])
def demultiplex_illumina(fastq_fnames, barcode_fnames, map_fname, output_fname,
                         verbose=True, qiime_opts={}):

    output_dir=os.path.dirname(output_fname)
    default_opts = {
        "i": ",".join(fastq_fnames),
        "b": ",".join(barcode_fnames),
        "m": map_fname,
        "o": output_dir
    }
    default_opts.update(qiime_opts)
    opts = dict_to_cmd_opts(default_opts)
    
    cmd = "qiime_cmd split_libraries_fastq.py "

    revcomp_map_fname = new_file(addtag(map_fname, "revcomp"),
                                 basedir=output_dir)
    revcomp_opts = default_opts.copy()
    recvomp_opts['m'] = revcomp_map_fname
    revcomp_opts = dict_to_cmd_opts(revcomp_opts)
    def _revcomp():
        from anadama.util import deserialize_map_file
        from Bio.Seq import Seq

        def _reverse(sample):
            seq = Seq(sample.BarcodeSequence).reverse_complement()
            sample.BarcodeSequence = str(seq)
            return sample

        with open(map_fname) as from_map:
            from_samples = deserialize_map_file(from_map)
            serialize_map_file(
                ( _reverse(s) for s in from_samples ),
                revcomp_map_fname
            )


    def run():
        strategies.backup(
            (CmdAction(cmd+opts, verbose=verbose),
             strategies.Group(
                 PythonAction(_revcomp),
                 CmdAction(cmd+revcomp_opts,verbose=verbose))),
        )


    return {
        "name": "demultiplex_illumina:"+output_fname,
        "actions": [run],
        "file_dep": list(fastq_fnames) + list(barcode_fnames) + [map_fname],
        "targets": [output_fname]
    }



@requires(binaries=['qiime_cmd', 'sequence_convert'])
def pick_otus_closed_ref(input_fname, output_dir, verbose=None, qiime_opts={}):
    """Workflow to perform OTU picking, generates a biom-formatted OTU
    table from demultiplexed 16S reads. This workflow (in general
    terms) wraps qiime's pick_closed_reference_otus.py, which itself
    wraps either uclust or usearch. Note that uclust and usearch
    require a fairly large memory footprint (1.5-2.0G in some cases).

    :param input_fname: String; File path to the input,
                        fasta-formatted 16S sequences
    :param output_dir: String; Path to the directory where the output OTU 
                       table will be saved as 'otu_table.biom'. Other 
                       qiime-specific logs will go there, too.
    :keyword verbose: Boolean: set to true to print the commands that are 
                      run as they are run
    :keyword qiime_opts: Dictionary; A dictionary of command line options to
                         be passed to the wrapped split_libraries.py script. 
                         No - or -- flags are necessary; the correct - or --t
                         flags are inferred based on the length of the option. 
                         For boolean options, use the key/value pattern 
                         of { "my-option": "" }.

    External dependencies:
      - Qiime 1.8.0: https://github.com/qiime/qiime-deploy
      - USEARCH: (only if using the usearch option) 
        http://www.drive5.com/usearch/

    Resource utilization:
      - RAM: >1.5 G

    """

    output_fname = new_file("otu_table.biom", basedir=output_dir)
    revcomp_fname = new_file(
        os.path.basename(rmext(input_fname))+"_revcomp.fna",
        basedir=os.path.dirname(input_fname))

    verbose = settings.workflows.verbose if verbose is None else verbose

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
        strategies.backup(
            (CmdAction(cmd.format(input_fname),verbose=verbose),
             strategies.Group(
                 CmdAction(revcomp_cmd),
                 CmdAction(cmd.format(revcomp_fname),verbose=verbose))),
            extra_conditions = [ 
                lambda ret, output_fname: os.stat(output_fname).st_size == 0
            ],
            output_fname=output_fname,
        )
             
    return {
        "name": "pick_otus_closed_ref:"+input_fname,
        "actions": [run],
        "targets": [output_fname],
        "file_dep": [input_fname]
    }


@requires(binaries=['qiime_cmd', 'sequence_convert'])
def pick_otus_open_ref(input_fname, output_dir, verbose=None, qiime_opts={}):
    """Workflow to perform open-reference OTU picking. Similar to
    closed-reference OTU picking, this workflow generates a
    biom-formatted OTU table from demultiplexed 16S reads. This
    workflow (in general terms) wraps qiime's
    pick_open_reference_otus.py, which itself wraps either uclust or
    usearch. Note that uclust and usearch require a fairly large
    memory footprint (1.5-2.0G in some cases).

    :param input_fname: String; File path to the input,
                        fasta-formatted 16S sequences
    :param output_dir: String; Path to the directory where the output OTU 
                       table will be saved as 'otu_table.biom'. Other 
                       qiime-specific logs will go there, too.
    :keyword verbose: Boolean: set to true to print the commands that are 
                      run as they are run
    :keyword qiime_opts: Dictionary; A dictionary of command line options to
                         be passed to the wrapped split_libraries.py script. 
                         No - or -- flags are necessary; the correct - or --t
                         flags are inferred based on the length of the option. 
                         For boolean options, use the key/value pattern 
                         of { "my-option": "" }.

    External dependencies:
      - Qiime 1.8.0: https://github.com/qiime/qiime-deploy
      - USEARCH: (only if using the usearch option) 
        http://www.drive5.com/usearch/

    Resource utilization:
      - RAM: >1.5 G

    """

    output_fname = new_file("otu_table.biom", basedir=output_dir)
    revcomp_fname = new_file(
        "revcomp.fna", basedir=os.path.dirname(input_fname))

    verbose = settings.workflows.verbose if verbose is None else verbose

    default_opts = {
        "reference_fp": settings.workflows.sixteen.otu_refseq
    }
    default_opts.update(qiime_opts)
    opts = dict_to_cmd_opts(default_opts)

    cmd = ("qiime_cmd pick_open_reference_otus.py"+
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
        strategies.backup(
            (CmdAction(cmd.format(input_fname),verbose=verbose),
             strategies.Group(
                 CmdAction(revcomp_cmd),
                 CmdAction(cmd.format(revcomp_fname),verbose=verbose))),
            extra_conditions = [ 
                lambda ret, output_fname: os.stat(output_fname).st_size == 0
            ],
            output_fname=output_fname
        )

    return {
        "name": "pick_otus_open_ref:"+input_fname,
        "actions": [run],
        "targets": [output_fname],
        "file_dep": [input_fname]
    }


@requires(binaries=['qiime_cmd'])
def merge_otu_tables(files_list, name):
    """Workflow to merge OTU tables into a single OTU table. Also accepts
    biom-formatted OTU tables. This workflow will skip otu tables with
    a file size of zero at runtime.

    :param files_list: List of strings; A list of file paths to the input 
                       OTU tables to be merged
    :param name: String; The file name of the merged OTU table

    External dependencies:
      - Qiime 1.8.0: https://github.com/qiime/qiime-deploy

    """
    
    def merge_filter(deps,targets):
        files = [file for file in deps
                 if os.path.exists(file) and os.stat(file).st_size > 0] 
        cmd = "qiime_cmd merge_otu_tables.py -i {filenames} -o {output}"
        cmd = cmd.format( filenames = ",".join(files), 
                          output    = name  )
        return CmdAction(cmd, verbose=True).execute()

    return {
        "name": "merge_otu_tables: "+name,
        "actions": [(merge_filter, [files_list])],
        "targets": [name],
        "file_dep": files_list
    }


@requires(binaries=['picrust_cmd', 'biom'],
          version_methods=["picrust_cmd print_picrust_config.py "
                           "| awk '/PICRUSt version/{print $NF;}'"])
def picrust(file, output_dir=None, verbose=True, **opts):
    """Workflow to predict metagenome functional content from 16S OTU tables.

    :param file: String; input OTU table.
    :keyword tab_in: Boolean; True if the input is a tabulated 
                     file (default:0)
    :keyword tab_out: Boolean; True if the output file is to be
                      tabulated (default:False)
    :keyword gg_version: String; the greengenes version to be used
                         (default:most recent version)
    :keyword t: String; option to use a different type of prediction
                   (default:KO)
    :keyword with_confidence: Boolean; Set to True to output confidence 
                              intervals (default:0)
    :keyword custom: String; specify a file containing a custom trait to 
                     predict metagenomes

    External Dependencies:
      - PICRUSt: Version 1.0.0, 
        http://picrust.github.io/picrust/install.html#install

    """
    norm_out = new_file(addtag(file, "normalized_otus"), basedir=output_dir)
    predict_out = new_file(addtag(file, "picrust"), basedir=output_dir)

    all_opts = { 'tab_in'          : 0,  'tab_out' : 0,
                 'gg_version'      : '', 't'       : '', 
                 'with_confidence' : 0,  'custom'  : '' }
    all_opts.update(opts)

    cmd1 = ("picrust_cmd normalize_by_copy_number.py "
            + "-i %s"
            + " -o " + norm_out)
    if all_opts['gg_version']:
        cmd1 += " -g " + all_opts['gg_version']
    if all_opts['tab_in']:
        cmd1 += " -f"

    cmd2 = ("picrust_cmd predict_metagenomes.py "
            + "-i %s"
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


    converted = addtag(file, "json")
    format_cmd = CmdAction('biom convert --table-type="OTU table"'
                           ' --header-key taxonomy --to-json'
                           ' -i {} -o {} '.format(file, converted),
                           verbose=verbose)
    def run(targets):
        # try to run without converting to json, if that fails,
        # convert first, then run on the json-converted biom file
        if os.stat(file).st_size < 1:
            for target in targets:
                open(target, "w").close()
            return True

        return strategies.backup(
            (strategies.Group(CmdAction(cmd1%(file), verbose=verbose),
                              CmdAction(cmd2%(file), verbose=verbose)),
             strategies.Group(format_cmd,
                              CmdAction(cmd1%(converted), verbose=verbose),
                              CmdAction(cmd2%(converted), verbose=verbose)))
        )
             

    return dict(name = "picrust:"+predict_out,
                actions = [run],
                file_dep = [file],
                targets = [predict_out, norm_out])

