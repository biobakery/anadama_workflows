import re
import sys
import logging
import optparse
import operator
from pprint import pformat
from functools import partial
from itertools import ifilter, imap

try:
    import bz2
    import gzip
except:
    pass

import pysam
from Bio import SeqIO, Seq, SeqRecord

HELP="""%prog [options] -f <format> [-t <format>] [<file> [<file> [ ...]]]

%prog - Convert sequence files from one format to another, printing
results to stdout

Available formats:
""" 


comparators = { '=': operator.eq, '==': operator.eq, 
                '!': operator.ne, '!=': operator.ne,
                '>': operator.gt, '>=': operator.ge,
                '<': operator.lt, '<=': operator.le }

opts_list = [
    optparse.make_option('-f', '--format', action="store", 
                         dest="from_format", type="string",
                         help="The file format to convert from"),
    optparse.make_option('-t', '--to', action="store", 
                         dest="to_format", type="string", default="fasta",
                         help="The file format to convert to"),
    optparse.make_option('-r', '--reverse_complement', action="store_true", 
                         dest="revcomp", 
                         help="Write the reverse complement sequence"),
    optparse.make_option('-l', '--logging', action="store", type="string",
                         dest="logging", default="INFO",
                         help="Logging verbosity, options are debug, info,"
                         " warning, and critical"),
    optparse.make_option('-n', '--filter_by_len', action='append', 
                         type="string", dest='lenfilters', 
                         help="Filter sequences by length. If a sequence has a"
                         " length that matches this condition, then it's kept."
                         " Otherwise, the sequence is discarded. recognized"
                         " conditions are {}".format(
                             comparators.keys())),
    optparse.make_option('-s', '--slice', action="store", type="string", 
                         dest="slice", default="",
                         help="Slice sequence from start:end")
]


def parse_comparison(s):
    match = re.match(r'^([><=!]{1,2})(\d+)$', s)
    if not match:
        raise ValueError("Unrecognized comparison string;"
                         " try something like >=60.")

    condition, int_ = match.group(1), int(match.group(2))
    if condition not in comparators:
        raise ValueError("Unrecognized condition."
                         " Try something like >=60.")

    comparator = comparators[condition]
    return lambda val: comparator(val, int_)


def zipcheck():
    ret = True
    for zipmod in ("bz2", "gzip"):
        if zipmod not in globals():
            logging.warn("%s module not available; "
                         "unable to handle %s encoded files"%(
                             zipmod, zipmod))
            ret = False

    return ret


def generate_filter(comparison_list):
    comparator_funcs = [ parse_comparison(str_) for str_ in comparison_list ]
    return lambda val: all( f(len(val)) for f in comparator_funcs )


def generate_slicer(slice_str):
    if slice_str:
        start, stop = map(int, slice_str.split(":"))
        return lambda val: val[start:stop]
    else:
        return lambda val: val


def my_open(f, *args, **kwargs):
    if f == '-':
        return sys.stdin
    elif f.endswith(".bz2"):
        return bz2.BZ2File(f, *args, **kwargs)
    elif f.endswith("gzip") or f.endswith("gz"):
        return gzip.GzipFile(f, *args, **kwargs)
    else:
        return open(f, *args, **kwargs)


def handle_samfile(file_str, filemode="r"):
    sam_file = pysam.Samfile(file_str, filemode, 
                             check_header=False, check_sq=False)
    for read in sam_file:
        seq = Seq.Seq(read.seq)
        qual = [ ord(x)-33 for x in read.qual ]
        if read.is_reverse:
            seq = seq.reverse_complement()
        yield SeqRecord.SeqRecord(
            seq, read.qname, "", "",
            letter_annotations={"phred_quality": qual}
        )
            

def handle_biopython(file_str, format=None):
    in_file = my_open(file_str)
    return SeqIO.parse(in_file, format)

formats = {
    None: handle_biopython,
    "sam": handle_samfile,
    "bam": partial(handle_samfile, filemode="rb"),
}
formats.update( (key, partial(handle_biopython, format=key)) 
                for key in SeqIO._FormatToWriter.keys() )

def handle_cli():
    global HELP
    HELP += pformat(formats.keys())
    parser = optparse.OptionParser(option_list=opts_list, 
                                   usage=HELP)
    return parser.parse_args()

def maybe_reverse_complement(seqs, do_reverse):
    if not do_reverse:
        for record in seqs:
            yield record
    else:
        for record in seqs:
            revcomp = record.reverse_complement()
            yield SeqIO.SeqRecord(
                revcomp.seq,
                letter_annotations=revcomp.letter_annotations,
                id=record.id,
                name=record.name,
                description=record.description
            )

def convert(*input_files, **opts):
    from_format = opts["format"]
    to_format = opts["to"]
    revcomp = opts["revcomp"]
    filter_func = opts['filter_']
    slice_ = opts['slice_']

    for in_file in input_files:
        logging.debug("Converting %s from %s to %s", 
                      in_file, from_format, to_format)
        sequences = formats[from_format](in_file)
        sequences = maybe_reverse_complement(sequences, revcomp)
        sequences = ifilter(filter_func, sequences)
        sequences = imap(slice_, sequences)
        for i, inseq in enumerate(sequences):
            SeqIO.write(inseq, sys.stdout, to_format)

            if logging.getLogger().isEnabledFor(logging.DEBUG):
                if i % 250 == 0 and i != 0:
                    logging.debug("Converted %d records", i)

def main():
    (opts, input_files) = handle_cli()
    logging.getLogger().setLevel(getattr(logging, opts.logging.upper()))
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s")

    if not opts.to_format:
        parser.print_usage()
        sys.exit(1)

    zipcheck()

    if not input_files:
        input_files = ['-']

    lenfilters = opts.lenfilters if opts.lenfilters else []

    try:
        convert(*input_files, 
                 format=opts.from_format, 
                 to=opts.to_format,
                 revcomp=opts.revcomp,
                 filter_=generate_filter(lenfilters),
                 slice_=generate_slicer(opts.slice))
    except IOError as e:
        if e.errno == 32:
            # That's the error for a broken pipe this usually happens
            # when someone piped to head or tail. I don't want to
            # raise an exception
            pass
        else:
            raise


if __name__ == '__main__':
    main()
