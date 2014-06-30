
import sys
import logging
import optparse
from pprint import pformat
from contextlib import nested
from functools import partial

import pysam
from Bio import SeqIO, Seq, SeqRecord

HELP="""%prog [options] -f <format> [-t <format>] [<file> [<file> [ ...]]]

%prog - Convert sequence files from one format to another, printing
results to stdout

Available formats:
""" 

opts_list = [
    optparse.make_option('-f', '--format', action="store", 
                         dest="from_format", type="string",
                         help="The file format to convert from"),
    optparse.make_option('-t', '--to', action="store", 
                         dest="to_format", type="string", default="fasta",
                         help="The file format to convert to"),
    optparse.make_option('-r', '--reverse_compliment', action="store_true", 
                         dest="revcomp", 
                         help="Write the reverse complement sequence"),
     optparse.make_option('-l', '--logging', action="store", type="string",
                         dest="logging", default="INFO",
                         help="Logging verbosity, options are debug, info, "+
                         "warning, and critical")
]


def my_open(f, *args, **kwargs):
    if f == '-':
        return sys.stdin
    else:
        return open(f, *args, **kwargs)



def handle_samfile(file_str, filemode="r"):
    sam_file = pysam.Samfile(file_str, filemode, 
                             check_header=False, check_sq=False)
    for read in sam_file:
        seq = Seq.Seq(read.seq)
        qual = [ ord(x)-32 for x in read.qual ]
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

    for in_file in input_files:
        logging.debug("Converting %s from %s to %s", 
                      in_file, from_format, to_format)
        sequences = formats[from_format](in_file)
        sequences = maybe_reverse_complement(sequences, revcomp)
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

    if not input_files:
        input_files = ['-']

    try:
        convert(*input_files, 
                format=opts.from_format, 
                to=opts.to_format,
                revcomp=opts.revcomp)
    except IOError as e:
        if e.errno == 32:
            pass
        else:
            raise


if __name__ == '__main__':
    main()
