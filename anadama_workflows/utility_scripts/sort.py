#!/usr/bin/env python

import sys
import optparse
import errno
from functools import partial
from pprint import pformat
from operator import attrgetter

from Bio import SeqIO

HELP="""%prog [options] -f <format> [-t <format>] [<file> [<file> [ ...]]]

%prog - Sort sequence files according to the sequence ID 
        (the thing that comes after the `>' symbol in FASTA files).

Available formats:
""" 

opts_list = [
    optparse.make_option('-f', '--format', action="store", 
                         dest="from_format", type="string",
                         help="The file format to convert from. Required!"),
    optparse.make_option('-t', '--to', action="store", 
                         dest="to_format", type="string", default="fasta",
                         help="The file format to convert to. Default fasta."),
]

id_ = attrgetter('id')

formats = SeqIO._FormatToWriter.keys()
HELP += pformat(formats)

def main():
    parser = optparse.OptionParser(option_list=opts_list, 
                                   usage=HELP)
    (opts, sequences) = parser.parse_args()

    if not opts.from_format:
        parser.print_usage()
        sys.exit(1)

    for file_ in sequences:
        with open(file_, 'r') as f_in:
            seqs = SeqIO.parse(f_in, opts.from_format)
            try:
                SeqIO.write(sorted(seqs, key=id_), sys.stdout, opts.to_format)
            except IOError as e:
                if e.errno == errno.EPIPE:
                    sys.exit(0)
            
            del seqs


if __name__ == '__main__':
    main()

