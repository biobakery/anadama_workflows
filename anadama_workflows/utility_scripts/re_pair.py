#!/usr/bin/env python

import sys
import optparse
from pprint import pformat
from itertools import chain, repeat
from collections import namedtuple

from Bio import SeqIO

HELP="""%prog [options] -f <format> [-t <format>] -b barcode.fa <file1.fa> <file2.fa> [...]

Available formats:
""" 

opts_list = [
    optparse.make_option('-f', '--format', action="store", 
                         dest="from_format", type="string",
                         help="The file format to convert from. Required!"),
    optparse.make_option('-t', '--to', action="store", 
                         dest="to_format", type="string", default="fasta",
                         help="The file format to convert to. Default fasta."),
    optparse.make_option('-b', '--barcode', action="store", 
                         dest="barcode_file", type="string", 
                         help="Barcode file to match headers"),
]

formats = SeqIO._FormatToWriter.keys()
HELP += pformat(formats)

nullobject = namedtuple("null", "id")(id=None)

class Cache(object):
    def __init__(self, filename, format):
        self.seqs = chain(SeqIO.parse(filename, format), repeat(nullobject))
        self.peek = self.seqs.next()

    def next(self):
        ret = self.peek
        self.peek = self.seqs.next()
        return ret
        

class Matcher(object):
    def __init__(self, seq_file_list, format):
        self.caches = [Cache(f, format) for f in seq_file_list]
        
    def match(self, sequence):
        for cache in self.caches:
            if sequence.id == cache.peek.id:
                return cache.next()

        raise Exception("Unable to find header %s in reads" %sequence.id)



def main():
    parser = optparse.OptionParser(option_list=opts_list, 
                                   usage=HELP)
    (opts, sequences) = parser.parse_args()

    if not opts.from_format or not opts.barcode_file:
        parser.print_usage()
        sys.exit(1)

    m = Matcher(sequences, opts.from_format)
    with open(opts.barcode_file) as bc_f:
        for sequence in SeqIO.parse(bc_f, opts.from_format):
            SeqIO.write(m.match(sequence), sys.stdout, opts.to_format)


if __name__ == '__main__':
    main()

