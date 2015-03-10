#!/usr/bin/env python

import sys
import optparse
import errno
from functools import partial
from pprint import pformat
from operator import attrgetter

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

id_ = attrgetter('id')

formats = SeqIO._FormatToWriter.keys()
HELP += pformat(formats)

class Cache(object):
    def __init__(self, filename, format):
        self.seqs = SeqIO.parse(filename, format)
        self._cached = self.seqs.next()

    def _read(self):
        if not self._cached:
            try:
                self._cached = self.seqs.next()
            except StopIteration:
                self._cached = None
        return self._cached

    def remove(self):
        self._cached = None


class Matcher(object):
    def __init__(self, seq_file_list, format):
        self.caches = [Cache(f, format) for f in seq_file_list]
        
    def match(self, sequence):
        for cache in self.caches:
            s = cache._read()
            if sequence.id == s.id:
                cache.remove()
                return s

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

