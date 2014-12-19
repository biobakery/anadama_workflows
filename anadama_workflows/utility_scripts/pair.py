#!/usr/bin/env python

import re
import logging
import optparse
from pprint import pformat
from itertools import izip_longest, ifilter, count, izip

from Bio import SeqIO
from Bio import BiopythonParserWarning

DEBUG = False

HELP="""%prog [options] -f <format> <read_1_file> <read_2_file>

%prog - Pair reads according to the sequence ID
        (the thing that comes after the `>' symbol in FASTA files).
        Unpaired reads are discarded.

Available formats:
""" 

opts_list = [
    optparse.make_option('-f', '--format', action="store", 
                         dest="from_format", type="string",
                         help="The file format to convert from. Required!"),
    optparse.make_option('-t', '--to', action="store", 
                         dest="to_format", type="string", default="fasta",
                         help="The file format to convert to. Default fasta."),
    optparse.make_option('-1', '--r1_out', action="store", 
                         dest="r1_out", type="string", default=None,
                         help=("Filename to save read1 output. Defaults to "
                               "the filename of read_1_file with `-paired' "
                               "added on the end")),
    optparse.make_option('-2', '--r2_out', action="store", 
                         dest="r2_out", type="string", default=None,
                         help=("Filename to save read2 output. Defaults to "
                               "the filename of read_2_file with `-paired' "
                               "added on the end")),
    optparse.make_option('-r', '--compare_by_regex', action="store", 
                         dest="compare_regex", type="string", 
                         default=r'^(\S+)\s?.*$',
                         help="Use this PCRE regex on the sequence ID to "
                         "match the sequence comparison key. Defaults to "
                         "everything before the first space character"),
    optparse.make_option('-l', '--logging', action="store", type="string",
                         dest="logging", default="INFO",
                         help="Logging verbosity, options are debug, info, "+
                         "warning, and critical"),
    optparse.make_option('-i', '--inner_join', action="store", 
                         dest="inner_join", default=None,
                         help="Assume the sequences are sorted and perform a "
                         "left or right inner join of one sequence set to the "
                         "other. Acceptable values are 'left' and 'right'. "
                         "Defaults to an equi-join, where any sequence on "
                         "either side is dropped if no pair is found. Using "
                         "either 'left' or 'right' reduces runtime and memory "
                         "resources"),
    optparse.make_option('-m', '--no_mangle', action="store_false", 
                         default=True, dest="mangle",
                         help="Output sequence labels with just the ID in "
                         "it. Defaults to true"),
]

def handle_cli():
    global HELP
    formats = SeqIO._FormatToWriter.keys()
    HELP += pformat(formats)
    parser = optparse.OptionParser(option_list=opts_list, 
                                   usage=HELP)
    return parser.parse_args()


def generate_fnames(opts, r1_fname, r2_fname):
    for in_, out in [(r1_fname, opts.r1_out), (r2_fname, opts.r2_out)]:
        if not out:
            out = in_+'-paired.'+opts.to_format
        logging.debug("Saving output to %s", out)
        yield out
            

read_cache1 = dict()
read_cache2 = dict()
def _pair_reads(keys, reads):
    global read_cache1, read_cache2
    r1, r2 = reads
    k1, k2 = keys
    # time optimization; if they already equal eachother, return them
    # immediately
    if k1 == k2:
        yield r1, r2

    if k1 in read_cache2:
        yield r1, read_cache2.pop(k1)
    else:
        if r1 and k1:
            read_cache1[k1] = r1

    if k2 in read_cache1:
        yield read_cache1.pop(k2), r2
    else:
        if r2 and k2:
            read_cache2[k2] = r2


def _pair_reads_cached(seqs, regex):
    for i, reads in enumerate(seqs):
        try:
            keys = extract_compare_key(reads, regex)
            maybe_matched = _pair_reads(keys, reads)
            for mated_pair in ifilter(None, maybe_matched):
                # now we know they're mated pairs, so return them
                yield mated_pair
        except BiopythonParserWarning as e:
            logging.warning(e)
            continue
        except ValueError as e:
            logging.warning('Exception on record %i', i)
            logging.exception(e)
            continue

    if DEBUG:
        global read_cache1, read_cache2
        logging.debug("Dropped %i orphaned reads:", 
                      len(read_cache1)+len(read_cache2))
        logging.debug(read_cache1.keys()+read_cache2.keys())


def _pair_reads_sorted(join_direction="left"):
    """Assume the sequences come pre-sorted and perform a right or left
    inner join, depending on the ``join_direction`` keyword argument."""

    def _join(forward_reads, reverse_reads, cmp_regex):
        counter = count(1)
        r1s, r2s = forward_reads, reverse_reads
        if join_direction.lower() == "left":
            def _discard(r1, r2, k1, k2):
                logging.debug("Dropped read: %s", k2)
                i, r2 = counter.next(), r2s.next()
                k2 = extract_compare_key([r2], cmp_regex, track=i).next()
                return i, r1, r2, k1, k2

        else:
            def _discard(r1, r2, k1, k2):
                logging.debug("Dropped read: %s", k1)
                i, r1 = counter.next(), r1s.next()
                k1 = extract_compare_key([r1], cmp_regex, track=i).next()
                return i, r1, r2, k1, k2

        while True:
            i, r1, r2 = counter.next(), r1s.next(), r2s.next()
            k1, k2 = extract_compare_key((r1, r2), cmp_regex, track=i)
            while k1 != k2:
                i, r1, r2, k1, k2 = _discard(r1, r2, k1, k2)
            yield r1, r2


    return _join


def extract_compare_key(reads, regex, track=0):
    for read in reads:
        if not read:
            yield None
            continue
        match = regex.match(read.id)
        if not match or not match.groups():
            logging.warning("Provided regex has no matches or gives no "
                            "groups against string %s on line %i", 
                            read.id, track)
        else:
            yield match.group(1)
        

def pair_reads(in1, in2, id_match_pattern, parse_format, join_direction=None):
    regex = re.compile(id_match_pattern)
    if join_direction:
        pairer = _pair_reads_sorted(join_direction=join_direction)
        return pairer(forward_reads = SeqIO.parse(in1, parse_format),
                      reverse_reads = SeqIO.parse(in2, parse_format),
                      cmp_regex     = regex)
    else:
        seqs = izip_longest(SeqIO.parse(in1, parse_format),
                            SeqIO.parse(in2, parse_format))
        return _pair_reads_cached(seqs, regex)


def _output(read_pairs, output_pair, output_format, only_id=True):
    out1, out2 = output_pair
    for i, (r1, r2) in enumerate(read_pairs):
        if only_id:
            r1.description = r1.id + " 1"
            r2.description = r2.id + " 2"
        SeqIO.write(r1, out1, output_format)
        SeqIO.write(r2, out2, output_format)
        if DEBUG:
            if i % 100 == 0 and i != 0:
                logging.debug("Wrote %d records", i)


def main():
    opts, (r1_fname, r2_fname) = handle_cli()

    global DEBUG
    level = getattr(logging, opts.logging.upper())
    logging.getLogger().setLevel(level)
    DEBUG = logging.getLogger().isEnabledFor(logging.DEBUG)
    logging.basicConfig(format="%(asctime)s %(levelname)s: %(message)s")

    r1out_fname, r2out_fname = generate_fnames(opts, r1_fname, r2_fname)

    with   open(r1out_fname, 'w') as r1out, \
           open(r2out_fname, 'w') as r2out, \
           open(r1_fname,    'r') as r1in,  \
           open(r2_fname,    'r') as r2in:
        paired = pair_reads(in1=r1in, in2=r2in,
                            id_match_pattern=opts.compare_regex,
                            parse_format=opts.to_format,
                            join_direction=opts.inner_join)
        _output(paired, (r1out, r2out),
                output_format=opts.to_format, only_id=opts.mangle)


    

if __name__ == '__main__':
    main()

