
import sys
import logging
import optparse
import tempfile
from pprint import pformat
from contextlib import nested

from Bio import SeqIO
from Bio import BiopythonParserWarning
HELP="""%prog [options] -F <format> --fasta-out <file> [--qual_out <file>]

%prog - Read in a sequence file from stdin, splitting sequence records
        to a pair of fasta and qual files, specified by the
        --fasta-out and the --qual-out options

Available formats:
""" 
HELP+=pformat(SeqIO._FormatToIterator.keys())

opts_list = [
    optparse.make_option('-F', '--format', action="store", 
                         dest="from_format", type="string",
                         help="The file format to convert from"),
    optparse.make_option('-f', '--fasta_out', action="store", 
                         dest="fasta_outfile", type="string",
                         help="File to which to write fasta records"),
     optparse.make_option('-q', '--qual_out', action="store", 
                         dest="qual_outfile", type="string", default=None,
                         help="File to which to write qual records"),
     optparse.make_option('-r', '--reverse_compliment', action="store_true", 
                          default=False,
                          help="Write the sequence reverse compliment. "+
                          "Defaults to True"),
     optparse.make_option('-l', '--logging', action="store", type="string",
                          dest="logging", default="INFO",
                          help="Logging verbosity, options are debug, info, "+
                          "warning, and critical")
]

def fa(record):
    return ">%s\n" %(record.id), "%s\n" %(str(record.seq))
        

def ql(record):
    return ">%s\n" %(record.id), "%s\n" %(
        ' '.join([ str(x) for x in record.letter_annotations['phred_quality'] ])
    )

def main():
    parser = optparse.OptionParser(option_list=opts_list, 
                                   usage=HELP)
    (opts, args) = parser.parse_args()
    logging.getLogger().setLevel(getattr(logging, opts.logging.upper()))
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s")

    if not opts.fasta_outfile:
        parser.print_usage()
        sys.exit(1)

    if not args:
        logging.debug("Writing input to temporary file")
        tmp_fp = tempfile.TemporaryFile()
        map(tmp_fp.write, sys.stdin)
        tmp_fp.seek(0)
        args = [tmp_fp]
    else:
        args = [ open(f) for f in args ]
        
            
    with nested(*args), open(opts.fasta_outfile, 'w') as fa_file:

        if opts.qual_outfile:
            # opening the qual_file here pains me
            qual_file = open(opts.qual_outfile, 'w')
            def _output(record):
                SeqIO.write(record, fa_file, "fasta")
                SeqIO.write(record, qual_file, "qual")
        else:
            def _output(record):
                SeqIO.write(record, fa_file, "fasta")

        if opts.reverse_compliment:
            def output(record):
                record =  SeqIO.SeqRecord(
                    record.reverse_complement().seq,
                    letter_annotations=record.letter_annotations,
                    id=record.id, 
                    name=record.name, 
                    description=record.description
                )
                return _output(record)
        else:
            output=_output

        for fp in args:
            for i, record in enumerate(SeqIO.parse(fp, opts.from_format)):
                tmp = record.letter_annotations.pop('phred_quality')[4:]
                record.seq = record.seq[4:]
                record.letter_annotations['phred_quality'] = tmp
                try:
                    output(record)
                except BiopythonParserWarning as e:
                    print >> sys.stderr, e
                if logging.getLogger().isEnabledFor(logging.DEBUG):
                    if i % 250 == 0 and i != 0:
                        logging.debug("Converted %d records", i)
            
                

if __name__ == '__main__':
    main()
