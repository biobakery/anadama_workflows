import sys
import optparse
import logging

import pysam
from Bio import SeqIO, Seq, SeqRecord

HELP = """ %prog <output_prefix>

%prog - Split from a sorted-by-name bam file (via stdin) three fastq
files: forward reads, reverse reads, and single-end reads.

<output_prefix> - file name prefix of the three resulting fastq files.
                  Resultant files are named output_prefix.r1.fastq,
                  output_prefix.r2.fastq, and output_prefix.single.fastq 
"""

opts_list = [
    optparse.make_option('-l', '--logging', action="store", type="string",
                         dest="logging", default="INFO",
                         help="Logging verbosity, options are debug, info,"
                         " warning, and critical"),
]

def handle_cli():
    parser = optparse.OptionParser(option_list=opts_list, usage=HELP)
    opts, args = parser.parse_args()
    if len(args) > 1:
        parser.print_usage()
        sys.exit(1)
    return opts, args[0]



def _write(bam_seq, f):
    seq = Seq.Seq(bam_seq.seq)
    qual = [ ord(x)-33 for x in bam_seq.qual ]
    if bam_seq.is_reverse:
        seq = seq.reverse_complement()

    record = SeqRecord.SeqRecord(
        seq, bam_seq.qname, "", "",
        letter_annotations={"phred_quality": qual}
    )
    SeqIO.write(record, f, "fastq")


def output(sequences, r1_file, r2_file, se_file):

    if logging.getLogger().isEnabledFor(logging.DEBUG):
        def log(i):
            if i % 1000 == 0 and i != 0:
                logging.debug("Converted %d records", i)
    else:
        def log(i):
            pass

    for i, seq in enumerate(sequences):
        if not seq.is_paired:
            _write(seq, se_file)
        elif seq.is_read1:
            _write(seq, r1_file)
        elif seq.is_read2:
            _write(seq, r2_file)
        log(i)


def main():
    opts, output_prefix = handle_cli()
    logging.getLogger().setLevel(getattr(logging, opts.logging.upper()))
    logging.basicConfig(
        format="%(asctime)s %(levelname)s: %(message)s")
    
    sequences = pysam.Samfile("-", 'rb', check_header=False, check_sq=False)
    with open(output_prefix+".r1.fastq", 'w') as r1_file, \
         open(output_prefix+".r2.fastq", 'w') as r2_file, \
         open(output_prefix+".single.fastq", 'w') as se_file:
        output(sequences, r1_file, r2_file, se_file)


if __name__ == '__main__':
    ret = main()
    sys.exit(ret)
