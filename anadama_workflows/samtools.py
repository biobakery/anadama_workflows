from anadama.util import dict_to_cmd_opts
from anadama.decorators import requires


@requires(binaries=["samtools"],
          version_methods=["samtools 2>&1 | awk '/Version/{print $2;}'"])
def sort(input_bam, output_prefix, 
         memory_level="768M", num_threads=1, 
         **kwargs):
    """Sort a bam file by sequence name with samtools.
    
    :param input_bam: String; file name of input bam file.

    :param output_prefix: String; file name of output sorted bam file
    without the .bam suffix

    :keyword memory_level: String; K/M/G human readable amount of ram
    to give each sorting thread

    :keyword num_threads: Int or string; number of threads to use when
    sorting.

    """


    output_file = output_prefix+".bam"

    opts = { 'n': "", '@': str(num_threads), 
             'm': memory_level }

    opts.update(kwargs)
    cmd = ("samtools sort "
           +" "+dict_to_cmd_opts(opts)
           +" "+input_bam
           +" "+output_prefix)

    return { "name": "samtools.sort: "+output_file,
             "file_dep": [input_bam],
             "actions": [cmd],
             "targets": [output_file] }


def to_paired_fastq(input_bam, output_prefix, **kwargs):
    """Sort a bam file by name, then split into three fastq files based on
    ``output_prefix``: forward reads, reverse reads, and singleton reads

    :param input_bam: String; file name of input bam file.

    :param output_prefix: String; file name of output fastq files
    without the ``.r1.fastq`` etc. suffix.

    :keyword **kwargs: All extra keyword arguments are passed to
    :py:func:`anadama_workflows.samtools.sort`

    """

    output_r1 = output_prefix+".r1.fastq"
    output_r2 = output_prefix+".r2.fastq"
    output_se = output_prefix+".single.fastq"

    opts = { 'o': "" }
    opts.update(kwargs)

    sort_cmd = sort(input_bam, output_prefix, **opts).next()['actions'][0]
    pe_split_cmd = "bam_pe_split "+output_prefix
    cmd = (sort_cmd+" | "+pe_split_cmd)

    return { "name": "samtools.to_paired_fastq %s..."%(output_r1),
             "file_dep": [input_bam],
             "targets": [output_r1, output_r2, output_se],
             "actions": [cmd] }
