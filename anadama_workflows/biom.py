from anadama.decorators import requires

@requires(binaries=['biom'],
          version_methods=["pip freeze | grep biom-format"])
def to_tsv(infile, outfile):
    """Convert a biom file to a tsv (also known as pcl) file using the
    biom package

    :param infile: file of the biom format

    :param outfile: file of the tsv format

    External dependencies
      - biom-format: http://biom-format.org/
    """

    cmd="biom convert -i " + infile + " -o " + \
        outfile + " -b --header-key taxonomy " + \
        "--output-metadata-id \"Consensus Lineage\" " + \
        "--table-type 'otu table'"

    return {
        "name": "biom_to_tsv: " + infile,
        "actions": [cmd],
        "file_dep": [infile],
        "targets": [outfile]
    }
biom_to_tsv = to_tsv


@requires(binaries=['biom'],
          version_methods=["pip freeze | grep biom-format"])
def add_metadata(infile, outfile, sample_metadata):
    """Add sample-level metadata to a biom file. Sample-level metadata
    should be in a format akin to
    http://qiime.org/tutorials/tutorial.html#mapping-file-tab-delimited-txt

    :param infile: String; name of the biom file to which metadata 
                   shall be added
    :param outfile: String; name of the resulting metadata-enriched biom file
    :param sample_metadata: String; name of the sample-level metadata 
                            tab-delimited text file. Sample attributes are
                            taken from this file. Note: the sample names in
                            the `sample_metadata` file must match the sample
                            names in the biom file.

    External dependencies
      - biom-format: http://biom-format.org/
    """

    return {
        "name": "biom_add_metadata: " + infile,
        "actions": [("biom add-metadata"
                     " -i "+infile+
                     " -o "+outfile+
                     " -m "+sample_metadata)],
        "file_dep": [infile],
        "targets": [outfile]
    }


@requires(binaries=['biom'],
          version_methods=["pip freeze | grep biom-format"])
def from_pcl(infile, outfile):
    """Convert a pcl file to biom format using the biom package

    :param infile: String; name of the input tsv or pcl file
    :param outfile: String; name of the resulting converted biom file

    External dependencies
      - biom-format: http://biom-format.org/

    """

    cmd="biom convert -i " + infile + " -o " + \
        outfile + \
        " --table-type 'otu table'"

    return {
        "name": "biom_from_pcl: " + infile,
        "actions": [cmd],
        "file_dep": [infile],
        "targets": [outfile]
    }
