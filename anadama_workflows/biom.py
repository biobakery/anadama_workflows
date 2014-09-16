from anadama.decorators import requires

@requires(binaries=['biom'])
def to_tsv(infile, outfile):
    """ Convert a biom file to a tsv using the biom package 

    param: infile: file of the biom format

    param: outfile: file of the tsv format

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

@requires(binaries=['biom'])
def add_metadata(infile, outfile, sample_metadata):
    return {
        "name": "biom_add_metadata: " + infile,
        "actions": [("biom add-metadata"
                     " -i "+infile+
                     " -o "+outfile+
                     " -m "+sample_metadata)],
        "file_dep": [infile],
        "targets": [outfile]
    }


@requires(binaries=['biom'])
def from_pcl(infile, outfile):
    """ Convert a biom file to a tsv using the biom package 

    param: infile: file of the biom format

    param: outfile: file of the tsv format

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
