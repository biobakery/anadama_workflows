import os, re

from anadama.decorators import requires

@requires(binaries=['biom'])
def biom_to_tsv(infile, outfile):
    """ Convert a biom file to a tsv using the biom package 

    param: infile: file of the biom format

    param: outfile: file of the tsv format

    External dependencies
    - biom-format: http://biom-format.org/
    """

    cmd="biom convert -i " + infile + " -o " + \
        outfile + " --to-tsv --header-key taxonomy " + \
        "--output-metadata-id \"Consensus Lineage\""

    return {
        "name": "biom_to_tsv: " + infile,
        "actions": [cmd],
        "file_dep": [infile],
        "targets": [outfile]
    }

@requires(binaries=['qiimeToMaaslin.py'])
def qiime_to_maaslin(in_datafile, outfile):
    """ Converts a tsv file from qiime to a maaslin format

    param: in_datafile: file of the qiime tsv format (from biom_to_tsv)

    param: outfile: file of the tsv format for input to maaslin

    External dependencies
    - QiimeToMaaslin: https://bitbucket.org/biobakery/qiimetomaaslin

    """

    cmd="qiimeToMaaslin.py < " + in_datafile + \
        " > " + outfile

    return {
        "name": "qiime_to_maaslin: " + in_datafile,
        "actions": [cmd],
        "file_dep": [in_datafile],
        "targets": [outfile]
    }

@requires(binaries=['merge_metadata.py'])
def merge_otu_metadata(otu_table,metadata_file,outfile):
    """ Use the merge python script from maaslin to merge
        the otu_table and the metadata into a single file.
        

    param: otu_table: file of the tsv format (not formatted by qiime)

    param: metadata_file: file of metadata in tsv format

    param: outfile: file that is a merge of the otu and metadata

    External dependencies
    - Maaslin: https://bitbucket.org/biobakery/maaslin
        
    """

    cmd="merge_metadata.py " + metadata_file + " < " + \
            otu_table + " > " + outfile

    return {
        "name": "merge_otu_metadata: " + otu_table,
        "actions": [cmd],
        "file_dep": [otu_table, metadata_file],
        "targets": [outfile]
    }

def create_maaslin_read_config(metadata_file, pcl_file, read_config_file):
    """ Creates a read config file for the maaslin run
        using the otu and metadata files
    """

    def create_config(metadata_file, pcl_file, read_config_file):
        file_handle=open(metadata_file,"r")
        header=file_handle.readline()
        file_handle.close()

        # sort the columns
        columns=header.split("\t")
        metadata_reference_column=sorted(columns)[-1]

        # identify the reference for the abundance column
        file_handle=open(pcl_file,"r")
        line=file_handle.readline()

        abundance_reference_column=""
        while line:
            if re.search("^"+metadata_reference_column,line):
                # this next line has the reference for the abundance
                ref_line=file_handle.readline()
                abundance_reference_column=ref_line.split("\t")[0]
                break                        
            line=file_handle.readline()
        file_handle.close()

        read_config_file_contents = "Matrix: Metadata\n" + \
            "Read_PCL_Rows: -" + metadata_reference_column  + "\n" + \
            "\n" + \
            "Matrix: Abundance\n" + \
            "Read_PCL_Rows: " + abundance_reference_column + "-"

        # write the read.config file
        file_handle=open(read_config_file,"w")
        file_handle.write(read_config_file_contents)
        file_handle.close()  


    return {
        "name": "create_maaslin_read_config: " + metadata_file,
        "actions": [(create_config,[metadata_file,pcl_file,read_config_file])],
        "file_dep": [metadata_file,pcl_file],
        "targets": [read_config_file]
    }

@requires(binaries=['transpose.py'])
def transpose(pcl_file, outfile):
    """ Transpose the merged pcl and metadata file
    
    param: pcl_file: file that is the merge of the otu and metadata    

    param: outfile: file that is the transpose of the pcl_file

    External dependencies
    - Maaslin: https://bitbucket.org/biobakery/maaslin
        
    """

    cmd="transpose.py < " + pcl_file + " > " + outfile

    return {
        "name": "transpose: " + pcl_file,
        "actions": [cmd],
        "file_dep": [pcl_file],
        "targets": [outfile]
    }

@requires(binaries=['Maaslin.R'])
def run_maaslin(read_config_file, maaslin_outfile, pcl_file):
    """ Run the maaslin software

    param: read_config_file: file generated based on pcl and metadata

    param: maaslin_outfile: one of the output files written by maaslin
    
    param: pcl_file: file generated to meet maaslin format requirements

    External dependencies
    - Maaslin: https://bitbucket.org/biobakery/maaslin
        
    """

    cmd = "Maaslin.R -i " + read_config_file + " " + \
        maaslin_outfile + " " + pcl_file
                
    return {
        "name": "run_maaslin: " + pcl_file,
        "actions": [cmd],
        "file_dep": [read_config_file, pcl_file],
        "targets": [maaslin_outfile]
    }


def maaslin(otu_table, metadata_file):
    """ Workflow to compute the significance of association in microbial community
        using a transform abundance or relative function
        table obtained from Qiime, HUMAnN or MetaPhlAn plus study metadata

    param: otu_table: file of the biom or tsv format

    param: metadata_file: file of metadata in tsv format

    External dependencies
    - Maaslin: https://bitbucket.org/biobakery/maaslin
        
    """
    
    # place the output in the same location as the input
    outdir = os.path.dirname(otu_table)
    
    # keep the same project naming convention
    project_name = os.path.splitext(os.path.basename(otu_table))[0] 

    # new files full path and project base name
    new_file_basename=os.path.join(outdir,project_name)

    # initial_targets[1] is used for the name by maaslin for output files
    # (ie maaslin_demo2 for maaslin_demo2.tsv")
    initial_targets = [new_file_basename + "_maaslin.pcl",
        new_file_basename + "_maaslin.tsv"]    
    
    # project specific read_config_file
    read_config_file = new_file_basename + "_maaslin.read.config"
    
    # the directory of the final_targets[0] is where the output will be written    
    final_targets = [new_file_basename + "_maaslin.txt", 
        new_file_basename + "_maaslin_log.txt"]
    
    # need to merge otu table and metadata
    if otu_table.endswith(".biom"):
        otu_table_tsv_format=new_file_basename + "_format.tsv"
        yield biom_to_tsv(otu_table,otu_table_tsv_format)
        otu_table_maaslin_format=new_file_basename + "_format_maaslin.tsv"
        yield qiime_to_maaslin(otu_table_tsv_format, otu_table_maaslin_format)
        yield merge_otu_metadata(otu_table_maaslin_format, metadata_file, initial_targets[0])
    else:
        yield merge_otu_metadata(otu_table, metadata_file,initial_targets[0])
    
    # create new project specific read.config file for maaslin input
    yield create_maaslin_read_config(metadata_file, initial_targets[0], read_config_file)

    # for command line need to transpose first
    yield transpose(initial_targets[0],initial_targets[1])
          
    yield run_maaslin(read_config_file, final_targets[0], initial_targets[1])  

