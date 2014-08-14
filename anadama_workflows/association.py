import os, re

from anadama.decorators import requires

@requires(binaries=['merge_metadata.py', 'transpose.py', 'Maaslin.R'])
def maaslin(otu_table, metadata_file):
    ''' Workflow to compute the significance of association in microbial community
        using a transform abundance or relative function
        table obtained from Qiime, HUMAnN or MetaPhlAn plus study metadata
    '''
    
    # place the output in the same location as the input
    outdir = os.path.dirname(otu_table)
    
    # keep the same project naming convention
    project_name = os.path.splitext(os.path.basename(otu_table))[0]
    
    # initial_targets[1] is used for the name by maaslin for output files
    # (ie maaslin_demo2 for maaslin_demo2.tsv")
    initial_targets = [outdir + "/" + project_name + "_maaslin.pcl",
        outdir + "/" + project_name + "_maaslin.tsv"]    
    
    # project specific read_config_file
    read_config_file = outdir + "/" + project_name + "_maaslin.read.config"
    
    # the directory of the final_targets[0] is where the output will be written    
    final_targets = [outdir + "/" + "demo.txt", 
        outdir + "/" + project_name + ".txt",  
        outdir + "/" + project_name + "_log.txt"]
    
    all_targets = final_targets + initial_targets + [ read_config_file ]

    # need to merge otu table and metadata
    merge_otu_and_metadata = "merge_metadata.py " + metadata_file + " < " + \
        otu_table + " > " + initial_targets[0]
    
    # determine ending column of otu vs metadata in pcl file used to create
    # new project specific read.config file
    def create_read_config_for_project(metadata_file, targets):
        # open the metadata file to read in the column ids
        pcl_file=targets[3]
        read_config_file=targets[5]

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

    # for command line need to transpose first
    transpose = "transpose.py < " + initial_targets[0] + " > " + initial_targets[1] 
            
    run_maaslin = "Maaslin.R -i " + read_config_file + " " + \
        final_targets[0] + " " + initial_targets[1]
                
    return {
        "name": "maaslin: " + all_targets[0],
        "actions": [merge_otu_and_metadata,  
            (create_read_config_for_project,[metadata_file]),
            transpose, run_maaslin],
        "file_dep": [otu_table, metadata_file],
        "targets": all_targets
    }
