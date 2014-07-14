import os

def maaslin(otu_table, metadata_file):
	''' Workflow to compute the significance of association in microbial community
		using a transform abundance or relative function
		table obtained from Qiime, HUMAnN or MetaPhlAn plus study metadata
	'''
	
	# default contents of the read config file
	read_config_file_contents = "Matrix: Metadata\n" + \
		"Read_PCL_Rows: -REPLACE\n" + \
		"\n" + \
		"Matrix: Abundance\n" + \
		"Read_PCL_Rows: Bacteria-"	
	
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
	
	# need to merge otu table and metadata
	merge_otu_and_metadata = "merge_metadata.py " + metadata_file + " < " + \
		otu_table + " > " + initial_targets[0]
	
	# create default read.config file
	create_default_read_config = "echo $\'" + read_config_file_contents + "\' > " + read_config_file
	
	# determine ending column of otu vs metadata in pcl file used to create
	# new project specific read.config file
	modify_read_config_for_project = "sed -i s/REPLACE/$(grep Bacteria " + initial_targets[0] + \
		" -B 1 | head -n 1 | awk '{ print $1 }')/ " + read_config_file
	
	# for command line need to transpose first
	transpose = "transpose.py < " + initial_targets[0] + " > " + initial_targets[1] 
			
	run_maaslin =  "Maaslin.R -i " + read_config_file + " " + \
		final_targets[0] + " " + initial_targets[1]
				
	return {
		"name": "maaslin: " + targets[0],
		"actions": [merge_otu_and_metadata, create_default_read_config, 
			modify_read_config_for_project, transpose, run_maaslin],
		"file_dep": [otu_table, metadata_file],
		"targets": final_targets
	}
