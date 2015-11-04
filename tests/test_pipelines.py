from anadama_workflows import pipelines

import subprocess
import tempfile
import os
import filecmp
import shutil
import sys

def eq(a, b):
    assert a == b

def deepeq(answers, result):
    assert (a==b for a, b in zip(answers, result))

def file_equal(file1,file2):
    assert filecmp.cmp(file1, file2, shallow=False)
    
def read_file(file):
    file_data=[]
    with open(file) as file_handle:
            file_data=[line for line in file_handle]
    return file_data

def read_file_precision(file,precision):
    """ Read in the file converting floats to the given precision """
    file_data=[]
    with open(file) as file_handle:
        for line in file_handle:
            formatted_data=[]
            for data in line.rstrip().split("\t"):
                try:
                    data_float=float(data)
                    data="{:.{digits}f}".format(data_float, digits=precision)
                except ValueError:
                    pass
                formatted_data.append(data)
            file_data.append("\t".join(formatted_data))
    return file_data

def file_equal_ignore_line_order(file1,file2):
    file1_data=read_file(file1)
    file2_data=read_file(file2)

    file1_data=sorted(file1_data)
    file2_data=sorted(file2_data)

    assert file1_data == file2_data

def file_almost_equal_ignore_line_order(file1,file2):
    """ Check if the files are almost equal to a precision for floats """
    precision=6
    
    file1_data=read_file_precision(file1,precision)
    file2_data=read_file_precision(file2,precision)

    file1_data=sorted(file1_data)
    file2_data=sorted(file2_data)

    assert file1_data == file2_data

def remove_temp_folder(folder):
    try:
        shutil.rmtree(folder)
    except EnvironmentError:
        print("Warning: Unable to remove temp folder: " +folder)

def run_anadama_pipeline_command(command,temp_directory):
    # get the current directory
    current_working_directory=os.getcwd()
    
    # move to the temp directory
    os.chdir(temp_directory)
    
    # add the anadama portion to the command
    command=["anadama","pipeline"]+command
    
    # run the pipeline
    try:
        print("Testing: "+" ".join(command))
        return_code=subprocess.call(command, stderr = subprocess.STDOUT)
    except (EnvironmentError, subprocess.CalledProcessError):
        print("Warning: Unable to execute anadama workflows pipeline in test.")
            
    if return_code != 0:
        print("Warning: Error running anadama workflows pipeline in test.")
    else:
        print("Done.")

    # move back to original directory
    os.chdir(current_working_directory)
    
    output_folder=os.path.join(temp_directory,"anadama_products")
    
    return output_folder

def list_subdirectories(folder):
    """ Return all of the subdirectories in a folder """
    folder=os.path.abspath(folder)
    subdirs=[]
    for file in os.listdir(folder):
        path_to_file=os.path.join(folder,file)
        if os.path.isdir(path_to_file):
            subdirs.append(path_to_file)
    return sorted(subdirs)

def pair_tsv_files(folder1,folder2,subdirectories=None):
    """ Yield pairs of tsv files from the two folders 
    Only search subdirectories if set """
    
    if not subdirectories is None:
        # use the first matching subdirectory
        folder1=list_subdirectories(folder1)[0]
        folder2=list_subdirectories(folder2)[0]
    
    for file in os.listdir(folder1):
        if file.endswith(".tsv"):
            file1=os.path.join(folder1,file)
            file2=os.path.join(folder2,file)
            yield file1, file2

def data_folder():
    """ Get the full path to the tests data folder """
    return os.path.join(os.path.dirname(os.path.abspath(__file__)),"data")

def test__regex_filter_pair():
    """ Test grouping single pair """
    a = ['1.fastq', '2.fastq']
    result = pipelines._regex_filter(a)
    yield deepeq, ([a[0]], [a[1]], []), result

def test__regex_filter_pair_diff_case():
    """ Test grouping single pair with different case"""
    a = ['r1.fastq', 'R2.fastq']
    result = pipelines._regex_filter(a)
    yield deepeq, ([a[0]], [a[1]], []), result

def test__regex_filter_six_files_one_pair():
    """ Test grouping six files representing one pair """
    a = ['sample_a.r1.fastq', 'sample_b.r2.fastq', 'sample_c.r1.fastq',
         'sample_b.r1.fastq', 'sample_c.r2.fastq', 'sample_a.r2.fastq']
    result = pipelines._regex_filter(a)
    yield deepeq, ([a[0], a[3], a[2]], [a[1], a[5], a[4]], []), result

def test__regex_filter_six_files_one_pair_diff_notation():
    """ Test grouping six files representing one pair different notation"""
    a = ['sample.a.r1.fastq', 'sample_b.R2.fastq', 'sample.c.r1.fastq',
         'sample-b.r1.fastq', 'sample_c.r2.fastq', 'sample-a.R2.fastq']
    result = pipelines._regex_filter(a)
    yield deepeq, ([a[0], a[3], a[2]], [a[1], a[5], a[4]], []), result

def test__to_merged_named():
    """ Test named merge tagging """
    yield eq, pipelines._to_merged("koji_R1.fastq", "foo"), "koji_foo.fastq"

def test__to_merged_unnamed():
    """ Test unnamed merge tagging """
    yield eq, pipelines._to_merged("koji_R1.fastq"), "koji_merged.fastq"

def test__to_merged_named_diff_notation():
    """ Test named merge tagging different notation """
    yield eq, pipelines._to_merged("kojiR1.fastq", "foo"), "kojifoo.fastq"

def test__to_merged_unamed_diff_default():
    """ Test unnamed merge tagging different default """
    yield eq, pipelines._to_merged("koji1.fastq"), "kojimerged.fastq"

def test__to_merged_stiched():
    """ Test stitch tagging """
    yield eq, pipelines._to_merged("2.fastq", "stitched"), "stitched.fastq"
 
def test_demultiplexed_usearch64_16S():
    """ Test the usearch64 bit 16S pipeline on a set of demultiplexed samples
    that are qiime fasta formatted """
    
    # create a temp directory to store output
    temp_directory=tempfile.mkdtemp(prefix="anadama_workflows_test_usearch64")
    
    input_files=os.path.join(data_folder(),"16S_demultiplexed","*.fasta")
    
    # run the anadama command
    command=["Usearch64_16S","-f","demuxed_fasta_files:glob:"+input_files]
    output_folder=run_anadama_pipeline_command(command, temp_directory)
        
    # test the output files are as expected
    # only check the main output tsv files
    benchmark_output_folder=os.path.join(data_folder(),"16S_demultiplexed","anadama_products_Usearch64_16S")
    for expected_result_file, actual_result_file in pair_tsv_files(benchmark_output_folder, output_folder):
        yield file_equal_ignore_line_order, actual_result_file, expected_result_file
    
    # remove the temp folder
    remove_temp_folder(temp_directory)
    
def test_raw_wgs():
    """ Test the wgs pipeline on a set of raw fastq samples """
    
    # create a temp directory to store output
    temp_directory=tempfile.mkdtemp(prefix="anadama_workflows_test_wgs")
    
    input_files=os.path.join(data_folder(),"wgs_raw","*.fastq")
    
    # run the anadama command
    command=["WGS","-f","raw_seq_files:glob:"+input_files,"-o","decontaminate.strategy:storage",
        "-o","humann.bypass-translated-search:"]
    output_folder=run_anadama_pipeline_command(command, temp_directory)
        
    # test the output files are as expected
    # only check the main output tsv files
    benchmark_output_folder=os.path.join(data_folder(),"wgs_raw","anadama_products_wgs_raw")
    for expected_result_file, actual_result_file in pair_tsv_files(benchmark_output_folder, output_folder,subdirectories=True):
        yield file_almost_equal_ignore_line_order, actual_result_file, expected_result_file
    
    # remove the temp folder
    remove_temp_folder(temp_directory)
    
