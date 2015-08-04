Workflows
#########

:Release: |version|
:Date: |today|

.. toctree::
   :maxdepth: 2

   association
   alignment
   biom
   general
   sixteen
   usearch
   visualization
   wgs


What's a workflow?
==================

Workflows are functions that, when called, return a dictionary with
attributes DoIt recognizes as a task. Workflows are intended to be an
atomic single step of work with (in general) a list of input files and
a list of output files. Workflow functions can be used within your
DoIt tasks.

Here's an example of using workflow functions within a DoIt Task:

.. code:: python

	  import anadama_workflows as workflows
	  fastq_files = glob.glob("*.fastq")

	  def task_metagenome_profile():
	      """perform metagenomic profiling on all fastq files"""

	      alignment_files = list()
	      for fastq_file in fastq_files:
	          yield workflows.metaphlan2([fastq_file])
		  alignment_file = fastq_file+".sam"
		  yield workflows.bowtie2_align([fastq_file],
		                                alignment_file)
	          alignment_files.append(alignment_file)
              yield workflows.humann([alignment_files], workdir='./humann')
		  
