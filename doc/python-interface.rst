#################
Python interfaces
#################


A word about dodo.py files
==========================

Spoiler alert: workflows are just python function that return
dictionaries understood by doit and AnADAMA. That means that using a
python interpreter or a python script to simply call the workflow
function won't actually produce the expected outputs. It's doit's or
AnADAMA's job to understand the dictionary created by the workflow
function and actually produce the expected products. How do we give
these dictionaries to AnADAMA or doit? Use a dodo.py file.

Although clever programming may lead you to believe otherwise, doit
and AnADAMA's core functionality is to read a python script, execute
the commands in the script to get some dictionaries, then finally use
the dictionaries to build the expected analysis products. The default
name of the script that AnADAMA or doit reads is ``dodo.py``.

You can write your own ``dodo.py`` script. Start by creating a text
file with a text editor of your choice (text editors are programs like
nano or Sublime text, NOT word processors like Microsoft Word, Google
Docs or Openoffice). The examples below show what goes into this text
file. For more information on writing dodo.py files, please see this
presentation_ or the doit docs_.

.. _presentation: http://rschwager-hsph.bitbucket.org/2014-07-11_lab-presentation/index.html#/3 

.. _docs: http://pydoit.org/



Pipelines
=========

Here's an example of using a pipeline in a dodo.py file:

.. code:: python

  import os
  from anadama_workflows import pipelines

  all_files = os.listdir('./raw_data_files/')
  
  def task_16S_analysis():
      """Use the 16S pipeline from anadama_workflows to get some OTU
      tables from some samples"""

      # Instantiate the pipeline
      my_pipeline = pipelines.SixteenSPipeline(
          raw_seq_files = all_files,
          sample_metadata = "./map.txt",
      )
      # have it build up the tasks necessary to finish processing
      my_pipeline.configure()

      # give the tasks over to DoIt
      yield my_pipeline.task_dicts


Workflows
=========

The example dodo.py file below uses the metaphlan2, bowtie2_align, and
humann workflows:


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

Then save the file and run by executing ``anadama run -f dodo.py``.
This is all python code; if you're new to Python, consider reading the
Python tutorial_ or the `code academy tutorial`_.

.. _tutorial: https://docs.python.org/2/tutorial/index.html

.. _`code academy tutorial`: https://www.codecademy.com/en/tracks/python
