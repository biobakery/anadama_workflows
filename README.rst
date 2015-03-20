AnADAMA workflows
#################

*Pre-built workflows and pipelines for use with* AnADAMA_.

.. _AnADAMA: https://bitbucket.org/biobakery/anadama

.. contents:: 

_________________________________________________________________________________


Overview
========

``anadama_workflows`` is simply a container for workflows and
pipelines ready for use with DoIt_ or AnADAMA_.

.. _DoIt: http://pydoit.org/


Installation
============

One liner::

  $ pip install -e 'git+https://bitbucket.org/biobakery/anadama_workflows.git@master#egg=anadama_workflows-0.0.1'


Documentation
============

Head on over to the docs_ site.

.. _docs: http://rschwager-hsph.bitbucket.org/documentation/anadama-workflows/index.html


How to get help
============

Submit an issue_.

.. _issue: https://bitbucket.org/biobakery/anadama_workflows/issues


Workflows
=========

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
		  

Pipelines
=========

Pipelines are collections of workflows. Workflows are reusable pieces,
but sometimes they're put together in predictable ways. For cases
where it's useful to get many data files of varying stages of analysis
all to a single end-point, use a pipeline.

Here's an example of using a pipeline within a DoIt task:

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


Adding your own workflows to anadama_workflows
==============================================

If you have a workflow that you think other people will use, submit a
pull request from your own fork of this repository.


Adding workflows
----------------

A workflow is just a python function that returns a dictionary, or a
generator that produces dictionaries. It's OK to nest generators
within generators.

The dictionaries produced by a workflow should conform to the
interface_ set forth by DoIt.

.. _interface: http://pydoit.org/tasks.html

The workflow should be documented by a python docstring.

For examples, don't be afraid to read the source!


Adding pipelines
----------------

A pipeline in the ``anadama_workflows`` sense of the term is a
subclass of `anadama.pipelines.Pipeline`_.

.. _anadama.pipelines.Pipeline: https://bitbucket.org/biobakery/anadama/src/a25a3953d962054fb5daef759807bba979ef2c56/anadama/pipelines.py?at=master

Read the ``anadama`` source for development guidelines and the
``anadama_workflows`` source for examples.