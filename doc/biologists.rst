########################
For Users and Biologists
########################

Introduction
============

AnADAMA is a tool that orchestrates analysis tasks.
``anadama_workflows`` is a collection of routines that plug in to
AnADAMA to perform specific tasks like running MetaPhlAn for taxonomic
profiling of WGS sequences. The routines in ``anadama_workflows``
often rely on third party tools being installed. Think of this
ecosystem like a puppet show: the puppeteer is AnADAMA, the brackets
and strings are in ``anadama_workflows``, the third party tools
are the puppets, and the director is the user (you!).

AnADAMA knows and works with workflows and
pipelines. ``anadama_workflows`` contains a healthy number of
workflows and also contains a handful of pipelines. AnADAMA also
exposes a number of interfaces to use workflows and pipelines. The
sections below detail how to use workflows and pipelines.


What's the difference between a workflow and a pipeline?
========================================================

A workflow is AnADAMA's smallest reproduceable unit of work. Workflows
make tasks, AnADAMA executes the tasks, moving the overall project or
analysis one step closer to completion. Workflows are like one-dish
recipies: the chef reads the recipe and uses tools like a pan and heat
to move a meal one step closer to completion. Just as a recipe isn't
the dish, a workflow is not a task.

A pipeline is AnADAMA's collection of workflows. Think of a pipeline
like a stack of recipies to make a breakfast: eggs benedict,
strawberry and green salad, coffee, and freshly pressed
juice. AnADAMA, like a good chef, knows what needs to happen
first. Also, if AnADAMA is interrupted, it will pick up where it left
off and only redo what is necessary.


Installation
============

Please refer to :ref:`installation-guide`.


How to use Pipelines
====================

AnADAMA provides three interfaces to pipelines: the directory skeleton
interface, the command-line interface, and the python interface.

Pipeline Directory Skeleton
___________________________

AnADAMA has a ``anadama skeleton`` command. This command creates a skeleton of
input directories for a pipeline and a working dodo.py file. Place
your input files into the appropriate directories, then run ``anadama
run`` or ``doit run``. Workflow options can be changed by editing the
files under the ``input/_options`` directory.

For more information on the ``skeleton`` command, please refer to
`AnADAMA's skeleton documentation <http://rschwager-hsph.bitbucket.org/documentation/anadama/your_own_pipeline.html#using-pipelines-via-the-directory-skeleton>`_.


Pipeline Command Line Interface
_______________________________

AnADAMA also lets you put all options and paramters on the command
line and execute everything with one (possibly big) command. Use the
``anadama pipeline`` command.

For more information on how to use the ``pipeline`` command, see
`AnADAMA's documentation <http://rschwager-hsph.bitbucket.org/documentation/anadama/your_own_pipeline.html#using-pipelines-via-the-command-line-interface>`_.


Pipeline Python Interface
_________________________

Pipelines are importable python classes.
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




How to use workflows
====================

Unlike pipelines, there's only one interface to individual
workflows: python. Start by creating your own dodo.py file with a text
editor of your choice (text editors are programs like nano or Sublime
text, NOT word processors like Microsoft Word, Google Docs or
Openoffice) and typing the following:

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
