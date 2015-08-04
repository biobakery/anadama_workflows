Pipelines
#########

:Release: |version|
:Date: |today|

.. automodule::
   :members: all


.. toctree::
   :maxdepth: 2

   pipelines.sixteen
   pipelines.vis
   pipelines.wgs


What's a pipeline?
==================

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

