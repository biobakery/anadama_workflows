Whole Genome Shotgun Pipelines
##############################

.. contents:: 
   :local:
.. currentmodule:: anadama_workflows.pipelines.wgs

.. automodule:: anadama_workflows.pipelines.wgs
   :members:
   :undoc-members:
   :special-members:


.. _wgsdemo:

Demo
####

Using the directory skeleton
____________________________

::

   source ~/anadama_dev/bin/activate # Do this if you have not done so already
   mkdir docdemo
   cd docdemo
   wget 'http://huttenhower.sph.harvard.edu/biobakery-shop/anadama/wgs_functional_test.tgz'
   tar -xvzf wgs_functional_test.tgz
   cd wgs_functional_test/
   anadama skeleton anadama_workflows.pipelines:WGSPipeline
   mv -iv *.fastq input/raw_seq_files/
   mv -iv map.txt input/sample_metadata/
   echo 'name: humann' > _skip.txt
   echo 'threads: 8' >> input/_options/decontaminate.txt
   echo 'nproc: 8' >> input/_options/metaphlan2.txt
   anadama run


Using the command line interface
________________________________

::

   mkdir docdemo
   cd docdemo
   wget 'http://huttenhower.sph.harvard.edu/biobakery-shop/anadama/wgs_functional_test.tgz'
   tar -xvzf wgs_functional_test.tgz
   cd wgs_functional_test/
   anadama pipeline anadama_workflows.pipelines:WGSPipeline \
       -f 'raw_seq_files: glob:*.fastq' \
       -o 'decontaminate.threads: 8' \
       -o 'metaphlan2.nproc: 8' \
       -k 'name: humann'

