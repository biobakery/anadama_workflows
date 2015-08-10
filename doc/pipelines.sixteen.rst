16S Pipelines
#############

.. contents:: 
   :local:
.. currentmodule:: anadama_workflows.pipelines.sixteen

.. automodule:: anadama_workflows.pipelines.sixteen
   :members:
   :undoc-members:
   :special-members:


.. _sixteendemo:


Demo
####

Via the directory skeleton
==========================

::

   mkdir docdemo 
   cd docdemo  
   wget 'http://huttenhower.sph.harvard.edu/biobakery-shop/anadama/16s_functional_test.tgz' 
   tar -xvzf 16s_functional_test.tgz 
   cd 16s_functional_test/ 
   anadama skeleton anadama_workflows.pipelines:SixteenSPipeline 
   mv -iv Fasting_Example.sff input/raw_seq_files/ 
   mv -iv map.txt input/sample_metadata/ 
   anadama run 


Via the command line interface
==============================

::

   mkdir docdemo 
   cd docdemo  
   wget 'http://huttenhower.sph.harvard.edu/biobakery-shop/anadama/16s_functional_test.tgz' 
   tar -xvzf 16s_functional_test.tgz 
   cd 16s_functional_test/ 
   anadama pipeline anadama_workflows.pipelines:SixteenSPipeline \
       -f 'raw_seq_files: Fasting_Example.sff' \
       -f'sample_metadata: map.txt' 
