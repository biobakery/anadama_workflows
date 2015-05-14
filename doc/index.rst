.. anadama_workflows documentation master file, created by
   sphinx-quickstart on Wed Jul 16 13:41:24 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to anadama_workflows's documentation!
=============================================

Contents:

.. toctree::
   :maxdepth: 2
   
   workflows
   pipelines


Installation Guide
==================

Manually
________

Start with the latest version of python2_, virtualenv_, and
setuptools_.

You'll also need a decent development environment. Refer to the
documentation for your distribution, but mostly it's something like
``apt-get install build-essential`` or ``sudo yum groupinstall
'Development Tools'``. Make sure you have ``git`` and ``hg``, too.

Sadly, you'll also need a JVM. Download and install the Oracle JVM
v1.7.

Start off with a virtualenv::

  virtualenv ~/anadama_env
  cd ~/anadama_env
  source bin/activate

Now we start filling out our ``src`` directory with anadama and
anadama_workflows::

  pip install -e 'git+https://bitbucket.org/biobakery/anadama.git@master#egg=anadama-0.0.1'

  pip install -e 'git+https://bitbucket.org/biobakery/anadama_workflows.git@master#egg=anadama_workflows-0.0.1'
  
  
That's all we need to drive things. Next, we'll install the
dependencies needed for AnADAMA's ``SixteenSPipeline`` and
``WGSPipeline``.


For these steps, let's define some bash functions to ease our
installation tasks::

  function link() { ln -sv $(readlink -f "$1") ~/anadama_env/bin/; }

  function download_unpack() { wget -O- "$1" | tar -xvzf - ; }

Let's start downloading some tools into our ``src`` directory. First
is bowtie2::

  cd src
  rm pip-delete-this-directory.txt
  # get this link by going to the bowtie project on sourceforge
  wget -O bowtie2.zip 'http://downloads.sourceforge.net/project/bowtie-bio/bowtie2/2.2.5/bowtie2-2.2.5-linux-x86_64.zip?r=http%3A%2F%2Fsourceforge.net%2Fprojects%2Fbowtie-bio%2Ffiles%2Fbowtie2%2F2.2.5%2F&ts=1431627195&use_mirror=softlayer-dal'
  
  unzip bowtie2.zip
  link bowtie2*/bowtie2
  link bowtie2*/bowtie2-build
  rm bowtie2.zip

Ok, now DIAMOND::

  download_unpack http://github.com/bbuchfink/diamond/releases/download/v0.7.9/diamond-linux64.tar.gz
  link diamond
  

Next up is metaphlan2::

  download_unpack https://bitbucket.org/biobakery/metaphlan2/get/tip.tgz
  link biobakery-metaphlan2*/metaphlan2.py
  link biobakery-metaphlan2*/db_v20

Now we do humann2::

  download_unpack https://bitbucket.org/biobakery/humann2/downloads/humann2_v0.1.9.tar.gz
  pip install -e humann2_v0.1.9
  cd humann2_v0.1.9
  python setup.py minpath
  
  # Now we need somewhere to put the databases
  cd -
  mkdir ../databases
  humann2_databases --download chocophlan full ../databases
  humann2_databases --download uniref diamond ../databases

And now we do KNEAD_data::

  download_unpack https://bitbucket.org/biobakery/kneaddata/get/v0.3.tar.gz
  pip install -e biobakery-kneaddata-*/

Next, we'll need Breadcrumbs::
  
  download_unpack https://bitbucket.org/biobakery/breadcrumbs/get/ed59079c2e5e.tgz
  pip install -e biobakery-breadcrumbs-*/


That's all for the WGS Pipeline.


Proceeding on to the installation dependencies ``SixteenSPipeline``.
First up is ea-utils::

  # Go to https://code.google.com/p/ea-utils/
  # look to the right side of the page
  # near the bottom
  # click on "downloads and binaries"
  # download ea-utils.1.1.2-806.tar.gz

  tar -xvzf ea-utils.1.1.2-806.tar.gz
  cd ea-utils.1.1.2-806
  make
  link fastq-join
  cd -

Now we install qiime. Qiime requires python packages (like numpy) that
have conflicting versions with some other packages we recently
installed (breadcrumbs). We'll get around this problem by using a
script called ``docent``. We first install docent, then we install
qiime::

  download_unpack https://bitbucket.org/biobakery/docent/get/HEAD.tgz
  pip install -e biobakery-docent-*/

  download_unpack 'https://github.com/biocore/qiime/archive/1.8.0.tar.gz'
  cd ~/anadama_env
  docent -v \
      -i 'numpy==1.7.1' -i' -e "~/anadama_env/src/qiime-1.8.0/"' \
      -j ~/anadama_env/src/biobakery-docent*/specs/qiime.json \
      -e ~/anadama_env/src/qiime-1.8.0/env

  # Need to download some databases
  cd ~/anadama_env/databases/
  download_unpack 'ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5_otus.tar.gz'


Finally, install picrust::

  cd ~/anadama_env/src
  download_unpack 'https://github.com/picrust/picrust/releases/download/1.0.0/picrust-1.0.0.tar.gz'
  cd ../
  docent -v -a '-p /usr/local/bin/python '\
      -i "numpy==1.5.1" -i "biom-format==1.3.1" -i "cogent==1.5.3" \
      -i '-e "~/anadama_env/src/picrust-1.0.0"'  \
      -j ~/anadama_env/src/biobakery-docent*/specs/picrust.json \
      -e ~/anadama_env/src/picrust-1.0.0/env

  # Finally, we install the last databases, but in an unusual spot
  mkdir -pv ~/anadama_env/src/picrust-1.0.0/picrust/data
  cd !$
  wget 'ftp://ftp.microbio.me/pub/picrust-references/picrust-1.0.0/16S_13_5_precalculated.tab.gz'
  wget 'ftp://ftp.microbio.me/pub/picrust-references/picrust-1.0.0/ko_13_5_precalculated.tab.gz'


.. _setuptools: https://pypi.python.org/pypi/setuptools
.. _python2: https://www.python.org/downloads/
.. _virtualenv: https://pypi.python.org/pypi/virtualenv


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

