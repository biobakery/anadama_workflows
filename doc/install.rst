##################
Installation Guide
##################

.. _installation-guide:

Manually
________

Currently, the best way to install ``anadama_workflows`` and all of
its dependencies is by hand. It's not daunting and shouldn't take
longer than 20 minutes.

The process generally proceeds as follows:

  * Install Python, virtualenv, and the anadama drivers.

  * Install the prerequisites for any pipeline or type of workflows
    you plan to use.

  * Update the ``anadama_workflows`` settings file



General Prerequisites
---------------------

Start with the latest version of python2_, virtualenv_, and
setuptools_ (click the links to go to their respective download
pages).

You'll also need a decent development environment. To install a
development environment, please refer to your distribution's
documentation. Likely, it's something akin to running ``apt-get
install build-essential`` or ``sudo yum groupinstall 'Development
Tools'``. Make sure you have ``git`` and ``hg`` (aka mercurial)
installed, too.

Sadly, you'll also need a JVM. Download and install the Oracle JVM
v1.7.

It's best to install the dependencies required by
``anadama_workflows`` in a virtualenv_.

Here's how we'll install a virtualenv::

  virtualenv ~/anadama_env
  cd ~/anadama_env
  source bin/activate

Next, we'll install the python programs and libraries that orchestrate
the workflows. That's just two commands::

  pip install -e 'git+https://bitbucket.org/biobakery/anadama.git@master#egg=anadama-0.0.1'

  pip install -e 'git+https://bitbucket.org/biobakery/anadama_workflows.git@master#egg=anadama_workflows-0.0.1'
  
  
That's all we need to drive things. Next, we'll install the
dependencies needed for AnADAMA's ``WGSPipeline`` and
``SixteenSPipeline``.


WGS Pipeline
------------

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

  download_unpack https://bitbucket.org/biobakery/kneaddata/get/9f2bc13440e5.tar.gz
  pip install -e biobakery-kneaddata-*/

  wget -O trimmomatic.zip 'http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip'
  unzip trimmomatic.zip

  # Go download jre1.7 from http://www.oracle.com/technetwork/java/javase/downloads/jre7-downloads-1880261.html
  tar -xvzf jre-7u80-linux-x64.gz
  link jre1.7.0_80/bin/java

  # And we'll also need a human reference database
  mkdir -pv ~/anadama_env/databases/bowtie2
  cd ~/anadama_env/databases/bowtie2
  download_unpack 'ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz'
  for file in *; do mv -iv $file $( echo "$file" | sed -e 's|GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index|humanGRCh38|g' ); done


Next, we'll need Breadcrumbs. Breadcrumbs requires python packages
(like numpy) that have conflicting versions with some other packages
you might install later (qiime). We'll get around this problem by using
a script called ``docent``. We first install docent, then we install
breadcrumbs::

  cd ~/anadama_env/src/
  download_unpack https://bitbucket.org/biobakery/docent/get/HEAD.tgz
  pip install -e biobakery-docent-*/
  
  download_unpack https://bitbucket.org/biobakery/breadcrumbs/get/ed59079c2e5e.tgz
  cd ~/anadama_env
  docent -e  ~/anadama_env/src/biobakery-breadcrumbs-ed59079c2e5e/env -v \
      -i "-e ~/anadama_env/src/biobakery-breadcrumbs-ed59079c2e5e" \
      -o bin/scriptConvertBetweenBIOMAndPCL.py \
      -o bin/scriptEnvToTable.py \
      -o bin/scriptManipulateTable.py \
      -o bin/scriptPcoa.py \
      -o bin/scriptPlotFeature.py


That's all for the WGS Pipeline.


16S Pipeline
------------

For these steps, let's define some bash functions to ease our
installation tasks::

  function link() { ln -sv $(readlink -f "$1") ~/anadama_env/bin/; }

  function download_unpack() { wget -O- "$1" | tar -xvzf - ; }


Start in the ``src`` directory of your virtualenv::

  cd src


Some software that we'll install requires python packages that
conflict with other 

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

And Biom-format::

  pip install biom-format==1.3.1


Next, we'll need Qiime. Qiime requires python packages (like numpy)
that have conflicting versions with some other packages you might
install (Breadrumbs). We'll get around this problem by using a script
called ``docent``.We first install docent, then we install qiime. If
you've already installed ``docent``, you can skip that step::

  # The docent install step. skip if you've already installed docent
  cd ~/anadama_env/src
  download_unpack https://bitbucket.org/biobakery/docent/get/HEAD.tgz
  pip install -e biobakery-docent-*/
  
  # regardless of installing docent, run this step to install qiime
  download_unpack 'https://github.com/biocore/qiime/archive/1.8.0.tar.gz'
  cd ~/anadama_env
  docent -v \
      -i 'numpy==1.7.1' -i' -e "~/anadama_env/src/qiime-1.8.0/"' \
      -j ~/anadama_env/src/biobakery-docent*/specs/qiime.json \
      -e ~/anadama_env/src/qiime-1.8.0/env

  # Need to download some databases
  mkdir -pv ~/anadama_env/databases/
  cd ~/anadama_env/databases/
  download_unpack 'ftp://greengenes.microbio.me/greengenes_release/gg_13_5/gg_13_5_otus.tar.gz'


Finally, install picrust::

  cd ~/anadama_env/src
  download_unpack 'https://github.com/picrust/picrust/releases/download/1.0.0/picrust-1.0.0.tar.gz'
  cd ../
  docent -v \
      -i "numpy==1.5.1" -i "biom-format==1.3.1" -i "cogent==1.5.3" \
      -i "-e ~/anadama_env/src/picrust-1.0.0"  \
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


Editing the anadama_workflows settings.py file
______________________________________________

Change the file locations to where you've installed them. Like so::

  # edit ~/anadama_env/src/anadama-workflows/anadama_workflows/settings.py

  class metaphlan2:
      bowtie2db = "/home/user/anadama_env/bin/db_v20/mpa_v20_m200"
      mpa_pkl   = "/home/user/anadama_env/bin/db_v20/mpa_v20_m200.pkl"
  class sixteen:
      otu_taxonomy = "/home/user/anadama_env/databases/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt"
      otu_refseq   = "/home/user/anadama_env/databases/gg_13_5_otus/rep_set/97_otus.fasta"

  class knead:
      reference_db = "/home/user/anadama_env/databases/bowtie2/humanGRCh38"
      trim_path = "/home/user/anadama_env/src/Trimmomatic-0.33/trimmomatic-0.33.jar"



Tricks for Debian 8
___________________

You'll need to install a few packages::

  sudo apt-get update
  sudo apt-get install build-essential \
      git mercurial \
      virtualenv python-pip python-dev python-pip \
      zlib1g-dev zlib1g unzip zip libbz2 libbz2-dev \
      libglpk-dev libglpk36 gfortran \
      swig \
      libfreetype6-dev libfreetype6 libpng12-0

Matplotlib can't find /usr/include/freetype2/ft2build.h, but it can
find /usr/include/ft2build.h, so link it up::

  sudo ln -sv /usr/include/freetype2/ft2build.h /usr/include/ft2build.h


Tricks for Ubuntu 14.04
_______________________

A few packages you must install::

  sudo apt-get update
  sudo apt-get upgrade
  sudo apt-get install build-essential \
      git mercurial python-virtualenv \
      zlib1g-dev unzip zip libbz2-dev \
      libglpk-dev libglpk36 gfortran \
      swig \
      libfreetype6-dev libfreetype6 libpng12-0

And, just like in debian 8, libfreetype2's headers can't be
found. Symlink them like so::

  sudo ln -sv /usr/include/freetype2/ft2build.h /usr/include/ft2build.h


