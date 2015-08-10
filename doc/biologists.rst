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
________________________________________________________

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

AnADAMA provides three interfaces to pipelines:

 * The directory skeleton (easiest for new users and those afraid of
   long shell commands).

 * The command line interface (more advanced, but only one command
   needed).

 * The python interface (advanced, intended for python-istas)


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
`AnADAMA's documentation <http://huttenhower.sph.harvard.edu/docs/anadama/your_own_pipeline.html#using-pipelines-via-the-command-line-interface>`_.


Pipeline Python Interface
_________________________

Pipelines are importable python classes. If writing python is your
thing, see :doc:`python-interface`.


Further pipeline documentation
______________________________

The pipelines are composed of many workflows, each having many
options.  The docs are your guide to understanding those
options. Below are the currently documented pipelines:

* :doc:`pipelines.sixteen`
* :doc:`pipelines.wgs`
* :doc:`pipelines.vis`


How to use workflows
====================

Unlike pipelines, there's only one interface to individual
workflows: python. If the default pipelines don't give you the
flexibility or features you want, see :doc:`python-interface`.


Demos
=====

*  16S Pipeline :ref:`sixteendemo`

*  WGS Pipeline :ref:`wgsdemo`

*  Python interface :doc:`python-interface`
