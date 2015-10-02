from setuptools import setup, find_packages

setup(
    name='anadama_workflows',
    version='0.0.1',
    description='Hutlab and associated workflows for AnADAMA',
    packages=find_packages(exclude=['ez_setup', 'tests', 'tests.*']),
    zip_safe=False,
    install_requires=[
        'biopython>=1.63',
        'pysam==0.7.8',
        # anadama should also install doit
        'anadama'
    ],
    dependency_links=[
        'git+https://bitbucket.org/biobakery/anadama.git@master#egg=anadama-0.0.1', 
    ],
    classifiers=[
        "Development Status :: 2 - Pre-Alpha"
    ],
    entry_points= {
        'console_scripts': [
            'sequence_convert = anadama_workflows.utility_scripts.convert:main',
            'fastq_split      = anadama_workflows.utility_scripts.seqsplit:main',
            'sequence_sort    = anadama_workflows.utility_scripts.sort:main',
            'sequence_pair    = anadama_workflows.utility_scripts.pair:main',
            'bam_pe_split     = anadama_workflows.utility_scripts.bam_pe_split:main',
            'sequence_re-pair = anadama_workflows.utility_scripts.re_pair:main',
            'uclust_otutable = anadama_workflows.utility_scripts.uclust_otutable:main'
        ],
        'anadama.pipeline': [
            '.SixteenS = anadama_workflows.pipelines.sixteen:SixteenSPipeline',
            '.Usearch16S = anadama_workflows.pipelines.usearch:Usearch16SPipeline',
            '.WGS = anadama_workflows.pipelines.wgs:WGSPipeline',
            '.RNA = anadama_workflows.pipelines.rna:RNAPipeline',
            '.Visualization = anadama_workflows.pipelines.vis:VisualizationPipeline',
        ],
    }
)
