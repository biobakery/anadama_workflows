class workflows:
    product_directory = "mibc_products"
    verbose           = True

    class metaphlan2:
        bowtie2db = "/vagrant/scripts/metaphlan2/db_v20/mpa_v20_m200"
        mpa_pkl   = "/vagrant/scripts/metaphlan2/db_v20/mpa_v20_m200.pkl"
    class sixteen:
        otu_taxonomy = "/vagrant/data/ccfa/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt"
        otu_refseq   = "/vagrant/data/ccfa/gg_13_5_otus/rep_set/97_otus.fasta"
    class alignment:
        kegg_bowtie2_db = "/vagrant/data/ccfa/ccfa/KEGG_reduced_bowtie2/KEGG_1"
                                             
    class humann2:
        install_dir = "/vagrant/scripts/humann2/"
        uniref_path = install_dir+"/uniref/"
        chocophlan_path = install_dir+"/chocophlan/"
        _pathway_dir = install_dir+"/pathways"
        pathways_databases = (_pathway_dir+"/metacyc_reactions.uniref,"+
                              _pathway_dir+"/metacyc_pathways")

    class knead:
        reference_db = "/vagrant/scripts/kneaddata/databases/bowtie2/humanGRCh38"
        trim_path = "/vagrant/scripts/Trimmomatic-0.33/trimmomatic-0.33.jar"

    class subread:
        index = "/vagrant/data/subread/GRCh38"
        annotations="/vagrant/data/subread/default_annotations.gtf"
