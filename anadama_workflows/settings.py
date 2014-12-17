class workflows:
    product_directory = "mibc_products"
    verbose           = True

    class metaphlan2:
        bowtie2db = "/aux/deploy2/var/lib/bowtie2db/mpa.200.ffn"
        mpa_pkl   = "/aux/deploy2/var/lib/bowtie2db/mpa.200.pkl"
    class sixteen:
        otu_taxonomy = "/aux/deploy2/var/lib/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt"
        otu_refseq   = "/aux/deploy2/var/lib/gg_13_5_otus/rep_set/97_otus.fasta"
    class alignment:
        kegg_bowtie2_db = "/aux/deploy2/var/lib/ccfa/KEGG_reduced_bowtie2/KEGG_1"
                                             
    class humann2:
        install_dir = "/aux/deploy2/var/lib/humann2/"
        uniref_path = install_dir+"/uniref/"
        chocophlan_path = install_dir+"/chocophlan/"
        _pathway_dir = install_dir+"/pathways"
        pathways_databases = (_pathway_dir+"/metacyc_reactions.uniref "+
                              _pathway_dir+"/metacyc_pathways")

    class knead:
        reference_db = "/aux/deploy2/var/lib/kneaddata/humanGRCh38"
