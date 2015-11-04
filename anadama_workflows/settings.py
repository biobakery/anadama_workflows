base = "/net/hutlab11/srv/export/hutlab11/share_root/"

class workflows:
    product_directory = "anadama_products"
    verbose           = True

    class metaphlan2:
        bowtie2db = base+"databases/metaphlan2/db_v20/mpa_v20_m200"
        mpa_pkl   = base+"databases/metaphlan2/db_v20/mpa_v20_m200.pkl"
    class sixteen:
        otu_taxonomy = base+"databases/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt"
        otu_refseq   = base+"databases/gg_13_5_otus/rep_set/97_otus.fasta"
    class alignment:
        kegg_bowtie2_db = base+"databases/ccfa/KEGG_reduced_bowtie2/KEGG_1"
    class knead:
        reference_db = [ base+"databases/bowtie2/humanGRCh38"]
    class picrust:
        copy_number = base+"tools/anadama/picrust/tmp/picrust/picrust/data/16S_13_5_precalculated.tab.gz"
    class subread:
        index = "/vagrant/data/subread/GRCh38"
        annotations="/vagrant/data/subread/default_annotations.gtf"
    class usearch:
        chimera_gold_standard = base+"databases/uchime/gold.fa"
        otu_db = base+"databases/usearch/gg_13_5_97-otus.udb"
